/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/microstructure/book.q
/ description: limit order book reconstruction, queue dynamics, and micro-price
/ ==============================================================================

\d .book

/ Initialize empty limit order book state
/ Represented as dictionary containing bids and asks tables, and order lookup map
initBook: {[sym]
  bids: ([price: `float$()] size: `float$(); orderCount: `long$());
  asks: ([price: `float$()] size: `float$(); orderCount: `long$());
  orders: ([orderId: `long$()] side: `symbol$(); price: `float$(); size: `float$());
  `sym`bids`asks`orders!(sym; bids; asks; orders)
 };

/ Insert limit order into book: updates resting queue and depth table
addOrder: {[book; orderId; side; price; size]
  if[size <= 0f; :book];
  
  / Record order in lookup map
  book[`orders]: book[`orders], (enlist orderId)!(enlist side; enlist price; enlist size);
  
  if[side == `buy;
    curr_size: $[price in key book`bids; book[`bids; price; `size]; 0f];
    curr_cnt: $[price in key book`bids; book[`bids; price; `orderCount]; 0];
    book[`bids]: book[`bids], (enlist price)!(enlist (curr_size + size); enlist (curr_cnt + 1))
  ;
    curr_size: $[price in key book`asks; book[`asks; price; `size]; 0f];
    curr_cnt: $[price in key book`asks; book[`asks; price; `orderCount]; 0];
    book[`asks]: book[`asks], (enlist price)!(enlist (curr_size + size); enlist (curr_cnt + 1))
  ];
  book
 };

/ Cancel resting order from book: reduces size or removes empty price level
cancelOrder: {[book; orderId; cancelSize]
  if[not orderId in key book`orders; :book];
  ord: book[`orders; orderId];
  side: ord`side;
  price: ord`price;
  curr_ord_size: ord`size;
  
  actual_cancel: min (cancelSize; curr_ord_size);
  new_ord_size: curr_ord_size - actual_cancel;
  
  if[new_ord_size <= 0f;
    book[`orders]: delete from book[`orders] where orderId = orderId
  ;
    book[`orders; orderId; `size]: new_ord_size
  ];
  
  tbl_name: $[side == `buy; `bids; `asks];
  if[price in key book[tbl_name];
    lvl_size: book[tbl_name; price; `size] - actual_cancel;
    lvl_cnt: book[tbl_name; price; `orderCount] - ($[new_ord_size <= 0f; 1; 0]);
    if[(lvl_size <= 0f) or (lvl_cnt <= 0);
      book[tbl_name]: delete from book[tbl_name] where price = price
    ;
      book[tbl_name; price; `size]: lvl_size;
      book[tbl_name; price; `orderCount]: lvl_cnt
    ]
  ];
  book
 };

/ Execute incoming market trade against resting book
executeTrade: {[book; side; tradeSize]
  matched_qty: 0f;
  remaining: tradeSize;
  
  tbl_name: $[side == `buy; `asks; `bids];
  prices: $[side == `buy; asc key book`asks; desc key book`bids];
  
  i: 0;
  while[(i < count prices) and (remaining > 0f);
    p: prices[i];
    avail: book[tbl_name; p; `size];
    fill_qty: min (remaining; avail);
    matched_qty: matched_qty + fill_qty;
    remaining: remaining - fill_qty;
    
    if[fill_qty >= avail;
      book[tbl_name]: delete from book[tbl_name] where price = p
    ;
      book[tbl_name; p; `size]: avail - fill_qty
    ];
    i: i + 1
  ];
  `book`executed`unfilled!(book; matched_qty; remaining)
 };

/ Extract top K price levels for bids (sorted descending) and asks (sorted ascending)
getDepth: {[book; k]
  top_bids: k # `price xdesc select price, size, orderCount from book`bids;
  top_asks: k # `price xasc select price, size, orderCount from book`asks;
  `bids`asks!(top_bids; top_asks)
 };

/ Best Bid and Offer (BBO) and basic spread metrics
getBBO: {[book]
  bid_prices: desc key book`bids;
  ask_prices: asc key book`asks;
  if[(count bid_prices == 0) or (count ask_prices == 0);
    :`bestBid`bestAsk`midPrice`spread!(0n; 0n; 0n; 0n)
  ];
  bb: first bid_prices;
  ba: first ask_prices;
  `bestBid`bestAsk`midPrice`spread!(bb; ba; 0.5 * (bb + ba); ba - bb)
 };

/ Volume-Weighted Mid-Price (VWMP): P_vwmp = (bb * ask_vol + ba * bid_vol) / (bid_vol + ask_vol)
getVWMP: {[book]
  bbo: getBBO[book];
  if[null bbo`bestBid; :0n];
  bb: bbo`bestBid;
  ba: bbo`bestAsk;
  vb: book[`bids; bb; `size];
  va: book[`asks; ba; `size];
  tot: vb + va;
  $[tot > 0f; ((bb * va) + (ba * vb)) % tot; bbo`midPrice]
 };

/ Stoikov micro-price formulation: adjustments based on order book queue imbalance
/ Imbalance I = (v_b - v_a) / (v_b + v_a)
/ Micro-price = I * ba + (1 - I) * bb adjusted for tick discretization
stoikovMicroPrice: {[book]
  bbo: getBBO[book];
  if[null bbo`bestBid; :0n];
  bb: bbo`bestBid;
  ba: bbo`bestAsk;
  vb: book[`bids; bb; `size];
  va: book[`asks; ba; `size];
  tot: vb + va;
  if[tot == 0f; :bbo`midPrice];
  imb: (vb - va) % tot;
  bbo[`midPrice] + (0.5 * imb * bbo`spread)
 };

/ Volume-Synchronized Probability of Toxicity (VPIN): Easley, Lopez de Prado, O'Hara (2012)
/ trades: table with columns `price`size
/ bucket_vol: volume threshold per synchronized bucket V
/ n_buckets: rolling window length of buckets N
computeVPIN: {[trades; bucket_vol; n_buckets]
  prices: trades`price;
  sizes: trades`size;
  n: count trades;
  if[n < 2; :`vpin`buckets!(0f; ())];
  
  / Price changes and standard deviation for Bulk Volume Classification (BVC)
  dp: 0.0, 1 _ (prices - prev prices);
  sigma: .stats.std[dp];
  if[sigma == 0f; sigma: 1e-6];
  
  / Probability of buy trade via normal CDF of standardized price movement
  p_buy: .math.normalCdf each dp % sigma;
  
  / Accumulate into constant volume buckets
  buy_vols: ();
  sell_vols: ();
  
  curr_b: 0f;
  curr_s: 0f;
  rem_bucket: bucket_vol;
  
  i: 0;
  while[i < n;
    trade_sz: sizes[i];
    pb: p_buy[i];
    
    while[trade_sz > 0f;
      alloc: min (trade_sz; rem_bucket);
      curr_b: curr_b + (alloc * pb);
      curr_s: curr_s + (alloc * (1f - pb));
      rem_bucket: rem_bucket - alloc;
      trade_sz: trade_sz - alloc;
      
      if[rem_bucket <= 1e-6;
        buy_vols: buy_vols, enlist curr_b;
        sell_vols: sell_vols, enlist curr_s;
        curr_b: 0f;
        curr_s: 0f;
        rem_bucket: bucket_vol
      ]
    ];
    i: i + 1
  ];
  
  tot_buckets: count buy_vols;
  if[tot_buckets < n_buckets;
    abs_imbalance: abs buy_vols - sell_vols;
    vpin_val: (sum abs_imbalance) % (tot_buckets * bucket_vol);
    :`vpin`bucketCount!(vpin_val; tot_buckets)
  ];
  
  / Rolling VPIN over last n_buckets
  recent_b: (neg n_buckets) # buy_vols;
  recent_s: (neg n_buckets) # sell_vols;
  abs_imb: abs recent_b - recent_s;
  vpin: (sum abs_imb) % (1f * n_buckets * bucket_vol);
  `vpin`bucketCount!(vpin; tot_buckets)
 };

/ Multi-Level Order Flow Imbalance (OFI): Cont, Kukanov, Stoikov (2014)
/ Computes net order flow across top K price levels between snapshots t-1 and t
/ prev_book, curr_book: dictionaries with `bids and `asks tables
computeOFI: {[prev_book; curr_book; k_levels]
  p_prev: getDepth[prev_book; k_levels];
  p_curr: getDepth[curr_book; k_levels];
  
  ofi_levels: ();
  lvl: 0;
  while[lvl < k_levels;
    / Bid side innovation
    pb_curr: p_curr[`bids; lvl; `price];
    vb_curr: p_curr[`bids; lvl; `size];
    pb_prev: p_prev[`bids; lvl; `price];
    vb_prev: p_prev[`bids; lvl; `size];
    
    delta_bid: $[pb_curr > pb_prev; vb_curr;
                 pb_curr == pb_prev; vb_curr - vb_prev;
                 neg vb_prev];
    
    / Ask side innovation
    pa_curr: p_curr[`asks; lvl; `price];
    va_curr: p_curr[`asks; lvl; `size];
    pa_prev: p_prev[`asks; lvl; `price];
    va_prev: p_prev[`asks; lvl; `size];
    
    delta_ask: $[pa_curr < pa_prev; va_curr;
                 pa_curr == pa_prev; va_curr - va_prev;
                 neg va_prev];
    
    ofi_k: delta_bid - delta_ask;
    ofi_levels: ofi_levels, enlist ofi_k;
    lvl: lvl + 1
  ];
  `ofiTotal`ofiByLevel!(sum ofi_levels; ofi_levels)
 };

/ Roll (1984) effective spread estimator: S = 2 * sqrt(-Cov(Delta p_t, Delta p_{t-1}))
rollSpread: {[prices]
  n: count prices;
  if[n < 3; :0f];
  dp: 1 _ (prices - prev prices);
  dp_t: 1 _ dp;
  dp_tm1: -1 _ dp;
  cov: (sum (dp_t - .stats.mean[dp_t]) * (dp_tm1 - .stats.mean[dp_tm1])) % (count dp_t - 1);
  $[cov < 0f; 2f * sqrt[neg cov]; 0f]
 };

/ Glosten-Milgrom (1985) information asymmetry spread decomposition
/ S = alpha * S + c * S where alpha is adverse selection share and c is order processing share
glostenMilgromDecomp: {[prices; trades_direction; horizon]
  / trades_direction: +1 for buyer-initiated, -1 for seller-initiated
  n: count prices;
  if[n <= horizon; :`adverseSelectionShare`orderHandlingShare!(0.5; 0.5)];
  
  dp_ahead: (horizon _ prices) - (neg horizon) _ prices;
  q_t: (neg horizon) _ trades_direction;
  
  / Regress forward price change on trade sign: Delta p_{t+k} = theta * q_t + epsilon
  cov_pq: (sum dp_ahead * q_t) % (count dp_ahead);
  var_q: (sum q_t * q_t) % (count q_t);
  theta: cov_pq % var_q;
  
  half_spread: 0.5 * .stats.mean[abs 1 _ (prices - prev prices)];
  alpha: $[half_spread > 0f; min (1.0; max (0.0; theta % half_spread)); 0.5];
  `adverseSelectionShare`orderHandlingShare!(alpha; 1f - alpha)
 };

\d .
