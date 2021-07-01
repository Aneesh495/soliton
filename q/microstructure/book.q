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

\d .
