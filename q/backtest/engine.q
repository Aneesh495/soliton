/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/backtest/engine.q
/ description: discrete-event limit order book backtester and TCA attribution
/ ==============================================================================

\d .backtest

/ Initialize backtest simulation environment
/ maker_fee: negative indicates rebate (e.g. -0.0002 for 2 bps rebate)
/ taker_fee: fee on aggressive crossing (e.g. 0.0005 for 5 bps fee)
/ latency_nanos: one-way network transit latency in nanoseconds
initEngine: {[maker_fee; taker_fee; latency_nanos]
  orders: ([orderId: `long$()] sym: `symbol$(); side: `symbol$(); price: `float$(); origSize: `float$(); remSize: `float$();
           status: `symbol$(); submitTime: `timestamp$(); ackTime: `timestamp$(); queueAhead: `float$());
  fills: ([] time: `timestamp$(); orderId: `long$(); price: `float$(); size: `float$(); fee: `float$(); isMaker: `boolean$());
  
  `makerFee`takerFee`latency`orders`fills`cash`pnl!(maker_fee; taker_fee; latency_nanos; orders; fills; 0f; 0f)
 };

/ Submit limit order: generates simulated ack with latency and tracks queue volume ahead
submitLimitOrder: {[engine; orderId; sym; side; price; size; submit_time; current_depth]
  ack_time: submit_time + engine`latency;
  
  / Estimate volume ahead at price level from current market depth
  vol_ahead: 0f;
  if[side == `buy;
    bids: current_depth`bids;
    if[price in bids`price; vol_ahead: bids[price; `size]]
  ;
    asks: current_depth`asks;
    if[price in asks`price; vol_ahead: asks[price; `size]]
  ];
  
  new_order: (enlist orderId)!(enlist sym; enlist side; enlist price; enlist size; enlist size;
                              enlist `pending; enlist submit_time; enlist ack_time; enlist vol_ahead);
  engine[`orders]: engine[`orders], new_order;
  engine
 };

/ Process market trade feed tick to match resting child orders with queue priority
processMarketTrade: {[engine; trade_time; sym; trade_price; trade_size; trade_side]
  active_orders: select from engine`orders where sym = sym, status in `pending`active;
  if[count active_orders == 0; :engine];
  
  order_ids: key active_orders;
  i: 0;
  while[i < count order_ids;
    oid: order_ids[i];
    ord: engine[`orders; oid];
    
    / Check if order has reached market after latency delay
    if[trade_time >= ord`ackTime;
      if[ord`status == `pending; engine[`orders; oid; `status]: `active];
      
      / Check matching conditions
      side: ord`side;
      price: ord`price;
      rem_sz: ord`remSize;
      q_ahead: ord`queueAhead;
      
      can_fill: false;
      fill_qty: 0f;
      
      if[(side == `buy) and (trade_price <= price);
        / Seller traded: if trade price < my price, full priority; if trade price == my price, deduct from queueAhead
        if[trade_price < price;
          fill_qty: min (rem_sz; trade_size);
          can_fill: true
        ;
          if[q_ahead <= trade_size;
            fill_qty: min (rem_sz; trade_size - q_ahead);
            engine[`orders; oid; `queueAhead]: 0f;
            can_fill: fill_qty > 0f
          ;
            engine[`orders; oid; `queueAhead]: q_ahead - trade_size
          ]
        ]
      ];
      
      if[(side == `sell) and (trade_price >= price);
        if[trade_price > price;
          fill_qty: min (rem_sz; trade_size);
          can_fill: true
        ;
          if[q_ahead <= trade_size;
            fill_qty: min (rem_sz; trade_size - q_ahead);
            engine[`orders; oid; `queueAhead]: 0f;
            can_fill: fill_qty > 0f
          ;
            engine[`orders; oid; `queueAhead]: q_ahead - trade_size
          ]
        ]
      ];
      
      if[can_fill;
        / Record fill
        fee_rate: engine`makerFee;
        fee_paid: fill_qty * price * fee_rate;
        new_fill: enlist `time`orderId`price`size`fee`isMaker!(trade_time; oid; price; fill_qty; fee_paid; true);
        engine[`fills]: engine[`fills], new_fill;
        
        new_rem: rem_sz - fill_qty;
        engine[`orders; oid; `remSize]: new_rem;
        if[new_rem <= 0f; engine[`orders; oid; `status]: `filled];
        
        / Update cash balance
        cash_delta: $[side == `buy; neg (fill_qty * price) + fee_paid; (fill_qty * price) - fee_paid];
        engine[`cash]: engine[`cash] + cash_delta
      ]
    ];
    i: i + 1
  ];
  engine
 };

/ Cancel existing resting order
cancelOrder: {[engine; orderId; cancel_time]
  if[orderId in key engine`orders;
    engine[`orders; orderId; `status]: `canceled
  ];
  engine
 };

/ Transaction Cost Analysis (TCA) and Implementation Shortfall attribution
/ benchmark_price: decision price P_0 when order was routed
/ market_trades: trade history over the execution window
computeTCA: {[engine; benchmark_price; market_trades]
  fills: engine`fills;
  if[count fills == 0;
    :`decisionPrice`executedShares`vwapExec`implementationShortfall`slippageBps!(benchmark_price; 0f; 0f; 0f; 0f)
  ];
  
  tot_shares: sum fills`size;
  tot_dollars: sum fills[`price] * fills`size;
  vwap_exec: tot_dollars % tot_shares;
  
  / Market VWAP over same period
  mkt_shares: sum market_trades`size;
  mkt_dollars: sum market_trades[`price] * market_trades`size;
  vwap_mkt: $[mkt_shares > 0f; mkt_dollars % mkt_shares; benchmark_price];
  
  / Implementation Shortfall: (VWAP_exec - P_0) * Total_Shares + Fees
  gross_is: (vwap_exec - benchmark_price) * tot_shares;
  tot_fees: sum fills`fee;
  net_is: gross_is + tot_fees;
  
  slippage_bps: (10000f * (vwap_exec - benchmark_price)) % benchmark_price;
  vwap_slippage_bps: (10000f * (vwap_exec - vwap_mkt)) % vwap_mkt;
  
  fill_rate: tot_shares % sum engine[`orders]`origSize;
  
  `decisionPrice`executedShares`vwapExec`vwapMarket`grossIS`netIS`totalFees`slippageBps`vwapSlippageBps`fillRate!(
    benchmark_price; tot_shares; vwap_exec; vwap_mkt; gross_is; net_is; tot_fees; slippage_bps; vwap_slippage_bps; fill_rate
  )
 };

\d .
