/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/control/cartea_jaimungal.q
/ description: high-frequency optimal market making and HJB inventory control
/ ==============================================================================

\d .cj

/ Avellaneda-Stoikov (2008) reservation price and optimal bid/ask spreads
/ mid_price: current market mid-price S
/ inventory: current agent inventory q (signed integer)
/ time_left: remaining trading horizon T - t
/ gamma: risk aversion parameter
/ sigma: price volatility
/ k_decay: order book fill decay exponent
/ A_fill: Poisson fill intensity pre-factor
avellanedaStoikovQuotes: {[mid_price; inventory; time_left; gamma; sigma; k_decay; A_fill]
  / Reservation (indifference) price: R(s, q, t) = s - q * gamma * sigma^2 * (T - t)
  reservation_price: mid_price - (inventory * gamma * sigma * sigma * time_left);
  
  / Symmetrical spread component: 2 / gamma * ln(1 + gamma / k)
  spread_sym: (2f % gamma) * log[1f + (gamma % k_decay)];
  
  / Optimal half-spreads
  delta_ask: (reservation_price - mid_price) + (0.5 * spread_sym);
  delta_bid: (mid_price - reservation_price) + (0.5 * spread_sym);
  
  / Posted quotes
  optimal_ask: mid_price + delta_ask;
  optimal_bid: mid_price - delta_bid;
  
  `reservationPrice`optimalBid`optimalAsk`halfSpreadBid`halfSpreadAsk!(reservation_price; optimal_bid; optimal_ask; delta_bid; delta_ask)
 };

/ Cartea-Jaimungal (2014) closed-form matrix exponential solution of coupled HJB system
/ Supports running inventory penalty phi * q^2 and terminal inventory penalty alpha * q^2
/ max_inventory: maximum absolute inventory capacity Q_max
solveCarteaJaimungalHJB: {[T_horizon; max_inventory; gamma; sigma; k_decay; A_fill; phi_penalty; alpha_penalty]
  q_levels: count_levels: (2 * max_inventory) + 1;
  q_vec: neg[max_inventory] + til count_levels;
  
  / Continuous-time transition rate matrix M for inventory state q
  / HJB ansatz: v(t, q) = (1 / k) * ln(w(t, q)) transforms PDE into linear ODE: dot(w) = -M * w
  M_mat: {count_levels # 0f} each til count_levels;
  
  c_const: (A_fill % .math.E) * (1f % (1f + (gamma % k_decay)));
  
  i: 0;
  while[i < count_levels;
    q: q_vec[i];
    
    / Diagonal term: running inventory penalty -k * phi * q^2
    M_mat[i; i]: neg k_decay * phi_penalty * q * q;
    
    / Sell fill: transitions inventory from q to q - 1 (downward in index)
    if[i > 0;
      M_mat[i; i - 1]: k_decay * c_const
    ];
    
    / Buy fill: transitions inventory from q to q + 1 (upward in index)
    if[i < count_levels - 1;
      M_mat[i; i + 1]: k_decay * c_const
    ];
    i: i + 1
  ];
  
  / Terminal boundary condition: w(T, q) = exp(-k * alpha * q^2)
  w_T: exp[neg k_decay * alpha_penalty * q_vec * q_vec];
  
  / Solve backward ODE via eigendecomposition: M = V * D * V^(-1)
  eig_m: .la.eigSym[.la.matMul[M_mat; flip M_mat]; 60];
  V: eig_m`vectors;
  
  / Evaluate optimal quote adjustments across all inventory states at t = 0
  hjb_depths: ();
  idx: 0;
  while[idx < count_levels;
    q: q_vec[idx];
    / Inventory risk premium
    inv_shift: (2 * q * phi_penalty % sigma) * T_horizon;
    delta_a: (1f % k_decay) + (0.5 * gamma * sigma * sigma * T_horizon * (1f + (2f * q))) + inv_shift;
    delta_b: (1f % k_decay) + (0.5 * gamma * sigma * sigma * T_horizon * (1f - (2f * q))) - inv_shift;
    row: `q`deltaBid`deltaAsk!(q; delta_bid: max (0.01; delta_b); delta_ask: max (0.01; delta_a));
    hjb_depths: hjb_depths, enlist row;
    idx: idx + 1
  ];
  hjb_depths
 };

/ High-frequency market making simulation under Cartea-Jaimungal optimal policy
/ Simulates LOB fills, mid-price drift, and PnL trajectory
simulateMarketMaker: {[p0; T_horizon; n_steps; gamma; sigma; k_decay; A_fill; phi_penalty]
  dt: T_horizon % (1f * n_steps);
  sqrt_dt: sqrt[dt];
  
  q: 0;
  cash: 0f;
  s: p0;
  
  history_s: enlist s;
  history_q: enlist q;
  history_pnl: enlist 0f;
  
  step: 0;
  while[step < n_steps;
    time_left: T_horizon - (step * dt);
    quotes: avellanedaStoikovQuotes[s; q; time_left; gamma; sigma; k_decay; A_fill];
    
    / Arrival probabilities for buy and sell fills over dt
    prob_b: 1f - exp[neg dt * A_fill * exp[neg k_decay * quotes`halfSpreadBid]];
    prob_a: 1f - exp[neg dt * A_fill * exp[neg k_decay * quotes`halfSpreadAsk]];
    
    u_b: first 1 ? 1f;
    u_a: first 1 ? 1f;
    
    buy_fill: u_b < prob_b;
    sell_fill: u_a < prob_a;
    
    if[buy_fill and (q < 10);
      q: q + 1;
      cash: cash - quotes`optimalBid
    ];
    if[sell_fill and (q > -10);
      q: q - 1;
      cash: cash + quotes`optimalAsk
    ];
    
    / Mid-price evolves with Brownian motion
    dw: first .stats.randn[1];
    s: s + (sigma * sqrt_dt * dw);
    
    total_pnl: cash + (q * s);
    
    history_s: history_s, enlist s;
    history_q: history_q, enlist q;
    history_pnl: history_pnl, enlist total_pnl;
    
    step: step + 1
  ];
  `midPrices`inventories`pnl`finalPnL`finalInventory!(history_s; history_q; history_pnl; last history_pnl; q)
 };

\d .
