/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/risk/surface.q
/ description: streaming Greeks surfaces, SVI parameterization, VaR and ES
/ ==============================================================================

\d .risk

/ Analytical Black-Scholes Greeks: Delta, Gamma, Vega, Theta, Rho, Vanna, Volga
analyticalGreeks: {[S0; K_strike; T; r; q; sigma; isCall]
  if[(T <= 0f) or (sigma <= 0f);
    intrinsic: $[isCall; max (0f; S0 - K_strike); max (0f; K_strike - S0)];
    delta_val: $[isCall; $[S0 > K_strike; 1f; 0f]; $[S0 < K_strike; -1f; 0f]];
    :`price`delta`gamma`vega`theta`rho`vanna`volga!(intrinsic; delta_val; 0f; 0f; 0f; 0f; 0f; 0f)
  ];
  
  sqrt_T: sqrt[T];
  d1: (log[S0 % K_strike] + ((r - q + (0.5 * sigma * sigma)) * T)) % (sigma * sqrt_T);
  d2: d1 - (sigma * sqrt_T);
  
  pdf_d1: .math.normalPdf[d1];
  cdf_d1: .math.normalCdf[d1];
  cdf_d2: .math.normalCdf[d2];
  
  exp_qT: exp[neg q * T];
  exp_rT: exp[neg r * T];
  
  price: $[isCall; (S0 * exp_qT * cdf_d1) - (K_strike * exp_rT * cdf_d2);
                   (K_strike * exp_rT * (1f - cdf_d2)) - (S0 * exp_qT * (1f - cdf_d1))];
  
  delta: $[isCall; exp_qT * cdf_d1; exp_qT * (cdf_d1 - 1f)];
  gamma: (exp_qT * pdf_d1) % (S0 * sigma * sqrt_T);
  vega: S0 * exp_qT * sqrt_T * pdf_d1;
  
  theta_term1: neg (S0 * exp_qT * sigma * pdf_d1) % (2f * sqrt_T);
  theta: $[isCall; theta_term1 - (r * K_strike * exp_rT * cdf_d2) + (q * S0 * exp_qT * cdf_d1);
                   theta_term1 + (r * K_strike * exp_rT * (1f - cdf_d2)) - (q * S0 * exp_qT * (1f - cdf_d1))];
  
  rho_val: $[isCall; K_strike * T * exp_rT * cdf_d2; neg K_strike * T * exp_rT * (1f - cdf_d2)];
  
  / Second-order cross sensitivities: Vanna = d(Vega)/dS = d(Delta)/d(sigma)
  vanna: neg exp_qT * pdf_d1 * (d2 % sigma);
  
  / Volga (Vomma) = d(Vega)/d(sigma) = Vega * d1 * d2 / sigma
  volga: vega * d1 * d2 % sigma;
  
  `price`delta`gamma`vega`theta`rho`vanna`volga!(price; delta; gamma; vega; theta; rho_val; vanna; volga)
 };

/ Construct 2D streaming Greeks surface over matrix grid of strikes and maturities
constructGreeksSurface: {[S0; strikes; maturities; r; q; vol_matrix; isCall]
  n_k: count strikes;
  n_t: count maturities;
  
  surface_rows: ();
  i: 0;
  while[i < n_t;
    T: maturities[i];
    j: 0;
    while[j < n_k;
      K_val: strikes[j];
      vol: vol_matrix[i; j];
      greeks: analyticalGreeks[S0; K_val; T; r; q; vol; isCall];
      row: `maturity`strike`impliedVol`price`delta`gamma`vega`vanna`volga!(T; K_val; vol; greeks`price; greeks`delta; greeks`gamma; greeks`vega; greeks`vanna; greeks`volga);
      surface_rows: surface_rows, enlist row;
      j: j + 1
    ];
    i: i + 1
  ];
  surface_rows
 };

\d .
