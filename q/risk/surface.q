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

/ Gatheral's Stochastic Volatility Inspired (SVI) raw total variance:
/ w(k) = a + b * (rho * (k - m) + sqrt((k - m)^2 + sigma^2))
/ k = ln(K / F) is log-moneyness
sviTotalVariance: {[k; a; b; rho; m; sigma]
  dk: k - m;
  a + (b * ((rho * dk) + sqrt[(dk * dk) + (sigma * sigma)]))
 };

/ Convert SVI total variance w(k) into annualized implied Black-Scholes volatility: sigma_BS = sqrt(w(k) / T)
sviImpliedVol: {[k; T; a; b; rho; m; sigma]
  w: sviTotalVariance[k; a; b; rho; m; sigma];
  sqrt[max (1e-6; w % T)]
 };

/ Verify Durrleman no-arbitrage bounds for SVI parameters
isSVIArbitrageFree: {[a; b; rho; m; sigma]
  c1: b >= 0f;
  c2: abs[rho] < 1f;
  c3: sigma > 0f;
  c4: (a + (b * sigma * sqrt[1f - (rho * rho)])) >= 0f;
  c5: (b * (1f + abs[rho])) < 4f; / Lee's moment formula bound for total variance slope
  c1 and c2 and c3 and c4 and c5
 };

/ Parametric Gaussian Value-at-Risk (VaR)
varParametric: {[returns; alpha_conf]
  mu: .stats.mean[returns];
  s: .stats.std[returns];
  z: .math.normalInv[1f - alpha_conf];
  neg (mu + (s * z))
 };

/ Cornish-Fisher expansion Value-at-Risk adjusting for skewness and excess kurtosis
varCornishFisher: {[returns; alpha_conf]
  mu: .stats.mean[returns];
  s: .stats.std[returns];
  sk: .stats.skew[returns];
  kt: .stats.kurt[returns];
  
  z: .math.normalInv[1f - alpha_conf];
  z2: z * z;
  z3: z2 * z;
  
  / Cornish-Fisher expansion quantile
  z_cf: z + ((z2 - 1f) * sk % 6f) + ((z3 - (3f * z)) * kt % 24f) - (((2f * z3) - (5f * z)) * sk * sk % 36f);
  neg (mu + (s * z_cf))
 };

/ Historical simulation Value-at-Risk
varHistorical: {[returns; alpha_conf]
  p_tail: 1f - alpha_conf;
  neg .stats.quantile[returns; p_tail]
 };

/ Coherent Expected Shortfall (CVaR): average loss beyond VaR threshold
expectedShortfall: {[returns; alpha_conf]
  var_thresh: varHistorical[returns; alpha_conf];
  losses: neg returns;
  tail_losses: losses where losses >= var_thresh;
  $[count tail_losses > 0; .stats.mean[tail_losses]; var_thresh]
 };

/ Euler Marginal Risk Contributions for a multi-asset portfolio with weights w and covariance Sigma
/ MRC_i = (Sigma * w)_i / sqrt(w^T * Sigma * w)
marginalRiskContribution: {[weights; Sigma]
  port_var: sum weights * .la.matVec[Sigma; weights];
  port_std: sqrt[max (1e-12; port_var)];
  cov_w: .la.matVec[Sigma; weights];
  cov_w % port_std
 };

\d .
