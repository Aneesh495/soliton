/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/microstructure/impact.q
/ description: price impact models, Hasbrouck SVAR, and transient decay kernels
/ ==============================================================================

\d .impact

/ Kyle (1985) linear price impact lambda: Delta p_t = lambda * Q_t + epsilon_t
/ prices: trade price series [N]
/ volumes: signed trade volume series [N] (+ for buy, - for sell)
fitKyleLambda: {[prices; volumes]
  n: count prices;
  if[n < 3; :`lambda`rSquared!(0f; 0f)];
  dp: 1 _ (prices - prev prices);
  q: -1 _ volumes;
  
  / OLS regression of dp on signed volume q
  mean_q: .stats.mean[q];
  mean_dp: .stats.mean[dp];
  cov_qp: sum (q - mean_q) * (dp - mean_dp);
  var_q: sum (q - mean_q) * (q - mean_q);
  var_dp: sum (dp - mean_dp) * (dp - mean_dp);
  
  lambda_val: $[var_q > 0f; cov_qp % var_q; 0f];
  r2: $[(var_q > 0f) and (var_dp > 0f); (cov_qp * cov_qp) % (var_q * var_dp); 0f];
  `lambda`rSquared`residuals!(lambda_val; r2; dp - (lambda_val * q))
 };

/ Hasbrouck (1991) Structural Vector Autoregression (SVAR) of returns r_t and trade signs x_t
/ Estimates bivariate VAR(p) to separate permanent price impact from temporary microstructure noise
fitHasbrouckVAR: {[returns; trade_signs; lags]
  n: count returns;
  p: lags;
  if[n <= (p + 2); :`coeffR`coeffX`permanentImpact!((); (); 0f)];
  
  / Build design matrix X with lagged returns and trade signs
  k: n - p;
  Y_r: (neg k) # returns;
  Y_x: (neg k) # trade_signs;
  
  design_cols: ();
  / Contemporaneous trade sign
  design_cols: design_cols, enlist Y_x;
  
  l: 1;
  while[l <= p;
    lag_r: (neg k) # (neg l) _ (l # 0f, returns);
    lag_x: (neg k) # (neg l) _ (l # 0f, trade_signs);
    design_cols: design_cols, enlist lag_r;
    design_cols: design_cols, enlist lag_x;
    l: l + 1
  ];
  
  X_mat: flip design_cols;
  Xt: design_cols;
  XtX: .la.matMul[Xt; X_mat];
  
  / Regularized normal equations: (X^T * X + lambda * I) * b = X^T * Y
  I: .la.eye[count XtX];
  i: 0;
  while[i < count XtX;
    XtX[i; i]: XtX[i; i] + 1e-6;
    i: i + 1
  ];
  
  Xty_r: .la.matVec[Xt; Y_r];
  beta_r: .la.solve[XtX; Xty_r];
  
  / Permanent impact is cumulative impulse response to unexpected trade shock
  perm_impact: sum beta_r;
  `contemporaneousCoeff`laggedCoeffs`permanentImpact!(beta_r[0]; 1 _ beta_r; perm_impact)
 };

/ Bouchaud power-law propagator kernel: G(t) = Gamma_0 / (1 + t / tau)^gamma
/ gamma typically around 0.5 (square-root decay of temporary impact)
propagatorPowerLaw: {[t; gamma0; tau; gamma]
  gamma0 % ((1f + (t % tau)) xexp gamma)
 };

/ Gatheral (2010) exponential transient impact decay kernel: G(t) = G_0 * exp(-rho * t)
propagatorExp: {[t; g0; rho]
  g0 * exp[neg rho * t]
 };

/ Simulate transient price impact trajectory under discrete trade execution schedule
/ v_schedule: executed shares at each step [M]
/ dt_step: time step between executions
/ kernel_type: `exp or `power
simulatePriceTrajectory: {[p0; v_schedule; dt_step; g0; decay_param; kernel_type]
  m: count v_schedule;
  prices: m # 0f;
  curr_p: p0;
  
  t: 0;
  while[t < m;
    impact_t: 0f;
    k: 0;
    while[k <= t;
      elapsed: (t - k) * dt_step;
      g_val: $[kernel_type == `exp; propagatorExp[elapsed; g0; decay_param];
               propagatorPowerLaw[elapsed; g0; 1.0; decay_param]];
      impact_t: impact_t + (g_val * v_schedule[k]);
      k: k + 1
    ];
    prices[t]: p0 + impact_t;
    t: t + 1
  ];
  prices
 };

/ Cross-Asset Price Impact Matrix Estimation
/ Relates price returns of asset i to order flow across all N correlated assets
/ Returns matrix Lambda [N x N] where Lambda_ij is cross-impact of asset j on asset i
fitCrossAssetImpact: {[returns_matrix; flows_matrix]
  / returns_matrix: [T x N], flows_matrix: [T x N]
  T: count returns_matrix;
  N: count first returns_matrix;
  
  Ft: flip flows_matrix;
  FtF: .la.matMul[Ft; flows_matrix];
  
  / Regularize FtF
  i: 0;
  while[i < N;
    FtF[i; i]: FtF[i; i] + 1e-5;
    i: i + 1
  ];
  
  / Lambda^T = (F^T * F)^(-1) * F^T * R
  FtR: .la.matMul[Ft; returns_matrix];
  Lambda_t: .la.solve[FtF; FtR];
  flip Lambda_t
 };

\d .
