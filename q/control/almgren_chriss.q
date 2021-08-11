/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/control/almgren_chriss.q
/ description: Almgren-Chriss risk-averse optimal execution and efficient frontiers
/ ==============================================================================

\d .ac

/ Compute Almgren-Chriss analytical optimal holding trajectory
/ total_shares: initial target position X_0
/ horizon: total execution window T
/ n_steps: number of discrete trading intervals N
/ sigma: asset price volatility
/ gamma: permanent market impact coefficient
/ eta: temporary market impact coefficient
/ lambda_risk: risk aversion parameter (lambda = 0 gives TWAP linear execution)
optimalTrajectory: {[total_shares; horizon; n_steps; sigma; gamma; eta; lambda_risk]
  tau: horizon % (1f * n_steps);
  
  if[lambda_risk <= 1e-12;
    / Risk-neutral limit: TWAP linear liquidation trajectory
    step_vol: total_shares % (1f * n_steps);
    holdings: total_shares - (step_vol * (til n_steps + 1));
    trades: n_steps # step_vol;
    :`holdings`trades`kappa`tau!(holdings; trades; 0f; tau)
  ];
  
  / Intrinsic urgency parameter kappa
  kappa_tilde_sq: (lambda_risk * sigma * sigma) % eta;
  / Exact discrete-time urgency parameter: 2 * (cosh(kappa * tau) - 1) = kappa_tilde_sq * tau^2
  cosh_val: 1f + (0.5 * kappa_tilde_sq * tau * tau);
  / Arccosh(y) = ln(y + sqrt(y^2 - 1))
  kappa_tau: log[cosh_val + sqrt[(cosh_val * cosh_val) - 1f]];
  kappa: kappa_tau % tau;
  
  / Holdings at discrete time grid points t_j = j * tau
  times: tau * (til n_steps + 1);
  sinh_kappa_T: 0.5 * (exp[kappa * horizon] - exp[neg kappa * horizon]);
  
  holdings: {
    rem_t: horizon - x;
    sinh_rem: 0.5 * (exp[kappa * rem_t] - exp[neg kappa * rem_t]);
    total_shares * (sinh_rem % sinh_kappa_T)
  } each times;
  / Force final holding to zero
  holdings[n_steps]: 0f;
  
  / Trades executed in each interval: n_j = x_{j-1} - x_j
  trades: (neg 1 _ holdings) - 1 _ holdings;
  `holdings`trades`kappa`tau!(holdings; trades; kappa; tau)
 };

/ Expected cost and variance of capture for Almgren-Chriss schedule
evaluateCostAndRisk: {[total_shares; trajectory; sigma; gamma; eta]
  trades: trajectory`trades;
  holdings: trajectory`holdings;
  tau: trajectory`tau;
  
  / Permanent impact cost: 0.5 * gamma * X^2
  perm_cost: 0.5 * gamma * total_shares * total_shares;
  
  / Temporary impact cost: eta * sum (n_j^2 / tau)
  temp_cost: (eta % tau) * sum trades * trades;
  expected_cost: perm_cost + temp_cost;
  
  / Variance of capture: sigma^2 * tau * sum(x_j^2)
  holding_nodes: -1 _ holdings;
  variance_cost: sigma * sigma * tau * sum holding_nodes * holding_nodes;
  std_cost: sqrt[variance_cost];
  
  `expectedCost`varianceCost`stdCost`permanentCost`temporaryCost!(expected_cost; variance_cost; std_cost; perm_cost; temp_cost)
 };

/ Compute efficient frontier of expected cost versus standard deviation of cost across lambdas
efficientFrontier: {[total_shares; horizon; n_steps; sigma; gamma; eta; lambda_grid]
  frontier: ();
  i: 0;
  while[i < count lambda_grid;
    lam: lambda_grid[i];
    traj: optimalTrajectory[total_shares; horizon; n_steps; sigma; gamma; eta; lam];
    eval_res: evaluateCostAndRisk[total_shares; traj; sigma; gamma; eta];
    row: `lambda`expectedCost`stdCost!(lam; eval_res`expectedCost; eval_res`stdCost);
    frontier: frontier, enlist row;
    i: i + 1
  ];
  frontier
 };

/ Nonlinear power-law temporary market impact optimization: h(v) = eta * |v|^alpha
/ Solves nonlinear calculus of variations trajectory via Nelder-Mead
nonlinearOptimalTrajectory: {[total_shares; horizon; n_steps; sigma; gamma; eta; power_alpha; lambda_risk]
  tau: horizon % (1f * n_steps);
  
  / Objective function to minimize: Total expected impact cost + lambda * Variance
  obj: {[mid_nodes; X0; tau; sigma; gamma; eta; alpha_p; lam]
    / Reconstruct full holding trajectory [X0, x_1, x_2, ..., x_{N-1}, 0]
    holdings: (enlist X0), mid_nodes, 0f;
    trades: (neg 1 _ holdings) - 1 _ holdings;
    
    / Penalize negative trades (strictly monotonic liquidation)
    penalty: 1e6 * sum max'[neg trades; 0f];
    
    / Nonlinear temporary impact cost: sum eta * |n_k / tau|^alpha * n_k
    rates: abs[trades] % tau;
    temp_cost: sum eta * (rates xexp alpha_p) * trades;
    
    / Variance cost: lambda * sigma^2 * tau * sum(x_k^2)
    var_cost: lam * sigma * sigma * tau * sum holdings * holdings;
    
    temp_cost + var_cost + penalty
  };
  
  / Initial guess: linear TWAP holdings
  init_mid: {total_shares * (1f - (x % (1f * n_steps)))} each 1 + til n_steps - 1;
  
  opt_res: .optim.nelderMead[obj[; total_shares; tau; sigma; gamma; eta; power_alpha; lambda_risk]; init_mid; 200; 1e-4];
  best_mid: opt_res`x;
  
  best_holdings: (enlist total_shares), best_mid, 0f;
  best_trades: (neg 1 _ best_holdings) - 1 _ best_holdings;
  `holdings`trades`tau!(best_holdings; best_trades; tau)
 };

/ Multi-Asset optimal liquidation with cross-asset asset covariance Sigma
/ and cross-asset temporary market impact matrix H_mat
/ Total shares vector X0 [D], Horizon T, N steps
multiAssetTrajectory: {[X0; horizon; n_steps; Sigma; H_mat; lambda_risk]
  d: count X0;
  tau: horizon % (1f * n_steps);
  
  / Urgency matrix M = lambda * H^(-1) * Sigma
  H_inv: .la.inv[H_mat];
  M_mat: lambda_risk * .la.matMul[H_inv; Sigma];
  
  / Eigendecomposition of symmetric positive-definite urgency matrix M
  eig_m: .la.eigSym[M_mat; 60];
  V: eig_m`vectors;
  lambdas_m: eig_m`values;
  
  / Eigen-coordinate initial positions: Y0 = V^T * X0
  Y0: .la.matVec[.la.transpose[V]; X0];
  
  / Solve 1D problem independently for each eigen-mode
  eigen_trajectories: ();
  i: 0;
  while[i < d;
    omega_i: sqrt[max (0f; lambdas_m[i])];
    
    if[omega_i <= 1e-8;
      / Linear liquidation
      traj_i: {Y0[i] * (1f - (x % horizon))} each tau * til n_steps + 1
    ;
      sinh_om_T: 0.5 * (exp[omega_i * horizon] - exp[neg omega_i * horizon]);
      traj_i: {
        rem_t: horizon - x;
        sinh_rem: 0.5 * (exp[omega_i * rem_t] - exp[neg omega_i * rem_t]);
        Y0[i] * (sinh_rem % sinh_om_T)
      } each tau * til n_steps + 1
    ];
    eigen_trajectories: eigen_trajectories, enlist traj_i;
    i: i + 1
  ];
  
  / Transform back to original asset coordinates: X(t) = V * Y(t)
  asset_holdings: ();
  step: 0;
  while[step <= n_steps;
    Y_step: {eigen_trajectories[x; step]} each til d;
    X_step: .la.matVec[V; Y_step];
    asset_holdings: asset_holdings, enlist X_step;
    step: step + 1
  ];
  
  `assetHoldings`tau!(asset_holdings; tau)
 };

/ Kyle-Obizhaeva (2016) execution under transient impact with exponential replenishment rate rho
/ Closed-form optimal block trades at t = 0 and t = T, with constant continuous rate between
kyleObizhaevaTrajectory: {[total_shares; horizon; rho_replenish; lambda_impact]
  / Block trade at t = 0: x0_block = X / (2 + rho * T)
  denom: 2f + (rho_replenish * horizon);
  block_0: total_shares % denom;
  block_T: block_0;
  continuous_rate: (rho_replenish * total_shares) % denom;
  
  expected_slippage: lambda_impact * block_0 * (1f + (0.5 * rho_replenish * horizon));
  `initialBlock`terminalBlock`continuousRate`expectedSlippage!(block_0; block_T; continuous_rate; expected_slippage)
 };

\d .
