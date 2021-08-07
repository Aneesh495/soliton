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

\d .
