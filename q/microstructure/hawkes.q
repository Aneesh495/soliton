/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/microstructure/hawkes.q
/ description: multidimensional mutually-exciting point processes and Ogata thinning
/ ==============================================================================

\d .hawkes

/ Initialize M-dimensional Hawkes process model parameters
/ mu: baseline intensity vector [M]
/ alpha: excitation magnitude matrix [M x M] (alpha[i; j] is excitation of i by j)
/ beta: exponential decay rate matrix [M x M]
initModel: {[mu; alpha; beta]
  m: count mu;
  / Compute branching ratio matrix: Gamma_ij = alpha_ij / beta_ij
  gamma: {alpha[x] % beta[x]} each til m;
  
  / Verify stationarity: spectral radius of gamma matrix must be strictly < 1
  eig_res: .la.eigSym[.la.matMul[gamma; flip gamma]; 50];
  spec_radius: sqrt max eig_res`values;
  is_stationary: spec_radius < 1.0;
  
  `dim`mu`alpha`beta`gamma`specRadius`isStationary!(m; mu; alpha; beta; gamma; spec_radius; is_stationary)
 };

/ Compute instantaneous conditional intensity vector lambda(t) given current history state
/ r_state: recursive state matrix R[i; j] representing sum_{t_k < t} exp(-beta_ij * (t - t_k))
evalIntensity: {[model; r_state]
  m: model`dim;
  mu: model`mu;
  alpha: model`alpha;
  
  / lambda_i = mu_i + sum_j (alpha_ij * R_ij)
  lambdas: m # 0f;
  i: 0;
  while[i < m;
    lambdas[i]: mu[i] + sum alpha[i] * r_state[i];
    i: i + 1
  ];
  lambdas
 };

/ Ogata's modified thinning algorithm for exact continuous-time Hawkes simulation
/ Generates exact event arrival times and event marks over interval [0, T_max]
simulateOgata: {[model; t_max; max_events]
  m: model`dim;
  mu: model`mu;
  alpha: model`alpha;
  beta: model`beta;
  
  / State matrices: R[i; j]
  r_state: {m # 0f} each til m;
  
  event_times: `float$();
  event_types: `long$();
  
  t: 0f;
  n_ev: 0;
  
  while[(t < t_max) and (n_ev < max_events);
    / Compute current intensity upper bound lambda_bar
    lambdas: evalIntensity[model; r_state];
    lambda_bar: sum lambdas;
    if[lambda_bar <= 0f; lambda_bar: sum mu];
    
    / Generate candidate arrival time from homogeneous Poisson process
    u1: max (1e-12; first 1 ? 1f);
    dt: neg log[u1] % lambda_bar;
    t_candidate: t + dt;
    if[t_candidate >= t_max; t: t_max];
    
    if[t_candidate < t_max;
      / Decay recursive state R to candidate time
      r_candidate: {r_state[x] * exp[neg beta[x] * dt]} each til m;
      lambdas_candidate: evalIntensity[model; r_candidate];
      lambda_tot_candidate: sum lambdas_candidate;
      
      / Rejection sampling test
      u2: first 1 ? 1f;
      if[u2 <= (lambda_tot_candidate % lambda_bar);
        / Accepted: select which dimension m fired
        cum_probs: (sums lambdas_candidate) % lambda_tot_candidate;
        u3: first 1 ? 1f;
        fired_type: first where cum_probs >= u3;
        if[null fired_type; fired_type: m - 1];
        
        / Record event
        event_times: event_times, enlist t_candidate;
        event_types: event_types, enlist fired_type;
        n_ev: n_ev + 1;
        
        / Update state with new jump in column fired_type
        j: fired_type;
        i: 0;
        while[i < m;
          r_candidate[i; j]: r_candidate[i; j] + 1f;
          i: i + 1
        ];
        r_state: r_candidate;
        t: t_candidate
      ;
        / Rejected: advance time and state
        r_state: r_candidate;
        t: t_candidate
      ]
    ]
  ];
  `times`types`count!(event_times; event_types; n_ev)
 };

/ Ozaki (1979) linear time O(N) recursive log-likelihood calculation
/ times: sorted arrival times [N]
/ types: integer event dimension tags [N] in {0, ..., M-1}
/ t_end: observation horizon T
logLikelihood: {[model; times; types; t_end]
  m: model`dim;
  mu: model`mu;
  alpha: model`alpha;
  beta: model`beta;
  n: count times;
  if[n == 0; :neg sum mu * t_end];
  
  / State R[i; j] for recursive decay
  r_state: {m # 0f} each til m;
  
  log_lik: 0f;
  last_t: 0f;
  
  k: 0;
  while[k < n;
    curr_t: times[k];
    curr_type: types[k];
    dt: curr_t - last_t;
    
    / Decay state to curr_t
    r_state: {r_state[x] * exp[neg beta[x] * dt]} each til m;
    
    / Evaluate intensity for dimension curr_type
    lambda_k: mu[curr_type] + sum alpha[curr_type] * r_state[curr_type];
    if[lambda_k <= 0f; lambda_k: 1e-12];
    log_lik: log_lik + log[lambda_k];
    
    / Increment jump in column curr_type
    i: 0;
    while[i < m;
      r_state[i; curr_type]: r_state[i; curr_type] + 1f;
      i: i + 1
    ];
    
    last_t: curr_t;
    k: k + 1
  ];
  
  / Subtract compensator integral: Lambda_i(T) = mu_i * T + sum_j (alpha_ij / beta_ij) * sum_k (1 - exp(-beta_ij * (T - t_k)))
  i: 0;
  while[i < m;
    log_lik: log_lik - (mu[i] * t_end);
    j: 0;
    while[j < m;
      b: beta[i; j];
      a: alpha[i; j];
      j_times: times where types = j;
      if[count j_times > 0;
        integrated_kernel: sum 1f - exp[neg b * (t_end - j_times)];
        log_lik: log_lik - ((a % b) * integrated_kernel)
      ];
      j: j + 1
    ];
    i: i + 1
  ];
  log_lik
 };

/ Maximum Likelihood Estimation (MLE) of 1D univariate Hawkes process parameters (mu, alpha, beta)
fitUnivariate: {[times; t_end]
  n: count times;
  types: n # 0;
  
  / Objective function to minimize: negative log-likelihood
  obj: {[params; times; types; t_end]
    / params: [log(mu), log(alpha), log(beta)]
    mu_val: exp params[0];
    alpha_val: exp params[1];
    beta_val: exp params[2];
    
    / Enforce stationarity penalty: alpha < beta
    penalty: $[alpha_val >= (0.999 * beta_val); 1e6 * (alpha_val - beta_val); 0f];
    
    m_test: initModel[enlist mu_val; enlist enlist alpha_val; enlist enlist beta_val];
    ll: logLikelihood[m_test; times; types; t_end];
    (neg ll) + penalty
  };
  
  / Initial parameter guess based on empirical event rate
  mean_rate: (1f * n) % t_end;
  init_params: (log[0.5 * mean_rate]; log[0.3 * mean_rate]; log[1.0]);
  
  opt_res: .optim.nelderMead[obj[; times; types; t_end]; init_params; 150; 1e-5];
  best_params: exp opt_res`x;
  
  mu_hat: best_params[0];
  alpha_hat: best_params[1];
  beta_hat: best_params[2];
  
  branching_ratio: alpha_hat % beta_hat;
  initModel[enlist mu_hat; enlist enlist alpha_hat; enlist enlist beta_hat]
 };

\d .
