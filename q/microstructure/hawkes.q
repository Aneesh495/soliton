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
  lambdas: n # 0f;
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

\d .
