/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/models/jump_diffusion.q
/ description: Merton and Kou jump-diffusion processes, characteristic functions
/ ==============================================================================

\d .jd

/ Merton (1976) log-normal jump-diffusion characteristic function
/ u: Fourier frequency argument
/ S0: initial spot price
/ T: maturity in years
/ r: risk-free rate
/ q: continuous dividend yield
/ sigma: diffusion volatility
/ lambda_jump: Poisson jump intensity lambda
/ mu_jump: mean log jump size
/ sigma_jump: volatility of log jump size
mertonCharFn: {[u; S0; T; r; q; sigma; lambda_jump; mu_jump; sigma_jump]
  / Expected percentage jump size: k_bar = exp(mu_j + 0.5 * sigma_j^2) - 1
  k_bar: exp[mu_jump + (0.5 * sigma_jump * sigma_jump)] - 1f;
  
  / Drift compensator: omega = r - q - lambda * k_bar - 0.5 * sigma^2
  omega: r - q - (lambda_jump * k_bar) - (0.5 * sigma * sigma);
  
  / Continuous diffusion component: i * u * omega * T - 0.5 * u^2 * sigma^2 * T
  diff_re: neg 0.5 * u * u * sigma * sigma * T;
  diff_im: u * (log[S0] + (omega * T));
  
  / Jump component: lambda * T * (exp(i * u * mu_j - 0.5 * u^2 * sigma_j^2) - 1)
  exp_jump_re: exp[neg 0.5 * u * u * sigma_jump * sigma_jump] * cos[u * mu_jump];
  exp_jump_im: exp[neg 0.5 * u * u * sigma_jump * sigma_jump] * sin[u * mu_jump];
  
  jump_re: lambda_jump * T * (exp_jump_re - 1f);
  jump_im: lambda_jump * T * exp_jump_im;
  
  total_re: diff_re + jump_re;
  total_im: diff_im + jump_im;
  
  / Complex exponentiation: exp(re + i * im)
  e_re: exp[total_re];
  (e_re * cos[total_im]; e_re * sin[total_im])
 };

/ Kou (1978) double-exponential jump-diffusion characteristic function
/ p_up: probability of positive jump
/ eta1: rate of upward jumps (> 1)
/ eta2: rate of downward jumps (> 0)
kouCharFn: {[u; S0; T; r; q; sigma; lambda_jump; p_up; eta1; eta2]
  / Expected jump size: E[exp(Y) - 1] = p * eta1 / (eta1 - 1) + (1 - p) * eta2 / (eta2 + 1) - 1
  k_bar: (p_up * eta1 % (eta1 - 1f)) + ((1f - p_up) * eta2 % (eta2 + 1f)) - 1f;
  omega: r - q - (lambda_jump * k_bar) - (0.5 * sigma * sigma);
  
  diff_re: neg 0.5 * u * u * sigma * sigma * T;
  diff_im: u * (log[S0] + (omega * T));
  
  / Jump expectation: E[exp(i * u * Y)] = p * eta1 / (eta1 - i*u) + (1-p) * eta2 / (eta2 + i*u)
  / First term: p * eta1 * (eta1 + i*u) / (eta1^2 + u^2)
  denom1: (eta1 * eta1) + (u * u);
  j1_re: p_up * eta1 * eta1 % denom1;
  j1_im: p_up * eta1 * u % denom1;
  
  / Second term: (1-p) * eta2 * (eta2 - i*u) / (eta2^2 + u^2)
  denom2: (eta2 * eta2) + (u * u);
  j2_re: (1f - p_up) * eta2 * eta2 % denom2;
  j2_im: neg (1f - p_up) * eta2 * u % denom2;
  
  jump_term_re: lambda_jump * T * (j1_re + j2_re - 1f);
  jump_term_im: lambda_jump * T * (j1_im + j2_im);
  
  total_re: diff_re + jump_term_re;
  total_im: diff_im + jump_term_im;
  
  e_re: exp[total_re];
  (e_re * cos[total_im]; e_re * sin[total_im])
 };

/ Exact Monte Carlo simulation of Merton jump-diffusion path
simulateMertonPaths: {[n_paths; n_steps; S0; T; r; q; sigma; lambda_jump; mu_jump; sigma_jump]
  dt: T % (1f * n_steps);
  sqrt_dt: sqrt[dt];
  k_bar: exp[mu_jump + (0.5 * sigma_jump * sigma_jump)] - 1f;
  drift: (r - q - (lambda_jump * k_bar) - (0.5 * sigma * sigma)) * dt;
  
  / Pre-allocate paths matrix [n_paths x (n_steps + 1)]
  paths: ();
  p: 0;
  while[p < n_paths;
    s_path: enlist S0;
    curr_s: S0;
    t: 0;
    while[t < n_steps;
      / Continuous diffusion shock
      z_diff: first .stats.randn[1];
      log_ret: drift + (sigma * sqrt_dt * z_diff);
      
      / Poisson jump count over dt
      pois_mean: lambda_jump * dt;
      u: first 1 ? 1f;
      n_jumps: $[u < pois_mean; 1; 0];
      
      if[n_jumps > 0;
        z_jump: first .stats.randn[1];
        jump_size: mu_jump + (sigma_jump * z_jump);
        log_ret: log_ret + jump_size
      ];
      
      curr_s: curr_s * exp[log_ret];
      s_path: s_path, enlist curr_s;
      t: t + 1
    ];
    paths: paths, enlist s_path;
    p: p + 1
  ];
  paths
 };

\d .
