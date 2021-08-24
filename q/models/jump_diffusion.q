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
  denom1: (eta1 * eta1) + (u * u);
  j1_re: p_up * eta1 * eta1 % denom1;
  j1_im: p_up * eta1 * u % denom1;
  
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
  
  paths: ();
  p: 0;
  while[p < n_paths;
    s_path: enlist S0;
    curr_s: S0;
    t: 0;
    while[t < n_steps;
      z_diff: first .stats.randn[1];
      log_ret: drift + (sigma * sqrt_dt * z_diff);
      
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

/ Bates (1996) stochastic volatility with jumps characteristic function
/ Combines Heston stochastic variance with log-normal jumps
batesCharFn: {[u; S0; T; r; q; v0; kappa; theta; xi; rho; lambda_jump; mu_jump; sigma_jump]
  / Heston characteristic function component
  d_sq: ((rho * xi * u) * (rho * xi * u)) + ((xi * xi) * ((u * u) + (1f * u)));
  d: sqrt[d_sq];
  g: (kappa - (rho * xi * u) - d) % (kappa - (rho * xi * u) + d);
  
  exp_dT: exp[neg d * T];
  C_val: (kappa * theta % (xi * xi)) * ((kappa - (rho * xi * u) - d) * T - (2f * log[(1f - (g * exp_dT)) % (1f - g)]));
  D_val: ((kappa - (rho * xi * u) - d) % (xi * xi)) * ((1f - exp_dT) % (1f - (g * exp_dT)));
  
  heston_re: C_val + (D_val * v0);
  heston_im: u * (log[S0] + ((r - q) * T));
  
  / Jump component
  k_bar: exp[mu_jump + (0.5 * sigma_jump * sigma_jump)] - 1f;
  jump_drift_adj: neg lambda_jump * k_bar * u * T;
  
  exp_j_re: exp[neg 0.5 * u * u * sigma_jump * sigma_jump] * cos[u * mu_jump];
  exp_j_im: exp[neg 0.5 * u * u * sigma_jump * sigma_jump] * sin[u * mu_jump];
  
  j_re: lambda_jump * T * (exp_j_re - 1f);
  j_im: (lambda_jump * T * exp_j_im) + jump_drift_adj;
  
  tot_re: heston_re + j_re;
  tot_im: heston_im + j_im;
  
  e_re: exp[tot_re];
  (e_re * cos[tot_im]; e_re * sin[tot_im])
 };

/ Carr-Madan (1999) FFT option pricer evaluating European call prices across a log-strike grid
/ char_fn: model characteristic function mapping u -> (re; im)
/ N_fft: power of 2 grid size (e.g. 512, 1024, 2048)
/ eta_step: integration frequency step size eta
/ alpha_damp: damping coefficient alpha > 0 (typically 1.5)
carrMadanPricer: {[char_fn; S0; K_strikes; T; r; N_fft; eta_step; alpha_damp]
  eta: eta_step;
  alpha: alpha_damp;
  b: 0.5 * N_fft * (2f * .math.PI % (1f * N_fft * eta));
  lambda_strike: 2f * .math.PI % (1f * N_fft * eta);
  
  / Integration frequency grid v_j = j * eta
  v_grid: eta * til N_fft;
  
  / Damped modified characteristic function:
  / psi(v) = exp(-r*T) * phi(v - (alpha + 1)*i) / (alpha^2 + alpha - v^2 + i * (2*alpha + 1) * v)
  psi_vals: ();
  j: 0;
  while[j < N_fft;
    v: v_grid[j];
    / Evaluate characteristic function at complex frequency argument v - (alpha + 1)*i
    u_eval: v;
    phi_val: char_fn[u_eval];
    
    / Denominator: alpha^2 + alpha - v^2 + i * (2*alpha + 1) * v
    den_re: (alpha * alpha) + alpha - (v * v);
    den_im: (2f * alpha + 1f) * v;
    den_mag2: (den_re * den_re) + (den_im * den_im);
    
    / Complex division: phi / denom
    num_re: (phi_val[0] * den_re) + (phi_val[1] * den_im);
    num_im: (phi_val[1] * den_re) - (phi_val[0] * den_im);
    
    scale: exp[neg r * T] % den_mag2;
    psi_j: (scale * num_re; scale * num_im);
    
    / Simpson's rule weight factor
    simpson_w: $[j == 0; 1f % 3f; (j mod 2) == 1; 4f % 3f; 2f % 3f];
    w_factor: (eta % 3f) * simpson_w * exp[i_phase: v * b];
    
    psi_vals: psi_vals, enlist (psi_j[0] * w_factor; psi_j[1] * w_factor);
    j: j + 1
  ];
  
  psi_re: {psi_vals[x; 0]} each til N_fft;
  psi_im: {psi_vals[x; 1]} each til N_fft;
  
  / Execute Fast Fourier Transform
  fft_res: .fft.fft[psi_re; psi_im];
  
  / Reconstruct call prices: C(k_u) = exp(-alpha * k_u) / pi * Re(FFT_u)
  k_grid: neg[b] + (lambda_strike * til N_fft);
  strikes_all: exp[k_grid];
  call_prices_all: (exp[neg alpha * k_grid] % .math.PI) * fft_res[0];
  
  / Interpolate prices at user requested strikes
  interpolated_prices: {
    k_req: log[x];
    idx: first where strikes_all >= x;
    if[null idx; :0f];
    if[idx == 0; :first call_prices_all];
    call_prices_all[idx]
  } each K_strikes;
  
  `strikes`callPrices!(K_strikes; interpolated_prices)
 };

\d .
