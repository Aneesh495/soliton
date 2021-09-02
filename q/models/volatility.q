/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/models/volatility.q
/ description: continuous-time stochastic volatility and Heston model calibration
/ ==============================================================================

\d .vol

/ Verify Feller condition for mean-reverting square-root CIR process: 2 * kappa * theta > xi^2
isFellerSatisfied: {[kappa; theta; xi]
  (2f * kappa * theta) > (xi * xi)
 };

/ Numerically stable "Little Heston Trap" formulation of Heston (1993) characteristic function
/ Avoids branch cut discontinuities in complex logarithm and square root (Albrecher et al. 2007)
/ u: frequency argument
/ S0: spot price, T: maturity, r: risk-free rate, q: dividend yield
/ v0: initial variance, kappa: mean reversion speed, theta: long term variance
/ xi: vol-of-vol, rho: spot-vol correlation
hestonCharFn: {[u; S0; T; r; q; v0; kappa; theta; xi; rho]
  / d = sqrt((rho * xi * u * i - kappa)^2 + xi^2 * (u * i + u^2))
  b: neg rho * xi * u;
  
  re_part: (kappa * kappa) - (b * b) + (xi * xi * u * u);
  im_part: (2f * kappa * b) + (xi * xi * u);
  
  mag: sqrt[(re_part * re_part) + (im_part * im_part)];
  d_re: sqrt[0.5 * (mag + re_part)];
  d_im: sign[im_part] * sqrt[0.5 * max (0f; mag - re_part)];
  
  num_re: kappa - d_re;
  num_im: neg (rho * xi * u) - d_im;
  
  den_re: kappa + d_re;
  den_im: neg (rho * xi * u) + d_im;
  den_mag2: (den_re * den_re) + (den_im * den_im);
  
  g_re: ((num_re * den_re) + (num_im * den_im)) % den_mag2;
  g_im: ((num_im * den_re) - (num_re * den_im)) % den_mag2;
  
  exp_dt_mag: exp[neg d_re * T];
  exp_dt_re: exp_dt_mag * cos[neg d_im * T];
  exp_dt_im: exp_dt_mag * sin[neg d_im * T];
  
  g_exp_re: (g_re * exp_dt_re) - (g_im * exp_dt_im);
  g_exp_im: (g_re * exp_dt_im) + (g_im * exp_dt_re);
  
  one_minus_g_exp_re: 1f - g_exp_re;
  one_minus_g_exp_im: neg g_exp_im;
  
  one_minus_g_re: 1f - g_re;
  one_minus_g_im: neg g_im;
  
  ratio_den: (one_minus_g_re * one_minus_g_re) + (one_minus_g_im * one_minus_g_im);
  ratio_re: ((one_minus_g_exp_re * one_minus_g_re) + (one_minus_g_exp_im * one_minus_g_im)) % ratio_den;
  ratio_im: ((one_minus_g_exp_im * one_minus_g_re) - (one_minus_g_exp_re * one_minus_g_im)) % ratio_den;
  
  log_ratio_re: 0.5 * log[(ratio_re * ratio_re) + (ratio_im * ratio_im)];
  log_ratio_im: atan2[ratio_im; ratio_re];
  
  c_term_re: ((kappa - d_re) * T) - (2f * log_ratio_re);
  c_term_im: (neg (rho * xi * u) - d_im) * T - (2f * log_ratio_im);
  scale_c: (kappa * theta) % (xi * xi);
  C_re: scale_c * c_term_re;
  C_im: scale_c * c_term_im;
  
  one_minus_exp_re: 1f - exp_dt_re;
  one_minus_exp_im: neg exp_dt_im;
  
  omg_den: (one_minus_g_exp_re * one_minus_g_exp_re) + (one_minus_g_exp_im * one_minus_g_exp_im);
  d_ratio_re: ((one_minus_exp_re * one_minus_g_exp_re) + (one_minus_exp_im * one_minus_g_exp_im)) % omg_den;
  d_ratio_im: ((one_minus_exp_im * one_minus_g_exp_re) - (one_minus_exp_re * one_minus_g_exp_im)) % omg_den;
  
  kd_re: kappa - d_re;
  kd_im: neg (rho * xi * u) - d_im;
  
  D_re: (((kd_re * d_ratio_re) - (kd_im * d_ratio_im)) % (xi * xi));
  D_im: (((kd_re * d_ratio_im) + (kd_im * d_ratio_re)) % (xi * xi));
  
  tot_re: C_re + (D_re * v0);
  tot_im: C_im + (D_im * v0) + (u * (log[S0] + ((r - q) * T)));
  
  e_re: exp[tot_re];
  (e_re * cos[tot_im]; e_re * sin[tot_im])
 };

sign: {[x] $[x < 0f; -1f; 1f]};

atan2: {[y; x]
  if[x > 0f; :atan[y % x]];
  if[x < 0f; :$[y >= 0f; atan[y % x] + .math.PI; atan[y % x] - .math.PI]];
  $[y > 0f; 0.5 * .math.PI; y < 0f; neg 0.5 * .math.PI; 0f]
 };

priceHestonCall: {[S0; K_strike; T; r; q; v0; kappa; theta; xi; rho]
  integrand1: {[u; S0; K_strike; T; r; q; v0; kappa; theta; xi; rho]
    phi: hestonCharFn[u - (1f * 0f); S0; T; r; q; v0; kappa; theta; xi; rho];
    re: phi[0]; im: phi[1];
    cos_uk: cos[u * log[K_strike]]; sin_uk: sin[u * log[K_strike]];
    num_re: (re * cos_uk) + (im * sin_uk);
    num_im: (im * cos_uk) - (re * sin_uk);
    num_im % (u * S0 * exp[(r - q) * T])
  };
  
  integrand2: {[u; S0; K_strike; T; r; q; v0; kappa; theta; xi; rho]
    phi: hestonCharFn[u; S0; T; r; q; v0; kappa; theta; xi; rho];
    re: phi[0]; im: phi[1];
    cos_uk: cos[u * log[K_strike]]; sin_uk: sin[u * log[K_strike]];
    num_re: (re * cos_uk) + (im * sin_uk);
    num_im: (im * cos_uk) - (re * sin_uk);
    num_im % u
  };
  
  int1: .math.quad[integrand1[; S0; K_strike; T; r; q; v0; kappa; theta; xi; rho]; 1e-4; 100f];
  int2: .math.quad[integrand2[; S0; K_strike; T; r; q; v0; kappa; theta; xi; rho]; 1e-4; 100f];
  
  p1: 0.5 + (int1 % .math.PI);
  p2: 0.5 + (int2 % .math.PI);
  
  call_price: (S0 * exp[neg q * T] * p1) - (K_strike * exp[neg r * T] * p2);
  max (0f; call_price)
 };

bsCall: {[S0; K_strike; T; r; sigma]
  if[sigma <= 0f; :max (0f; S0 - (K_strike * exp[neg r * T]))];
  d1: (log[S0 % K_strike] + ((r + (0.5 * sigma * sigma)) * T)) % (sigma * sqrt[T]);
  d2: d1 - (sigma * sqrt[T]);
  (S0 * .math.normalCdf[d1]) - (K_strike * exp[neg r * T] * .math.normalCdf[d2])
 };

impliedVol: {[market_price; S0; K_strike; T; r]
  intrinsic: max (0f; S0 - (K_strike * exp[neg r * T]));
  if[market_price <= intrinsic; :0.001];
  
  obj: {[vol; S0; K_strike; T; r; target]
    bsCall[S0; K_strike; T; r; vol] - target
  };
  
  .math.brentRoot[obj[; S0; K_strike; T; r; market_price]; 0.001; 4.0; 1e-5]
 };

/ Davies-Harte algorithm for exact Fractional Gaussian Noise (fGn) path synthesis
/ n: number of increments; H: Hurst parameter in (0, 1) (rough volatility if H < 0.5)
simulateFGN: {[n; H]
  two_H: 2f * H;
  / Autocovariance sequence of fGn: gamma(k) = 0.5 * (|k+1|^{2H} - 2*|k|^{2H} + |k-1|^{2H})
  k_vec: til n;
  gamma_k: {
    k: 1f * x;
    if[k == 0f; :1f];
    0.5 * ((abs[k + 1f] xexp two_H) - (2f * (abs[k] xexp two_H)) + (abs[k - 1f] xexp two_H))
  } each k_vec;
  
  / Circulant embedding vector of length 2 * (n - 1)
  c_len: 2 * (n - 1);
  c_vec: gamma_k, (reverse 1 _ -1 _ gamma_k);
  
  / Circulant eigenvalues via FFT
  f_res: .fft.fft[c_vec; c_len # 0f];
  eigenvals: f_res[0];
  / Ensure numerical positivity of eigenvalues
  eigenvals: max'[eigenvals; 0f];
  sqrt_eig: sqrt[eigenvals];
  
  / Independent standard normal variates
  z1: .stats.randn[c_len];
  z2: .stats.randn[c_len];
  
  w_re: sqrt_eig * z1;
  w_im: sqrt_eig * z2;
  
  / Inverse FFT to obtain correlated fGn path
  inv_res: .fft.ifft[w_re; w_im];
  n # inv_res[0]
 };

/ Simulate rough Bergomi (rBergomi) instantaneous variance paths:
/ v_t = xi_0(t) * exp(eta * W_t^H - 0.5 * eta^2 * t^{2H})
simulateRoughBergomi: {[n_steps; T; H; eta; xi0]
  dt: T % (1f * n_steps);
  fgn: simulateFGN[n_steps; H];
  / Cumulative fractional Brownian motion path
  fbm_increments: (dt xexp H) * fgn;
  fbm_path: 0f, sums fbm_increments;
  
  times: dt * til n_steps + 1;
  variance_path: {
    t: times[x];
    drift_adj: 0.5 * eta * eta * (t xexp (2f * H));
    xi0 * exp[(eta * fbm_path[x]) - drift_adj]
  } each til n_steps + 1;
  `times`variance`H`eta!(times; variance_path; H; eta)
 };

/ Realized Variance (RV): sum of squared intraday log returns
realizedVariance: {[prices]
  dp: 1 _ (log[prices] - prev log[prices]);
  sum dp * dp
 };

/ Realized Bipower Variation (BV): jump-robust estimator of integrated variance
/ BV = (pi / 2) * (n / (n - 1)) * sum(|r_i| * |r_{i-1}|)
realizedBipowerVar: {[prices]
  n: count prices;
  if[n < 3; :0f];
  r: 1 _ (log[prices] - prev log[prices]);
  abs_r: abs r;
  scale: (0.5 * .math.PI) * (1f * (n - 1)) % (1f * (n - 2));
  scale * sum (1 _ abs_r) * -1 _ abs_r
 };

/ Jump variation separation via Lee-Mykland test: J = max(RV - BV, 0)
jumpVariation: {[prices]
  rv: realizedVariance[prices];
  bv: realizedBipowerVar[prices];
  max (0f; rv - bv)
 };

/ Two-Scale Realized Volatility (TSRV) for filtering market microstructure noise
/ TSRV = RV_slow - (n_slow / n_fast) * RV_fast
twoScaleRV: {[prices; k_subsample]
  n: count prices;
  if[n <= (2 * k_subsample); :realizedVariance[prices]];
  
  / Fast scale RV (all ticks)
  rv_fast: realizedVariance[prices];
  
  / Slow scale RV averaged across k_subsample interleaved grids
  slow_rvs: ();
  sub: 0;
  while[sub < k_subsample;
    sub_prices: prices[sub + (k_subsample * til `long$floor (n - sub) % k_subsample)];
    slow_rvs: slow_rvs, enlist realizedVariance[sub_prices];
    sub: sub + 1
  ];
  rv_slow_avg: .stats.mean[slow_rvs];
  
  ratio: (1f * (n - k_subsample + 1)) % (1f * k_subsample * n);
  tsrv_val: rv_slow_avg - (ratio * rv_fast);
  max (0f; tsrv_val)
 };

\d .
