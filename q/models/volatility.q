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
  a: kappa - (rho * xi * u * 0f);
  b: neg rho * xi * u;
  
  / (b*i + kappa)^2 + xi^2 * (u^2 + i*u)
  / = kappa^2 - b^2 + 2*i*kappa*b + xi^2 * u^2 + i * xi^2 * u
  re_part: (kappa * kappa) - (b * b) + (xi * xi * u * u);
  im_part: (2f * kappa * b) + (xi * xi * u);
  
  / Complex square root: sqrt(re + i*im)
  mag: sqrt[(re_part * re_part) + (im_part * im_part)];
  d_re: sqrt[0.5 * (mag + re_part)];
  d_im: sign[im_part] * sqrt[0.5 * max (0f; mag - re_part)];
  
  / g = (kappa - rho*xi*u*i - d) / (kappa - rho*xi*u*i + d)
  num_re: kappa - d_re;
  num_im: neg (rho * xi * u) - d_im;
  
  den_re: kappa + d_re;
  den_im: neg (rho * xi * u) + d_im;
  den_mag2: (den_re * den_re) + (den_im * den_im);
  
  g_re: ((num_re * den_re) + (num_im * den_im)) % den_mag2;
  g_im: ((num_im * den_re) - (num_re * den_im)) % den_mag2;
  
  / exp(-d * T)
  exp_dt_mag: exp[neg d_re * T];
  exp_dt_re: exp_dt_mag * cos[neg d_im * T];
  exp_dt_im: exp_dt_mag * sin[neg d_im * T];
  
  / 1 - g * exp(-d*T)
  g_exp_re: (g_re * exp_dt_re) - (g_im * exp_dt_im);
  g_exp_im: (g_re * exp_dt_im) + (g_im * exp_dt_re);
  
  one_minus_g_exp_re: 1f - g_exp_re;
  one_minus_g_exp_im: neg g_exp_im;
  
  / 1 - g
  one_minus_g_re: 1f - g_re;
  one_minus_g_im: neg g_im;
  
  / Complex log of (1 - g*exp(-d*T)) / (1 - g)
  ratio_den: (one_minus_g_re * one_minus_g_re) + (one_minus_g_im * one_minus_g_im);
  ratio_re: ((one_minus_g_exp_re * one_minus_g_re) + (one_minus_g_exp_im * one_minus_g_im)) % ratio_den;
  ratio_im: ((one_minus_g_exp_im * one_minus_g_re) - (one_minus_g_exp_re * one_minus_g_im)) % ratio_den;
  
  log_ratio_re: 0.5 * log[(ratio_re * ratio_re) + (ratio_im * ratio_im)];
  log_ratio_im: atan2[ratio_im; ratio_re];
  
  / C(u) = (kappa * theta / xi^2) * [(kappa - rho*xi*u*i - d)*T - 2 * log_ratio]
  c_term_re: ((kappa - d_re) * T) - (2f * log_ratio_re);
  c_term_im: (neg (rho * xi * u) - d_im) * T - (2f * log_ratio_im);
  scale_c: (kappa * theta) % (xi * xi);
  C_re: scale_c * c_term_re;
  C_im: scale_c * c_term_im;
  
  / D(u) = (kappa - rho*xi*u*i - d) / xi^2 * (1 - exp(-d*T)) / (1 - g*exp(-d*T))
  one_minus_exp_re: 1f - exp_dt_re;
  one_minus_exp_im: neg exp_dt_im;
  
  omg_den: (one_minus_g_exp_re * one_minus_g_exp_re) + (one_minus_g_exp_im * one_minus_g_exp_im);
  d_ratio_re: ((one_minus_exp_re * one_minus_g_exp_re) + (one_minus_exp_im * one_minus_g_exp_im)) % omg_den;
  d_ratio_im: ((one_minus_exp_im * one_minus_g_exp_re) - (one_minus_exp_re * one_minus_g_exp_im)) % omg_den;
  
  kd_re: kappa - d_re;
  kd_im: neg (rho * xi * u) - d_im;
  
  D_re: (((kd_re * d_ratio_re) - (kd_im * d_ratio_im)) % (xi * xi));
  D_im: (((kd_re * d_ratio_im) + (kd_im * d_ratio_re)) % (xi * xi));
  
  / Exponent: C(u) + D(u)*v0 + i * u * (ln(S0) + (r - q)*T)
  tot_re: C_re + (D_re * v0);
  tot_im: C_im + (D_im * v0) + (u * (log[S0] + ((r - q) * T)));
  
  e_re: exp[tot_re];
  (e_re * cos[tot_im]; e_re * sin[tot_im])
 };

/ Sign helper for atan2
sign: {[x] $[x < 0f; -1f; 1f]};

/ Two-argument arctangent atan2(y, x)
atan2: {[y; x]
  if[x > 0f; :atan[y % x]];
  if[x < 0f; :$[y >= 0f; atan[y % x] + .math.PI; atan[y % x] - .math.PI]];
  $[y > 0f; 0.5 * .math.PI; y < 0f; neg 0.5 * .math.PI; 0f]
 };

/ Heston European call option price via Gauss-Kronrod numerical quadrature
priceHestonCall: {[S0; K_strike; T; r; q; v0; kappa; theta; xi; rho]
  / Inversion integral: C = S0 * P1 - K * exp(-r*T) * P2
  integrand1: {[u; S0; K_strike; T; r; q; v0; kappa; theta; xi; rho]
    phi: hestonCharFn[u - (1f * 0f); S0; T; r; q; v0; kappa; theta; xi; rho];
    / Real part of exp(-i * u * ln(K)) * phi(u - i) / (i * u * phi(-i))
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
  
  / Truncated integration over [1e-4, 100]
  int1: .math.quad[integrand1[; S0; K_strike; T; r; q; v0; kappa; theta; xi; rho]; 1e-4; 100f];
  int2: .math.quad[integrand2[; S0; K_strike; T; r; q; v0; kappa; theta; xi; rho]; 1e-4; 100f];
  
  p1: 0.5 + (int1 % .math.PI);
  p2: 0.5 + (int2 % .math.PI);
  
  call_price: (S0 * exp[neg q * T] * p1) - (K_strike * exp[neg r * T] * p2);
  max (0f; call_price)
 };

/ Black-Scholes European call pricing for implied volatility inversion
bsCall: {[S0; K_strike; T; r; sigma]
  if[sigma <= 0f; :max (0f; S0 - (K_strike * exp[neg r * T]))];
  d1: (log[S0 % K_strike] + ((r + (0.5 * sigma * sigma)) * T)) % (sigma * sqrt[T]);
  d2: d1 - (sigma * sqrt[T]);
  (S0 * .math.normalCdf[d1]) - (K_strike * exp[neg r * T] * .math.normalCdf[d2])
 };

/ Invert market call price to Black-Scholes implied volatility via Brent root solver
impliedVol: {[market_price; S0; K_strike; T; r]
  intrinsic: max (0f; S0 - (K_strike * exp[neg r * T]));
  if[market_price <= intrinsic; :0.001];
  
  obj: {[vol; S0; K_strike; T; r; target]
    bsCall[S0; K_strike; T; r; vol] - target
  };
  
  .math.brentRoot[obj[; S0; K_strike; T; r; market_price]; 0.001; 4.0; 1e-5]
 };

\d .
