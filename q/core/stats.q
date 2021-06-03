/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/core/stats.q
/ description: probability distributions, empirical cumulants, and copula theory
/ ==============================================================================

\d .stats

/ Mean of vector
mean: {[x]
  (sum x) % count x
 };

/ Unbiased sample variance s^2
var: {[x]
  n: count x;
  if[n < 2; :0f];
  m: mean[x];
  d: x - m;
  (sum d * d) % (n - 1)
 };

/ Sample standard deviation s
std: {[x]
  sqrt var[x]
 };

/ Sample skewness (Fisher-Pearson standardized 3rd moment)
skew: {[x]
  n: count x;
  if[n < 3; :0f];
  m: mean[x];
  s: std[x];
  if[s == 0f; :0f];
  z: (x - m) % s;
  factor: (1f * n) % ((1f * (n - 1)) * (1f * (n - 2)));
  factor * sum z * z * z
 };

/ Sample excess kurtosis (Fisher standardized 4th moment minus 3)
kurt: {[x]
  n: count x;
  if[n < 4; :0f];
  m: mean[x];
  s: std[x];
  if[s == 0f; :0f];
  z: (x - m) % s;
  c1: (1f * n * (n + 1)) % ((1f * (n - 1)) * (1f * (n - 2)) * (1f * (n - 3)));
  c2: (3f * (1f * (n - 1)) * (1f * (n - 1))) % ((1f * (n - 2)) * (1f * (n - 3)));
  (c1 * sum z * z * z * z) - c2
 };

/ Jarque-Bera test for normality: JB = (n/6) * (S^2 + (K^2 / 4))
/ Asymptotically distributed as Chi-square with 2 degrees of freedom
jarqueBera: {[x]
  n: count x;
  s: skew[x];
  k: kurt[x];
  jb_stat: ((1f * n) % 6f) * ((s * s) + (0.25 * k * k));
  / p-value from Chi-square distribution with 2 df: P(X >= jb) = exp(-0.5 * jb)
  p_val: exp[-0.5 * jb_stat];
  `stat`pValue`isNormal!(jb_stat; p_val; p_val > 0.05)
 };

/ Student-t probability density function with nu degrees of freedom
studentPdf: {[x; nu]
  c: exp[.math.logGamma[0.5 * (nu + 1f)] - .math.logGamma[0.5 * nu]];
  denom: sqrt[nu * .math.PI] * ((1f + ((x * x) % nu)) xexp (0.5 * (nu + 1f)));
  c % denom
 };

/ Regularized incomplete beta function I_x(a, b) via continued fraction
betaIncomplete: {[x; a; b]
  if[(x < 0f) or (x > 1f); 'domain];
  if[(x == 0f) or (x == 1f); :x];
  
  / Pre-factor: x^a * (1 - x)^b / (a * Beta(a, b))
  lbeta: .math.logGamma[a] + .math.logGamma[b] - .math.logGamma[a + b];
  front: exp[(a * log[x]) + (b * log[1f - x]) - lbeta];
  
  / Continued fraction expansion of Lentz
  if[x < ((a + 1f) % (a + b + 2f));
    / Direct continued fraction
    c: 1f; d: 1f % (1f - ((a + b) * x % (a + 1f)));
    h: d;
    m: 1;
    while[m <= 100;
      m2: 2 * m;
      / Even step
      aa: (1f * m * (b - m) * x) % ((a + (1f * m2) - 1f) * (a + (1f * m2)));
      d: 1f + (aa * d); if[abs[d] < 1e-30; d: 1e-30];
      c: 1f + (aa % c); if[abs[c] < 1e-30; c: 1e-30];
      d: 1f % d; h: h * (d * c);
      / Odd step
      aa: neg (1f * (a + m) * (a + b + m) * x) % ((a + (1f * m2)) * (a + (1f * m2) + 1f));
      d: 1f + (aa * d); if[abs[d] < 1e-30; d: 1e-30];
      c: 1f + (aa % c); if[abs[c] < 1e-30; c: 1e-30];
      d: 1f % d; del: d * c; h: h * del;
      if[abs[del - 1f] < 1e-8; m: 999];
      m: m + 1
    ];
    :(front % a) * h
  ];
  1f - betaIncomplete[1f - x; b; a]
 };

/ Student-t cumulative distribution function
studentCdf: {[x; nu]
  if[x == 0f; :0.5];
  val: nu % (nu + (x * x));
  ib: 0.5 * betaIncomplete[val; 0.5 * nu; 0.5];
  $[x > 0f; 1f - ib; ib]
 };

/ Chi-squared probability density function with k degrees of freedom
chi2Pdf: {[x; k]
  if[x <= 0f; :0f];
  c: 1f % ((exp[0.5 * k * log[2f]]) * exp[.math.logGamma[0.5 * k]]);
  c * (x xexp ((0.5 * k) - 1f)) * exp[-0.5 * x]
 };

/ Chi-squared cumulative distribution function
chi2Cdf: {[x; k]
  if[x <= 0f; :0f];
  .math.gammaP[0.5 * k; 0.5 * x]
 };

/ Gamma probability density function with shape alpha and scale theta
gammaPdf: {[x; alpha; theta]
  if[x <= 0f; :0f];
  denom: (theta xexp alpha) * exp[.math.logGamma[alpha]];
  ((x xexp (alpha - 1f)) * exp[neg[x] % theta]) % denom
 };

/ Gamma cumulative distribution function
gammaCdf: {[x; alpha; theta]
  if[x <= 0f; :0f];
  .math.gammaP[alpha; x % theta]
 };

/ Pareto Type I probability density function: f(x) = (alpha * x_m^alpha) / x^(alpha + 1)
paretoPdf: {[x; xm; alpha]
  if[x < xm; :0f];
  (alpha * (xm xexp alpha)) % (x xexp (alpha + 1f))
 };

/ Pareto Type I cumulative distribution function: F(x) = 1 - (xm / x)^alpha
paretoCdf: {[x; xm; alpha]
  if[x < xm; :0f];
  1f - ((xm % x) xexp alpha)
 };

/ Generalized Error Distribution (GED) density with parameter beta (leptokurtic if beta < 2)
gedPdf: {[x; beta]
  lambda: sqrt[(exp[.math.logGamma[1f % beta] - .math.logGamma[3f % beta]])];
  c: beta % (2f * lambda * exp[.math.logGamma[1f % beta]]);
  c * exp[neg (abs[x % lambda] xexp beta)]
 };

/ Bivariate Gaussian Copula density: c(u, v; rho)
copulaGaussianPdf: {[u; v; rho]
  x: .math.normalInv[u];
  y: .math.normalInv[v];
  denom: sqrt[1f - (rho * rho)];
  num: neg ((rho * rho * ((x * x) + (y * y))) - (2f * rho * x * y)) % (2f * (1f - (rho * rho)));
  (1f % denom) * exp[num]
 };

/ Clayton Copula CDF: C(u, v; theta) = (u^(-theta) + v^(-theta) - 1)^(-1/theta)
copulaClaytonCdf: {[u; v; theta]
  if[theta <= 0f; 'domain];
  term: (u xexp neg theta) + (v xexp neg theta) - 1f;
  max (0f; term) xexp (neg 1f % theta)
 };

/ Gumbel Copula CDF: C(u, v; theta) = exp(-((-ln u)^theta + (-ln v)^theta)^(1/theta))
copulaGumbelCdf: {[u; v; theta]
  if[theta < 1f; 'domain];
  term: ((neg log[u]) xexp theta) + ((neg log[v]) xexp theta);
  exp[neg (term xexp (1f % theta))]
 };

/ Empirical Kendall's tau rank correlation between vectors x and y
kendallTau: {[x; y]
  n: count x;
  concordant: 0;
  discordant: 0;
  i: 0;
  while[i < n - 1;
    dx: x[i] - (x[(i + 1) + til n - 1 - i]);
    dy: y[i] - (y[(i + 1) + til n - 1 - i]);
    prod: dx * dy;
    concordant: concordant + sum prod > 0f;
    discordant: discordant + sum prod < 0f;
    i: i + 1
  ];
  total_pairs: 0.5 * (1f * n * (n - 1));
  (1f * (concordant - discordant)) % total_pairs
 };

\d .
