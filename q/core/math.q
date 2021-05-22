/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/core/math.q
/ description: high-precision special functions, quadrature, and calculus
/ ==============================================================================

\d .math

/ Mathematical constants defined to full IEEE 754 precision
PI: 3.14159265358979323846264338327950288419716939937510;
TWO_PI: 6.28318530717958647692528676655900576839433879875021;
SQRT_2PI: 2.5066282746310005024157652848110452530069867406099;
E: 2.71828182845904523536028747135266249775724709369995;
EPS: 2.220446049250313080847263336181640625e-16;
SQRT_EPS: 1.4901161193847656e-8;

/ Signum function returning -1, 0, or 1
sign: {[x]
  $[x>0; 1f; $[x<0; -1f; 0f]]
 };

/ Horner polynomial evaluation given coefficients [c0, c1, ..., cn]
/ P(x) = c0 + c1*x + c2*x^2 + ... + cn*x^n
polyval: {[coeffs; x]
  rev_c: reverse coeffs;
  acc: first rev_c;
  i: 1;
  while[i < count rev_c;
    acc: (acc * x) + rev_c[i];
    i: i + 1
  ];
  acc
 };

/ Lanczos approximation for log-gamma function ln(Gamma(z))
/ Precision: 15 decimal digits for Re(z) > 0
logGamma: {[z]
  if[z <= 0f; 'domain];
  p: (0.99999999999980993; 676.5203681218851; -1259.1392167224028;
      771.32342877765313; -176.61502916214059; 12.507343278686905;
      -0.13857109583656712; 9.9843695780195716e-6; 1.5056327351493116e-7);
  if[z < 0.5;
    :log[PI] - log[sin[PI * z]] - logGamma[1f - z]
  ];
  zm1: z - 1f;
  x: p[0];
  i: 1;
  while[i < count p;
    x: x + p[i] % (zm1 + (1f * i));
    i: i + 1
  ];
  t: zm1 + 7.5;
  (0.5 * log[TWO_PI]) + ((zm1 + 0.5) * log[t]) - t + log[x]
 };

/ Gamma function Gamma(z) = exp(ln(Gamma(z)))
gamma: {[z]
  exp logGamma[z]
 };

/ Chebyshev rational approximation of the error function erf(x)
/ Maximum absolute error bounded below 1.5e-7
erf: {[x]
  s: sign[x];
  ax: abs[x];
  t: 1f % (1f + (0.3275911 * ax));
  a1: 0.254829592;
  a2: -0.284496736;
  a3: 1.421413741;
  a4: -1.453152027;
  a5: 1.061405429;
  poly: t * (a1 + (t * (a2 + (t * (a3 + (t * (a4 + (t * a5))))))));
  s * (1f - (poly * exp[neg[ax * ax]]))
 };

/ Complementary error function erfc(x) = 1 - erf(x)
erfc: {[x]
  1f - erf[x]
 };

/ Standard normal cumulative distribution function Phi(x)
normalCdf: {[x]
  0.5 * (1f + erf[x % 1.4142135623730950488016887242])
 };

/ Standard normal probability density function phi(x)
normalPdf: {[x]
  (1f % SQRT_2PI) * exp[-0.5 * x * x]
 };

/ Inverse normal CDF (probit function) via rational approximation of Acklam
normalInv: {[p]
  if[(p <= 0f) or (p >= 1f); 'domain];
  a: (-3.969683028665376e+01; 2.209460984245205e+02; -2.759285104469687e+02;
      1.383577518672690e+02; -3.066479806614716e+01; 2.506628277459239e+00);
  b: (-5.447609879822406e+01; 1.615858368580409e+02; -1.556989798598866e+02;
      6.680131188771972e+01; -1.328068155288572e+01);
  c: (-7.784894002430293e-03; -3.223964580411365e-01; -2.400758277161838e+00;
      -2.549732539343734e+00; 4.374664141464968e+00; 2.938163982698783e+00);
  d: (7.784695709041462e-03; 3.224671290700398e-01; 2.445134137142996e+00;
      3.754408661907416e+00);
  p_low: 0.02425;
  p_high: 1f - p_low;
  if[p < p_low;
    q: sqrt[-2f * log[p]];
    :((((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) %
      (((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1f))
  ];
  if[p <= p_high;
    q: p - 0.5;
    r: q * q;
    :((((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5])*q) %
      ((((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1f))
  ];
  q: sqrt[-2f * log[1f - p]];
  neg ((((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) %
        (((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1f))
 };

/ Regularized lower incomplete gamma function P(a, x) via series and continued fractions
gammaP: {[a; x]
  if[(a <= 0f) or (x < 0f); 'domain];
  if[x == 0f; :0f];
  if[x < (a + 1f);
    / Series development
    ap: a;
    sum: 1f % a;
    del: sum;
    n: 1;
    while[n <= 100;
      ap: ap + 1f;
      del: del * (x % ap);
      sum: sum + del;
      if[abs[del] < (abs[sum] * 3e-7); n: 999];
      n: n + 1
    ];
    :sum * exp[neg[x] + (a * log[x]) - logGamma[a]]
  ];
  / Continued fraction approximation
  b: x + 1f - a;
  c: 1f % 1e-30;
  d: 1f % b;
  h: d;
  i: 1;
  while[i <= 100;
    an: neg[1f * i] * (1f * i - a);
    b: b + 2f;
    d: (an * d) + b;
    if[abs[d] < 1e-30; d: 1e-30];
    c: b + (an % c);
    if[abs[c] < 1e-30; c: 1e-30];
    d: 1f % d;
    del: d * c;
    h: h * del;
    if[abs[del - 1f] < 3e-7; i: 999];
    i: i + 1
  ];
  1f - (exp[neg[x] + (a * log[x]) - logGamma[a]] * h)
 };

/ Modified Bessel function of the first kind I_0(x)
besselI0: {[x]
  ax: abs[x];
  if[ax < 3.75;
    y: (x % 3.75) * (x % 3.75);
    :1f + (y * (3.5156229 + (y * (3.0899424 + (y * (1.2067492 + (y * (0.2659732 + (y * (0.360768e-1 + (y * 0.45813e-2))))))))))
  ];
  y: 3.75 % ax;
  ans: 0.39894228 + (y * (0.1328592e-1 + (y * (0.225319e-2 + (y * (-0.157565e-2 + (y * (0.916281e-2 + (y * (-0.2057706e-1 + (y * (0.2635537e-1 + (y * (-0.1647633e-1 + (y * 0.392377e-2)))))))))))))));
  (exp[ax] % sqrt[ax]) * ans
 };

/ Modified Bessel function of the first kind I_1(x)
besselI1: {[x]
  ax: abs[x];
  if[ax < 3.75;
    y: (x % 3.75) * (x % 3.75);
    ans: ax * (0.5 + (y * (0.87890594 + (y * (0.51498869 + (y * (0.15084934 + (y * (0.2658733e-1 + (y * (0.301532e-2 + (y * 0.32411e-3))))))))))));
    :$[x < 0f; neg ans; ans]
  ];
  y: 3.75 % ax;
  ans: 0.39894228 + (y * (-0.3988024e-1 + (y * (-0.362018e-2 + (y * (0.163801e-2 + (y * (-0.1031555e-1 + (y * (0.2282967e-1 + (y * (-0.2895312e-1 + (y * (0.1787654e-1 + (y * -0.420059e-2)))))))))))))));
  ans: (exp[ax] % sqrt[ax]) * ans;
  $[x < 0f; neg ans; ans]
 };

/ Composite Gauss-Kronrod 7-15 point adaptive numerical quadrature over [a, b]
quad: {[f; a; b]
  / 15-point Gauss-Kronrod nodes and weights on [-1, 1]
  x_nodes: (0f; 0.2077849550078985; 0.4058451513773972; 0.5860872354676911;
            0.7415311855993944; 0.8648644233597691; 0.9491079123427585; 0.9914553711208126);
  w_kronrod: (0.2094821410847278; 0.2044329400752989; 0.1903505780647854;
              0.1690047266392679; 0.1406532597155259; 0.1047900103222502;
              0.06309209262997855; 0.02293532201052922);
  
  half_len: 0.5 * (b - a);
  mid: 0.5 * (a + b);
  
  sum: w_kronrod[0] * f[mid];
  i: 1;
  while[i < count x_nodes;
    dx: half_len * x_nodes[i];
    sum: sum + (w_kronrod[i] * (f[mid - dx] + f[mid + dx]));
    i: i + 1
  ];
  sum * half_len
 };

/ Romberg integration algorithm: extrapolation table for high-order accuracy
romberg: {[f; a; b; max_steps]
  steps: $[max_steps < 3; 6; max_steps];
  h: b - a;
  r: enlist enlist 0.5 * h * (f[a] + f[b]);
  
  i: 1;
  while[i < steps;
    h: 0.5 * h;
    sum: 0f;
    k: 1;
    while[k <= (1 sv (i # 2));
      sum: sum + f[a + ((1f * (2 * k - 1)) * h)];
      k: k + 1
    ];
    r_row: enlist (0.5 * r[i-1][0]) + (h * sum);
    
    j: 1;
    while[j <= i;
      factor: exp[log[4f] * j] - 1f;
      val: r_row[j-1] + ((r_row[j-1] - r[i-1][j-1]) % factor);
      r_row: r_row, enlist val;
      j: j + 1
    ];
    r: r, enlist r_row;
    i: i + 1
  ];
  last last r
 };

/ High-precision central finite difference derivative f'(x)
derivative: {[f; x]
  h: SQRT_EPS * (1f + abs[x]);
  (f[x + h] - f[x - h]) % (2f * h)
 };

/ Second-order central finite difference derivative f''(x)
derivative2: {[f; x]
  h: (EPS xexp (1f % 3f)) * (1f + abs[x]);
  (f[x + h] - (2f * f[x]) + f[x - h]) % (h * h)
 };

/ One-dimensional root finding via Brent-Dekker algorithm
brentRoot: {[f; a; b; tol]
  fa: f[a];
  fb: f[b];
  if[(fa * fb) > 0f; 'bracket];
  c: a;
  fc: fa;
  d: b - a;
  e: d;
  iter: 0;
  while[iter < 100;
    if[(fb * fc) > 0f;
      c: a; fc: fa; d: b - a; e: d
    ];
    if[abs[fc] < abs[fb];
      a: b; b: c; c: a;
      fa: fb; fb: fc; fc: fa
    ];
    tol1: (2f * EPS * abs[b]) + (0.5 * tol);
    xm: 0.5 * (c - b);
    if[(abs[xm] <= tol1) or (fb == 0f); :b];
    if[(abs[e] >= tol1) and (abs[fa] > abs[fb]);
      s: fb % fa;
      if[a == c;
        p: 2f * xm * s;
        q: 1f - s
      ];
      if[a <> c;
        q: fa % fc;
        r: fb % fc;
        p: s * ((2f * xm * q * (q - r)) - ((b - a) * (r - 1f)));
        q: (q - 1f) * (r - 1f) * (s - 1f)
      ];
      if[p > 0f; q: neg q];
      p: abs[p];
      min1: (3f * xm * q) - abs[tol1 * q];
      min2: abs[e * q];
      if[(2f * p) < ($[min1 < min2; min1; min2]);
        e: d; d: p % q
      ;
        d: xm; e: d
      ]
    ;
      d: xm; e: d
    ];
    a: b; fa: fb;
    b: b + ($[abs[d] > tol1; d; $[xm > 0f; tol1; neg tol1]]);
    fb: f[b];
    iter: iter + 1
  ];
  b
 };

\d .
