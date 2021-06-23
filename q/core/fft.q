/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/core/fft.q
/ description: Cooley-Tukey radix-2 FFT, Fractional FFT, and spectral filtering
/ ==============================================================================

\d .fft

/ Complex addition: (xr + yr; xi + yi)
cAdd: {[x; y]
  (x[0] + y[0]; x[1] + y[1])
 };

/ Complex subtraction: (xr - yr; xi - yi)
cSub: {[x; y]
  (x[0] - y[0]; x[1] - y[1])
 };

/ Complex multiplication: (xr*yr - xi*yi; xr*yi + xi*yr)
cMul: {[x; y]
  ((x[0] * y[0]) - (x[1] * y[1]); (x[0] * y[1]) + (x[1] * y[0]))
 };

/ Complex exponential: exp(a + i*b) = exp(a) * (cos(b) + i*sin(b))
cExp: {[x]
  ea: exp[x[0]];
  (ea * cos[x[1]]; ea * sin[x[1]])
 };

/ Bit-reversal permutation index generator for N = 2^k
bitReverse: {[n]
  bits: `long$(log[1f * n] % log[2f]);
  rev: n # 0;
  i: 0;
  while[i < n;
    val: i;
    r: 0;
    b: 0;
    while[b < bits;
      r: (r * 2) + (val mod 2);
      val: val div 2;
      b: b + 1
    ];
    rev[i]: r;
    i: i + 1
  ];
  rev
 };

/ Radix-2 Cooley-Tukey Decimation-in-Time Forward FFT
/ Input: real vector xr and imaginary vector xi of length N = 2^k
/ Returns: (X_real; X_imag)
fft: {[xr; xi]
  n: count xr;
  if[n == 1; :(xr; xi)];
  
  / Bit-reversal permutation
  rev_idx: bitReverse[n];
  re: xr[rev_idx];
  im: xi[rev_idx];
  
  s: 1;
  while[s <= `long$(log[1f * n] % log[2f]);
    m: 1 sv (s # 2);        / Sub-transform length 2^s
    m2: m div 2;            / Half-length
    theta: neg .math.TWO_PI % (1f * m);
    
    k: 0;
    while[k < n;
      j: 0;
      while[j < m2;
        w_re: cos[j * theta];
        w_im: sin[j * theta];
        
        u_re: re[k + j];
        u_im: im[k + j];
        
        t_re: (w_re * re[k + j + m2]) - (w_im * im[k + j + m2]);
        t_im: (w_re * im[k + j + m2]) + (w_im * re[k + j + m2]);
        
        re[k + j]: u_re + t_re;
        im[k + j]: u_im + t_im;
        
        re[k + j + m2]: u_re - t_re;
        im[k + j + m2]: u_im - t_im;
        
        j: j + 1
      ];
      k: k + m
    ];
    s: s + 1
  ];
  (re; im)
 };

/ Inverse Fast Fourier Transform with 1/N scaling
ifft: {[xr; xi]
  n: count xr;
  / Conjugate input, take forward FFT, conjugate and divide by n
  conj_im: neg xi;
  f_res: fft[xr; conj_im];
  (f_res[0] % (1f * n); neg f_res[1] % (1f * n))
 };

/ Fast convolution of two real sequences via FFT: h = f * g
convolve: {[f; g]
  n: count f;
  m: count g;
  total_len: n + m - 1;
  / Pad to next power of 2
  pow2: 1;
  while[pow2 < total_len; pow2: pow2 * 2];
  
  f_pad: f, (pow2 - n) # 0f;
  g_pad: g, (pow2 - m) # 0f;
  zeros: pow2 # 0f;
  
  F: fft[f_pad; zeros];
  G: fft[g_pad; zeros];
  
  / Complex pointwise multiplication: H = F * G
  H_re: (F[0] * G[0]) - (F[1] * G[1]);
  H_im: (F[0] * G[1]) + (F[1] * G[0]);
  
  res: ifft[H_re; H_im];
  total_len # res[0]
 };

/ Fractional FFT (Chirp-Z transform) evaluating Fourier sums at arbitrary frequency step alpha
/ X_k = sum_{j=0}^{N-1} x_j * exp(-2 * pi * i * j * k * alpha)
/ Allows decoupling frequency step size from grid point count
frfft: {[xr; xi; alpha]
  n: count xr;
  m: n;
  total_len: n + m - 1;
  pow2: 1;
  while[pow2 < total_len; pow2: pow2 * 2];
  
  / Chirp sequence y_k = exp(-i * pi * alpha * k^2)
  k_vec: til n;
  phase_y: neg .math.PI * alpha * k_vec * k_vec;
  y_re: cos[phase_y];
  y_im: sin[phase_y];
  
  / Modulate input: x_mod = x * y
  x_re: (xr * y_re) - (xi * y_im);
  x_im: (xr * y_im) + (xi * y_re);
  
  / Chirp kernel z_k = exp(i * pi * alpha * k^2)
  phase_z: .math.PI * alpha * k_vec * k_vec;
  z_re: cos[phase_z];
  z_im: sin[phase_z];
  
  / Construct padded vectors for circular convolution
  u_re: x_re, (pow2 - n) # 0f;
  u_im: x_im, (pow2 - n) # 0f;
  
  v_re: z_re, ((pow2 - (2 * n) + 1) # 0f), reverse 1 _ z_re;
  v_im: z_im, ((pow2 - (2 * n) + 1) # 0f), reverse 1 _ z_im;
  
  U: fft[u_re; u_im];
  V: fft[v_re; v_im];
  
  W_re: (U[0] * V[0]) - (U[1] * V[1]);
  W_im: (U[0] * V[1]) + (U[1] * V[0]);
  
  conv_res: ifft[W_re; W_im];
  
  / Final modulation by y_k
  c_re: n # conv_res[0];
  c_im: n # conv_res[1];
  
  out_re: (c_re * y_re) - (c_im * y_im);
  out_im: (c_re * y_im) + (c_im * y_re);
  (out_re; out_im)
 };

/ Hann smoothing window: w_k = 0.5 * (1 - cos(2 * pi * k / (N - 1)))
windowHann: {[n]
  0.5 * (1f - cos[.math.TWO_PI * (til n) % (1f * (n - 1))])
 };

/ Power Spectral Density estimation via Welch periodogram method
psd: {[x; window]
  n: count x;
  xw: x * window;
  f_res: fft[xw; n # 0f];
  power: (f_res[0] * f_res[0]) + (f_res[1] * f_res[1]);
  scale: sum window * window;
  power % scale
 };

\d .
