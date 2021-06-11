/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/core/optim.q
/ description: vectorized unconstrained and constrained nonlinear optimization
/ ==============================================================================

\d .optim

/ Numerical gradient vector grad(f)(x) via central finite differences
grad: {[f; x]
  n: count x;
  h_base: 1.4901161193847656e-8;
  g: n # 0f;
  i: 0;
  while[i < n;
    h: h_base * max (1f; abs[x[i]]);
    xp: x; xm: x;
    xp[i]: xp[i] + h;
    xm[i]: xm[i] - h;
    g[i]: (f[xp] - f[xm]) % (2f * h);
    i: i + 1
  ];
  g
 };

/ Numerical Hessian matrix H(f)(x) via forward-backward second finite differences
hessian: {[f; x]
  n: count x;
  h_base: 1.0e-5;
  H: {n # 0f} each til n;
  f0: f[x];
  
  i: 0;
  while[i < n;
    hi: h_base * max (1f; abs[x[i]]);
    x_pi: x; x_pi[i]: x_pi[i] + hi;
    x_mi: x; x_mi[i]: x_mi[i] - hi;
    f_pi: f[x_pi];
    f_mi: f[x_mi];
    H[i; i]: (f_pi - (2f * f0) + f_mi) % (hi * hi);
    
    j: i + 1;
    while[j < n;
      hj: h_base * max (1f; abs[x[j]]);
      x_pij: x_pi; x_pij[j]: x_pij[j] + hj;
      f_pij: f[x_pij];
      
      x_mij: x_mi; x_mij[j]: x_mij[j] - hj;
      f_mij: f[x_mij];
      
      d2: (f_pij - f_pi - f[x] + f_mi) % (hi * hj);
      H[i; j]: d2;
      H[j; i]: d2;
      j: j + 1
    ];
    i: i + 1
  ];
  H
 };

/ Backtracking line search satisfying the Armijo condition
/ f(x + alpha * p) <= f(x) + c1 * alpha * grad(f)^T * p
lineSearchArmijo: {[f; x; p; g0; f0]
  c1: 1e-4;
  rho: 0.5;
  alpha: 1f;
  slope: .la.dot[g0; p];
  
  / Ensure descent direction
  if[slope >= 0f; :0f];
  
  iter: 0;
  while[iter < 30;
    x_trial: x + (alpha * p);
    f_trial: f[x_trial];
    if[f_trial <= (f0 + (c1 * alpha * slope));
      :alpha
    ];
    alpha: alpha * rho;
    iter: iter + 1
  ];
  alpha
 };

/ BFGS quasi-Newton algorithm for unconstrained nonlinear minimization
/ f: objective function mapping vector -> scalar
/ x0: initial parameter guess
/ max_iter: iteration cap
/ tol: gradient tolerance stopping threshold
bfgs: {[f; x0; max_iter; tol]
  n: count x0;
  x: x0;
  H_inv: .la.eye[n];
  g: grad[f; x];
  fx: f[x];
  
  max_it: $[max_iter < 10; 100; max_iter];
  tolerance: $[tol <= 0f; 1e-6; tol];
  
  iter: 0;
  converged: false;
  
  while[(iter < max_it) and not converged;
    gnorm: .la.norm2[g];
    if[gnorm < tolerance;
      converged: true
    ;
      / Search direction: p = -H_inv * g
      p: neg .la.matVec[H_inv; g];
      
      / Perform Armijo line search
      alpha: lineSearchArmijo[f; x; p; g; fx];
      if[alpha <= 1e-12;
        / Line search failure: reset H_inv to identity
        H_inv: .la.eye[n];
        p: neg g;
        alpha: lineSearchArmijo[f; x; p; g; fx];
        if[alpha <= 1e-12; converged: true]
      ];
      
      if[not converged;
        s: alpha * p;
        x_next: x + s;
        fx_next: f[x_next];
        g_next: grad[f; x_next];
        y: g_next - g;
        
        ys: .la.dot[y; s];
        if[ys > 1e-12;
          / BFGS rank-2 update:
          / H_inv = (I - rho * s * y^T) * H_inv * (I - rho * y * s^T) + rho * s * s^T
          rho_val: 1f % ys;
          I: .la.eye[n];
          V1: I - (rho_val * .la.outer[s; y]);
          V2: I - (rho_val * .la.outer[y; s]);
          H_inv: .la.matMul[.la.matMul[V1; H_inv]; V2] + (rho_val * .la.outer[s; s])
        ];
        
        x: x_next;
        fx: fx_next;
        g: g_next;
        iter: iter + 1
      ]
    ]
  ];
  `x`fval`grad`iterations`converged!(x; fx; g; iter; converged)
 };

\d .
