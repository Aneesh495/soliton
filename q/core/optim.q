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

/ Nelder-Mead downhill simplex algorithm for derivative-free optimization
nelderMead: {[f; x0; max_iter; tol]
  n: count x0;
  alpha: 1.0;
  gamma: 2.0;
  rho: 0.5;
  sigma: 0.5;
  
  / Initialize simplex vertices with n + 1 points
  simplex: enlist x0;
  i: 0;
  while[i < n;
    vertex: x0;
    step: $[abs[x0[i]] > 0f; 0.05 * x0[i]; 0.00025];
    vertex[i]: vertex[i] + step;
    simplex: simplex, enlist vertex;
    i: i + 1
  ];
  
  fvals: f each simplex;
  iter: 0;
  max_it: $[max_iter < 10; 200; max_iter];
  
  while[iter < max_it;
    / Sort vertices by function value
    order: iasc fvals;
    simplex: simplex[order];
    fvals: fvals[order];
    
    / Check termination tolerance
    diff: max abs fvals - first fvals;
    if[diff < tol; iter: 9999];
    
    / Centroid of best n vertices
    best_n: n # simplex;
    centroid: (sum best_n) % (1f * n);
    
    / Reflection step
    x_worst: last simplex;
    f_worst: last fvals;
    x_refl: centroid + (alpha * (centroid - x_worst));
    f_refl: f[x_refl];
    
    if[(f_refl >= fvals[0]) and (f_refl < fvals[n - 1]);
      simplex[n]: x_refl;
      fvals[n]: f_refl
    ;
      if[f_refl < fvals[0];
        / Expansion step
        x_exp: centroid + (gamma * (x_refl - centroid));
        f_exp: f[x_exp];
        if[f_exp < f_refl;
          simplex[n]: x_exp; fvals[n]: f_exp
        ;
          simplex[n]: x_refl; fvals[n]: f_refl
        ]
      ;
        / Contraction step
        if[f_refl < f_worst;
          / Outside contraction
          x_cont: centroid + (rho * (x_refl - centroid));
          f_cont: f[x_cont];
          if[f_cont <= f_refl;
            simplex[n]: x_cont; fvals[n]: f_cont
          ;
            / Shrink simplex
            k: 1;
            while[k <= n;
              simplex[k]: simplex[0] + (sigma * (simplex[k] - simplex[0]));
              fvals[k]: f[simplex[k]];
              k: k + 1
            ]
          ]
        ;
          / Inside contraction
          x_cont: centroid - (rho * (centroid - x_worst));
          f_cont: f[x_cont];
          if[f_cont < f_worst;
            simplex[n]: x_cont; fvals[n]: f_cont
          ;
            / Shrink simplex
            k: 1;
            while[k <= n;
              simplex[k]: simplex[0] + (sigma * (simplex[k] - simplex[0]));
              fvals[k]: f[simplex[k]];
              k: k + 1
            ]
          ]
        ]
      ]
    ];
    iter: iter + 1
  ];
  `x`fval`iterations!(first simplex; first fvals; iter)
 };

/ Numerical Jacobian of vector residual function r(x)
jacobian: {[res_fn; x]
  r0: res_fn[x];
  m: count r0;
  n: count x;
  h_base: 1.4901161193847656e-8;
  J_cols: ();
  
  j: 0;
  while[j < n;
    hj: h_base * max (1f; abs[x[j]]);
    xp: x; xm: x;
    xp[j]: xp[j] + hj;
    xm[j]: xm[j] - hj;
    col: (res_fn[xp] - res_fn[xm]) % (2f * hj);
    J_cols: J_cols, enlist col;
    j: j + 1
  ];
  flip J_cols
 };

/ Levenberg-Marquardt nonlinear least squares algorithm
/ Solves min 0.5 * ||r(x)||^2 via adaptive damping (J^T * J + lambda * I) * delta = -J^T * r
levenbergMarquardt: {[res_fn; x0; max_iter; tol]
  x: x0;
  r: res_fn[x];
  cost: 0.5 * sum r * r;
  lambda: 1e-3;
  
  iter: 0;
  max_it: $[max_iter < 10; 100; max_iter];
  converged: false;
  
  while[(iter < max_it) and not converged;
    J: jacobian[res_fn; x];
    Jt: flip J;
    JtJ: .la.matMul[Jt; J];
    grad_cost: .la.matVec[Jt; r];
    
    if[.la.norm2[grad_cost] < tol;
      converged: true
    ;
      n: count x;
      / Damped normal equations: (J^T * J + lambda * diag(J^T * J)) * delta = -grad
      A: JtJ;
      d: .la.diag[JtJ];
      i: 0;
      while[i < n;
        A[i; i]: A[i; i] + (lambda * max (1e-6; d[i]));
        i: i + 1
      ];
      
      delta: .la.solve[A; neg grad_cost];
      x_trial: x + delta;
      r_trial: res_fn[x_trial];
      cost_trial: 0.5 * sum r_trial * r_trial;
      
      if[cost_trial < cost;
        / Successful step: decrease damping
        x: x_trial;
        r: r_trial;
        cost: cost_trial;
        lambda: max (1e-7; lambda * 0.1);
        if[abs[cost - cost_trial] < (tol * cost); converged: true]
      ;
        / Unsuccessful step: increase damping
        lambda: min (1e7; lambda * 10f)
      ];
      iter: iter + 1
    ]
  ];
  `x`cost`iterations`converged!(x; cost; iter; converged)
 };

/ Projected gradient descent under box constraints lower <= x <= upper
projectedGradientDescent: {[f; x0; lower; upper; max_iter; tol]
  x: min'[upper; max'[lower; x0]];
  fx: f[x];
  alpha: 0.1;
  iter: 0;
  max_it: $[max_iter < 10; 100; max_iter];
  
  while[iter < max_it;
    g: grad[f; x];
    x_trial: min'[upper; max'[lower; x - (alpha * g)]];
    fx_trial: f[x_trial];
    
    if[fx_trial < fx;
      if[abs[fx - fx_trial] < tol; iter: 9999];
      x: x_trial;
      fx: fx_trial;
      alpha: min (1.0; alpha * 1.1)
    ;
      alpha: alpha * 0.5
    ];
    iter: iter + 1
  ];
  `x`fval`iterations!(x; fx; iter)
 };

\d .
