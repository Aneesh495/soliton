/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/core/linalg.q
/ description: vectorized linear algebra and matrix factorizations in pure Q
/ ==============================================================================

\d .la

/ Inner product of two vectors
dot: {[x; y]
  sum x * y
 };

/ Matrix-vector multiplication A * x
matVec: {[A; x]
  dot[; x] each A
 };

/ Matrix multiplication A * B
matMul: {[A; B]
  Bt: flip B;
  {dot[x;] each y}[; Bt] each A
 };

/ Transpose of a matrix
transpose: {[A]
  flip A
 };

/ Identity matrix of dimension n
eye: {[n]
  {x = y}'[til n;] each til n
 };

/ Vector outer product x * y^T
outer: {[x; y]
  x */: y
 };

/ Matrix trace: sum of diagonal elements
trace: {[A]
  n: min (count A; count first A);
  sum {A[x; x]} each til n
 };

/ Matrix diagonal extraction
diag: {[A]
  n: min (count A; count first A);
  {A[x; x]} each til n
 };

/ Construct diagonal matrix from vector
diagMat: {[v]
  n: count v;
  {?[x = y; v[x]; 0f]}'[til n;] each til n
 };

/ Frobenius norm of a matrix
normF: {[A]
  sqrt sum raze A * A
 };

/ Euclidean norm of a vector
norm2: {[x]
  sqrt sum x * x
 };

/ LU decomposition with partial row pivoting: P * A = L * U
/ Returns dictionary: `L`U`p`sign
lu: {[A]
  n: count A;
  U: A;
  L: eye[n];
  p: til n;
  sign: 1f;
  
  k: 0;
  while[k < n;
    / Find pivot element in column k
    col_vals: abs {U[x; k]} each k + til n - k;
    max_idx: k + first where col_vals = max col_vals;
    
    if[max_idx <> k;
      / Swap rows in U
      tmp: U[k]; U[k]: U[max_idx]; U[max_idx]: tmp;
      / Swap rows in P
      ptmp: p[k]; p[k]: p[max_idx]; p[max_idx]: ptmp;
      / Swap elements in L
      if[k > 0;
        ltmp: L[k; til k]; L[k; til k]: L[max_idx; til k]; L[max_idx; til k]: ltmp
      ];
      sign: neg sign
    ];
    
    pivot: U[k; k];
    if[abs[pivot] > 1e-15;
      i: k + 1;
      while[i < n;
        factor: U[i; k] % pivot;
        L[i; k]: factor;
        U[i; k + til n - k]: U[i; k + til n - k] - (factor * U[k; k + til n - k]);
        i: i + 1
      ]
    ];
    k: k + 1
  ];
  `L`U`p`sign!(L; U; p; sign)
 };

/ Forward substitution solving L * y = b for unit lower triangular matrix L
forwardSub: {[L; b]
  n: count b;
  y: n # 0f;
  i: 0;
  while[i < n;
    sum_prev: sum L[i; til i] * y[til i];
    y[i]: b[i] - sum_prev;
    i: i + 1
  ];
  y
 };

/ Back substitution solving U * x = y for upper triangular matrix U
backSub: {[U; y]
  n: count y;
  x: n # 0f;
  i: n - 1;
  while[i >= 0;
    sum_next: sum U[i; (i + 1) + til n - 1 - i] * x[(i + 1) + til n - 1 - i];
    denom: U[i; i];
    x[i]: $[abs[denom] > 1e-15; (y[i] - sum_next) % denom; 0f];
    i: i - 1
  ];
  x
 };

/ Solve linear system A * x = b via LU decomposition with partial pivoting
solve: {[A; b]
  decomp: lu[A];
  bp: b[decomp`p];
  y: forwardSub[decomp`L; bp];
  backSub[decomp`U; y]
 };

/ Determinant of a square matrix det(A)
det: {[A]
  decomp: lu[A];
  decomp[`sign] * prd diag[decomp`U]
 };

/ Matrix inversion A^(-1) via LU column-by-column solution
inv: {[A]
  n: count A;
  I: eye[n];
  flip solve[A;] each I
 };

/ Cholesky decomposition of symmetric positive-definite matrix A: A = L * L^T
/ Returns lower triangular matrix L
cholesky: {[A]
  n: count A;
  L: {n # 0f} each til n;
  
  i: 0;
  while[i < n;
    j: 0;
    while[j <= i;
      s: sum L[i; til j] * L[j; til j];
      if[i == j;
        val: A[i; i] - s;
        if[val <= 0f; 'not_pos_def];
        L[i; j]: sqrt[val]
      ;
        denom: L[j; j];
        L[i; j]: $[denom > 0f; (A[i; j] - s) % denom; 0f]
      ];
      j: j + 1
    ];
    i: i + 1
  ];
  L
 };

/ Solve symmetric positive-definite system A * x = b via Cholesky factor L
solveChol: {[L; b]
  y: forwardSub[L; b];
  backSub[flip L; y]
 };

/ Log-determinant of symmetric positive-definite matrix: ln(det(A)) = 2 * sum(ln(L_ii))
logDetChol: {[L]
  2f * sum log diag[L]
 };

/ QR decomposition via Householder reflectors: A = Q * R
/ Returns dictionary: `Q`R
qr: {[A]
  m: count A;
  n: count first A;
  R: A;
  Q: eye[m];
  
  k: 0;
  while[k < min (m - 1; n);
    x: {R[x; k]} each k + til m - k;
    norm_x: norm2[x];
    alpha: $[x[0] >= 0f; neg norm_x; norm_x];
    
    u: x;
    u[0]: u[0] - alpha;
    norm_u: norm2[u];
    
    if[norm_u > 1e-15;
      v: u % norm_u;
      
      / Update submatrix of R: R = R - 2 * v * (v^T * R)
      sub_R: {R[x; k + til n - k]} each k + til m - k;
      vt_R: dot[v;] each flip sub_R;
      v_vt_R: outer[v; vt_R];
      
      r_idx: 0;
      while[r_idx < count v;
        R[k + r_idx; k + til n - k]: sub_R[r_idx] - (2f * v_vt_R[r_idx]);
        r_idx: r_idx + 1
      ];
      
      / Accumulate orthogonal matrix Q: Q = Q - 2 * (Q * v) * v^T
      sub_Q: {Q[x; k + til m - k]} each til m;
      Q_v: dot[; v] each sub_Q;
      
      q_idx: 0;
      while[q_idx < m;
        Q[q_idx; k + til m - k]: sub_Q[q_idx] - (2f * Q_v[q_idx] * v);
        q_idx: q_idx + 1
      ]
    ];
    k: k + 1
  ];
  `Q`R!(Q; R)
 };

/ Symmetric eigenvalue problem via QR algorithm with Rayleigh quotient shifts
/ Returns dictionary: `values`vectors
eigSym: {[A; max_iter]
  n: count A;
  V: eye[n];
  T: A;
  iter: 0;
  max_it: $[max_iter < 10; 80; max_iter];
  
  while[iter < max_it;
    / Rayleigh shift from bottom 2x2 corner
    d1: T[n - 1; n - 1];
    d0: T[n - 2; n - 2];
    b: T[n - 1; n - 2];
    tr: d0 + d1;
    dt: (d0 * d1) - (b * b);
    disc: sqrt[max (0f; (0.25 * tr * tr) - dt)];
    lambda1: (0.5 * tr) + disc;
    lambda2: (0.5 * tr) - disc;
    shift: $[abs[lambda1 - d1] < abs[lambda2 - d1]; lambda1; lambda2];
    
    / Shifted QR decomposition: T - shift * I = Q * R
    T_shifted: T;
    i: 0;
    while[i < n;
      T_shifted[i; i]: T_shifted[i; i] - shift;
      i: i + 1
    ];
    qr_res: qr[T_shifted];
    Q: qr_res`Q;
    R: qr_res`R;
    
    / Reconstruct shifted matrix: T = R * Q + shift * I
    T_next: matMul[R; Q];
    i: 0;
    while[i < n;
      T_next[i; i]: T_next[i; i] + shift;
      i: i + 1
    ];
    T: T_next;
    V: matMul[V; Q];
    
    / Check off-diagonal convergence
    off_diag: 0f;
    r: 0;
    while[r < n;
      c: 0;
      while[c < r;
        off_diag: off_diag + abs[T[r; c]];
        c: c + 1
      ];
      r: r + 1
    ];
    if[off_diag < 1e-12; iter: 9999];
    iter: iter + 1
  ];
  
  vals: diag[T];
  / Sort eigenvalues and eigenvectors in descending order
  order: idesc vals;
  `values`vectors!(vals[order]; flip {flip[V][x]} each order)
 };

/ Singular Value Decomposition: A = U * S * V^T
/ Computed via symmetric eigen-decomposition of A^T * A and A * A^T
/ Returns dictionary: `U`S`V
svd: {[A]
  m: count A;
  n: count first A;
  At: flip A;
  AtA: matMul[At; A];
  
  eig_v: eigSym[AtA; 100];
  v_vals: eig_v`values;
  V: eig_v`vectors;
  
  / Singular values are square roots of eigenvalues of A^T * A
  s_vals: sqrt[max'[v_vals; 0f]];
  
  / Compute left singular vectors U_k = (1 / s_k) * A * V_k
  U_cols: ();
  k: 0;
  while[k < count s_vals;
    s: s_vals[k];
    v_col: {V[x; k]} each til n;
    u_col: $[s > 1e-14; matVec[A; v_col] % s; m # 0f];
    U_cols: U_cols, enlist u_col;
    k: k + 1
  ];
  U: flip U_cols;
  `U`S`V!(U; s_vals; V)
 };

/ Moore-Penrose pseudoinverse via SVD with thresholding
pinv: {[A]
  s_res: svd[A];
  U: s_res`U;
  S: s_res`S;
  V: s_res`V;
  
  max_s: max S;
  tol: 1e-14 * max_s;
  
  S_inv: {?[x > tol; 1f % x; 0f]} each S;
  diag_S_inv: diagMat[S_inv];
  
  / A^+ = V * S^+ * U^T
  matMul[matMul[V; diag_S_inv]; flip U]
 };

/ Levinson-Durbin algorithm for symmetric Toeplitz system T * x = y
/ Time complexity: O(n^2) instead of O(n^3)
/ r: autocovariance vector [r0, r1, ..., r_{n-1}] defining Toeplitz rows
/ y: right hand side vector
levinsonDurbin: {[r; y]
  n: count y;
  if[n == 1; :enlist (y[0] % r[0])];
  
  a: enlist (neg r[1] % r[0]);
  e: r[0] * (1f - (a[0] * a[0]));
  
  f: enlist (y[0] % r[0]);
  
  k: 1;
  while[k < n;
    / Compute reflection coefficient
    q: r[k + 1] + sum r[1 + til k] * reverse a;
    refl: neg q % e;
    
    a_new: a + (refl * reverse a);
    a: a_new, enlist refl;
    e: e * (1f - (refl * refl));
    
    / Update solution vector f
    gamma: (y[k] - sum r[1 + til k] * reverse f) % e;
    f_new: f + (gamma * reverse a[til k]);
    f: f_new, enlist gamma;
    
    k: k + 1
  ];
  f
 };

\d .
