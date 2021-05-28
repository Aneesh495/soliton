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

\d .
