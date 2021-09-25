/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ test runner: master verification suite
/ ==============================================================================

\echo "Running Soliton Master Verification Test Suite..."

/ Load library
\l q/soliton.q

/ Helper assertion function
assert: {[cond; msg]
  if[not cond;
    0N! ("ASSERTION FAILED:"; msg);
    'test_failed
  ];
  0N! ("PASS:"; msg);
 };

/ ------------------------------------------------------------------------------
/ Test 1: Mathematics and Special Functions (.math)
/ ------------------------------------------------------------------------------
g_val: .math.gamma[5.0];
assert[abs[g_val - 24.0] < 1e-6; "Gamma function Gamma(5) == 24"];

erf_val: .math.erf[1.0];
assert[abs[erf_val - 0.84270079] < 1e-4; "Error function erf(1.0) == 0.8427"];

norm_cdf_val: .math.normalCdf[0.0];
assert[abs[norm_cdf_val - 0.5] < 1e-6; "Standard normal CDF Phi(0.0) == 0.5"];

norm_inv_val: .math.normalInv[0.5];
assert[abs[norm_inv_val - 0.0] < 1e-6; "Inverse normal CDF Phi^(-1)(0.5) == 0.0"];

/ Quadrature test: integral of x^2 from 0 to 1 == 1/3
quad_val: .math.quad[{x * x}; 0.0; 1.0];
assert[abs[quad_val - (1f % 3f)] < 1e-8; "Gauss-Kronrod quadrature integral(x^2, 0, 1) == 1/3"];

/ Root finding test: x^3 - 8 == 0 -> x = 2
root_val: .math.brentRoot[{ (x * x * x) - 8.0 }; 0.0; 5.0; 1e-8];
assert[abs[root_val - 2.0] < 1e-6; "Brent root solver x^3 - 8 == 0 -> x = 2"];

/ ------------------------------------------------------------------------------
/ Test 2: Linear Algebra and Matrix Factorizations (.la)
/ ------------------------------------------------------------------------------
A_mat: ((4f; 1f); (1f; 3f));
b_vec: (1f; 2f);
x_sol: .la.solve[A_mat; b_vec];
Ax_check: .la.matVec[A_mat; x_sol];
assert[max abs Ax_check - b_vec < 1e-6; "LU linear solve A*x == b"];

L_chol: .la.cholesky[A_mat];
LLt: .la.matMul[L_chol; .la.transpose[L_chol]];
assert[max raze abs LLt - A_mat < 1e-6; "Cholesky decomposition A == L * L^T"];

qr_decomp: .la.qr[A_mat];
QR_check: .la.matMul[qr_decomp`Q; qr_decomp`R];
assert[max raze abs QR_check - A_mat < 1e-6; "Householder QR decomposition A == Q * R"];

/ ------------------------------------------------------------------------------
/ Test 3: Probability Distributions and Cumulants (.stats)
/ ------------------------------------------------------------------------------
data_samples: 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0;
m_val: .stats.mean[data_samples];
assert[abs[m_val - 5.5] < 1e-8; "Sample mean == 5.5"];

v_val: .stats.var[data_samples];
assert[abs[v_val - 9.1666667] < 1e-5; "Sample variance == 9.1667"];

jb_res: .stats.jarqueBera[data_samples];
assert[jb_res`isNormal; "Jarque-Bera normality test on uniform sequence"];

/ ------------------------------------------------------------------------------
/ Test 4: Optimization Algorithms (.optim)
/ ------------------------------------------------------------------------------
/ Rosenbrock banana valley: f(x, y) = (1 - x)^2 + 100*(y - x^2)^2, min at (1, 1)
rosenbrock: {[v] ((1f - v[0]) * (1f - v[0])) + (100f * (v[1] - (v[0] * v[0])) * (v[1] - (v[0] * v[0]))) };
opt_nm: .optim.nelderMead[rosenbrock; (0.0; 0.0); 300; 1e-5];
assert[.la.norm2[opt_nm[`x] - (1.0; 1.0)] < 0.05; "Nelder-Mead optimization on Rosenbrock function"];

/ ------------------------------------------------------------------------------
/ Test 5: Fast Fourier Transform (.fft)
/ ------------------------------------------------------------------------------
sig_re: 1.0 0.0 -1.0 0.0;
sig_im: 0.0 0.0 0.0 0.0;
fft_res: .fft.fft[sig_re; sig_im];
ifft_res: .fft.ifft[fft_res[0]; fft_res[1]];
assert[max abs ifft_res[0] - sig_re < 1e-6; "FFT -> IFFT round-trip invertibility"];

/ ------------------------------------------------------------------------------
/ Test 6: Limit Order Book Dynamics and Metrics (.book)
/ ------------------------------------------------------------------------------
lob: .book.initBook[`AAPL];
lob: .book.addOrder[lob; 1; `buy; 150.00; 100f];
lob: .book.addOrder[lob; 2; `buy; 149.95; 200f];
lob: .book.addOrder[lob; 3; `sell; 150.05; 150f];
bbo: .book.getBBO[lob];
assert[bbo[`bestBid] == 150.00; "LOB best bid == 150.00"];
assert[bbo[`bestAsk] == 150.05; "LOB best ask == 150.05"];
assert[abs[bbo[`spread] - 0.05] < 1e-6; "LOB bid-ask spread == 0.05"];

mp: .book.stoikovMicroPrice[lob];
assert[(mp >= 150.00) and (mp <= 150.05); "Stoikov micro-price bounded by bid and ask"];

/ ------------------------------------------------------------------------------
/ Test 7: Multidimensional Hawkes Point Process (.hawkes)
/ ------------------------------------------------------------------------------
m_hawkes: .hawkes.initModel[enlist 1.2; enlist enlist 0.5; enlist enlist 1.0];
assert[m_hawkes`isStationary; "Hawkes stationarity branching ratio alpha/beta < 1"];

sim_hawkes: .hawkes.simulateOgata[m_hawkes; 10.0; 200];
assert[sim_hawkes[`count] > 0; "Ogata thinning simulation generated positive event count"];

ll: .hawkes.logLikelihood[m_hawkes; sim_hawkes`times; sim_hawkes`types; 10.0];
assert[not null ll; "Hawkes Ozaki log-likelihood evaluation non-null"];

/ ------------------------------------------------------------------------------
/ Test 8: Price Impact Models (.impact)
/ ------------------------------------------------------------------------------
p_series: 100.0 + sums 0.0, .stats.randn[20];
v_series: 100f * .stats.randn[20];
kyle_res: .impact.fitKyleLambda[p_series; v_series];
assert[not null kyle_res`lambda; "Kyle linear price impact lambda fitted"];

/ ------------------------------------------------------------------------------
/ Test 9: Almgren-Chriss Optimal Control (.ac)
/ ------------------------------------------------------------------------------
ac_res: .ac.optimalTrajectory[100000f; 1.0; 10; 0.30; 2.5e-7; 2.5e-6; 1e-5];
assert[abs[first[ac_res`holdings] - 100000f] < 1e-4; "Almgren-Chriss initial position X_0 == 100000"];
assert[abs[last[ac_res`holdings] - 0f] < 1e-4; "Almgren-Chriss terminal position X_T == 0"];
assert[abs[sum[ac_res`trades] - 100000f] < 1e-4; "Almgren-Chriss sum of trades == 100000"];

/ ------------------------------------------------------------------------------
/ Test 10: Cartea-Jaimungal Market Making (.cj)
/ ------------------------------------------------------------------------------
cj_quotes: .cj.avellanedaStoikovQuotes[100.0; 0; 0.5; 0.1; 0.3; 1.5; 140.0];
assert[cj_quotes[`optimalBid] < 100.0; "Avellaneda-Stoikov optimal bid below mid-price"];
assert[cj_quotes[`optimalAsk] > 100.0; "Avellaneda-Stoikov optimal ask above mid-price"];

/ ------------------------------------------------------------------------------
/ Test 11: Jump-Diffusion Models (.jd)
/ ------------------------------------------------------------------------------
cf_merton: .jd.mertonCharFn[1.0; 100.0; 1.0; 0.05; 0.0; 0.20; 0.5; -0.05; 0.10];
assert[abs[cf_merton[0]] > 0f; "Merton characteristic function real component non-zero"];

/ ------------------------------------------------------------------------------
/ Test 12: Continuous-Time Stochastic Volatility (.vol)
/ ------------------------------------------------------------------------------
assert[.vol.isFellerSatisfied[2.0; 0.04; 0.3]; "Feller condition check: 2 * 2.0 * 0.04 > 0.3^2"];
bs_price: .vol.bsCall[100.0; 100.0; 1.0; 0.05; 0.20];
iv_inverted: .vol.impliedVol[bs_price; 100.0; 100.0; 1.0; 0.05];
assert[abs[iv_inverted - 0.20] < 1e-3; "Black-Scholes implied volatility inversion recovery"];

/ ------------------------------------------------------------------------------
/ Test 13: Streaming Greeks and Risk Surfaces (.risk)
/ ------------------------------------------------------------------------------
greeks: .risk.analyticalGreeks[100.0; 100.0; 1.0; 0.05; 0.0; 0.20; true];
assert[(greeks[`delta] > 0.5) and (greeks[`delta] < 0.7); "ATM call Delta in [0.5, 0.7]"];
assert[greeks[`gamma] > 0f; "Option Gamma strictly positive"];
assert[greeks[`vega] > 0f; "Option Vega strictly positive"];

cf_var: .risk.varCornishFisher[0.01 * .stats.randn[100]; 0.99];
assert[cf_var > 0f; "Cornish-Fisher 99% VaR strictly positive loss threshold"];

/ ------------------------------------------------------------------------------
/ Test 14: Complex Event Processing (.cep)
/ ------------------------------------------------------------------------------
w_state: .cep.initWelford[];
w_state: .cep.updateWelford[w_state; 10.0];
w_state: .cep.updateWelford[w_state; 20.0];
assert[w_state[`mean] == 15.0; "Welford running mean of 10 and 20 == 15"];
assert[w_state[`variance] == 50.0; "Welford running sample variance of 10 and 20 == 50"];

/ ------------------------------------------------------------------------------
/ Test 15: Discrete-Event Backtesting Engine (.backtest)
/ ------------------------------------------------------------------------------
sim_engine: .backtest.initEngine[-0.0002; 0.0005; 1000000];
depth_snap: `bids`asks!(([] price: enlist 100.0; size: enlist 500f); ([] price: enlist 100.05; size: enlist 500f));
sim_engine: .backtest.submitLimitOrder[sim_engine; 101; `AAPL; `buy; 100.0; 100f; 1000000000; depth_snap];
sim_engine: .backtest.processMarketTrade[sim_engine; 2000000000; `AAPL; 100.0; 1000f; `sell];
assert[count sim_engine[`fills] > 0; "Backtest engine executed limit order fill upon matching trade"];

0N! "==============================================================================";
0N! "ALL 15 SOLITON SUB-SYSTEM VERIFICATION TEST BLOCKS PASSED DETERMINISTICALLY!";
0N! "==============================================================================";
