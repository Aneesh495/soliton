/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ module: q/stream/cep.q
/ description: complex event processing (CEP), sliding windows, and change-point tests
/ ==============================================================================

\d .cep

/ Initialize online Welford running mean and variance accumulator
initWelford: {
  `count`mean`M2!(0; 0f; 0f)
 };

/ Update Welford accumulator with new stream scalar sample x
updateWelford: {[state; x]
  c: state[`count] + 1;
  delta: x - state`mean;
  new_mean: state[`mean] + (delta % (1f * c));
  delta2: x - new_mean;
  new_M2: state[`M2] + (delta * delta2);
  var_val: $[c > 1; new_M2 % (1f * (c - 1)); 0f];
  `count`mean`M2`variance`stdDev!(c; new_mean; new_M2; var_val; sqrt[var_val])
 };

/ Exponential Moving Average (EMA) update: S_t = alpha * x_t + (1 - alpha) * S_{t-1}
updateEMA: {[prev_ema; x; alpha]
  if[null prev_ema; :x];
  (alpha * x) + ((1f - alpha) * prev_ema)
 };

/ Holt-Winters double exponential smoothing with trend: Level l_t and Trend b_t
initHoltWinters: {[x0; trend0; alpha; beta]
  `level`trend`alpha`beta!(x0; trend0; alpha; beta)
 };

updateHoltWinters: {[state; x]
  l_prev: state`level;
  b_prev: state`trend;
  a: state`alpha;
  b: state`beta;
  
  new_l: (a * x) + ((1f - a) * (l_prev + b_prev));
  new_b: (b * (new_l - l_prev)) + ((1f - b) * b_prev);
  forecast_1step: new_l + new_b;
  `level`trend`alpha`beta`forecast!(new_l; new_b; a; b; forecast_1step)
 };

/ Tumbling temporal window aggregation: buckets time-series into non-overlapping bins of length dt
tumblingWindow: {[tbl; time_col; dt_nanos; agg_dict]
  / Floor timestamps to interval start: dt * floor(time / dt)
  t_vec: tbl[time_col];
  bins: dt_nanos * `long$floor t_vec % dt_nanos;
  tbl_binned: update windowBin: bins from tbl;
  
  / Group by windowBin and evaluate aggregations
  grouped: ?[tbl_binned; (); (enlist `windowBin)!enlist `windowBin; agg_dict];
  `windowBin xasc grouped
 };

/ Sliding temporal window: extracts events within [t - window_nanos, t] for each event
slidingWindowSlice: {[tbl; time_col; curr_time; window_nanos]
  start_time: curr_time - window_nanos;
  ?[tbl; enlist (within; time_col; (start_time; curr_time)); 0b; ()]
 };

/ Cumulative Sum (CUSUM) change-point detection for shift in mean of process
/ Detects sudden regime change in expected return or order flow
initCUSUM: {[target_mean; drift_slack; threshold_h]
  `target`slack`threshold`sPlus`sMinus`alarm!(target_mean; drift_slack; threshold_h; 0f; 0f; false)
 };

updateCUSUM: {[state; x]
  target: state`target;
  slack: state`slack;
  h: state`threshold;
  
  new_sp: max (0f; state[`sPlus] + x - target - slack);
  new_sm: max (0f; state[`sMinus] - (x - target) - slack);
  
  alarm_triggered: (new_sp > h) or (new_sm > h);
  `target`slack`threshold`sPlus`sMinus`alarm!(target; slack; h; new_sp; new_sm; alarm_triggered)
 };

/ Page-Hinkley statistical test for detecting changes in variance or volatility
initPageHinkley: {[delta_tolerance; threshold_h]
  `delta`threshold`cumSum`minCumSum`maxCumSum`alarm!(delta_tolerance; threshold_h; 0f; 0f; 0f; false)
 };

updatePageHinkley: {[state; x; mean_ref]
  delta: state`delta;
  h: state`threshold;
  
  new_cumsum: state[`cumSum] + x - mean_ref - delta;
  new_min: min (state`minCumSum; new_cumsum);
  
  ph_stat: new_cumsum - new_min;
  alarm_triggered: ph_stat > h;
  `delta`threshold`cumSum`minCumSum`alarm!(delta; h; new_cumsum; new_min; alarm_triggered)
 };

\d .
