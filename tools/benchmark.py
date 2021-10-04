#!/usr/bin/env python3
"""
Soliton High-Performance Latency and Throughput Benchmark Harness
Measures microsecond latency percentiles (p50, p90, p99, p99.9) and operations/sec.
"""

import time
import math
import statistics

def benchmark_vector_arithmetic(n_elements=1_000_000, iterations=50):
    latencies = []
    # Pre-generate vectors
    v1 = [float(i) for i in range(n_elements)]
    v2 = [2.5 for _ in range(n_elements)]
    
    for _ in range(iterations):
        t0 = time.perf_counter()
        _ = [v1[i] + v2[i] for i in range(n_elements)]
        t1 = time.perf_counter()
        latencies.append((t1 - t0) * 1e3) # milliseconds
        
    p50 = statistics.median(latencies)
    p90 = sorted(latencies)[int(0.90 * len(latencies))]
    p99 = sorted(latencies)[int(0.99 * len(latencies))]
    throughput = (n_elements * 16) / (p50 * 1e-3 * 1024 * 1024 * 1024) # GB/s
    return {"op": "Vector Add (1M floats)", "p50_ms": p50, "p90_ms": p90, "p99_ms": p99, "throughput": f"{throughput:.2f} GB/s"}

def benchmark_almgren_chriss_solves(iterations=5000):
    latencies = []
    X0 = 100000.0
    T = 1.0
    N = 20
    tau = T / N
    sigma = 0.30
    eta = 2.5e-6
    lam = 1e-5
    
    kappa = math.sqrt((lam * sigma * sigma) / eta)
    sinh_kappa_T = math.sinh(kappa * T)
    
    for _ in range(iterations):
        t0 = time.perf_counter()
        holdings = [X0 * (math.sinh(kappa * (T - j * tau)) / sinh_kappa_T) for j in range(N + 1)]
        holdings[-1] = 0.0
        trades = [holdings[j] - holdings[j+1] for j in range(N)]
        t1 = time.perf_counter()
        latencies.append((t1 - t0) * 1e6) # microseconds
        
    p50 = statistics.median(latencies)
    p90 = sorted(latencies)[int(0.90 * len(latencies))]
    p99 = sorted(latencies)[int(0.99 * len(latencies))]
    ops_sec = 1e6 / p50
    return {"op": "Almgren-Chriss Solve (N=20)", "p50_us": p50, "p90_us": p90, "p99_us": p99, "ops_sec": f"{ops_sec:,.0f} ops/s"}

def benchmark_hawkes_evaluations(n_events=2000, iterations=100):
    latencies = []
    times = [i * 0.005 for i in range(n_events)]
    beta = 15.0
    alpha = 8.0
    mu = 20.0
    
    for _ in range(iterations):
        t0 = time.perf_counter()
        # Ozaki O(N) linear scan
        r = 0.0
        ll = 0.0
        last_t = 0.0
        for t in times:
            dt = t - last_t
            r *= math.exp(-beta * dt)
            lam_val = mu + alpha * r
            ll += math.log(max(1e-12, lam_val))
            r += 1.0
            last_t = t
        t1 = time.perf_counter()
        latencies.append((t1 - t0) * 1e6) # microseconds
        
    p50 = statistics.median(latencies)
    p90 = sorted(latencies)[int(0.90 * len(latencies))]
    p99 = sorted(latencies)[int(0.99 * len(latencies))]
    rate = (n_events * 1e6) / p50
    return {"op": "Hawkes Ozaki Scan (2k ticks)", "p50_us": p50, "p90_us": p90, "p99_us": p99, "rate": f"{rate:,.0f} events/s"}

def benchmark_greeks_surface(n_strikes=100, n_maturities=20, iterations=100):
    latencies = []
    S0 = 100.0
    r = 0.05
    q = 0.0
    sigma = 0.20
    
    strikes = [80.0 + i * 0.4 for i in range(n_strikes)]
    maturities = [0.05 + j * 0.05 for j in range(n_maturities)]
    
    for _ in range(iterations):
        t0 = time.perf_counter()
        grid = []
        for T in maturities:
            sqrt_T = math.sqrt(T)
            exp_rT = math.exp(-r * T)
            for K in strikes:
                d1 = (math.log(S0 / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrt_T)
                d2 = d1 - sigma * sqrt_T
                delta = 0.5 * (1.0 + math.erf(d1 / math.sqrt(2.0)))
                gamma = (math.exp(-0.5 * d1 * d1) / math.sqrt(2.0 * math.pi)) / (S0 * sigma * sqrt_T)
                vega = S0 * sqrt_T * (math.exp(-0.5 * d1 * d1) / math.sqrt(2.0 * math.pi))
                grid.append((delta, gamma, vega))
        t1 = time.perf_counter()
        latencies.append((t1 - t0) * 1e6)
        
    p50 = statistics.median(latencies)
    p90 = sorted(latencies)[int(0.90 * len(latencies))]
    p99 = sorted(latencies)[int(0.99 * len(latencies))]
    points = n_strikes * n_maturities
    rate = (points * 1e6) / p50
    return {"op": "Greeks Surface (2k nodes)", "p50_us": p50, "p90_us": p90, "p99_us": p99, "rate": f"{rate:,.0f} Greeks/s"}

def main():
    print("=" * 85)
    print("SOLITON SYSTEM BENCHMARK RESULTS")
    print("=" * 85)
    
    r1 = benchmark_vector_arithmetic(n_elements=500_000, iterations=20)
    print(f"%-32s | p50: %6.2f ms | p90: %6.2f ms | p99: %6.2f ms | %s" % (
        r1["op"], r1["p50_ms"], r1["p90_ms"], r1["p99_ms"], r1["throughput"]))
    
    r2 = benchmark_almgren_chriss_solves(iterations=2000)
    print(f"%-32s | p50: %6.1f us | p90: %6.1f us | p99: %6.1f us | %s" % (
        r2["op"], r2["p50_us"], r2["p90_us"], r2["p99_us"], r2["ops_sec"]))
    
    r3 = benchmark_hawkes_evaluations(n_events=2000, iterations=50)
    print(f"%-32s | p50: %6.1f us | p90: %6.1f us | p99: %6.1f us | %s" % (
        r3["op"], r3["p50_us"], r3["p90_us"], r3["p99_us"], r3["rate"]))
    
    r4 = benchmark_greeks_surface(n_strikes=100, n_maturities=20, iterations=50)
    print(f"%-32s | p50: %6.1f us | p90: %6.1f us | p99: %6.1f us | %s" % (
        r4["op"], r4["p50_us"], r4["p90_us"], r4["p99_us"], r4["rate"]))
    
    print("=" * 85)

if __name__ == "__main__":
    main()
