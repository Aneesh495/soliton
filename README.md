# soliton

> Non-equilibrium market microstructure and stochastic point-process laboratory in Q: multidimensional Hawkes self/cross-excitation, jump-diffusion calibration, Kyle-Obizhaeva transient impact, limit order book hydrodynamics, optimal stochastic control, and real-time streaming risk surfaces.

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)]()
[![Language: Q](https://img.shields.io/badge/language-Q%20%2F%20kdb+-blue.svg)]()
[![Engine: C11](https://img.shields.io/badge/engine-C11%20zero--dependency-orange.svg)]()
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)]()

---

## Overview

`soliton` is an autonomous, mathematically rigorous computational laboratory designed for continuous-time stochastic processes, high-frequency limit order book mechanics, and vectorized array analytics. Built from first principles, `soliton` brings together:

1. **A Standalone Native C Array Engine**: An ultra-low latency, zero-dependency C11 execution engine implementing the core `K` object model, reference-counted vector memory, higher-order adverbs (`each`, `over`, `scan`, `prior`, `each-right`, `each-left`), qSQL relational queries, and the complete kdb+ v3 binary IPC wire protocol. Any engineer can clone the repository and run the entire suite in seconds with zero proprietary licensing barriers.
2. **Multidimensional Point-Process Dynamics**: High-dimensional Hawkes processes with mutual cross-excitation, Ozaki linear-time $O(N)$ recursive likelihood evaluation, Ogata continuous-time thinning simulation, and Papangelou compensator residual analysis.
3. **Stochastic Optimal Control**: Almgren-Chriss risk-averse execution frontiers, nonlinear power-law temporary impact trajectories, multi-asset portfolio liquidation with correlated inventories, and Cartea-Jaimungal Hamilton-Jacobi-Bellman (HJB) optimal market making.
4. **Continuous-Time Volatility & Jump Mechanics**: Merton and Kou jump-diffusion processes, Bates stochastic volatility with jumps, Carr-Madan Fast Fourier Transform (FFT) option pricing, and rough Bergomi fractional Brownian motion simulation via Davies-Harte circulant embedding.
5. **Streaming Microstructure & Event Processing**: Depth reconstruction for limit order books (LOBs), Volume-Synchronized Probability of Toxicity (VPIN), multi-level Order Flow Imbalance (OFI), Stoikov micro-price estimation, CUSUM / Page-Hinkley change-point detection, and microsecond discrete-event backtesting with queue priority modeling.

---

## Architecture

```mermaid
flowchart TD
    subgraph Data Layer
        TP["Streaming Tickerplant (.tp)"] -->|IPC Protocol| RDB["In-Memory Real-Time Database (.rdb)"]
        TP -->|WAL Sync| WAL[("Write-Ahead Log (WAL)")]
        RDB -->|EOD Flush| HDB[("Partitioned Historical DB (.hdb)")]
    end

    subgraph Analytical Core
        LOB["Order Book Mechanics (.book)"] --> VPIN["Toxicity & Micro-Price"]
        HWK["Hawkes Point Processes (.hawkes)"] --> INT["Conditional Intensity & Compensators"]
        IMP["Impact Dynamics (.impact)"] --> SVAR["Hasbrouck SVAR & Gatheral Kernels"]
    end

    subgraph Stochastic Control & Pricing
        AC["Optimal Execution (.ac)"] --> TRAJ["Almgren-Chriss Trajectories"]
        CJ["Optimal Market Making (.cj)"] --> HJB["Cartea-Jaimungal Viscosity Solutions"]
        JD["Jump-Diffusion & Volatility (.jd, .vol)"] --> FFT_MOD["Carr-Madan FFT & rough Bergomi"]
    end

    subgraph Streaming Risk & Backtesting
        CEP["Complex Event Processing (.cep)"] --> CUSUM["Online Change-Point Detection"]
        RSK["Risk Surfaces (.risk)"] --> VAR_ES["Greeks, Cornish-Fisher VaR & ES"]
        BKT["Discrete-Event Backtester (.backtest)"] --> TCA["Implementation Shortfall & TCA"]
    end

    Data Layer --> Analytical Core
    Analytical Core --> Stochastic Control & Pricing
    Stochastic Control & Pricing --> Streaming Risk & Backtesting
```

---

## Component Breakdown

| Namespace | Module | Key Capabilities |
| :--- | :--- | :--- |
| `.math` | `q/core/math.q` | Lanczos log-gamma, error functions, incomplete gamma, Bessel functions, Gauss-Kronrod quadrature, Brent root solver |
| `.la` | `q/core/linalg.q` | Vector dot/outer products, LU with partial pivoting, Cholesky, Householder QR, symmetric eigensolver, SVD, Levinson-Durbin Toeplitz inversion |
| `.stats` | `q/core/stats.q` | Probability distributions (Normal, Student-t, Chi2, Gamma, Pareto, GED), Jarque-Bera test, copulas, Silverman/Scott KDE |
| `.optim` | `q/core/optim.q` | BFGS quasi-Newton optimizer with Armijo line search, Nelder-Mead simplex, Levenberg-Marquardt, projected gradient descent |
| `.fft` | `q/core/fft.q` | Radix-2 Cooley-Tukey FFT, inverse FFT, Fractional FFT (Chirp-Z transform), Welch power spectral density |
| `.book` | `q/microstructure/book.q` | L2/L3 order book matching, Stoikov micro-price, VPIN toxicity, multi-level OFI, Roll spread estimator |
| `.hawkes` | `q/microstructure/hawkes.q` | Multidimensional Hawkes process, Ogata thinning, Ozaki $O(N)$ log-likelihood, MLE, time-change residual test |
| `.impact` | `q/microstructure/impact.q` | Kyle lambda, Hasbrouck structural VAR, Bouchaud power-law propagator, Gatheral transient impact |
| `.ac` | `q/control/almgren_chriss.q` | Almgren-Chriss mean-variance execution, nonlinear impact trajectories, multi-asset portfolio liquidation |
| `.cj` | `q/control/cartea_jaimungal.q` | Avellaneda-Stoikov reservation prices, Cartea-Jaimungal HJB market making, continuous inventory control |
| `.jd` | `q/models/jump_diffusion.q` | Merton log-normal and Kou double-exponential jump-diffusion, Bates stochastic volatility, Carr-Madan FFT pricing |
| `.vol` | `q/models/volatility.q` | Heston characteristic function, implied volatility inversion, rough Bergomi fGn simulation, realized variance and TSRV |
| `.risk` | `q/risk/surface.q` | First/second-order Greeks (Delta, Gamma, Vega, Vanna, Volga), SVI variance parameterization, Cornish-Fisher VaR, Expected Shortfall |
| `.tp` / `.rdb` | `q/stream/tickerplant.q` | Low-latency pub-sub multicast, transaction logging, rolling intraday tables, end-of-day partitioning |
| `.cep` | `q/stream/cep.q` | Streaming Welford variance, Holt-Winters double exponential smoothing, tumbling/sliding windows, CUSUM and Page-Hinkley |
| `.backtest` | `q/backtest/engine.q` | Microsecond discrete-event simulation, transit latency delay, queue priority modeling, Implementation Shortfall TCA |

---

## Benchmarks

All performance benchmarks evaluated on Apple Silicon (arm64) using Clang `-O3`:

| Operation | Scale / Grid Size | Median Latency (p50) | Tail Latency (p99) | Throughput / Rate |
| :--- | :--- | :--- | :--- | :--- |
| **Vector Addition** | 1,000,000 Float64 | 13.38 ms | 16.43 ms | 0.56 GB/s |
| **Almgren-Chriss Solve** | $N = 20$ discrete steps | 2.50 us | 3.30 us | 406,836 solves/s |
| **Hawkes Ozaki Scan** | 2,000 tick events | 298.6 us | 455.2 us | 6,697,363 events/s |
| **Greeks Surface Evaluation** | 2,000 strike/maturity nodes | 880.6 us | 1225.4 us | 2,271,114 Greeks/s |

---

## Quickstart

### 1. Build the Native Engine
```bash
git clone https://github.com/Aneesh495/soliton.git
cd soliton
make all
```

### 2. Run the Verification Suite
```bash
make test
```

### 3. Run the End-to-End Demonstration
```bash
make demo
```

### 4. Interactive REPL Mode
```bash
./bin/soliton
q) x: 10 20 30 40 50
q) y: 1 2 3 4 5
q) x + y
11.000000 22.000000 33.000000 44.000000 55.000000
q) \
```

### 5. Launch the IPC Tickerplant Daemon
```bash
./bin/soliton -p 5001
```
In another terminal, stream synthetic market events into the running tickerplant:
```bash
python3 tools/synthetic_feed.py --stream
```

---

## Documentation

* [Mathematical Specification](docs/MATHEMATICAL_SPEC.md): Formal derivations, stochastic differential equations, characteristic function proofs, and control formulations.
* [System Architecture](docs/ARCHITECTURE.md): Deep dive into memory allocators, AST evaluation, right-to-left reduction, and IPC wire protocol frames.

---

## License

MIT License. Copyright (c) 2026 Aneesh Krishna.
