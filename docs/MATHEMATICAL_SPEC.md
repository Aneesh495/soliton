# Soliton: Mathematical and Algorithmic Specification

## 1. Multidimensional Mutually-Exciting Hawkes Processes

A multivariate point process $N(t) = (N_1(t), \dots, N_M(t))^T$ is defined by its stochastic conditional intensity vector $\lambda(t) = (\lambda_1(t), \dots, \lambda_M(t))^T$:

$$\lambda_m(t) = \mu_m + \sum_{j=1}^M \int_0^t \alpha_{mj} e^{-\beta_{mj}(t - s)} dN_j(s)$$

where:
* $\mu_m > 0$ denotes baseline exogenous event intensity.
* $\alpha_{mj} \ge 0$ measures the excitation impact of an event on channel $j$ upon channel $m$.
* $\beta_{mj} > 0$ defines the exponential decay rate of memory.

### 1.1 Branching Ratio Matrix and Subcritical Stability
The overall branching matrix $\Gamma \in \mathbb{R}^{M \times M}$ has entries:

$$\Gamma_{mj} = \int_0^\infty \alpha_{mj} e^{-\beta_{mj} s} ds = \frac{\alpha_{mj}}{\beta_{mj}}$$

By the spectral theorem for non-negative matrices, the process is stationary and subcritical if and only if the spectral radius satisfies:

$$\rho(\Gamma) = \max_i |\lambda_i(\Gamma)| < 1$$

### 1.2 Ozaki Linear-Time $O(N)$ Likelihood Recursion
The log-likelihood of observation over $[0, T]$ with arrival events $\{(t_k, m_k)\}_{k=1}^N$ is given by:

$$\ln L(\theta) = \sum_{k=1}^N \ln \lambda_{m_k}(t_k) - \sum_{m=1}^M \int_0^T \lambda_m(s) ds$$

Using the Ozaki state variable $R_{mj}(t_k) = \sum_{t_l < t_k, m_l = j} e^{-\beta_{mj}(t_k - t_l)}$, the state obeys the constant-time step update:

$$R_{mj}(t_k) = e^{-\beta_{mj}(t_k - t_{k-1})} \left( R_{mj}(t_{k-1}) + \mathbf{1}_{\{m_{k-1} = j\}} \right)$$

This reduces log-likelihood evaluation from $O(N^2)$ to $O(M^2 N)$.

---

## 2. Stochastic Optimal Execution (Almgren-Chriss Framework)

An agent liquidates an initial position $X_0$ shares across horizon $T$ over $N$ discrete intervals $\tau = T / N$ at times $t_k = k \tau$.

### 2.1 Price Dynamics and Market Impact
Mid-price follows arithmetic Brownian motion with permanent impact:

$$S_k = S_{k-1} + \sigma \tau^{1/2} \xi_k - \tau \gamma \left(\frac{n_k}{\tau}\right)$$

Execution price for child order $n_k = x_{k-1} - x_k$ incorporates temporary impact:

$$\tilde{S}_k = S_{k-1} - \eta \left(\frac{n_k}{\tau}\right)$$

### 2.2 Mean-Variance Objective Functional
The objective minimizes expected capture shortfall plus variance penalty:

$$\min_{\{x_k\}} E[x] + \lambda V[x]$$

$$E[x] = \frac{1}{2} \gamma X_0^2 + \frac{\eta}{\tau} \sum_{k=1}^N (x_{k-1} - x_k)^2, \quad V[x] = \sigma^2 \tau \sum_{k=1}^N x_k^2$$

The Euler-Lagrange finite difference equation yields the hyperbolic trajectory:

$$x_j = \frac{\sinh(\kappa (T - t_j))}{\sinh(\kappa T)} X_0$$

where $\kappa$ is the unique solution to $\cosh(\kappa \tau) = 1 + \frac{1}{2} \frac{\lambda \sigma^2}{\eta} \tau^2$.

---

## 3. High-Frequency Market Making (Cartea-Jaimungal HJB)

Let the agent quote bid and ask depths $(\delta_t^b, \delta_t^a)$ relative to mid-price $S_t$.

### 3.1 Intensity of Order Fill
Incoming order executions follow controlled Poisson processes:

$$\lambda^b(\delta_t^b) = A e^{-k \delta_t^b}, \quad \lambda^a(\delta_t^a) = A e^{-k \delta_t^a}$$

### 3.2 Value Function and Indifference Quotes
Under inventory running penalty $\phi q^2$, the Hamilton-Jacobi-Bellman equation is:

$$\partial_t V + \frac{1}{2} \sigma^2 \partial_{SS} V - \phi q^2 + \max_{\delta^a} \lambda^a(\delta^a) [V(t, S, q - 1) + S + \delta^a - V] + \max_{\delta^b} \lambda^b(\delta^b) [V(t, S, q + 1) - S + \delta^b - V] = 0$$

Separation ansatz $V(t, S, q) = q S + h(t, q)$ yields the Avellaneda-Stoikov reservation price:

$$R(S, q, t) = S - q \gamma \sigma^2 (T - t)$$

$$\delta^{a*}(q) = R(S, q, t) - S + \frac{1}{\gamma} \ln\left(1 + \frac{\gamma}{k}\right)$$

$$\delta^{b*}(q) = S - R(S, q, t) + \frac{1}{\gamma} \ln\left(1 + \frac{\gamma}{k}\right)$$

---

## 4. Continuous-Time Jump-Diffusion & Stochastic Volatility

### 4.1 Merton (1976) Log-Normal Jump-Diffusion
Spot price satisfies:

$$\frac{dS_t}{S_{t^-}} = (r - q - \lambda \bar{k}) dt + \sigma dW_t + (J - 1) dN_t$$

where $\ln J \sim \mathcal{N}(\mu_J, \sigma_J^2)$ and $\bar{k} = e^{\mu_J + \frac{1}{2}\sigma_J^2} - 1$.

The characteristic function is given in closed form:

$$\phi_{\text{Merton}}(u) = \exp\left( i u \left(\ln S_0 + (r - q - \lambda \bar{k} - \frac{1}{2}\sigma^2)T\right) - \frac{1}{2} u^2 \sigma^2 T + \lambda T \left(e^{i u \mu_J - \frac{1}{2} u^2 \sigma_J^2} - 1\right) \right)$$

### 4.2 Heston (1993) Characteristic Function
Under CIR variance dynamics $dv_t = \kappa(\theta - v_t) dt + \xi \sqrt{v_t} dW_t^v$:

$$\phi(u) = \exp\left( C(u, T) + D(u, T) v_0 + i u (\ln S_0 + (r - q)T) \right)$$

Numerical stability is guaranteed via the Albrecher formulation:

$$D(u, T) = \frac{\kappa - \rho \xi u i - d}{\xi^2} \left( \frac{1 - e^{-d T}}{1 - g e^{-d T}} \right)$$

$$d = \sqrt{(\rho \xi u i - \kappa)^2 + \xi^2 (u i + u^2)}$$

$$g = \frac{\kappa - \rho \xi u i - d}{\kappa - \rho \xi u i + d}$$

---

## 5. Microstructure Toxicity: VPIN

Volume-Synchronized Probability of Toxicity (Easley, Lopez de Prado, O'Hara 2012) divides continuous volume into constant packets of size $V$:

$$\text{VPIN} = \frac{\sum_{\tau=1}^N |V_\tau^B - V_\tau^S|}{N \cdot V}$$

where buy volume $V_\tau^B$ and sell volume $V_\tau^S$ are classified via the Bulk Volume Classification (BVC) normal CDF kernel over standard returns.
