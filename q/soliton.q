/ ==============================================================================
/ soliton: non-equilibrium statistical mechanics and stochastic point processes
/ master library loader
/ ==============================================================================

\echo "Initializing soliton stochastic array laboratory..."

/ Load core mathematical foundations
\l q/core/math.q
\l q/core/linalg.q
\l q/core/stats.q
\l q/core/optim.q
\l q/core/fft.q

/ Load market microstructure and point processes
\l q/microstructure/book.q
\l q/microstructure/hawkes.q
\l q/microstructure/impact.q

/ Load stochastic optimal control
\l q/control/almgren_chriss.q
\l q/control/cartea_jaimungal.q

/ Load continuous-time stochastic models
\l q/models/jump_diffusion.q
\l q/models/volatility.q

/ Load risk and streaming analytics
\l q/risk/surface.q
\l q/stream/tickerplant.q
\l q/stream/cep.q
\l q/backtest/engine.q

\echo "Soliton laboratory loaded successfully. Namespaces active: .math .la .stats .optim .fft .book .hawkes .impact .ac .cj .jd .vol .risk .tp .rdb .cep .backtest"
