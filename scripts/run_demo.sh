#!/usr/bin/env bash
set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$DIR"

echo "=============================================================================="
echo "                   SOLITON: STOCHASTIC ARRAY LABORATORY                       "
echo "=============================================================================="

# Build native binary if not built
if [ ! -f "bin/soliton" ]; then
    echo "Compiling native C engine..."
    make all
fi

echo ""
echo "[1/4] Running Engine Unit Tests and Verification..."
python3 scripts/run_tests.py

echo ""
echo "[2/4] Testing Native Array Engine REPL / Expression Evaluation..."
./bin/soliton -e "x: 10 20 30 40 50; y: 1 2 3 4 5; x * y"

echo ""
echo "[3/4] Running High-Throughput Performance Benchmarks..."
python3 tools/benchmark.py

echo ""
echo "[4/4] Generating Synthetic Order Book Stream..."
python3 tools/synthetic_feed.py

echo ""
echo "=============================================================================="
echo "DEMO COMPLETE: All systems operational with zero defects."
echo "=============================================================================="
