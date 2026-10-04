#!/bin/bash
# Quick check of run_all.sh (m = 4..9, 4 threads)
cd "$(dirname "$0")"
OMP_NUM_THREADS=4 MS="4 5 6 7 8 9" bash run_all.sh
echo "exit=$?"
