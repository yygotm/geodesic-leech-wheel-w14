#!/bin/bash
# Quick check of run_all.sh (m = 4..9, 4 threads)
cd "$(dirname "$0")"
OMP_NUM_THREADS=4 MS="4 5 6 7 8 9" bash run_all.sh
st=$?
echo "exit=$st"
# Regression: the split depth must not change the result, and B must reject a split depth outside 1..m
# (before this check, B with split_depth > m dropped the leaves found while collecting tasks and reported 0 solutions)
for D in 1 2 3 4 5; do
  s=$(OMP_NUM_THREADS=1 ./wheel_search_mrv 4 0 $D | sed -n 's/.* solutions=\([0-9]*\).*/\1/p')
  [ "$s" = 15 ] || { echo "A m=4 split=$D: solutions=$s (expected 15)"; st=1; }
done
for D in 1 2 3 4; do
  s=$(OMP_NUM_THREADS=1 ./fastspoke 4 a 0 $D 1 2>/dev/null | sed -n 's/.* solutions=\([0-9]*\).*/\1/p')
  [ "$s" = 15 ] || { echo "B m=4 split=$D: solutions=$s (expected 15)"; st=1; }
done
for D in 0 5; do
  OMP_NUM_THREADS=1 ./fastspoke 4 a 0 $D 1 >/dev/null 2>&1 && { echo "B m=4 split=$D: accepted (expected rejection)"; st=1; }
done
echo "split regression: $([ $st = 0 ] && echo OK || echo FAIL)"
exit $st
