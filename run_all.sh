#!/bin/bash
# Reproduce the exhaustive search with two independent programs and compare them.
#   sh run_all.sh          m = 4..12 (a few minutes on 16 cores)
#   FULL=1 sh run_all.sh   also m = 13, i.e. W_14 (about 35 minutes in all on 16 cores)
# Requirements: gcc/g++ with OpenMP, python3. Exit status 0 iff every comparison and every check passes.
cd "$(dirname "$0")"
set -u
gcc -O3 -march=native -fopenmp -o wheel_search_mrv wheel_search_mrv.c || exit 2
g++ -O3 -march=native -fopenmp -std=c++17 -o fastspoke fastspoke.cpp || exit 2
mkdir -p results
python3 check_certificates.py wheel_certificates_W5_W13.json | tee results/published_check.txt
grep -q "False" results/published_check.txt && { echo "published labelings: FAIL"; exit 1; }
MS="${MS:-4 5 6 7 8 9 10 11 12}"
[ "${FULL:-0}" = 1 ] && MS="$MS 13"
status=0
for m in $MS; do
  pl=0; [ "$m" -ge 12 ] && pl=1   # print the leaves (LEAF lines) for m >= 12
  ./wheel_search_mrv "$m" 100000000 4 "$pl" > "results/A_m$m.txt"; ea=$?
  ./fastspoke "$m" a 1 3 1 > "results/B_m$m.txt" 2> "results/B_m$m.err"; eb=$?
  # B aborts with a message containing BUG if a recomputed geodesic weight disagrees
  if [ "$ea" != 0 ] || [ "$eb" != 0 ] || grep -qi "bug" "results/B_m$m.err"; then
    echo "W_$((m+1)) (m=$m): program failed (exit A=$ea B=$eb)"; status=1; continue
  fi
  la=$(sed -n 's/.*spoke_arrangements=\([0-9]*\).*/\1/p' "results/A_m$m.txt")
  sa=$(sed -n 's/.* solutions=\([0-9]*\).*/\1/p' "results/A_m$m.txt")
  lb=$(sed -n 's/.*spoke_leaves=\([0-9]*\).*/\1/p' "results/B_m$m.txt")
  sb=$(sed -n 's/.* solutions=\([0-9]*\).*/\1/p' "results/B_m$m.txt")
  va=$(python3 verify_sols.py < "results/A_m$m.txt")
  vb=$(python3 verify_sols.py < "results/B_m$m.txt")
  # the two programs must print the same set of normalized solutions
  ca=$(grep '^SOL' "results/A_m$m.txt" | sort | md5sum | cut -c1-12)
  cb=$(grep '^SOL' "results/B_m$m.txt" | sort | md5sum | cut -c1-12)
  ok=OK
  { [ "$la" = "$lb" ] && [ "$sa" = "$sb" ] && [ "$ca" = "$cb" ] && echo "$va" | grep -q "bad=0" && echo "$vb" | grep -q "bad=0"; } || { ok=MISMATCH; status=1; }
  echo "W_$((m+1)) (m=$m): leaves A=$la B=$lb  solutions A=$sa B=$sb  solution-set A=$ca B=$cb  A:[$va] B:[$vb]  $ok"
done
exit $status
