#!/bin/sh
# W_15 (m = 14): program A (sum bound and mod-3 test, prints leaves), then program B (nosum), all threads.
# Needs only gcc/g++ with OpenMP. About 5.5 hours in all on 16 cores. With Docker, for example:
#   docker run -d --name w15 -v "$PWD":/w -w /w gcc:14 sh run_m14.sh
# Outputs go to results/m14/; log.txt records the start and end times and the exit status of each program.
cd "$(dirname "$0")"
mkdir -p results/m14
gcc -O3 -march=native -fopenmp -o wheel_search_mrv wheel_search_mrv.c || exit 2
g++ -O3 -march=native -fopenmp -std=c++17 -o fastspoke fastspoke.cpp || exit 2
echo "nproc=$(nproc) A start $(date '+%Y/%m/%d %H:%M:%S %Z')" > results/m14/log.txt
./wheel_search_mrv 14 100000000 4 1 > results/m14/A_m14.txt; echo "A exit=$? end $(date '+%Y/%m/%d %H:%M:%S %Z')" >> results/m14/log.txt
echo "B start $(date '+%Y/%m/%d %H:%M:%S %Z')" >> results/m14/log.txt
./fastspoke 14 a 1 3 1 > results/m14/B_m14.txt 2> results/m14/B_m14.err; echo "B exit=$? end $(date '+%Y/%m/%d %H:%M:%S %Z')" >> results/m14/log.txt
echo done >> results/m14/log.txt
