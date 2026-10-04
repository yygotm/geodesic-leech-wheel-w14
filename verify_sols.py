"""Check each SOL line of a search output by enumerating all shortest paths (method A, no formula).

Usage: python verify_sols.py < output_file
"""
import re
import sys
from collections import Counter

from check_certificates import all_geodesic_weights, wheel_edges

ok = bad = 0
for line in sys.stdin:
    mt = re.match(r"SOL m=(\d+) spokes=([\d,]+) rims=([\d,]+)", line)
    if not mt:
        continue
    m = int(mt.group(1))
    s = [int(x) for x in mt.group(2).split(",")]
    r = [int(x) for x in mt.group(3).split(",")]
    t = m * (m + 3) // 2
    nv, w = wheel_edges(s, r)
    if Counter(all_geodesic_weights(nv, w)) == Counter(range(1, t + 1)):
        ok += 1
    else:
        bad += 1
        print("BAD", line.strip())
print(f"verified ok={ok} bad={bad}")
