"""Check the published labelings of W_5..W_13.

The multiset of geodesic weights is built in two ways and compared with {1..t_gp}:
  A: enumerate every unweighted shortest path between every pair of vertices (no formula)
  B: the formula of Lemma 1 (spokes, rim edges, both paths of each distance-2 pair,
     spoke sums of the pairs at distance >= 3)
"""
import json
import sys
from collections import Counter, deque
from itertools import combinations
from pathlib import Path


def wheel_edges(spokes, rims):
    m = len(spokes)
    hub = m
    w = {}
    for i in range(m):
        w[frozenset((hub, i))] = spokes[i]
        w[frozenset((i, (i + 1) % m))] = rims[i]
    return m + 1, w


def all_geodesic_weights(nv, w):
    adj = {v: [] for v in range(nv)}
    for e in w:
        a, b = tuple(e)
        adj[a].append(b)
        adj[b].append(a)
    out = []
    for s, t in combinations(range(nv), 2):
        dist = {s: 0}
        q = deque([s])
        while q:
            u = q.popleft()
            for x in adj[u]:
                if x not in dist:
                    dist[x] = dist[u] + 1
                    q.append(x)
        # walk back from t along decreasing dist to enumerate every shortest path
        def paths(v):
            if v == s:
                yield 0
                return
            for x in adj[v]:
                if dist.get(x) == dist[v] - 1:
                    for p in paths(x):
                        yield p + w[frozenset((x, v))]
        out.extend(paths(t))
    return out


def formula_weights(spokes, rims):
    m = len(spokes)
    out = list(spokes) + list(rims)
    out += [rims[i] + rims[(i + 1) % m] for i in range(m)]
    for i, j in combinations(range(m), 2):
        if (j - i) % m not in (1, m - 1):
            out.append(spokes[i] + spokes[j])
    return out


def main():
    data = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    for c in data["certificates"]:
        n, m, t = c["n"], c["m"], c["t_gp"]
        s, r = c["spokes"], c["rims"]
        nv, w = wheel_edges(s, r)
        target = Counter(range(1, t + 1))
        a = all_geodesic_weights(nv, w)
        b = formula_weights(s, r)
        ok_a = Counter(a) == target
        ok_b = Counter(b) == target
        # for m=4 each diagonal pair has two rim paths; the four sums r_i+r_{i+1} count both, so B still holds
        print(f"W_{n}: m={m} t_gp={t} m(m+3)/2={m*(m+3)//2} "
              f"#A={len(a)} A_ok={ok_a} #B={len(b)} B_ok={ok_b}")


if __name__ == "__main__":
    main()
