// Program B (reviewer session's independent implementation): a faster version of its earlier
// program rimfirst.cpp (not included here), fast enough for m = 11..13.
//
// Differences from program A (wheel_search_mrv.c):
//   - Geodesics come from BFS. The spoke pairs whose sums count are read off the length-2 spoke
//     geodesics (no formula is used).
//   - Fixed placement order (default 0,2,4,...,1,3,5,..., so non-adjacent partners appear early). No MRV.
//   - Pruning: forward checking of the candidate sets of all unplaced positions (empty set, size of
//     the union), and a sum bound
//       (m-2)·Σs = T(T+1)/2 − ΣL,  ΣL ≥ sum of the 2m smallest currently unused values,
//       Σs ≥ placed sum + sum of the (number unplaced) smallest values in the union of candidates
//     (the coefficients m−2 and 3 are taken from how often each edge occurs in the geodesics).
//   - At a leaf, the spoke geodesic weights are recomputed from the BFS paths and compared with `used`
//     (abort on mismatch); the rims are then searched by a generic per-path DFS (as in that earlier program).
//
// Usage: fastspoke m [order(a=alternating/s=sequential)] [print(0/1)] [split_depth] [nosum(0/1)]
#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
typedef unsigned __int128 u128;

static int m, T, CS, CR;
static long long TOT;
static bool doprint = true;
static bool nosum = false;  // true: drop the sum bound and the mod test entirely (every leaf goes to the rim DFS)
static vector<vector<int>> pathEdges;
static vector<int> rimEndsAt[16];          // rim geodesics completed by rim edge i (0..m-1)
static vector<int> spokePaths;             // indices of the spoke geodesics
static bool pairs_[16][16];                // is {p,q} a length-2 spoke geodesic (= a pair whose sum counts)
static int order_[16];
static u128 FULLM;

static inline u128 B(int v) { return ((u128)1) << v; }
static inline int pc(u128 x) { return __builtin_popcountll((uint64_t)x) + __builtin_popcountll((uint64_t)(x >> 64)); }
static inline int lo1(u128 x) { uint64_t l = (uint64_t)x; return l ? __builtin_ctzll(l) : 64 + __builtin_ctzll((uint64_t)(x >> 64)); }
static inline u128 upto(int hi) { return hi >= 127 ? ~(u128)0 >> 1 : (B(hi + 1) - 1); }  // bits 0..hi
static long long lowk(u128 x, int k) {  // sum of the k smallest elements of x (a huge value if x has fewer than k)
    long long s = 0;
    for (int c = 0; c < k; c++) { if (!x) return 1LL << 40; int v = lo1(x); s += v; x &= x - 1; }
    return s;
}

struct St {
    int val[16];  // 0 = unplaced
    u128 used;
    int depth;
    long long ssum;
};
struct Res {
    long long nodes = 0, leaves = 0, sols = 0, rimseq = 0;
    vector<string> out;
};

// ---- rims (generic) ----
static void rimDfs(const St &s, int k, u128 used, long long rsum, long long target, int *r, Res &R) {
    if (k == m) {
        if ((!nosum && rsum != target) || used != FULLM || !(r[1] < r[m - 1])) return;
        R.sols++;
        if (doprint) {
            string o = "SOL m=" + to_string(m) + " spokes=";
            for (int i = 0; i < m; i++) o += to_string(s.val[i]) + (i + 1 < m ? "," : "");
            o += " rims=";
            for (int i = 0; i < m; i++) o += to_string(r[i]) + (i + 1 < m ? "," : "");
            R.out.push_back(o);
        }
        return;
    }
    int lo = k == 0 ? 1 : r[0] + 1;
    for (int v = lo; v <= T; v++) {
        if (used & B(v)) continue;
        int rest = m - k - 1;
        long long base = k == 0 ? v + 1 : r[0] + 1;
        if (rsum + v + rest * base + rest * (rest - 1) / 2 > target) break;
        r[k] = v;
        u128 u = used;
        bool ok = true;
        for (int p : rimEndsAt[k]) {
            int sm = 0;
            for (int e : pathEdges[p]) sm += r[e - m];
            if (sm > T || (u & B(sm))) { ok = false; break; }
            u |= B(sm);
        }
        if (ok) rimDfs(s, k + 1, u, rsum + v, target, r, R);
    }
}

static void leaf(const St &s, Res &R) {
    R.leaves++;
    // recompute the spoke geodesic weights from the paths and compare
    u128 chk = 0;
    for (int p : spokePaths) {
        int sm = 0;
        for (int e : pathEdges[p]) sm += s.val[e];
        if (sm < 1 || sm > T || (chk & B(sm))) { fprintf(stderr, "BUG: spoke path collision at leaf\n"); exit(3); }
        chk |= B(sm);
    }
    if (chk != s.used) { fprintf(stderr, "BUG: used mismatch at leaf\n"); exit(3); }
    if (m >= 12) {
        // print the leaf (an admissible spoke sequence), with ΣL mod 3
        long long sL = 0;
        for (int v = 1; v <= T; v++) if (!(s.used & B(v))) sL += v;
        string o = "LEAF spokes=";
        for (int i = 0; i < m; i++) o += to_string(s.val[i]) + (i + 1 < m ? "," : "");
        o += " sum_s=" + to_string(s.ssum) + " sumL=" + to_string(sL) + " sumL_mod3=" + to_string(sL % 3);
        #pragma omp critical
        fprintf(stderr, "%s\n", o.c_str());
    }
    long long rem = TOT - (long long)CS * s.ssum;
    if (!nosum && (rem <= 0 || rem % CR)) return;
    R.rimseq++;
    int r[16] = {0};
    rimDfs(s, 0, s.used, 0, nosum ? (1LL << 40) : rem / CR, r, R);
}

// ---- spokes ----
// candidate set of position p
static u128 domain(const St &s, int p) {
    int s0 = s.val[0];
    u128 forb = s.used;
    int mx = 0;
    for (int q = 0; q < m; q++) {
        if (!s.val[q] || !pairs_[p][q]) continue;
        forb |= s.used >> s.val[q];
        mx = max(mx, s.val[q]);
    }
    int hi = T - mx;
    if (p == 1 && s.val[m - 1]) hi = min(hi, s.val[m - 1] - 1);  // reflection: s_1 < s_{m-1}
    int lo = s0 + 1;
    if (p == m - 1 && s.val[1]) lo = max(lo, s.val[1] + 1);
    if (hi < lo) return 0;
    return upto(hi) & ~upto(lo - 1) & ~forb;
}

// forward checking and the sum bound; true if the state survives
static bool feasible(const St &s) {
    int un = m - s.depth;
    if (un == 0) return true;
    u128 uni = 0;
    for (int i = s.depth; i < m; i++) {
        u128 d = domain(s, order_[i]);
        if (!d) return false;
        uni |= d;
    }
    if (pc(uni) < un) return false;
    if (nosum) return true;
    long long low2m = lowk(FULLM & ~s.used, 2 * m);
    long long smax = (TOT - low2m) / CS;  // floor
    if (s.ssum + lowk(uni, un) > smax) return false;
    return true;
}

static void dfs(St &s, Res &R, int stopDepth, vector<St> *collect) {
    R.nodes++;
    if (collect && s.depth == stopDepth) { collect->push_back(s); return; }
    if (s.depth == m) { leaf(s, R); return; }
    int p = order_[s.depth];
    u128 d = domain(s, p);
    u128 partners = 0;
    for (int q = 0; q < m; q++) if (s.val[q] && pairs_[p][q]) partners |= B(s.val[q]);
    u128 saved = s.used;
    while (d) {
        int x = lo1(d);
        d &= d - 1;
        s.val[p] = x;
        s.used = saved | B(x) | (partners << x);
        s.depth++;
        s.ssum += x;
        if (feasible(s)) dfs(s, R, stopDepth, collect);
        s.ssum -= x;
        s.depth--;
        s.val[p] = 0;
    }
    s.used = saved;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s m [a|s] [print] [split]\n", argv[0]); return 2; }
    m = atoi(argv[1]);
    char ord = argc > 2 ? argv[2][0] : 'a';
    if (argc > 3) doprint = atoi(argv[3]);
    int split = argc > 4 ? atoi(argv[4]) : 3;
    if (argc > 5) nosum = atoi(argv[5]);
    T = m * (m + 3) / 2;
    if (m < 4 || m > 15 || T > 126) return 2;
    if (split < 1 || split > m) { fprintf(stderr, "split_depth must be in 1..m (leaves above the split depth would not be counted)\n"); return 2; }
    FULLM = upto(T) & ~(u128)1;
    int hub = m;
    vector<vector<pair<int, int>>> adj(m + 1);
    for (int i = 0; i < m; i++) {
        adj[hub].push_back({i, i}); adj[i].push_back({hub, i});
        int j = (i + 1) % m;
        adj[i].push_back({j, m + i}); adj[j].push_back({i, m + i});
    }
    for (int s = 0; s <= m; s++) {
        vector<int> d(m + 1, -1); d[s] = 0;
        deque<int> q{s};
        while (!q.empty()) { int u = q.front(); q.pop_front(); for (auto [x, e] : adj[u]) if (d[x] < 0) { d[x] = d[u] + 1; q.push_back(x); } }
        for (int v = s + 1; v <= m; v++) {
            vector<int> cur;
            function<void(int)> go = [&](int u) {
                if (u == s) { pathEdges.push_back(cur); return; }
                for (auto [x, e] : adj[u]) if (d[x] == d[u] - 1) { cur.push_back(e); go(x); cur.pop_back(); }
            };
            go(v);
        }
    }
    if ((int)pathEdges.size() != T) { fprintf(stderr, "path count\n"); return 1; }
    vector<int> cnt(2 * m, 0);
    for (int p = 0; p < T; p++) {
        auto &pe = pathEdges[p];
        bool aS = false, aR = false;
        for (int e : pe) { cnt[e]++; (e < m ? aS : aR) = true; }
        if (aS && aR) { fprintf(stderr, "mixed\n"); return 1; }
        if (aS) {
            spokePaths.push_back(p);
            if (pe.size() == 2) pairs_[pe[0]][pe[1]] = pairs_[pe[1]][pe[0]] = true;
            else if (pe.size() != 1) { fprintf(stderr, "long spoke path\n"); return 1; }
        } else {
            rimEndsAt[*max_element(pe.begin(), pe.end()) - m].push_back(p);
        }
    }
    CS = cnt[0]; CR = cnt[m];
    for (int e = 0; e < 2 * m; e++) if (cnt[e] != (e < m ? CS : CR)) { fprintf(stderr, "uneven\n"); return 1; }
    TOT = (long long)T * (T + 1) / 2;
    int k = 0;
    if (ord == 's') for (int i = 0; i < m; i++) order_[k++] = i;
    else { for (int i = 0; i < m; i += 2) order_[k++] = i; for (int i = 1; i < m; i += 2) order_[k++] = i; }
    if (order_[0] != 0) return 1;
    fprintf(stderr, "m=%d T=%d CS=%d CR=%d order=", m, T, CS, CR);
    for (int i = 0; i < m; i++) fprintf(stderr, "%d ", order_[i]);
    fprintf(stderr, "\n");

    double t0 = omp_get_wtime();
    vector<St> tasks;
    Res R0;
    for (int a = 1; a <= T; a++) {
        St s; memset(&s, 0, sizeof s);
        s.val[0] = a; s.used = B(a); s.depth = 1; s.ssum = a;
        if (feasible(s)) dfs(s, R0, split, &tasks);
    }
    long long nodes = 0, leaves = 0, sols = 0, rimseq = 0;
    vector<string> all;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+ : nodes, leaves, sols, rimseq)
    for (long long i = 0; i < (long long)tasks.size(); i++) {
        Res R;
        St s = tasks[i];
        dfs(s, R, -1, nullptr);
        nodes += R.nodes; leaves += R.leaves; sols += R.sols; rimseq += R.rimseq;
        {
            static long long done = 0;
            long long d;
            #pragma omp atomic capture
            d = ++done;
            if (d % 2000 == 0 || R.sols)
                fprintf(stderr, "progress %lld/%zu tasks, %.0fs, sols_in_task=%lld\n", d, tasks.size(), omp_get_wtime() - t0, R.sols);
        }
        if (doprint && !R.out.empty()) {
            #pragma omp critical
            all.insert(all.end(), R.out.begin(), R.out.end());
        }
    }
    sort(all.begin(), all.end());
    for (auto &o : all) puts(o.c_str());
    printf("m=%d t=%d order=%c tasks=%zu nodes=%lld spoke_leaves=%lld rim_calls=%lld solutions=%lld time=%.2fs threads=%d\n",
           m, T, ord, tasks.size(), nodes, leaves, rimseq, sols, omp_get_wtime() - t0, omp_get_max_threads());
    return 0;
}
