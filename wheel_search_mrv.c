/*
 * Program A: exhaustive search for geodesic Leech labelings of the wheel W_{m+1},
 * placing spokes in dynamic order (the unplaced position with the fewest candidates first).
 *
 * Normalization: the spoke at position 0 is the minimum, s_1 < s_{m-1};
 * rims: r_0 = min L, r_1 < r_{m-1}.
 * At each node, for every unplaced position p, build the candidate set
 *   A_p = { y : s_0 < y ≤ T - max(placed values not adjacent to p), y ∉ used,
 *           y + s_q ∉ used (q placed and not adjacent to p) }
 * and prune if some A_p is empty, or if the union of the A_p has fewer elements than
 * the number of unplaced positions. The position with the fewest candidates is placed next.
 *
 * Usage:  wheel_search_mrv m [max_print] [split_depth] [print_leaves]
 * Build:  gcc -O3 -march=native -fopenmp -o wheel_search_mrv wheel_search_mrv.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>

typedef unsigned __int128 u128;

static int M, T, SMAX;
static u128 FULL;
static long long max_print = 10;

static inline u128 bit(int v) { return ((u128)1) << v; }
static inline int popc(u128 x) {
    return __builtin_popcountll((uint64_t)x) + __builtin_popcountll((uint64_t)(x >> 64));
}
static inline int lowbit(u128 x) {
    uint64_t lo = (uint64_t)x;
    if (lo) return __builtin_ctzll(lo);
    return 64 + __builtin_ctzll((uint64_t)(x >> 64));
}
static inline u128 range_mask(int lo, int hi) {  /* bits lo..hi (hi ≤ 126) */
    if (hi < lo) return 0;
    return (bit(hi + 1) - 1) & ~(bit(lo) - 1);
}
static inline int adjacent(int p, int q) {
    int d = (p - q + M) % M;
    return d == 1 || d == M - 1;
}

typedef struct {
    int val[16];        /* 0 = unplaced */
    u128 used;          /* S ∪ P */
    int nplaced;
    int ssum;
} State;

typedef struct {
    long long nodes, leaves, mod3, rimok, sols;
} Stats;

static long long printed = 0;

static int rim_dfs(const State *st, u128 rem, int k, int prev, int *r) {
    if (k == M) {
        int closing = prev + r[0];
        if (closing <= T && rem == bit(closing) && r[1] < r[M - 1]) {
            #pragma omp critical
            {
                if (printed++ < max_print) {
                    printf("SOL m=%d spokes=", M);
                    for (int i = 0; i < M; i++) printf("%d%c", st->val[i], i + 1 < M ? ',' : ' ');
                    printf("rims=");
                    for (int i = 0; i < M; i++) printf("%d%c", r[i], i + 1 < M ? ',' : '\n');
                    fflush(stdout);
                }
            }
            return 1;
        }
        return 0;
    }
    int found = 0;
    u128 cand = rem;
    while (cand) {
        int y = lowbit(cand);
        cand &= cand - 1;
        int sum = prev + y;
        if (sum > T) break;
        u128 b = bit(y) | bit(sum);
        if ((rem & b) != b) continue;
        r[k] = y;
        found += rim_dfs(st, rem & ~b, k + 1, y, r);
    }
    return found;
}

static int print_leaves = 0;
static void leaf(const State *st, Stats *S) {
    S->leaves++;
    u128 L = FULL & ~st->used;
    int sum = 0;
    for (u128 x = L; x; x &= x - 1) sum += lowbit(x);
    if (print_leaves) {
        #pragma omp critical
        {
            printf("LEAF m=%d spokes=", M);
            for (int i = 0; i < M; i++) printf("%d%c", st->val[i], i + 1 < M ? ',' : ' ');
            printf("L=");
            int first = 1;
            for (u128 x = L; x; x &= x - 1) { printf(first ? "%d" : ",%d", lowbit(x)); first = 0; }
            printf(" sumL=%d sumL_mod3=%d\n", sum, sum % 3);
            fflush(stdout);
        }
    }
    if (sum % 3) return;
    S->mod3++;
    int r[16];
    r[0] = lowbit(L);
    int f = rim_dfs(st, L & ~bit(r[0]), 1, r[0], r);
    if (f) { S->rimok++; S->sols += f; }
}

/* Build the candidate set of every unplaced position; return the one with the fewest (-1: prune) */
static int choose(const State *st, u128 *candp) {
    u128 sh[16];
    for (int q = 0; q < M; q++) if (st->val[q]) sh[q] = st->used >> st->val[q];
    int best = -1, bestc = 1 << 30;
    u128 uni = 0;
    int nun = 0;
    int s0 = st->val[0];
    for (int p = 0; p < M; p++) {
        if (st->val[p]) continue;
        nun++;
        u128 forb = st->used;
        int mx = 0;
        for (int q = 0; q < M; q++) {
            if (!st->val[q] || adjacent(p, q)) continue;
            forb |= sh[q];
            if (st->val[q] > mx) mx = st->val[q];
        }
        int lo = s0 + 1;
        /* reflection normalization: of positions 1 and M-1, s_1 < s_{M-1} whichever is placed first */
        if (p == M - 1 && st->val[1]) { if (st->val[1] + 1 > lo) lo = st->val[1] + 1; }
        if (p == 1 && st->val[M - 1]) { /* s_1 < s_{M-1} */
            int hi1 = st->val[M - 1] - 1;
            u128 a = range_mask(lo, hi1 < T - mx ? hi1 : T - mx) & ~forb;
            int c = popc(a);
            if (c == 0) return -1;
            uni |= a;
            if (c < bestc) { bestc = c; best = p; *candp = a; }
            continue;
        }
        u128 a = range_mask(lo, T - mx) & ~forb;
        int c = popc(a);
        if (c == 0) return -1;
        uni |= a;
        if (c < bestc) { bestc = c; best = p; *candp = a; }
    }
    if (popc(uni) < nun) return -1;
    return best;
}

static void dfs(State *st, Stats *S) {
    S->nodes++;
    if (st->nplaced == M) { leaf(st, S); return; }
    u128 cand;
    int p = choose(st, &cand);
    if (p < 0) return;
    /* placed values not adjacent to p (their sums with the new value join P) */
    u128 partners = 0;
    for (int q = 0; q < M; q++) if (st->val[q] && !adjacent(p, q)) partners |= bit(st->val[q]);
    int rest = M - st->nplaced - 1;
    int restmin = rest * (st->val[0] + 1) + rest * (rest - 1) / 2;
    u128 saved = st->used;
    while (cand) {
        int x = lowbit(cand);
        cand &= cand - 1;
#ifndef NO_SUMBOUND
        if (st->ssum + x + restmin > SMAX) break;
#endif
        st->val[p] = x;
        st->used = saved | bit(x) | (partners << x);
        st->nplaced++;
        st->ssum += x;
        dfs(st, S);
        st->ssum -= x;
        st->nplaced--;
        st->val[p] = 0;
    }
    st->used = saved;
}

/* For parallelization, collect the states at depth D */
static State *tasks;
static long long ntasks, cap;
static void collect(State *st, int D) {
    if (st->nplaced == D || st->nplaced == M) {
        if (ntasks == cap) { cap = cap ? cap * 2 : 1024; tasks = realloc(tasks, cap * sizeof(State)); }
        tasks[ntasks++] = *st;
        return;
    }
    u128 cand;
    int p = choose(st, &cand);
    if (p < 0) return;
    u128 partners = 0;
    for (int q = 0; q < M; q++) if (st->val[q] && !adjacent(p, q)) partners |= bit(st->val[q]);
    u128 saved = st->used;
    while (cand) {
        int x = lowbit(cand);
        cand &= cand - 1;
        st->val[p] = x;
        st->used = saved | bit(x) | (partners << x);
        st->nplaced++;
        st->ssum += x;
        collect(st, D);
        st->ssum -= x;
        st->nplaced--;
        st->val[p] = 0;
    }
    st->used = saved;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s m [max_print] [split_depth]\n", argv[0]); return 2; }
    M = atoi(argv[1]);
    if (argc > 2) max_print = atoll(argv[2]);
    int D = argc > 3 ? atoi(argv[3]) : 3;
    print_leaves = argc > 4 ? atoi(argv[4]) : 0;
    T = M * (M + 3) / 2;
    if (M < 4 || M > 15 || T > 126) { fprintf(stderr, "m out of range\n"); return 2; }
    FULL = range_mask(1, T);
    SMAX = (T * (T + 1) / 2 - 3 * (M * (M + 1) / 2)) / (M - 2);

    double t0 = omp_get_wtime();
    for (int a = 1; a <= T; a++) {
        State st;
        memset(&st, 0, sizeof st);
        st.val[0] = a; st.used = bit(a); st.nplaced = 1; st.ssum = a;
        collect(&st, D);
    }
    long long nodes = 0, leaves = 0, mod3 = 0, rimok = 0, sols = 0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:nodes,leaves,mod3,rimok,sols)
    for (long long i = 0; i < ntasks; i++) {
        Stats S = {0};
        State st = tasks[i];
        dfs(&st, &S);
        nodes += S.nodes; leaves += S.leaves; mod3 += S.mod3; rimok += S.rimok; sols += S.sols;
    }
    double t1 = omp_get_wtime();
    printf("m=%d n=%d t=%d tasks=%lld nodes=%lld spoke_arrangements=%lld mod3_pass=%lld "
           "spoke_with_rim=%lld solutions=%lld time=%.2fs threads=%d\n",
           M, M + 1, T, ntasks, nodes, leaves, mod3, rimok, sols, t1 - t0, omp_get_max_threads());
    return 0;
}
