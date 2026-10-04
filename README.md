# The wheels W_14 and W_15 have no geodesic Leech labeling

Prepared 2026/10/4 by Hirotaka Shimizu with the AI assistant Claude (see Credits); W_15 added 2026/10/5.
Status: computational result, reproduced by two separately written programs. Not yet formally verified
(no proof-checker certificate), and not yet seen by the author of the paper below.

## Claim

The wheel W_14 (a hub joined to every vertex of a 13-cycle) is **not** geodesic Leech.
Neither is W_15 (a hub joined to every vertex of a 14-cycle).

If correct, this settles the first two open cases left in

> Junyeop Yim, *A finiteness theorem for geodesic Leech wheels*, arXiv:2609.02544v1 (2026/9/2),

which shows that W_5, ..., W_13 are geodesic Leech and that a geodesic Leech wheel W_n (n ≥ 5) has n ≤ 40,
leaving W_14, ..., W_40 open. The author's repository
(github.com/junyeobe0315/geodesic-leech-wheels, README) reports a simulated-annealing search for W_14
(8 runs × 18 h, about 1.47×10^11 steps) that found no labeling; that search was not exhaustive.

## Definitions

A labeling assigns a positive integer to every edge. The weight of a path is the sum of its labels.
A geodesic is a shortest path in the unlabeled graph; t_gp(G) is the number of geodesics
(over unordered vertex pairs, every shortest path counted once).
G is geodesic Leech if some labeling makes the multiset of geodesic weights equal to {1, 2, ..., t_gp(G)}.

W_n: hub h, rim vertices 0, ..., m−1 with m = n − 1, spokes s_i = label of h–i, rims r_i = label of i–(i+1) (indices mod m).

## Reduction (m ≥ 4)

**Lemma 1 (geodesics of a wheel).** The geodesics of W_n are
- the m spokes and the m rim edges (distance 1);
- for each rim pair at cyclic distance 2, i.e. {i, i+2}: the rim path i, i+1, i+2 and the spoke path i, h, i+2;
- for each rim pair at cyclic distance ≥ 3: the spoke path i, h, j only.

For m = 4 the pairs {0,2}, {1,3} have two rim paths each; these are exactly the four paths r_i + r_{i+1},
so the description below is unchanged. In all cases t_gp = 2m + m + m(m−3)/2 = m(m+3)/2 (= 104 for W_14).

Hence the multiset of geodesic weights is

    S ∪ P ∪ R ∪ Q,   S = {s_i},  P = {s_i + s_j : i, j not adjacent on the rim},
                      R = {r_i},  Q = {r_i + r_{i+1}}.

**Lemma 2 (decoupling).** S ∪ P depends only on the cyclic sequence of spokes, and R ∪ Q only on the cyclic
sequence of rims. No geodesic uses both a spoke and a rim edge. So W_n is geodesic Leech iff there are a cyclic
spoke sequence and a cyclic rim sequence with S ∪ P ∪ R ∪ Q = {1, ..., T} as multisets, T = m(m+3)/2
(so all T values are distinct).
The rotation/reflection of one sequence relative to the other is irrelevant. Each sequence can therefore be
normalized independently under the dihedral group:
- spokes: s_0 = min, s_1 < s_{m−1};
- rims: r_0 = min L (the smallest leftover value must be a rim label, since every element of Q exceeds two
  elements of R), r_1 < r_{m−1}.

Both programs rely on Lemma 2 and on the spoke normalization. A check that does not assume them
(no normalization, raw counts) was run only for m = 4, 5, 6 (below).

Both lemmas were re-derived by the reviewer session and checked by machine: program B below builds the geodesics
by BFS (it does not use the formula), and confirms at run time that there are T of them, that none mixes spokes
and rims, and that each spoke lies on m−2 geodesics and each rim edge on 3.
The un-normalized count of solutions equals (2m)^2 × the normalized count for m = 4, 5, 6 (960 = 15·64,
14600 = 146·100, 22896 = 159·144), which confirms the two independent dihedral normalizations.

**Search.** Enumerate every normalized cyclic spoke sequence whose values S ∪ P are pairwise distinct and lie in
[1, T] (a "leaf"). For each leaf, let L = {1..T} \ (S ∪ P) (2m values) and decide by depth-first search whether
L = R ∪ Q for a normalized cyclic rim sequence.

## Results

Counts of normalized (spoke sequence, rim sequence) pairs. A and B are the two programs below.
"Leaves" for B are all admissible normalized spoke sequences; A's table entries for m ≤ 12 come from runs with
its sum bound switched on. The two counts agreeing is an observation: in these cases the sum bound removed no leaf.
For m = 10..13, A was also run without the sum bound and gave the same leaf counts.

| W_n | m | T | leaves (A = B) | solutions (A = B) |
|---|---:|---:|---:|---:|
| W_5 | 4 | 14 | 376 | 15 |
| W_6 | 5 | 20 | 3572 | 146 |
| W_7 | 6 | 27 | 28447 | 159 |
| W_8 | 7 | 35 | 131758 | 477 |
| W_9 | 8 | 44 | 279822 | 1653 |
| W_10 | 9 | 54 | 214668 | 319 |
| W_11 | 10 | 65 | 44021 | 27 |
| W_12 | 11 | 77 | 2773 | 5 |
| W_13 | 12 | 90 | 137 | 1 |
| **W_14** | **13** | **104** | **2** | **0** |
| **W_15** | **14** | **119** | **0** | **0** |

- For m ≤ 12 the two programs print the same set of solutions, and every printed solution was checked
  by recomputing all geodesics by BFS (`verify_sols.py`).
- The published labelings of W_5, ..., W_13 (Table 1 of the paper; `wheel_certificates_W5_W13.json` from the
  author's repository) all occur among the solutions.
- For W_13 there is exactly one normalized pair, the published one:
  spokes 2,42,3,43,46,21,25,10,36,16,30,44; rims 1,34,50,20,9,8,6,77,4,7,15,74.
  Note that a normalized pair is a class under *independent* dihedral actions on spokes and rims. Under the
  automorphism group of the wheel (one dihedral action on both) this one pair gives 2m = 24 classes,
  and (2m)^2 = 576 labelings in all.
- **W_14:** both programs find exactly two leaves and no solution. The two leaves are:

      spokes 1,50,2,44,4,47,8,38,13,28,18,33,51
      spokes 1,47,4,51,33,18,28,13,38,8,44,2,50

  Both leave the same 26 values for the rims,

      L = 7,11,16,23,24,25,27,43,67,70,73,74,76,81,84,86,87,90,92,93,96,99,100,102,103,104,

  with Σ L = 1753 ≡ 1 (mod 3). Since Σ L = Σ R + Σ Q = 3 Σr, no rim sequence fits (program A rejects both by this
  test; program B, which does not use it, rejects both by exhaustive rim search).
- **W_15:** neither program finds any leaf: no normalized spoke sequence of length 14 has its 91 values
  S ∪ P pairwise distinct in [1, 119]. So the rims are never reached.
  - A (with the sum bound): 6.77×10^11 nodes, 0 leaves, 0 solutions.
  - B (nosum): 2.56×10^11 nodes, 0 leaves, 0 solutions. This count does not depend on the sum bound or the
    mod-3 test.

  Both programs ran in Docker (image `gcc:14`, 16 threads), A from 2026/10/4 19:22:49 to 22:09:22 JST
  and B from 22:09:22 to 2026/10/5 00:52:02 JST, both with exit status 0 (`run_m14.sh`, outputs in
  `results/m14/`). A reports 9993 s, which includes 14 minutes (19:33–19:47) when the container was paused
  so that a timing comparison could use the whole machine. B reports 9761 s.

## Programs

**A: `wheel_search_mrv.c`** (search session). C, 128-bit value sets, OpenMP. Places spokes in dynamic order
(the position with the fewest candidates first), with forward checking over all unplaced positions.
Optional pruning, each a necessary condition:
- Σs ≤ SMAX from the identity (m−2)Σs + 3Σr = T(T+1)/2 and Σr ≥ m(m+1)/2 (disable with `-DNO_SUMBOUND`);
- at a leaf, Σ L ≡ 0 (mod 3), since Σ L = 3Σr.

m = 13: 8.44×10^10 nodes, 1895 s on 16 threads (with the sum bound); 8.44×10^10 nodes, 2278 s on 12 threads
(without it). Both runs: 2 leaves, 0 passing the mod-3 test, 0 solutions.

**B: `fastspoke.cpp`** (reviewer session; written separately, after that session had reviewed A, with a different design).
C++17. Geodesics from BFS; the pairs whose spoke sums count are read off the length-2 spoke geodesics.
Fixed placement order 0, 2, 4, ..., 1, 3, 5, ...; forward checking. With the last argument `1` ("nosum") it uses
**neither the sum bound nor the mod-3 test**: every leaf goes to a generic rim DFS over the geodesics.
At each leaf it recomputes the spoke geodesic weights from the BFS paths and aborts on any mismatch.

m = 13 with nosum: 3.12×10^10 nodes, 2 leaves, 2 rim searches, 0 solutions, 1715 s on 12 threads.

So the conclusion for W_14 does not depend on the sum bound or the mod-3 test.

The times above are from our first runs (other jobs may have shared the machine). In the run of `run_all.sh`
recorded below (16 threads, the two programs one after the other), m = 13 took 870 s for A (with the sum bound)
and 799 s for B (nosum), with the same node, leaf and solution counts.

## How to verify

Requirements: gcc and g++ with OpenMP, python3, and GNU coreutils (`md5sum`), e.g. Linux or WSL.

    sh run_all.sh            # m = 4..12 with both programs, a few minutes on 16 cores
    FULL=1 sh run_all.sh     # also m = 13 (W_14), about 35 minutes in all on 16 cores

W_15 (m = 14) is run by a separate script, about 5.5 hours on 16 cores; it does not compare the outputs:

    sh run_m14.sh            # or: docker run -d --name w15 -v "$PWD":/w -w /w gcc:14 sh run_m14.sh

For each m `run_all.sh` compares the leaf counts, the solution counts and the sets of printed solutions of A and B,
and checks every solution with `verify_sols.py`. It also checks the published labelings (`check_certificates.py`).
Exit status 0 iff everything agrees. Outputs go to `results/`.

For m ≥ 12 both programs also list their leaves (A: `LEAF` lines in `results/A_m*.txt`; B: `LEAF` lines in
`results/B_m*.err`), so the leaves themselves can be compared, not only their number.

Our run of `FULL=1 sh run_all.sh`: 2026/10/4 03:58:41–04:31:23 JST (32 min 42 s), WSL2 Ubuntu, 16 logical cores,
exit status 0; every m = 4..13 reports OK (`results/run_all_full.txt`).
That run was made before the source comments were translated (see below), so the first lines of
`results/run_all_full.txt` show the old labels `A一致`/`B一致` ("A agrees"/"B agrees") instead of `A_ok`/`B_ok`.
The files `results/*_m4..12*` are from a later `sh run_all.sh` (m = 4..12, exit status 0) after the translation;
`results/*_m13*` are from the FULL run.

## Files

| File | Role |
|---|---|
| `wheel_search_mrv.c` | program A |
| `fastspoke.cpp` | program B |
| `check_certificates.py` | checks the published W_5..W_13 labelings two ways (BFS geodesics, and the formula of Lemma 1) |
| `verify_sols.py` | checks printed solutions by BFS geodesics |
| `wheel_certificates_W5_W13.json` | the author's published labelings, copied unchanged from `data/` of github.com/junyeobe0315/geodesic-leech-wheels (CC BY 4.0) |
| `run_all.sh` | reproduction script |
| `run_m14.sh` | runs A and B for W_15 (m = 14) |
| `results/` | outputs of our runs (`results/m14/` for W_15) |
| `test_small.sh` | quick check of `run_all.sh` (m = 4..9, 4 threads), plus a check that the split depth does not change the result for m = 4 |
| `LICENSE` | MIT license for the code |

Source comments were translated from Japanese to English on 2026/10/4; the compiled assembly of both programs
is unchanged by the translation (checked with gcc/g++ -O3 -S).

## Limitations

- This is an exhaustive computer search, cross-checked by two programs written in separate sessions by the same
  kind of model (Claude). It is not a proof-checked certificate (for example, no SAT/DRAT proof).
- W_16, ..., W_40 remain open. The leaf counts drop quickly with m (137 for W_13, 2 for W_14, 0 for W_15),
  but the number of search nodes of program A grows by a factor of about 8 per step.
- For W_15, A was run only with the sum bound; the count without it comes from B alone.
- Program B used to accept a split depth outside 1..m; with a split depth above m it dropped the leaves found while
  collecting tasks and reported 0 solutions (found by a GPT review on 2026-10-05). It now rejects such a value.
  None of the recorded runs is affected: all of them use split depth 3 for B (and 4 for A, which handles any depth).

## Credits

Hirotaka Shimizu (GitHub: yygotm). The computations, the programs and the write-up were produced with the AI
assistant Claude (Anthropic), in two separate sessions: one wrote program A; the other reviewed it and then wrote
program B with a different design.

## License

- Code (`wheel_search_mrv.c`, `fastspoke.cpp`, `check_certificates.py`, `verify_sols.py`, `run_all.sh`,
  `test_small.sh`): MIT, see `LICENSE`.
- `wheel_certificates_W5_W13.json`: © Junyeop Yim, CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/),
  from https://github.com/junyeobe0315/geodesic-leech-wheels, unchanged.
- Our outputs in `results/` and this README: CC0 1.0 Universal (public domain dedication),
  https://creativecommons.org/publicdomain/zero/1.0/
