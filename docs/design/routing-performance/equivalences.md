# What may change without changing any result

Six families of changes were examined for exactness: each statement below says what may change, under which
assumptions, and the smallest example that shows the nearest wrong rule fails. Proofs are given in full for the
statements a future change would rely on: the [sink stop](#sink-stop-in-the-flow-solver), the
[reverse cut of a max-flow pair](#the-reverse-pair), the
[class-2 test and program restriction](#failure-classes-and-baseline-derived-programs), the
[prefix rule and tunnel-setup reuse](#reuse-of-the-failure-free-te-solve-and-tunnel-setup) and the
[swap of totals](#the-swap-of-totals). Terms are defined in the [README](README.md).

**What "same results" means.** A change keeps every public output of every call it changes, as `docs/semantics.md`
states them: every field of every returned object (distances, the full equal-cost DAG, the nearest-target DAG, the
canonical path, the KSP order, programs and plans with entry ids, change flags, placement amounts and their per-edge
split, termination and tunnel statuses, flow values, costs, curves, the source-side minimum cut, per-source and
per-tunnel rows, optional outputs, evaluator rows in pair order); the ledger and cache state after the call; identities
and sharing; errors (class, full text, the order of checks and of failing pairs); the state after a raising call; and
the threading contract. Names that start with an underscore are not contract. Equal distances or equal total delivery
do not prove an equal distribution.

**How the statements were checked.** Each question was answered by three independent analyses with the same inputs.
Each wrote its own proof and its own oracle (a program sharing no code with Core) and planted the near-miss variants;
a statement stands when two of three accept its proof with no open objection, the oracle finds no mismatch over the
named domain and the near miss is caught. The equivalence checks include independent models against Core and scratch
C++ probe builds (for example six variants of the flow solver for the sink stop). Correctness evidence is distinct from
a timed landing candidate: only the two proposals and the changes named in [candidates](candidates.md) as "built" were
timed, and no speed is inferred from the counts below, which say how many cases each analysis ran.

## Early-stopping searches

Every proposed early stop can be made exact, and each has a near miss that is wrong. With the correct rule every
target tied at the stopping distance, every tied route, every status and every error, with its validation order, is
unchanged.

| Statement | Assumptions and domain | Near miss that fails | Checked |
| --- | --- | --- | --- |
| A* for a targeted `spf` with a lower bound π (distances to the target set on a graph containing every edge the search may use) returns `target_cost`, the pruned DAG and the canonical path under every tie rule if it keeps popping every key ≤ d* (the nearest targets' distance) and still records tight predecessors of an already-popped node. A targeted `spf` returns no `costs`, so nothing beyond d* is owed | masks, residual gates, several roots, anycast, both directions; g + π fits int64 (costs sum < 2^62) | stopping at the first popped target: root 0, target 1, edges 0→1 cost 2, 0→2 cost 1, 2→1 cost 1; both routes cost 2, so Core's DAG holds all three edges; with π exact every key is 2, the (key, node) order pops target 1 before relay 2, and the stop returns only 0→1 | oracles of 147.7 M, 311,296 and 187,500 + 1,200 queries |
| On Core's bucket queue A* needs three changes: roots enter at π(r), not 0; a push at the current key joins the active bucket; the ring covers the largest reduced cost | the bucket queue | the third condition alone: a push at the current key waits a whole cycle (6,530 of 8,000 cases in one analysis) | a faithful model of the queue |
| Bidirectional search with the strict stop "top_f + top_b > µ" gives the full pruned DAG; with ≥ it is exact only with a rebuild from the detected sums (every pair of labels whose sum equals µ) | live roots and targets count as settled from the start | ≥ with a rebuild from settled meeting edges: 2,701,705 differences of 27.7 M | three oracles |
| KSP spurs: A* spurs with one reverse search per `ksp` call and strict cost bounds (the ceiling, the pending candidate at the k-th place) return Core's list, order and ties | every tie rule | a bound with ≥ instead of > loses the k-th place tie | 13.5 M calls, 43,392, 24,576 + 821 |
| The canonical-suffix shortcut (answer a spur by the best allowed first edge plus the static canonical suffix) is exact only under `FEWER_HOPS` and only after checking that the first edge's head and the whole suffix are allowed for this spur | `FEWER_HOPS` | under `HIGHER_*` it returns the wrong path: three nodes, 0→1 (capacity 4), two members 1→2 (capacities 5 then 10), costs 1; both walks have bottleneck 4, the smaller EdgeId sequence [0, 1] wins, while the static suffix of node 1 is the capacity-10 member; without the head check it differs on 2,160 of 72,000 cases | three oracles |
| T1's CSPF: one static potential for every attempt of a setup, plus a memory of dead nodes kept per request value and reset per setup, returns the same tunnels, order, statuses and holds. A node dead under a lower request is dead under a higher one, so the lower value's set may be unioned in | only exhausted searches enter the memory | a CSPF stopped at the ceiling merges `INSUFFICIENT_BANDWIDTH` and `COST_CEILING_EXCEEDED`; one memory shared in both directions is wrong | 19.8 M setups, 47,200, 415 + 189 on Core |

Keeping π across separate `spf` calls needs a new input or a store, and Core's only store, `ShortestPathCache`,
publishes `hits`, `misses`, `entries` and `bytes`: a contract decision, worth raising only for a workload with repeated
targeted `spf` toward one target set. Inside one `ksp` call (the reverse search CostStretch already runs) and one T1
policy the reuse is internal.

## Sink stop in the flow solver

Core's exact min-cost flow solver runs, per cost tier, one Dijkstra over reduced costs that scans the whole reachable
residual network, then raises every node's potential by its distance and pushes all shortest augmenting paths of the
tier. **The search may stop as soon as the sink is settled**, raising each settled node by its distance `d` and every
other node by the sink's distance `D`: every public result of every solve from zero stays identical, step for step
(tiers, Dinic levels, augmenting paths in arc order, flows, cost, marginal curve, minimum cut, termination). Domain:
every caller that solves from zero (`solve_min_cost_flow`, `solve_subset`, `max_flow_subset` with zero costs, the
private entry), any caps and source caps, ties at `D`, zero-capacity and parallel edges, several sources and sinks,
capped solves ending inside a tier, the final failing search. The private potentials and counters may differ.

Example: arcs s→t cost 1, s→a cost 5, a→t cost 1, capacity 1 each. Tier 1 pops s, then t at `D = 1`, with a still
pending at 5. The full search sets π(a) = 5; the early stop sets π(a) = min(5, 1) = 1. After the push on s→t, tier 2's
reduced costs are s→a 0, a→t 5 under the full update and s→a 4, a→t 1 under the early stop, both non-negative, and both
give the sink distance 5 over the same path s→a→t. Flows, total cost 7 and the curve [(1, 1), (6, 1)] are equal.

**Proof.** Notation: `N_f` the residual network of flow `f` (forward arcs with residual capacity, reverse arcs with the
flow, virtual arcs `s* → source`, `destination → t*` of cost 0); real costs are ≥ 1. For a potential π,
`c_π(x, y) = c(x, y) + π(x) − π(y)`. `R`: vertices reachable from `s*` over positive arcs; `dist(x)`: least cost of an
`s* → x` walk; `SP`: arcs on least-cost `s* → t*` walks; `W`: their vertices. Invariant I(π): `π(s*) = 0` and
`c_π(a) ≥ 0` on every positive arc with tail in `R`; I*(π): the same on every positive arc; *tight*:
`π(t*) = dist(t*)`.

- *Lemma 1.* Under I, `π(x) ≤ dist(x)` on `R`. Under I and tight, `π = dist` on `W`, every arc of `SP` has `c_π = 0`,
  and the admissible `s* → t*` walks (positive arcs, `c_π = 0`) are exactly the least-cost ones. (A least-cost walk has
  reduced cost `dist(t*) − π(t*) = 0`, a sum of non-negative terms.)
- *Lemma 2: a Dinic phase depends only on `N_f` and the slot order.* If π and π′ both satisfy I and tight, the level
  build and the DFS behave identically: the forward BFS levels only vertices of `W`, where the test "positive, `c_π =
  0`, head marked" holds iff the arc is in `SP`, the same set for both potentials; the DFS tests arcs with both ends in
  `W`, where `c_π(u, v) = c(u, v) + dist(u) − dist(v) = c_π′(u, v)`. Residuals change identically push by push. After a
  phase I still holds (a push creates reverses of arcs with `c_π = 0`; `R` cannot grow), and the tier ends at the same
  phase for both.
- *Lemma 3: the tier update.* The full search sets `π′ = dist` on `R`. The early stop pops `t*` with key
  `D = dist(t*) − π(t*)`; every vertex with `d < D` was popped before, popped vertices have exact `d ≤ D`, and
  `π′ = π + min(d, D)` with `d = ∞` off `R`. For a positive arc `(x, y)` with `x ∈ R`, `d(y) ≤ d(x) + c_π(x, y)` and
  `min(·, D)` is monotone and 1-Lipschitz, so `c_π′ ≥ 0`; for `x ∉ R`, `x` gains `D` ≥ what `y` gains. At zero flow
  (π = 0) I* holds, so the early stop keeps I* and tightness (`π′(t*) = dist(t*)`, `π′(s*) = 0`) at every step.
- *Theorem.* By induction over tier updates, phases and pushes, the two solvers hold the same residual network at every
  step; the reached flag of each search, each tier's cost `dist(t*)` and units, each phase (Lemma 2), the stop at the
  cap and the final failing search (`t*` is never popped, so it never stops early) are functions of that network. Hence
  flow value, curve, edge flows, the minimum cut (closure of the sources in the final residual), termination and the
  phase and push counts coincide. ∎

**Capture and repair.** Edge sensitivity's repairs and outputs depend on the captured potentials only through their
validity (`c_ω ≥ 0` on every positive arc of the captured residual): each repair's searches, phases and pushes are then
functions of the residual network alone. So the capture may store the early-stopped potentials of the last successful
tier as they are, with no extra full search. The final, failing search must not apply the old "update only the reached
nodes" rule: that breaks validity on arcs from unreached into reached vertices (invalid in 10,188 of 2.07 M capture
cases on 3 nodes). Either skip that update or raise the unreached by the largest reached distance.

**Widths.** With `C` the cost sum (< 2^62): in a solve from zero both forms keep potentials in `[0, C]` and keys below
`2C`, so int64. Early-stop capture potentials are at most `C`. Repairs still need 128 bits: the repair bound shrinks
from (22n+15)C to (14n+15)C (n the vertex count) and still exceeds 2^63 for large nC.

**Near miss**, smallest: three nodes, edges `0 → 1` and `0 → 2` (capacity 1, cost 1), source 0, destinations {1, 2}.
Tier 1 pops `t*` at `D = 1` with node 1 unsettled. Leaving unsettled potentials unchanged gives `0 → 1` reduced cost 1
and `1 → t*` cost −1: tier 1 moves one unit and tier 2 reports the same cost again, curve `[(1, 1), (1, 1)]` instead
of `[(1, 2)]`. Checked: one analysis compared release and early-stop solvers on 10.78 M exhaustive and 63,200 random
cases, all three analyses proved it, and every planted variant was caught.

**Related exact rules in the same solver.**

- *Zero-flow certificate.* When every failed edge carries zero baseline flow and no endpoint fails, delivered amount and
  total cost equal the baseline's (unrestricted CentralizedTE or the same fixed cap, OPTIMIZED sources, fresh ledger),
  and only those. Termination: deleting an unused capacity-0 edge 0→1 changes `EXHAUSTED` to `NO_ROUTE` at zero
  delivery. Edge flows: 8 nodes, source 0, sink 7, edges 0→1, 0→2, 0→3, 1→4, 1→5, 1→6, 2→5, 2→6, 3→4, 4→7, 5→7, 6→7,
  each capacity 1 and cost 1; deleting 1→4, which carries no baseline flow, keeps value 3 and cost 9 but changes the
  edge flows (in that edge order) from [1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1] to [1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1].
  Curve: edges `0->1 (1, 1)`, `0->2 (1, 3)`, `1->2 (1, 1)`, `1->3 (1, 3)`, `2->3 (1, 1)` (capacity, cost), 0 to 3: curve
  `[(3,1), (5,1)]` against `[(4,2)]` without the unused `1->2`, same delivered 2 and cost 8.
- *Exact certificates instead.* The cut: one closure search over the baseline residual without the failed arcs,
  seeded by the case's roster (seeding with all sources was wrong in 46 of 1585 cases). The curve: equal when no failed
  edge carries flow at the end of any baseline tier (sufficient). 0 wrong in 79,172 cases of one analysis.
- *Cost-free deletion loss*: `loss = f − min(f, cost-free max flow u→v)` for `loss()` only; `resolve()` keeps its
  repairs.
- *Cut steps instead of bisection* for the largest prescribed total: equality proved, at most `|V| − |D| + 1` solves for
  prescribed allocations; no step bound is proved for T2's sizing.

## Pop order, FIFO queues, bit lanes and max-flow pair symmetry

**No public result of Core's searches depends on the order in which nodes at the same distance leave the queue.** With
positive costs and any min-key queue, a non-stale pop of `x` has the true distance, and at that moment `x`'s
predecessor list holds every tight edge into `x` exactly once; with targets, the search stops at `d*` with every tied
target final. Every reader of the lists sorts them or takes a set function before returning (DAG emission sorts each
list, `of_edges` sorts, `extract` sorts nodes and edges, a Collapse run is one scan's group). So distances, the DAG
arrays, the canonical path under every tie rule, KSP lists, ECMP and WCMP splits and remainders, route programs and
records, the cache, and errors are the same under every tie order, and the per-bucket NodeId sort can be deleted with
no replacement order. Checked: one oracle against the release build on 35,442,319 queries and seven pop orders of a
probe build, 0 differences.

- *Example.* Nodes 1 and 2 at distance 1 from root 0, each with a cost-1 edge to node 3. If node 2 leaves first, node
  3's list is (1→3, 2→3) instead of (2→3, 1→3); the returned DAG arrays are byte-identical because emission sorts.
- *FIFO* is exact when every edge of the graph has one cost (`max_edge_cost == edge_cost_gcd`, an O(1) test): the queue
  then holds keys within `{D, D + c}` in order, so it is a min-key queue. *Near miss*, FIFO with two costs on a targeted
  search: roots {0, 1}, edges 0→2 cost 2 and 1→3 cost 1, targets {2, 3}; FIFO pops 0, 1, then node 2 before node 3 and
  returns `target_cost` 2 and DAG [0→2]; the right answer is 1 and [1→3]. On full searches FIFO stays exact
  (label-correcting) but loses Dijkstra's bound; a check made only of full searches cannot catch the near miss.
- *Bit lanes* (one sweep for up to 64 roots or masks) are exact per lane, but each lane still needs its own emission
  of Core's lists: batching saves the sweep, not the output. Retiring a lane at the first target bit inside the bucket,
  rather than after the bucket, is the near miss.
- *The header comment.* `include/netgraph/core/search.hpp` says both queues "pop in (distance, node) order and give
  identical results". No public API returns a pop order and the user documentation says results never depend on the
  queue, but the first half becomes false under sort removal, FIFO or a search that pops by distance plus a potential.
  Any such patch rewrites that sentence first, for example to "Both pop the nearest pending distance first (the order
  among equal distances is unspecified) and give identical results".

### The reverse pair

For max-flow pairs on a graph where every edge `e = (u, v)` has a partner `p(e) = (v, u)` with the same capacity and
cost, the Eligibility excludes an edge exactly when it excludes its partner, the configuration is unrestricted
`CentralizedTE` with `OPTIMIZED` sources, and both directions get the same offer:

- **Value, `total_cost`, `marginal_cost_curve` and `termination` of (t, s) equal those of (s, t).** Reversing every
  edge and renaming it by `p` maps the graph onto itself; a feasible S → T flow `f` maps to the feasible T → S flow
  `f′(p(e)) = f(e)` of the same value and cost, and back. So the maximum value and the minimum cost `C(F)` of every
  flow value are the same both ways; the curve is `C`'s slope sequence (tiers in strictly increasing cost, each with the
  amount sent at that marginal cost), and termination follows from the value and the offer. ∎
- **Per-edge flows and `min_cut` do not mirror.** Path 0 − 1 − 2, capacities 1: (0, 2) reports `[0 -> 1]`, (2, 0)
  reports `[2 -> 1]`, while the mirror of `[0 -> 1]` is `[1 -> 0]`. Edge flows differ when the optimum is not unique: 5
  nodes, edges (tail, head, capacity, cost) (0, 4, 2, 2), (0, 2, 3, 3), (0, 3, 1, 3), (1, 2, 2, 2), (1, 4, 2, 3),
  (2, 1, 2, 2), (2, 0, 3, 3), (3, 0, 1, 3), (4, 0, 2, 2), (4, 1, 2, 3), pair (3, 1): value 1 and cost 8 both ways, but
  edge flows [1, 0, 0, 0, 0, 0, 0, 1, 0, 1] for (3, 1) against [0, 0, 1, 1, 0, 0, 1, 0, 0, 0] for (1, 3).
- **When the cut is requested and the solve reaches the maximum on its view, the reverse cut is derivable.** Core
  reports `min_cut` only then (`docs/semantics.md`, "`marginal_cost_curve` and `min_cut`"): if the offer stops the
  solve below that maximum, `min_cut` is `None` in both directions; an unrequested cut is also `None`; with an empty
  roster, a requested cut is `[]` (below). Witness: two nodes, edges 0→1 and 1→0, each capacity 10 and cost 1, offer 5,
  cut requested: both directions deliver 5 at cost 5 and report `None`, while the formula below would give `[]`. For a
  solve that reached the maximum, `min_cut(t, s)` is the set of partners of the edges leaving `S_max`, the largest
  source side of a minimum (s, t) cut (the nodes that cannot reach `t` in the final residual; that set is the same for
  every maximum flow). Under the mirror, a residual arc `a → b` of the (t, s) flow is a residual arc `b → a` of the
  mirrored (s, t) flow, so the set `X` that `t` reaches in the (t, s) residual is the set of nodes that reach `t` in the
  (s, t) residual, whose complement is `S_max`. An edge leaves `X` exactly when its partner enters `X`, that is, leaves
  `S_max`. So the reverse cut needs the (s, t) final residual and one reverse reachability pass, not its public
  outputs. ∎
- **Guard: an empty roster.** When no allowed walk joins the endpoints, Core returns `NO_ROUTE` and the cut `[]`. Nodes
  0, 1, 2, link 0–1 of capacity 0, node 2 isolated: Core returns `NO_ROUTE` and `[]` both ways; the bare closure gives
  `[0]`. With the guard, 5,400 of 5,400 random symmetric cases agree.

Checked on NSFNET: all 43,848 reverse cuts derived by the closure equal the separately solved reverse results, and the
ten non-cut exported fields agree in all 43,848 (s, t)/(t, s) comparisons; the plainly reversed cut differs in 20,103.
**The pair reduction is exact only for the eleven `PlacementResult` fields NetGraph's export reads** (`requested`,
`edge_drops`, `offered`, `delivered`, `unadmitted`, `source_drops`, `termination`, `total_cost`, `min_cut`,
`marginal_cost_curve`, `program`). Under Core's result contract it is not: the reverse solve's `by_tunnel` rows (always
returned) and, on request, its edge flows and program paths differ, and NetGraph's API keeps the whole
`PlacementResult`. So it is an opt-in call or output mode returning only those fields, a contract decision. Its
derived cut keeps Core's rule: the closure only for a solve that reached the maximum with the cut requested, `None`
otherwise, `[]` for an empty roster.

## Failure classes and baseline-derived programs

Within one workflow step NetGraph places the same demands under many failure masks, each removing elements from the
step's baseline mask. For a ShortestPaths demand the shortest-path structure of a failure case can be derived from the
baseline's without a search whenever every surviving source of the baseline roster keeps its distance; that test is
cheap and exact. For NetGraph's default ECMP options the route program is then the baseline program restricted to the
surviving elements. The placement is always recomputed on the case's own ledger.

**Definitions.** Positive integer costs. Baseline mask M0, case mask M with allowed(M) ⊆ allowed(M0) (deletions
only); F: the nodes and edges allowed in M0 and not in M (an edge is allowed only with both endpoints). T0, T′: the
destinations allowed in M0, M. `d`, `d′`: distances toward T0, T′ over allowed edges, so `d′ ≥ d`. Tight edges
(`d(u) = c + d(v)`) form Core's support. Roster R: allowed sources with `0 < d < ∞`. D_S: the nodes and edges of the
support reachable from R; every support out-edge of a D_S node is in D_S. *What a ShortestPaths call reads*:
construction keeps the support and the roster; every later read starts at the roster; preparation reads the service
view only for HIGHER_AVAILABLE/HEADROOM scores and RESIDUAL_SNAPSHOT weights, and admission reads capacity only on the
plan's forwarding edges. So program and placement are a function of (reachable support, roster, configuration, static
graph data, service view on the reachable support); `residual`, `saturated_boundary_edges` and the records'
`destinations` read more.

**Class 1** (no node or edge of D_S is in F): the case's support reachable from its roster is D_S and its roster is R.
A D_S node has a tight walk to T0 inside D_S that avoids F, so `d′ = d` on D_S; a tight edge from a D_S node stays
tight, and a non-tight one stays non-tight; deletions add no source. ∎

**Class 2.** Let A be the D_S nodes with no walk to T′ inside D_S minus F (failed nodes included). If no allowed
roster source is in A, the case roster is R ∩ allowed(M) and the case support reachable from it is P: the D_S edges not
in F whose head is not in A, reachable from that roster through such edges.

- (i) For an allowed D_S node x, `d′(x) = d(x)` iff x ∉ A. If x ∉ A, a surviving tight walk costs d(x). Conversely
  take a shortest case walk x = v0 … vk, vk ∈ T′, of cost `d′(x) = d(x)`; for each i,
  `d(x) ≤ prefix_i + d(v_i) ≤ prefix_i + d′(v_i) ≤ prefix_i + suffix_i = d(x)`, so `d(v_i) = suffix_i`: every edge is
  tight, the walk lies in D_S and avoids F.
- (ii) For x ∉ A: a case-tight `e = (x, y)` has `c + d′(y) = d(x) ≤ c + d(y) ≤ c + d′(y)`, so e is baseline-tight
  with `d′(y) = d(y)`, so y ∉ A; conversely an edge of D_S \ F into y ∉ A is case-tight.
- (iii) A failed source leaves the roster; (iv) induction along the walk from the case roster. ∎

(i) also shows that A is exactly the D_S nodes whose distance grew, plus failed ones, so "every allowed root keeps its
distance" and "every surviving source still reaches T′ inside D_S minus F" are the same test. Counting failed sources
as well is stricter, never wrong.

**The program.** With `multipath=True`, `parallel_edges=True` and EQUAL or static WEIGHTED, no route choice prunes the
DAG and EQUAL reads no residual, so the case's members are the baseline members restricted to P and pruned to walks
from the surviving roster to allowed destinations, with the same frozen weights; roster order, members per router in
EdgeId order and entry ids by roster position are unchanged. For other options the route choice and plan are rebuilt
on P, still without a search.

**The test.** Per demand keep `cnt(x)`, the number of D_S out-edges of x (every tight edge; parallel members are
separate edges). Mark every failed D_S node; kill every failed D_S edge. kill(e), once per edge: unless its tail is
marked, decrement `cnt(tail)`; at 0, mark the tail unless it is an allowed destination. mark(x): kill every D_S in-edge
of x. Class 2 iff no allowed roster source is marked. Correct: a node with a surviving walk keeps a live first edge
into an unmarked node (induction from the destination), so it is never marked; a node without one has every out-edge
failed or into a node without one (marked first, by induction on d since costs are positive), so its counter reaches
0. Cost O(|F ∩ D_S| + Σ over marked x of indeg(x)) plus restoring the touched counters; each node is queued at most
once. The counters are int32 (a degree is below 2^31); no distance and no sum is read, so no width can overflow. The
test reads only immutable inputs (graph, capture mask, case mask, endpoints, configuration), so it holds for any
capture lifetime.

**Near misses, as numbers a test can use.**

| Case | Wrong rule | Wrong answer | Right answer |
| --- | --- | --- | --- |
| nodes 0 and 1, two parallel edges 0→1 of cost 1 and capacity 10, `parallel_edges=False`; the baseline chooses edge 0, the case fails edge 0 | restrict the baseline program | delivers 0 | a fresh evaluation chooses edge 1 and delivers 10 |
| edges 0→1, 0→2, 1→3, 2→3 of cost 1 and 0→3 of cost 3, capacity 10 each; demand 0→3, offer 100; the case fails 0→1, and every distance to 3 stays [2, 1, 1, 0] | copy the baseline's placed amounts | 20 | 10 |
| edges 0→1, 0→2, 1→3, 2→3, each capacity 10 and cost 1; an ordered ECMP list places 0→3, then 1→3, each at offer 10; the case fails 0→2, which leaves the second demand's support unchanged | copy the second demand's baseline result | the baseline's 5 | 0 |
| edges 0→1, 1→2, capacity 10, cost 1; demand 0→2, offer 5; the case excludes node 0, so the roster is empty | replay the restricted program | `NO_ROUTE`, 0, `routes_changed` and `weights_changed` True | `NO_ROUTE`, 0, both flags False, as a fresh evaluation reports |
| class 3 (some source lost every tight path) | run the affected-set walk over D_S only | wrong in 24 % of one analysis's class-3 cases | run it over the whole baseline support, then repair distances locally |
| `multipath=False`: edges 0→1, 0→2, 1→3, 2→3, capacity 10, cost 1; demand 0→3, offer 100; the case fails 0→1 | restrict the baseline program | delivers 0 | the canonical path moves to 0→2→3 and delivers 10 |

Checked: one analysis ran 2.25 million failure cases (every four-node topology of its family with every single and
double failure, random five-node, fabric and backbone shapes, a base mask, 16 options) with seven planted variants
caught; the others 123,435 cases in five suites and 20,898 exhaustive plus 301 random. Not oracle-run by anyone:
simultaneous lists, OPTIMIZED end to end, multi-placement epochs. The opt-in prototype passed an oracle of 2,720,512
rows. On DC-BB every one of the 4,884 failure demand-cases is class 2, none class 1 or class 3; among the 1,942 cases in
which some surviving node's distance changed, the changed region has a median of 17 and a maximum of 256 nodes.

## Admission sweeps and ordered publication

LOSSLESS admission searches for the largest total T that fits (every member's share of T stays within its residual
view); Core evaluates one estimate, then doubles and bisects, each step a full sweep over the plan.

- **The search's base is sound.** Each member's share never falls when T grows (round robin and highest averages both
  give each extra unit to exactly one member), so the totals that fit are `0..T_max` and every exact search returns
  Core's total and per-edge split. Near miss: largest remainders (not house-monotone) in place of highest averages
  breaks monotonicity on 1,252 and the prefix on 58 of 3,000 weighted plans. A replacement divisor method must keep the
  same prefix order and the lower-index tie rule.
- **A window from the estimate.** A fixed window around Core's floating-point estimate is wrong: two parallel weighted
  members with weights (w, 1) and views (2w, 0) give estimate 0 and answer w (checked at w = 1, 2^50 + 123 and 2^60 +
  123). A window from the plan's structure, `|a_m(T) − T·φ_m| ≤ E_m` (a_m(T): member m's allocation at total T; φ_m: its
  exact fractional share per unit of T; E_m: the rounding error the splits upstream of m can add), is proved in exact
  arithmetic, but it can be as wide as the whole range; on DC-BB it is 19,000-76,000 units wide while Core's estimate
  misses by about 780, so it alone saves 0-20 % of the sweeps.
- **Several totals per sweep (lanes)** are exact if each lane has its own inflow and output (a shared-inflow plant is
  caught). Core's internal 62-probe budget must be restated: at extreme widths a plain 4- or 8-lane search evaluates 106
  or 154 totals where Core evaluates 62. The checked rule: choose q ≤ lanes with
  `u + q + bit_width(ceil(d/(q+1)) − 1) ≤ 62`, else q = 1, where u is the number of totals already evaluated and
  d = hi − lo the open interval; 252,656 budgeted searches agree with Core.
- **Clearing only the written slots** is exact if the written set includes a failed probe's partial writes, the failed
  target check's writes and what earlier calls left in the thread-local probe scratch. Built and timed, it lost
  ([candidates](candidates.md)).
- **Deriving a total's allocation** from a complete fitting one is exact as a local recomputation of the changed
  positions; adding two allocations is wrong (`RR(1) + RR(1) ≠ RR(2)`), and so is deriving from a failed probe.
- **T_max without search** is exact on forests and in-trees; every local formula tried fails on a merge of 2-7 nodes,
  and DC-BB's plans merge. Smallest witnesses: the estimate made exact, `floor(min_m v_m / φ_m)` (v_m: member m's
  view), on 2 nodes with parallel members 0→1 of views 1 and 0 gives 0 where T_max is 1; the forest pass on any DAG,
  3 nodes, roots 1 and 2, 1→2 and 2→3 of view 1, gives 2 where T_max is 1; the same on a split followed by a merge,
  4 nodes, 0→1, 0→2, 2→1, 1→3 of view 1, gives 2 where T_max is 1.
- **Ordered lists without publication** (run the ordered loop, skip publication) return every `PlacementResult` field
  and the same validation and computation errors in the same order, but the caller then sees no records or FlowIds, no
  installed `policy.program`, no epochs or ledger change, and no publication-stage failures: a new opt-in mode. On DC-BB
  even removing all identified publication work is 8.88 % of Core time, under the bar ([candidates](candidates.md)).

### The swap of totals

The [ordered-publication edit](ordered-publication.md) rests on this statement: after an ordered list's loop, the
working totals equal what publication would derive, so swapping them into the ledger after the last fallible step
keeps every result and state.

- *The working vector.* After the loop, `W = L0 + Σ_i (new_i − old_i)` over the actions that publish (`Place(0)` neither
  publishes nor exchanges), with `old_i` the ledger's stored entry of the action's policy and `new_i` its pending entry.
  The old path adds exactly `Σ (value − old)` over the same entries (the same objects, unchanged between the loop and
  the commit). A policy appears at most once in a list (`duplicate policy in list` refuses a repeat before any
  computation), so no entry is counted twice: `W = L0 + Δ` edge by edge, and the same for holds. By induction over the
  list this is the sum of the final entries, so the swapped ledger keeps `carried`, `reserved`, the derived `residual`
  and `unreserved`, and every later scalar, ordered or simultaneous call.
- *Checks and their order.* The new check pass visits every edge of every old and new aggregate with its final value,
  the same edge set and values the old path's first-touch scratch held, in first-touch order: old edges first, new-only
  edges after; a repeat changes no first failure.
- *Atomicity.* The swap comes after the last allocation and is `noexcept`; nothing after it throws (the remaining
  record, identity and epoch updates do not throw). A failure anywhere earlier leaves the ledger untouched and discards
  `W`. The Python result holders are created before the native transaction.
- *Widths.* No new sum: `W[e] − L[e]` and `L[e] + x` stay inside int64 because every per-edge total is in `[0, B_e]`,
  `B_e ≤ c_e < 2^62`; with `M = 2^62 − 1` both totals are in `[0, M]`, a difference in `[−M, M]`, any such sum in
  `[−M, 2M]`.
- *Edge cases.* An empty list returns before publication; `Place(0)`; zero-delivery records; `replace(0)`; T1 with
  every tunnel down; the first hold on a ledger with no reservations; a mode switch (the next simultaneous or scalar
  call reads the swapped ledger, equal to the old one). ∎

Removing the feasibility pass as well would keep every result but turn a computation bug into wrong ledger state, so one
check pass stays. Checked: 231,913 fixed call histories, 16,000 random batches and an independent generator of 4,000
histories (26,067 operations, 10,899 successful ordered lists, 4,341 states with holds) give byte-identical outputs;
its negative control (forgetting the released old flows) reports 610 errors on 200 histories.

## Reuse of the failure-free TE solve and tunnel setup

A failure case of a controller-placed workload places the same demands as the failure-free baseline on a graph with
some edges removed. The rule "if no removed edge was used by any step of the baseline's solve, the case's solve is the
baseline's" is false as stated and true with one more clause.

**Setting.** Core's solve: per tier one full Dijkstra from `s*` on reduced costs, then Dinic phases. For tier k let
`N_k` be the arcs with positive residual at its start, `D_k(x)` the cost of a cheapest `s* → x` path in `N_k`, `S_k`
the arcs of `N_k` on a cheapest `s* → t*` path; the recorded set SP is the edges with an arc in some `S_k`. The case
removes a set F of edges (edges of failed nodes included); the views of other edges are equal.

**Theorem (one demand).** If F ∩ SP = ∅, no source or destination fails and the roster is kept, the case runs the same
tiers, levels, augmenting paths and amounts. Proof by induction over tiers with "equal flows, zero flow on F": the
case's `N′_k` is `N_k` without F's arcs. (1) `N′_k ⊆ N_k` gives `D′ ≥ D`; `S_k ⊆ N′_k` gives `D′ = D` on every vertex
of a cheapest path, so `S′_k = S_k` and the tier costs are equal; when the baseline's search fails, the case's fails
too. (2) Within a tier, every arc that is admissible on an admissible `s* → t*` path at any time is in `S_k` or
reverses an arc pushed earlier in the tier, so the level graph avoids F and is the same in both runs. (3) An arc passes
the DFS test only if it lies on such a path, so pointers, paths and amounts agree. (4) Nothing is pushed on F. ∎

Outputs then equal by the theorem and the roster: value, curve, edge flows, path decomposition, records, program
paths and weights, `total_cost`, termination, flags, `by_source`, `by_tunnel`, `edge_flows`, `residual`,
`edge_reservations`, errors. Not equal: `min_cut` (recompute it by the closure certificate of the
[sink-stop section](#sink-stop-in-the-flow-solver)) and `saturated_boundary_edges` (one forward closure), which read
the case's mask; and a reused program must carry the case's mask as its replay basis, because two programs with every
public field equal replay differently later (`NO_ACTIVE_WEIGHT` against `NO_ROUTE`).

**The roster clause cannot be dropped.** One edge 0→1 of capacity 0, offer 5: no tier reaches `t*`, nothing is
recorded; failing the edge empties the roster and termination turns `EXHAUSTED` into `NO_ROUTE`.

**Prefix rule (ordered lists).** On a fresh ledger demand j solves on capacity minus the earlier demands' flows. If
demands 1..j−1 have the baseline's flows, demand j's view is the baseline's on every surviving edge and the theorem
applies with its own recorded set. So a demand is reusable while every earlier demand of the case was certified; the
rest of the list is solved fresh on the ledger holding the certified prefix's records (created in the same order:
same sequence numbers, policy ids, epochs) and published as usual, which keeps every field, record order, epoch and
error. FlowIds cannot be copied between case ledgers; ordinary publication gives the right ones. Testing a demand
alone is wrong: it fails on 148 of 1,082 `backbone_clos` demand-cases that pass their own test. A finer sound rule
also certifies a later demand whose view changed only outside its recorded set, provided no edge's view rises from 0.

**Tunnel setup (T1).** With the setup view equal on the surviving edges (a fresh ledger; in a list, every earlier T1
setup certified), the case's setup equals the baseline's (every tunnel's path, cost, destination, status, request, hold
and entry id; program, plan, records, holds) if (a) no head or destination fails, (b) the roster is kept, (c) no failed
edge is on an UP tunnel's path, and (d) every `COST_CEILING_EXCEEDED` tunnel's CSPF path survives and, under
`CostStretch` with such a tunnel, the reference cost c0 is unchanged. Proof by induction over attempts: (b) fixes the
heads, the requests and their order; an equal view gives the same filter and scores, so the case's walks are the
baseline's minus those through F (with `parallel_edges=False`, a pair whose kept member failed keeps its next member,
which never precedes the baseline's choice in the tie order), so the minimum is the baseline's path for UP and exceeded
tunnels, with the same cost and, by (d), the same verdict. ∎

| Clause dropped | Witness | Effect |
| --- | --- | --- |
| roster | one head behind a capacity-0 edge | roster empty, `EXHAUSTED` → `NO_ROUTE` |
| (d), ceiling-rejected path | `AbsoluteCost(3)`, the only filtered path costs 5, fail it | `COST_CEILING_EXCEEDED` turns into `INSUFFICIENT_BANDWIDTH` |
| (d), c0 | `CostStretch(2)`, c0 = 1 from a capacity-0 edge, fail it | c0 rises to 3 and the tunnel comes up (delivered 0 → 5) |
| earlier tunnels' holds | reuse one tunnel ignoring them | wrong in 18,552 of 1,752,362 cases |
| T2 under the same rule | one head, links A (10, cost 1) and B (10, cost 5), two tunnels, design 2: both tunnels on A in both runs | failing B lowers the plan's largest common size from 10 to 5 |

How often it applies (counts; 0 wrong against Core in every run):

| Input | Unit | Reusable by the prefix rule |
| --- | --- | --- |
| `backbone_clos`, 591 failure cases | 3,546 demand-cases | 934 (26.3 %); 140 cases fully |
| synthetic backbone, 49 failure cases | 6,125 demand-cases | 2,473 (40.4 %); 6 cases fully |
| synthetic backbone with tunnels, 49 failure cases | 6,125 T1 setups | 2,666 (43.5 %) under the T1 rule; 9 cases fully |

The certificates nest strictly: every case the trajectory rule admits, the curve certificate admits, and every case
that one admits, the zero-flow certificate admits; the zero-flow certificate certifies delivered and cost, the curve
certificate adds the curve, the trajectory rule everything but the two closures. Checked: one analysis ran 26,907,167
exhaustive and 33,774,662 random model cases and 541,543 random Core cases for the one-demand theorem, and 1,217,916
random Core cases for the T1 rule, 0 wrong.
