# Policy fixture reference

## In plain words

`tests/data/policy_fixtures.json` holds 59 small scenarios (graphs, ledgers, policies and their
calls) together with the result every call must produce. Both replay tests read it:
`tests/cpp/flow_policy_fixture_tests.cpp` runs the scenarios on the C++ engine and
`tests/py/test_flow_policy_fixtures.py` runs them through the Python binding; each compares
every result, trace, ledger state and policy state exactly.

The expected values were **not** computed by this library. They come from an independent
reference: two plain-Python prototypes written from the design contract before the C++ engine
existed (prototype 1, the amounts, search and kernels; prototype 2, the placement
transaction and the policy engine), with every fixture's expected values derived by hand in
the fixture text. This directory is a frozen copy of that reference and of the exporter that
ran it to write the JSON file, so a clean checkout can reproduce the data byte for byte without
the library and without the untracked design tree the code came from.

The data is frozen: it is never regenerated from the build. The command below proves that it
still equals what the reference computes; a change of the data is a deliberate data change,
reviewed as such. Two such changes have been made: when edges of capacity 0 became ordinary
routing edges ("Zero-capacity links stay in routing" below), and when the sources of one
demand were pooled and the origin arrived ("Pooled sources and the origin" below).

## The historical `source` string inside the data

The JSON's `source` field reads "prototype 2 of the design tree dev/design/2026-10-streamline
(untracked by design): prototypes/p2_transaction in integer units (APPROACH 9), exported by
gates/1c/export_fixture_data.py". It names the exporter that first wrote the file, under its
earlier name `tests/data/p2_fixtures.json`, in the design tree's records of sub-PR 1c. The
exporter copied here (`export_policy_fixtures.py`) replaced it with the same derivation and
the same bytes when the file was renamed, and deliberately kept the string, because it is part
of the content: changing it would change the bytes of the data. It is therefore kept
unchanged, and this document is where it is explained. "APPROACH 9" is the section of the design
contract on placement transactions; the contract is not part of the repository.

## Reproduce and check

From the repository root, with any Python 3.11 or later (only the standard library is used;
`-I` keeps the environment and the user's site out, `-S` additionally proves the build is not
needed, since without site-packages the library cannot be imported):

    venv/bin/python -I -B tests/reference/policy/export_policy_fixtures.py --check
    # prints "policy fixture data: identical" and exits 0; "DIFFERS" and exit 1 otherwise

    venv/bin/python -I -B -S tests/reference/policy/export_policy_fixtures.py --out /tmp/policy_fixtures.json
    cmp /tmp/policy_fixtures.json tests/data/policy_fixtures.json

Without `--out` the exporter writes `tests/data/policy_fixtures.json` relative to this
directory (three levels up). It runs in under a second. Before running the fixtures it checks
that every module of the reference was imported from this directory (an installed package of
the same name, such as `graph` or `units`, would otherwise shadow it) and that no module of the
library (`netgraph_core`, `_netgraph_core`) was imported.

The data file: 741,400 bytes, sha256
`e353ef06eb62f2ab3d60b33d586a6ccf46851c5d22f89dce67f5621e06859d26` (59 fixtures, 636
operations).

This directory is excluded from the source distribution (`sdist.exclude` in `pyproject.toml`);
the wheel holds only the `netgraph_core` package. No test imports it.

## Contents and origin

The copies are the complete set of modules the exporter imports, found by running it and
listing the imported modules (lazy imports inside functions included: every fixture runs).
Their bodies are unchanged except for the two rule changes at the end of this file; they
cite the design contract's sections and notes, which are not in the repository. `p2_transaction/_p1.py` appends the sibling `p1_amount/` to `sys.path`, so
the two directory names are part of the layout. `run_fixtures.py` is used only for its
`run_one`; its own `main` would write snapshots beside it and is not used.

Original paths are relative to the design tree `dev/design/2026-10-streamline/` (untracked);
the copies were taken on 2026-10-09.

| Copy (this directory) | Original | sha256 of the original | sha256 of the copy |
| --- | --- | --- | --- |
| `export_policy_fixtures.py` | `equivalence/export_policy_fixtures.py` | `0aea76a050b2d83e7398559af6a7d484a5a82a48b35bac9a3686f5fdf73211d2` | `43ff1a092de611e813507376080ac98d6eb19ee08348b72982ae1796a61cde4f` (relocated; pooled sources, below) |
| `p2_transaction/_p1.py` | `prototypes/p2_transaction/_p1.py` | `3ec5468e1ee8357faa9d85e6d748c718de63baee22e064c204324e42de42ec98` | same |
| `p2_transaction/fixtures.py` | `prototypes/p2_transaction/fixtures.py` | `ea05c85ce851fd7b7c6ff34b6a9eedbfd7b280187d0bf7e7965e626e4fe52071` | `64f616349a76887741879afac04c9d09ad82bd962b0c424ff557503c0451d8a6` (zero capacity; pooled sources, below) |
| `p2_transaction/policy.py` | `prototypes/p2_transaction/policy.py` | `b367b71c6052cb70c9b107ec7fdd15f491349cfd6fdfcb12a18f24ab69a4d27b` | `aeab8cc02c67309d1cb55dd00b926b98b04886b30f228be234df1f89fbe22bb7` (zero capacity; pooled sources, below) |
| `p2_transaction/recorder.py` | `prototypes/p2_transaction/recorder.py` | `fb4a8f20ccc6ecef4bc6d934226df33d84ee4dbf70c59aad9f33a097b9d0eb8a` | same |
| `p2_transaction/run_fixtures.py` | `prototypes/p2_transaction/run_fixtures.py` | `c20d655a2324c42201e9b743c51ef16bb47942dda56782e2d16b135a9c4db2fa` | same |
| `p1_amount/apportion.py` | `prototypes/p1_amount/apportion.py` | `c1e2f5082a99db1d0e7de5cf4e713ad5790aca3d85060583d536f5bc7fdcd4b1` | same |
| `p1_amount/certify.py` | `prototypes/p1_amount/certify.py` | `e4f0cce8a68fd512026e6e1c84c558b370d0741564b8191b8893de5cf173ad95` | `8d27afb711fd54f7063e75fe60137a5e37dbeea7dafb01d5788f6922d7f9c5b7` (zero capacity, below) |
| `p1_amount/flowcore.py` | `prototypes/p1_amount/flowcore.py` | `72b3d8c91a077017cf759d4c834576b1babd513024e9e227650d273d8b08da2c` | `63a6cf7047cb2fe91d5992f23e5d9bada1b89e70a4cff595be0a24ddabb681df` (zero capacity, below) |
| `p1_amount/graph.py` | `prototypes/p1_amount/graph.py` | `44051d4a496a423c73ed473b6d4da045c32e30b40701bec0ecb19b8d9e5822a8` | same |
| `p1_amount/kernels.py` | `prototypes/p1_amount/kernels.py` | `6e8db4d6de8c275d0335a9fddad28b7f84daa7af04d3ab872e4078cc4f709729` | `b01f6f24846e39679308fa80ea0c95804fcb653f0aea7c999662423dfab0afc5` (zero capacity, below) |
| `p1_amount/ledger.py` | `prototypes/p1_amount/ledger.py` | `19a8ab41f7ae26ece8b8a0f151c55a287f0cdf22a779d8c6d2f9d8e5dc43fcf6` | `fb55b11fa7c416024dd8c6faef63f08f831f77bf20585e78f3847b515aa3a6c2` (pooled sources, below) |
| `p1_amount/search.py` | `prototypes/p1_amount/search.py` | `3ab124b99fa68e25d467cfbb22a4e9900ade86e10001269df84f3c99ad31e04a` | `c79bfdc6d2a11197364f9bf7f3f41c8594843c72742f3285cbf685a4c117560c` (zero capacity, below) |
| `p1_amount/snap.py` | `prototypes/p1_amount/snap.py` | `0d62ff46d3c804bdbc80fae3c2c59a25aa0b34c013a086e7e81b7f722bdc0552` | same |
| `p1_amount/units.py` | `prototypes/p1_amount/units.py` | `38b37a22ae56f592a63554683e860b44036cce678a7080f406a69c281a620cfe` | same |

The rest of the two prototype directories (their notes, tests, snapshots, the MSD study
`msd.py` and the historical grid prototypes) is not part of the reference and was not copied.

## Relocation changes

At the relocation only the exporter changed, and only where it named the design tree's layout:

1. The usage line names the new path; a paragraph says the file is a relocated copy and points
   here (the rest of the docstring describes the original layout and is kept).
2. Paths: `DESIGN = HERE.parent` is gone; `REPO = HERE.parents[2]` (was `DESIGN.parents[2]`);
   `P2 = HERE / "p2_transaction"` (was `DESIGN / "prototypes" / "p2_transaction"`). `OUT` is
   unchanged in meaning (`tests/data/policy_fixtures.json` of the repository).
3. The origin check: it compared five modules with the prototypes directory
   (`origin.is_relative_to(DESIGN / "prototypes")`); it now compares every module of the
   closure (`CLOSURE`, all fourteen) with its exact file in this directory, and also fails if
   the library was imported.
4. The comment on removing the script's directory from `sys.path` says the directory held
   another `recorder.py` in the original layout; the removal is kept.

The output, including the `source` string, is unchanged: the check above compares the bytes.

## Zero-capacity links stay in routing

The library's rule changed: an edge of capacity 0 is no longer treated as missing. Searches,
static supports, rosters and Stretch baselines read no capacity (only an explicit residual
filters), an `EqualBalanced` route keeps such an edge as a member that receives its round-robin
share and carries nothing, a `Weighted` route gives it weight 0 (no member), a priority read
without a residual reads the capacity as the free capacity (`Headroom` 1, or 0 on an edge of
capacity 0), and a minimum cut lists every edge leaving the sources' residual closure, whatever its
free capacity. The reference was changed in exactly these places, by hand from that rule:

1. `p1_amount/search.py`: `member_graph` no longer skips edges below 1 unit; `priority_fn`
   reads the capacity as the residual when none is given, and `Headroom` is 0 on an edge of
   capacity 0; the module docstring says so.
2. `p1_amount/kernels.py`: `make_plan` takes every route edge as a member for `MaxFlow` and
   `EqualBalanced` (the `Weighted` rules are unchanged); the module docstring says so.
3. `p2_transaction/policy.py`: a `Pinned` bundle's usable route keeps the allowed edges whatever
   their capacity.
4. `p1_amount/flowcore.py` and `p1_amount/certify.py`: the minimum cut (and its certificate) is
   every active edge leaving the residual closure of the sources, an edge of view 0 too: it adds
   0 to the cut's value and is listed because every edge leaving the closure is. No fixture's
   data changes by this rule.
5. `p2_transaction/fixtures.py`: the two fixtures whose own hand-derived checks stated the old
   rule. `admission_zero_and_no_route`: `CostOnly` (MaxFlow) on the diamond whose four links have
   capacity 0 now has source 0 in its roster and, on `place(1)`, creates one zero-delivery record
   whose bundle is the whole static support (both routes cost 2), route cost 2, `exhausted`
   (previously an empty roster and `no_route`). `termination_reduction`: its empty-roster case
   now uses a graph whose only edge leads away from the destination (`1 -> 0`, capacity 1),
   since a zero-capacity edge would put the source in the roster; its three results are unchanged
   (only the ledger's residual reads that edge's 1 unit).

Every other fixture is unchanged: the exported data differs from the earlier file only in these
two fixtures' operations (`admission_zero_and_no_route` operations 6 and 7,
`termination_reduction` operations 20 to 22) and that graph.

## Pooled sources and the origin

The library's rule changed again: the sources of one demand are pooled and how traffic enters is
the demand's origin. A CostOnly policy with `EqualBalanced` or `Weighted` prepares, at its
epoch's first positive call, each roster source's route on the residual at that call's entry and
forwards over their union: a router several sources pass splits their traffic once, and each
member's amount (and, where a link drops, what it carries) is given back to the sources by
highest averages over what each holds; Lossless admits the largest total whose round robin over
the roster forwards with no drop and gives no share to a source without a member. CostOnly
`MaxFlow` and Global are controllers: one solve over the union of the routes (CostOnly) or every
allowed edge (Global), with the optimized origin (one virtual source) by default or the equal
origin (source i's arc capped at its round-robin share; the largest total that fits, by
bisection). Results report the demand's accounting and one row per input source instead of the
per-source outcomes. The reference was changed in exactly these places, by hand from that rule:

1. `p2_transaction/policy.py`: `FlowPolicyConfig` takes an optional placement and selection,
   `Outputs` a `min_cut` flag, and `PlacementResult` the new fields (`delivery_delta`,
   `call_drops`, `call_source_loss`, `requested`, `offered`, `admitted`, `delivered`,
   `unadmitted`, `source_loss`, epoch `drops`, `by_source`, `total_cost`) in place of
   `local_outcomes`; `ha` (highest averages with zero weights), `Program`, `make_program` and
   `forward` (the pooled program and its attribution); `_prepare`, `_source_flow`,
   `_place_pooled` and `_place_controller` (optimized origin only) replace the per-source
   CostOnly stage and `_place_global`; validation follows the new table (Global takes no placement or selection; the
   pooled single-path rule); the epoch also resets the offered total, the prepared program and the
   curve; `Pinned` reports `NO_ACTIVE_WEIGHT` and its Lossy offered total and source loss.
   The fixture reference has no origin selector: split CostOnly uses EQUAL, and controllers use
   OPTIMIZED. Equal-origin controllers are checked separately by `tests/py/origin_oracle.py`
   (`tests/py/test_origin.py`) and the native controller tests (`GlobalEqual.*`,
   `CostOnlyMaxFlow.*` in `tests/cpp/pooled_origin_tests.cpp`).
   CapacityAware keeps its per-source transaction and its reduction of per-source outcomes, now
   only in the trace.
2. `export_policy_fixtures.py`: a configuration's `selection` may be `None`, and the outputs are
   written as four flags (`min_cut` last); a record has no `placement` (below).
3. `p2_transaction/fixtures.py`: `glob()` builds Global without a placement or selection;
   `lossy_diamond` reads this call's drops as `call_drops`; `termination`'s empty roster reads
   its zero `by_source` row; `global_capped_underfill` requests `min_cut` by its flag. The
   multi-source fixtures whose own hand-derived checks stated the per-source rule were derived
   again by hand from the new rule: `sources_lossy_per_source`
   (router X splits 8 once; X->3 carries 3 = H(3, [4, 4]) = [2, 1]: deliveries 2 and 1, drops 2
   and 3), `sources_attribution` and `sources_inf_every_kernel` (Lossless: MaxFlow 3 under the
   optimized origin, Weighted and EqualBalanced 2 = [1, 1] under the equal origin; Lossy
   unchanged), `sources_eb_per_source` (8 = [4, 4], not 4 + 10), `sources_1_9_request_4` (MaxFlow
   [1, 3], 4, complete; Weighted and EqualBalanced [1, 1], 2, exhausted),
   `sources_partition_1_1_vs_2_0` ([2, 0] for both partitions), `sources_split_cumulative`
   (MaxFlow supplies [1, 0, 0] then [2, 0, 0]; Lossy offers [1, 0] then [0, 1]),
   `epochs_place_remove_place` (the supplies, read from the flow core's trace) and
   `sources_maxflow_vs_global` (its prose: one controller on the union of the shortest walks;
   value unchanged, 1). The other fixtures change only in the result's and the trace's shape.
4. A record no longer carries its kernel (`FlowRecord.placement` is gone from the library: the
   kernel is the policy's configuration). `p1_amount/ledger.py` drops the field from `FlowRecord`
   and `Draft` and from the record's JSON, `p2_transaction/policy.py` no longer passes it, and the
   exporter no longer writes it. The data lost exactly the 759 record `placement` keys; nothing
   else changed.
