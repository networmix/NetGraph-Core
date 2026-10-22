"""Prescribed source quotas, captured-weight replay and independent integer witnesses."""

import copy
import gc
import itertools
import json
import math
import random
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
from helpers import build, config, dag_of, exact, ledger_state, plain

import netgraph_core as n

CASES = json.loads(
    (Path(__file__).parents[1] / "data/weighted_sources.json").read_text()
)
M = 2**62 - 1
OUT = n.OutputOptions(
    edge_flows=True,
    residual=True,
    edge_reservations=True,
    saturated_boundary_edges=True,
    min_cut=True,
)


def weighted(weights, **kwargs):
    return config(
        n.CentralizedTE(
            source_allocation=n.SourceAllocation.WEIGHTED,
            source_weights=weights,
            **kwargs,
        )
    )


def ha(total, weights):
    """Independent exact rational priorities; no library allocation helpers."""
    if not sum(weights):
        return [0] * len(weights)
    a = [total * w // sum(weights) for w in weights]
    for _ in range(total - sum(a)):
        i = max(range(len(a)), key=lambda j: (Fraction(weights[j], a[j] + 1), -j))
        a[i] += 1
    return a


def sequential(total, weights):
    a = [0] * len(weights)
    for _ in range(total if any(weights) else 0):
        i = max(range(len(a)), key=lambda j: (Fraction(weights[j], a[j] + 1), -j))
        a[i] += 1
    return a


def subset(got, want):
    """Strict recursive schema/value comparison of the promised fixture fields."""
    if isinstance(want, dict):
        assert isinstance(got, dict)
        assert want.keys() <= got.keys()
        for key, value in want.items():
            subset(got[key], value)
    elif isinstance(want, list):
        assert type(got) is list and len(got) == len(want)
        for a, b in zip(got, want, strict=True):
            subset(a, b)
    else:
        assert type(got) is type(want) and got == want


def decode(value):
    if not isinstance(value, list):
        return value
    tag, *args = value
    if tag == "list":
        return [decode(x) for x in args]
    if tag == "range":
        return range(*args)
    if tag == "array":
        return np.array(args[1], dtype=args[0])
    if tag == "numpy":
        return getattr(np, args[0])(args[1])
    if tag in ("bool", "str"):
        return args[0]
    raise AssertionError(value)


def config_case(case):
    kw = dict(case["kwargs"])
    for name in ("source_weights", "edge_weights"):
        if name in kw:
            kw[name] = decode(kw[name])
    for name, enum in (
        ("source_allocation", n.SourceAllocation),
        ("allocation", n.RouteAllocation),
    ):
        if name in kw:
            kw[name] = getattr(enum, kw[name])
    if "tunnel_count" in kw:
        kw["tunnel_count"] = n.TunnelCount(int(kw["tunnel_count"]))
    try:
        model = getattr(n, case["model"])(**kw)
    except (TypeError, ValueError) as error:
        return f"{type(error).__name__}: {error}"
    g = build(case.get("num_nodes", 4), [(0, 3, 10, 1), (1, 3, 10, 1)])
    try:
        n.FlowPolicy(
            n.FlowLedger(g),
            [0, 1],
            3,
            config(
                model, admission=getattr(n.Admission, case.get("admission", "LOSSLESS"))
            ),
            cache=n.ShortestPathCache(g) if case.get("cache") else None,
        )
    except ValueError as error:
        return str(error)
    return (
        f"source_allocation={model.source_allocation.name} support=None tunnel_count=None"
        f" prefer_headroom=False source_weights={list(model.source_weights)}"
    )


def run_case(case):
    data = case["input"]
    if case["oracle"] == "model_config":
        return [config_case(c) for c in data["cases"]]
    if case["oracle"] == "weighted_equal":
        g = build(4, data["edges"])
        out = {}
        for name, cfg in (
            ("weighted", weighted([data["weight"], data["weight"], 0, 0])),
            (
                "equal",
                config(n.CentralizedTE(source_allocation=n.SourceAllocation.EQUAL)),
            ),
        ):
            r = n.evaluate_demand(
                g,
                data["sources"],
                data["destinations"],
                n.EvaluationOptions(cfg, offer=data["offer"], outputs=OUT),
            )
            out[name] = dict(
                delivered=r.delivered,
                cost=r.total_cost,
                edge_flows=r.edge_flows.tolist(),
                by_source=[[s.source, s.delivered] for s in r.by_source],
            )
        return out
    g = build(len(data["source_weights"]), data["edges"])
    mask = g.eligibility(
        excluded_nodes=data.get("excluded_nodes", []),
        excluded_edges=data.get("excluded_edges", []),
    )
    cfg = weighted(data["source_weights"])
    ledger = n.FlowLedger(g)
    p = n.FlowPolicy(
        ledger, data["sources"], data["destinations"], cfg, eligibility=mask
    )
    offer = data.get("offer", math.inf)
    before = ledger_state(ledger)
    evaluated = p.evaluate(offer, outputs=OUT)
    assert ledger_state(ledger) == before and p.program is None
    r = p.replace(offer, outputs=OUT)
    if r.delivered:
        assert r.routes_changed and r.weights_changed and not r.reservations_changed
    assert plain(evaluated) == plain(r)
    assert r.plan is None and r.min_cut is None and r.marginal_cost_curve is None
    assert r.offered is None and r.admitted == r.delivered
    assert r.source_drops == 0 and r.edge_drops == []
    assert r.edge_reservations.tolist() == [0] * g.num_edges()
    assert sum(s.delivered for s in r.by_source) == r.delivered
    assert r.total_cost == sum(
        int(c) * int(f) for c, f in zip(g.edge_costs(), r.edge_flows, strict=True)
    )
    assert r.requested == offer
    assert r.unadmitted == (None if offer == math.inf else offer - r.delivered)
    assert ledger.carried().tolist() == r.edge_flows.tolist()
    assert ledger.residual().tolist() == [
        int(c) - int(f)
        for c, f in zip(g.edge_capacity_units(), r.edge_flows, strict=True)
    ]
    captured = r.program.source_weights
    out = dict(
        roster=[s for s, _ in captured],
        program_source_weights=[list(x) for x in captured],
        delivered=r.delivered,
        termination=r.termination.name,
        by_source=[[s.source, s.delivered] for s in r.by_source],
        cost=r.total_cost,
        edge_flows=r.edge_flows.tolist(),
        tunnels=[
            [t.entry_id, t.head, list(t.edges), t.weight] for t in r.program.tunnels
        ],
    )
    if "replays" in data:
        out["replays"] = []
        for x in data["replays"]:
            runs = {}
            for admission in (n.Admission.LOSSY, n.Admission.LOSSLESS):
                q = n.FlowPolicy(
                    n.FlowLedger(g),
                    data["sources"],
                    data["destinations"],
                    config(n.Pinned(r.program), admission=admission),
                    eligibility=mask,
                )
                runs[admission] = q.replace(x, outputs=OUT)
            a, b = runs[n.Admission.LOSSY], runs[n.Admission.LOSSLESS]
            rows = {s.source: s for s in a.by_source}
            out["replays"].append(
                dict(
                    offer=x,
                    edge_flows=a.edge_flows.tolist(),
                    source_drops=[rows[s].source_drops for s, _ in captured],
                    delivered=[rows[s].delivered for s, _ in captured],
                    demand_source_drops=a.source_drops,
                    lossy_termination=a.termination.name,
                    lossless_termination=b.termination.name,
                )
            )
        q = n.FlowPolicy(
            n.FlowLedger(g),
            data["sources"],
            data["destinations"],
            config(n.Pinned(r.program)),
            eligibility=mask,
        )
        out["lossless_bound"] = q.replace(math.inf).delivered
    return out


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_contract_cases(case):
    subset(run_case(case), case["expected"])


def test_oracles_and_comparator_reject_wrong_rules():
    assert sequential(2, [2, 1]) == [2, 0]
    assert sequential(50, [100, 1, 1, 1]) == [50, 0, 0, 0]
    for t in range(25):
        for w in itertools.product(range(4), repeat=3):
            assert ha(t, w) == sequential(t, w)
    assert ha(15, [2, 1]) == [10, 5] != ha(15, [11, 5])

    # Every promised fixture leaf: wrong value, wrong type, missing and null must fail.
    def mutants(value):
        yield None
        if isinstance(value, dict):
            for key in value:
                bad = copy.deepcopy(value)
                del bad[key]
                yield bad
                for replacement in mutants(value[key]):
                    bad = copy.deepcopy(value)
                    bad[key] = replacement
                    yield bad
        elif isinstance(value, list):
            yield value + [None]
            for i, leaf in enumerate(value):
                for replacement in mutants(leaf):
                    bad = copy.deepcopy(value)
                    bad[i] = replacement
                    yield bad
        elif type(value) is int:
            yield value + 1
            yield float(value)
        elif isinstance(value, str):
            yield value + "!"

    for case in CASES:
        expected = case["expected"]
        subset(expected, expected)
        for changed in mutants(expected):
            with pytest.raises((AssertionError, TypeError)):
                subset(changed, expected)


@pytest.mark.parametrize("ordered", [False, True])
def test_refused_traffic_is_not_retried_after_relief(ordered):
    g = build(4, [(0, 2, 10, 1), (1, 2, 10, 1), (2, 3, 10, 1)])
    ledger = n.FlowLedger(g)
    other = n.FlowPolicy(ledger, 2, 3, config(n.ShortestPaths()))
    p = n.FlowPolicy(ledger, [0, 1], 3, weighted([2, 1, 0, 0]))
    other.replace(8)

    def call(v):
        return (
            n.apply_actions(ledger, [n.Place(p, v)], outputs=OUT)[0]
            if ordered
            else p.place(v, outputs=OUT)
        )

    initial = call(5)
    assert initial.delivered == 2
    check_live_rows(initial, [[0, 0, "UP", 2, 2, 0]], 3)
    other.replace(5)
    r = call(1)
    assert (r.call_delivered, r.delivered, r.total_cost, r.requested, r.unadmitted) == (
        1,
        3,
        6,
        6,
        3,
    )
    assert [s.delivered for s in r.by_source] == [2, 1]
    assert r.edge_flows.tolist() == [2, 1, 3]
    check_live_rows(r, [[0, 0, "UP", 2, 2, 0], [1, 1, "UP", 1, 1, 0]], 0)


def test_exact_replay_paths_and_all_zero_termination():
    g = build(3, [(0, 2, 10, 1), (1, 2, 10, 1)])
    live = n.FlowPolicy(n.FlowLedger(g), [0, 1], 2, weighted([0, 0, 0]))
    saved = live.replace(6)
    assert saved.termination == n.Termination.NO_ACTIVE_WEIGHT
    for support in (None, dag_of(g, [0, 1])):
        q = n.FlowPolicy(
            n.FlowLedger(g), [0, 1], 2, weighted([0, 0, 0], support=support)
        )
        unlimited = q.evaluate(math.inf)
        assert unlimited.termination == n.Termination.NO_ACTIVE_WEIGHT
        assert unlimited.requested == unlimited.call_undelivered == math.inf
        assert unlimited.unadmitted is None and unlimited.delivered == 0
    assert live.forward(6).termination == n.Termination.EXHAUSTED
    for admission in (n.Admission.LOSSY, n.Admission.LOSSLESS):
        q = n.FlowPolicy(
            n.FlowLedger(g),
            [0, 1],
            2,
            config(n.Pinned(saved.program), admission=admission),
        )
        with pytest.raises(ValueError, match=exact("Pinned forwards with replace")):
            q.forward(6)
        for method in (q.replace, q.evaluate, q.place):
            q.remove()
            r = method(6)
            assert r.termination == (
                n.Termination.COMPLETE
                if admission == n.Admission.LOSSY
                else n.Termination.EXHAUSTED
            )
            assert r.source_drops == (6 if admission == n.Admission.LOSSY else 0)
            assert all(s.source_drops == s.delivered == 0 for s in r.by_source)
            assert method(0).termination == n.Termination.COMPLETE


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("reverse_creation", [False, True])
def test_simultaneous_ties_follow_creation_not_action_order(reverse, reverse_creation):
    g = build(2, [(0, 1, 3, 1)])
    p = n.FlowPolicy(n.FlowLedger(g), 0, 1, weighted([1, 0]))
    program = p.replace(3).program
    ledger = n.FlowLedger(g)
    first = n.FlowPolicy(
        ledger, 0, 1, config(n.Pinned(program), admission=n.Admission.LOSSY)
    )
    second = n.FlowPolicy(
        ledger, 0, 1, config(n.Pinned(program), admission=n.Admission.LOSSY)
    )
    a, b = (second, first) if reverse_creation else (first, second)
    owners = [b, a] if reverse else [a, b]
    r = n.apply_actions(
        ledger, [n.Replace(q, 3) for q in owners], mode=n.BatchMode.SIMULTANEOUS
    )
    assert [x.delivered for x in r] == [2 if q is first else 1 for q in owners]
    assert first.delivered() == 2 and second.delivered() == 1
    first.remove()
    second.remove()
    r = n.apply_actions(ledger, [n.Replace(q, 3) for q in owners])
    assert [x.delivered for x in r] == [3, 0]


@pytest.mark.parametrize("count", [1, 2, 3, 4, 5, 16, 17, 40])
@pytest.mark.parametrize("unequal", [False, True])
def test_cumulative_captured_weights_cross_sixty_four_bits(count, unequal):
    weights = [M - i if unequal else M for i in range(count)]
    g = build(count + 1, [(i, count, 5, 1) for i in range(count)])
    p = n.FlowPolicy(
        n.FlowLedger(g), list(range(count)), count, weighted(weights + [0])
    )
    program = p.replace(count).program
    q = n.FlowPolicy(
        n.FlowLedger(g),
        list(range(count)),
        count,
        config(n.Pinned(program), admission=n.Admission.LOSSY),
    )
    total = 0
    for x in [M, M, M, 7, M, 5, 17, 0, 1]:
        try:
            r = q.place(x)
        except RuntimeError as error:
            raise AssertionError("a valid cumulative offer must succeed") from error
        total += x
        assert [s.offered for s in r.by_source] == ha(total, weights)
        assert r.offered == r.requested == total
        assert r.delivered + r.source_drops + sum(v for _, v in r.edge_drops) == total


def test_owned_weights_and_all_evaluator_readers():
    data = next(c["input"] for c in CASES if c["id"] == "W-witness")
    g = build(6, data["edges"])
    weights = np.array(data["source_weights"], dtype=np.int64)
    cfg = weighted(weights)
    weights[:] = 0
    assert cfg.routing.source_weights == (2, 1, 0, 0, 0, 0)
    with pytest.raises(AttributeError):
        cfg.routing.source_weights = (0,) * 6
    opts = n.EvaluationOptions(cfg, offer=15, outputs=OUT)
    want = n.evaluate_demand(g, [1, 0, 1], 5, opts)
    assert [s.delivered for s in want.by_source] == [10, 5]
    for workers in (1, 2, 4):
        results = n.evaluate_pairs(
            g, [([0, 1], [5]), ([1, 0, 1], [5])], opts, workers=workers
        )
        assert all(plain(r) == plain(want) for r in results)
        failures = n.evaluate_edge_failures(g, [0, 1], 5, opts, workers=workers)
        for entry in failures.entries:
            deleted = g.eligibility(excluded_edges=[entry.edge])
            r = n.evaluate_demand(
                g, [0, 1], 5, n.EvaluationOptions(cfg, eligibility=deleted, offer=15)
            )
            assert entry.delivery_reduction == want.delivered - r.delivered
    program = want.program
    del opts, want, weights, cfg, results, failures
    gc.collect()
    q = n.FlowPolicy(n.FlowLedger(g), [0, 1], 5, config(n.Pinned(program)))
    assert [s.delivered for s in q.replace(15).by_source] == [10, 5]


def test_edge_failures_recompute_unused_zero_capacity_roster_edge():
    g = build(3, [(0, 2, 0, 1), (1, 2, 10, 1)])
    opts = n.EvaluationOptions(weighted([2, 1, 0]))
    r = n.evaluate_edge_failures(g, [0, 1], 2, opts, workers=2)
    assert r.baseline.delivered == 0
    assert [(e.edge, e.delivery_reduction) for e in r.entries] == [(0, -10), (1, 0)]


def test_small_integer_objective_with_masks_and_support():
    rng = random.Random(20261022)
    for _ in range(45):
        edges = [
            (0, 2, rng.randrange(3), 1),
            (0, 3, rng.randrange(3), 7),
            (1, 2, rng.randrange(3), 3),
            (1, 3, rng.randrange(3), 1),
            (2, 3, rng.randrange(3), 1),
        ]
        g = build(4, edges)
        weights = [rng.randrange(4), rng.randrange(4)]
        allowed = [rng.randrange(4) != 0 for _ in edges]
        # Capacity-blind reachability; support is independently selected from the DAG.
        roster = [
            s
            for s in (0, 1)
            if any(
                allowed[e]
                and a == s
                and (
                    b == 3
                    or any(
                        allowed[j] and c == b and d == 3
                        for j, (c, d, *_x) in enumerate(edges)
                    )
                )
                for e, (a, b, *_x) in enumerate(edges)
            )
        ]
        effective = [weights[s] for s in roster]
        best = {0: 0}
        for f in itertools.product(
            *(
                range(c + 1) if allowed[i] else [0]
                for i, (_a, _b, c, _k) in enumerate(edges)
            )
        ):
            bal = [0] * 4
            for x, (a, b, _c, _k) in zip(f, edges, strict=True):
                bal[a] += x
                bal[b] -= x
            t = bal[0] + bal[1]
            shares = dict(zip(roster, sequential(t, effective), strict=True))
            if bal[2] or any(bal[s] != shares.get(s, 0) for s in (0, 1)):
                continue
            cost = sum(x * edge[3] for x, edge in zip(f, edges, strict=True))
            best[t] = min(best.get(t, cost), cost)
        top = max(best)
        assert sorted(best) == list(range(top + 1))
        # Map input positions to sorted EdgeId for support/masks and exact cost.
        support = dag_of(
            g, [e for e, pos in enumerate(g.input_positions()) if allowed[int(pos)]]
        )
        cfg = weighted(weights + [0, 0], support=support)
        for offer in (0, 1, 5, math.inf):
            r = n.evaluate_demand(
                g,
                [0, 1],
                3,
                n.EvaluationOptions(
                    cfg, offer=None if offer == math.inf else offer, outputs=OUT
                ),
            )
            t = top if offer == math.inf else min(top, offer)
            assert (r.delivered, r.total_cost) == (t, best[t])
            assert [s.delivered for s in r.by_source] == [
                dict(zip(roster, ha(t, effective), strict=True)).get(s, 0)
                for s in (0, 1)
            ]


def test_source_vector_and_model_refusals_before_zero_offer():
    g = build(3, [(0, 2, 2, 1), (1, 2, 2, 1)])
    for model, message in (
        (
            n.CentralizedTE(
                source_allocation=n.SourceAllocation.WEIGHTED,
                source_weights=[1, 0, 0],
                tunnel_count=n.TunnelCount(1),
            ),
            "source_allocation=WEIGHTED requires CentralizedTE without tunnel_count",
        ),
        (
            n.ShortestPaths(source_allocation=n.SourceAllocation.WEIGHTED),
            "source_allocation=WEIGHTED requires CentralizedTE without tunnel_count",
        ),
        (
            n.CentralizedTE(source_allocation=n.SourceAllocation.WEIGHTED),
            "source_allocation=WEIGHTED requires source_weights",
        ),
    ):
        with pytest.raises(ValueError, match=exact(message)):
            n.evaluate_demand(g, 0, 2, n.EvaluationOptions(config(model), offer=0))
    for weights, message in (
        ([1, 0], "source_weights: expected 3 entries"),
        ([-1, 0], "source_weights: expected 3 entries"),
        ([0, M + 1], "source_weights: expected 3 entries"),
        ([1, 0, 0, M + 1], "source_weights: expected 3 entries"),
        ([1, 0, -1], "source_weights[2] must be an integer in [0, 2^62)"),
        ([1, 0, M + 1], "source_weights[2] must be an integer in [0, 2^62)"),
    ):
        with pytest.raises(ValueError, match=exact(message)):
            n.evaluate_demand(g, 0, 2, n.EvaluationOptions(weighted(weights), offer=0))
        with pytest.raises(ValueError, match=exact(message)):
            n.FlowPolicy(n.FlowLedger(g), 0, 2, weighted(weights))
    assert (
        n.evaluate_demand(
            g, 0, 2, n.EvaluationOptions(weighted([M, 0, 0]), offer=0)
        ).delivered
        == 0
    )


def test_support_black_holes_zero_calls_and_list_refusals():
    g = build(3, [(0, 2, 10, 1), (1, 2, 10, 1)])
    ledger = n.FlowLedger(g)
    p = n.FlowPolicy(ledger, [0, 1], 2, weighted([2, 1, 0]))
    zero = p.place(0)
    assert zero.program is None and zero.termination == n.Termination.COMPLETE
    assert not (
        zero.routes_changed or zero.weights_changed or zero.reservations_changed
    )
    program = p.replace(3).program
    mask = g.eligibility(excluded_nodes=[0])
    repaired = n.FlowPolicy(
        n.FlowLedger(g),
        [0, 1],
        2,
        config(n.Pinned(program), admission=n.Admission.LOSSY),
        eligibility=mask,
    )
    stale = n.FlowPolicy(
        n.FlowLedger(g),
        [0, 1],
        2,
        config(n.Pinned(program, black_holes=mask), admission=n.Admission.LOSSY),
    )
    a, b = repaired.replace(3), stale.replace(3)
    assert [r.delivered for r in a.by_source] == [0, 3]
    assert [r.delivered for r in b.by_source] == [0, 1]
    assert b.source_drops + sum(x for _, x in b.edge_drops) == 2
    assert program.source_weights == ((0, 2), (1, 1))
    saved = ledger_state(ledger)
    with pytest.raises(
        ValueError,
        match=exact(
            "SIMULTANEOUS Replace requires ShortestPaths EQUAL or WEIGHTED, or Pinned"
        ),
    ):
        n.apply_actions(ledger, [n.Replace(p, 3)], mode=n.BatchMode.SIMULTANEOUS)
    with pytest.raises(ValueError, match=exact("SIMULTANEOUS does not accept Place")):
        n.apply_actions(ledger, [n.Place(p, 0)], mode=n.BatchMode.SIMULTANEOUS)
    assert ledger_state(ledger) == saved
    r = p.replace(0)
    assert r.termination == n.Termination.COMPLETE and r.program.source_weights == (
        (0, 2),
        (1, 1),
    )
    supported = weighted([2, 1, 0], support=dag_of(g, [1]))
    r = n.evaluate_demand(g, [0, 1], 2, n.EvaluationOptions(supported, outputs=OUT))
    assert r.edge_flows.tolist() == [0, 10] and [s.delivered for s in r.by_source] == [
        0,
        10,
    ]


def test_wide_cost_and_zero_effect_participants_across_program_forms():
    g = build(2, [(0, 1, M, M)])
    r = n.evaluate_demand(g, 0, 1, n.EvaluationOptions(weighted([M, 0]), outputs=OUT))
    assert type(r.total_cost) is int and r.total_cost == M * M
    data = CASES[0]["input"]
    g = build(6, data["edges"])
    empty = n.evaluate_demand(
        g, [0, 1], 5, n.EvaluationOptions(weighted([0] * 6))
    ).program
    positive = n.evaluate_demand(
        g, [0, 1], 5, n.EvaluationOptions(weighted([2, 1, 0, 0, 0, 0]))
    ).program
    dag = dag_of(g, range(g.num_edges()))
    models = [
        n.Pinned(positive),
        n.ShortestPaths(),
        n.Pinned(n.ForwardingProgram.from_dag(g, [0, 1], 5, dag)),
        n.Pinned(
            n.ForwardingProgram.from_routes(
                g, [0, 1], 5, [n.RouteDefinition(0, dag), n.RouteDefinition(1, dag)]
            )
        ),
    ]
    for model in models:
        for offer in (0, 7):
            ledger = n.FlowLedger(g)
            p = n.FlowPolicy(
                ledger, [0, 1], 5, config(model, admission=n.Admission.LOSSY)
            )
            zero = n.FlowPolicy(
                ledger, [0, 1], 5, config(n.Pinned(empty), admission=n.Admission.LOSSY)
            )
            base = n.apply_actions(
                ledger, [n.Replace(p, 19)], mode=n.BatchMode.SIMULTANEOUS, outputs=OUT
            )[0]
            p.remove()
            got = n.apply_actions(
                ledger,
                [n.Replace(zero, offer), n.Replace(p, 19)],
                mode=n.BatchMode.SIMULTANEOUS,
                outputs=OUT,
            )[1]
            assert plain(got) == plain(base)


def test_transit_through_another_source_and_installed_flags():
    g = build(3, [(0, 1, 5, 1), (1, 2, 10, 1)])
    ledger = n.FlowLedger(g)
    p = n.FlowPolicy(ledger, [0, 1], 2, weighted([1, 1, 0]))
    first = p.replace(math.inf, outputs=OUT)
    assert first.delivered == 10 and first.total_cost == 15
    assert [r.delivered for r in first.by_source] == [5, 5]
    assert first.edge_flows.tolist() == [5, 10]
    assert [r.delivered for r in p.records()] == [5, 5]
    assert (
        first.routes_changed
        and first.weights_changed
        and not first.reservations_changed
    )
    r = p.forward(8)
    assert r.delivered == 8
    assert not (r.routes_changed or r.weights_changed or r.reservations_changed)
    zero = p.place(0)
    assert not (
        zero.routes_changed or zero.weights_changed or zero.reservations_changed
    )


def check_live_rows(result, rows, undelivered):
    actual = dict(
        by_tunnel=[
            [
                r.entry_id,
                r.head,
                r.status.name,
                r.assigned,
                r.delivered,
                r.total_edge_drops,
            ]
            for r in result.by_tunnel
        ],
        call_undelivered=result.call_undelivered,
    )
    subset(actual, dict(by_tunnel=rows, call_undelivered=undelivered))


@pytest.mark.parametrize("offer", [3, 7])
def test_live_tunnel_rows_and_undelivered_across_readers(offer):
    g = build(2, [(0, 1, 5, 1)])
    cfg = weighted([1, 0])
    for operation in ("place", "replace", "evaluate"):
        for ordered in (False, True):
            ledger = n.FlowLedger(g)
            p = n.FlowPolicy(ledger, 0, 1, cfg)
            if ordered and operation != "evaluate":
                action = (
                    n.Place(p, offer) if operation == "place" else n.Replace(p, offer)
                )
                r = n.apply_actions(ledger, [action], outputs=OUT)[0]
            else:
                r = getattr(p, operation)(offer, outputs=OUT)
            check_single_edge_result(r, offer)
    options = n.EvaluationOptions(cfg, offer=offer, outputs=OUT)
    check_single_edge_result(n.evaluate_demand(g, 0, 1, options), offer)
    for workers in (1, 2, 4):
        for r in n.evaluate_pairs(g, [([0], [1])] * 5, options, workers=workers):
            check_single_edge_result(r, offer)
        check_single_edge_result(
            n.evaluate_edge_failures(g, 0, 1, options, workers=workers).baseline, offer
        )


def test_live_tunnel_row_schema_rejects_missing_null_and_wrong_distribution():
    expected = dict(
        by_tunnel=[[0, 0, "UP", 3, 3, 0], [1, 0, "UP", 8, 8, 0], [2, 1, "UP", 5, 5, 0]],
        call_undelivered=0,
    )
    subset(expected, expected)
    bad_values = [
        dict(),
        dict(by_tunnel=None, call_undelivered=0),
        dict(by_tunnel=[], call_undelivered=0),
    ]
    for key in expected:
        bad = copy.deepcopy(expected)
        del bad[key]
        bad_values.append(bad)
        bad = copy.deepcopy(expected)
        bad[key] = None
        bad_values.append(bad)
    bad = copy.deepcopy(expected)
    bad["call_undelivered"] = 1
    bad_values.append(bad)
    for i in range(3):
        for j in range(6):
            for replacement in (
                None,
                "FAILED" if j == 2 else expected["by_tunnel"][i][j] + 1,
            ):
                bad = copy.deepcopy(expected)
                bad["by_tunnel"][i][j] = replacement
                bad_values.append(bad)
    for j in (3, 4, 5):
        bad = copy.deepcopy(expected)
        bad["by_tunnel"][0][j] += 1
        bad["by_tunnel"][1][j] -= 1
        bad_values.append(bad)
    for bad in bad_values:
        with pytest.raises((AssertionError, TypeError)):
            subset(bad, expected)
    data = next(c["input"] for c in CASES if c["id"] == "W-witness")
    g = build(6, data["edges"])
    r = n.evaluate_demand(
        g,
        [0, 1],
        5,
        n.EvaluationOptions(weighted(data["source_weights"]), offer=16, outputs=OUT),
    )
    check_live_rows(r, expected["by_tunnel"], 0)


@pytest.mark.parametrize("admission", [n.Admission.LOSSY, n.Admission.LOSSLESS])
def test_remaining_zero_captured_weights_with_a_program_path(admission):
    g = build(3, [(0, 2, 10, 1), (1, 2, 10, 1)])
    program = n.evaluate_demand(
        g, [0, 1], 2, n.EvaluationOptions(weighted([2, 0, 0]), offer=5)
    ).program
    assert len(program.tunnels) == 1 and program.source_weights == ((0, 2), (1, 0))
    mask = g.eligibility(excluded_nodes=[0])
    p = n.FlowPolicy(
        n.FlowLedger(g),
        [0, 1],
        2,
        config(n.Pinned(program), admission=admission),
        eligibility=mask,
    )
    for method in (p.replace, p.evaluate, p.place):
        p.remove()
        r = method(4, outputs=OUT)
        assert r.termination == n.Termination.NO_ACTIVE_WEIGHT
        assert r.delivered == 0 and r.call_undelivered == 4
        assert r.source_drops == (4 if admission == n.Admission.LOSSY else 0)
        assert r.unadmitted == (0 if admission == n.Admission.LOSSY else 4)
        assert all(s.delivered == s.source_drops == 0 for s in r.by_source)


def check_single_edge_result(r, offer):
    carried = min(offer, 5)
    check_live_rows(r, [[0, 0, "UP", carried, carried, 0]], offer - carried)
    assert r.delivery_delta == r.call_delivered == r.delivered == carried
    assert r.call_edge_drops == r.edge_drops == []
    assert r.call_source_drops == r.source_drops == 0
    assert r.requested == offer and r.admitted == carried
    assert r.unadmitted == offer - carried and r.offered is None
    assert r.termination == (
        n.Termination.COMPLETE if offer <= 5 else n.Termination.EXHAUSTED
    )
    assert r.routes_changed and r.weights_changed and not r.reservations_changed
    assert r.program is not None and r.program.source_weights == ((0, 1),)
    assert [
        (t.entry_id, t.head, list(t.edges), t.weight, t.path_cost)
        for t in r.program.tunnels
    ] == [(0, 0, [0], carried, 1)]
    assert r.plan is r.min_cut is r.marginal_cost_curve is None
    assert r.total_cost == carried
    assert [
        (
            s.source,
            s.offered,
            s.admitted,
            s.delivered,
            s.source_drops,
            s.total_edge_drops,
        )
        for s in r.by_source
    ] == [(0, None, carried, carried, 0, 0)]
    assert r.edge_flows.tolist() == [carried]
    assert r.residual.tolist() == [5 - carried]
    assert r.edge_reservations.tolist() == [0]
    assert r.saturated_boundary_edges == ([0] if carried == 5 else [])
