"""Native progress, exclusion, publication and lease witnesses, with parent deadlines."""

import sysconfig

import pytest
from concurrency_support import run_python
from private_hooks import needs_hooks

pytestmark = needs_hooks
GIL_ONLY = pytest.mark.skipif(
    bool(sysconfig.get_config_var("Py_GIL_DISABLED")), reason="GIL-proof-only"
)
SETUP = """
network = graph()
ledger = ngc.FlowLedger(network)
p = policy(ledger)
Probe = new_probe
"""


def run(source):
    run_python(SETUP + source)


def test_conversions_run_before_native_locks():
    run("""
class Index:
    def __index__(self):
        assert ledger.residual().tolist() == [10]
        assert p.program is None
        return 3
probe = Probe()
r = probe.run(p.place, Index())
assert r.delivered == 3
before = probe.counts()
with pytest.raises(ValueError, match='offer must be'):
    probe.run(p.forward, -1)
assert probe.counts() == before
q = policy(ledger)
cache = ngc.ShortestPathCache(network)
options = ngc.EvaluationOptions(ngc.FlowPolicyConfig(routing=ngc.ShortestPaths()), cache=cache)
invalid = {
    'scalar': lambda: p.replace(-1),
    'list': lambda: ngc.apply_actions(ledger, [ngc.Forward(p, 1), ngc.Replace(q, -1)]),
    'reset': lambda: ledger.reset([-1]),
    'demand': lambda: ngc.evaluate_demand(network, [], 1, options),
    'pairs': lambda: ngc.evaluate_pairs(network, [([], 1)], options, workers=2),
    'edge_failures': lambda: ngc.evaluate_edge_failures(network, [], 1, options, workers=2),
    'graph': lambda: ngc.Graph(2, [0], [1], [1], [0], 1.0),
    'mask': lambda: network.eligibility(excluded_edges=[99]),
    'path': lambda: ngc.Path(network, []),
    'dag': lambda: ngc.PredecessorDAG.from_csr(network, [0,1,1], [1], [0]),
    'program': lambda: ngc.ForwardingProgram.from_routes(network, [], 1, []),
    'spf': lambda: ngc.spf(network, [], ngc.SpfOptions(ngc.Direction.FORWARD)),
    'to_units': lambda: ngc.to_units([float('nan')], 1.0),
    'as_float': lambda: ngc.as_float([0.5], 1.0),
}
for name, call in invalid.items():
    failure = Probe()
    with pytest.raises((TypeError, ValueError)):
        failure.run(call)
    assert all(value == 0 for value in failure.counts().values()), (name, failure.counts())
assert ledger.carried().tolist() == [3] and cache.stats().entries == 0
# Every immutable entry check precedes native construction, for Path and DAG.
from helpers import dag_of
chain = ngc.Graph(4, [0,1,2], [1,2,3], [10]*3, [1]*3, 1.0)
foreign = ngc.Graph(4, [0,1,2], [1,2,3], [10]*3, [1]*3, 1.0)
whole = ngc.Path(chain, [[0],[1],[2]])
split = ngc.Graph(3, [0,0], [1,1], [10,10], [1,1], 1.0)
cases = [
    (chain, 0, 3, ngc.RouteDefinition(0, ngc.Path(foreign, [[0]])), 'route belongs to another graph'),
    (chain, 0, 3, ngc.RouteDefinition(1, whole), 'head 1 is not a source'),
    (chain, 0, 3, ngc.RouteDefinition(0, ngc.Path(chain, [[1],[2]])), 'the path starts at 1, not at head 0'),
    (chain, 0, 3, ngc.RouteDefinition(0, ngc.Path(chain, [[0]])), 'the path ends at 1, not at a destination'),
    (chain, 0, [1,3], ngc.RouteDefinition(0, whole), 'the path meets destination 1 before its end'),
    (split, 0, 1, ngc.RouteDefinition(0, ngc.Path(split, [[0,1]])), 'a tunnel has one edge per hop'),
    (chain, 0, 3, ngc.RouteDefinition(0, whole, weight=-1), 'weight must be an integer in [0, 2^62)'),
    (chain, 0, 3, ngc.RouteDefinition(0, dag_of(foreign, [0,1,2])), 'route belongs to another graph'),
    (chain, 0, 3, ngc.RouteDefinition(1, dag_of(chain, [0,1,2])), 'head 1 is not a source'),
    (chain, 0, 3, ngc.RouteDefinition(0, dag_of(chain, [1,2])), 'no walk from head 0 to a destination'),
    (chain, 0, 3, ngc.RouteDefinition(0, dag_of(chain, [0,1,2]), weight=-1), 'weight must be an integer in [0, 2^62)'),
]
for route_graph, sources, destinations, route, error in cases:
    probe = Probe()
    with pytest.raises(ValueError) as caught:
        probe.run(ngc.ForwardingProgram.from_routes, route_graph, sources, destinations, [route])
    assert str(caught.value) == 'routes[0]: ' + error
    assert probe.counts()['native_begin'] == 0, (error, probe.counts())
for route in (whole, dag_of(chain, [0,1,2])):
    probe = Probe()
    result = probe.run(ngc.ForwardingProgram.from_routes, chain, 0, 3, [ngc.RouteDefinition(0, route)])
    assert result.sources == (0,) and result.destinations == (3,)
    assert probe.counts()['native_begin'] == probe.counts()['native_complete'] == 1
    assert probe.counts()['route_definition_check'] == 1, probe.counts()
    assert probe.counts()['route_walk'] == (0 if isinstance(route, ngc.Path) else 1), probe.counts()
# The DAG's walk/refusal precedes its weight, and all of route 0 precedes route 1.
for routes, error in [
    ([ngc.RouteDefinition(0, dag_of(chain, [1,2]), weight=-1)], 'routes[0]: no walk from head 0 to a destination'),
    ([ngc.RouteDefinition(0, dag_of(chain, [0,1,2]), weight=-1), ngc.RouteDefinition(1, whole)], 'routes[0]: weight must be an integer in [0, 2^62)'),
]:
    probe = Probe()
    with pytest.raises(ValueError) as caught:
        probe.run(ngc.ForwardingProgram.from_routes, chain, 0, 3, routes)
    assert str(caught.value) == error
    assert probe.counts()['native_begin'] == 0
""")


def test_initializing_policy_is_not_exposed():
    run("""
candidate = ngc.FlowPolicy.__new__(ngc.FlowPolicy)
config = ngc.FlowPolicyConfig(routing=ngc.ShortestPaths())
holder = Probe(pause='ledger_acquired')
task = start(lambda: candidate.__init__(ledger, 0, 1, config), holder)
assert holder.wait('ledger_acquired')
for read in (lambda: candidate.program, lambda: candidate.config,
             candidate.delivered, candidate.records, lambda: candidate.place(1),
             lambda: ngc.Place(candidate, 1)):
    with pytest.raises(TypeError, match='FlowPolicy is not constructed'):
        read()
holder.release()
finish(task)
assert candidate.place(3).delivered == 3
""")


def test_inputs_are_copied_before_native_work():
    run(
        """
service = np.array([7], dtype=np.int64)
probe = Probe(pause='native_begin')
task = start(lambda: ledger.reset(service), probe)
assert probe.wait('native_begin')
service[0] = 1
probe.release()
finish(task)
assert ledger.residual().tolist() == [7]
""".replace("task = start", "del p\ntask = start")
    )


def test_options_snapshot_at_entry():
    run("""
options = ngc.OutputOptions(residual=True, edge_flows=True)
probe = Probe(pause='ledger_acquired')
task = start(lambda: p.place(3, options), probe)
assert probe.wait('ledger_acquired')
with pytest.raises(AttributeError):
    options.residual = False
probe.release()
r = finish(task)
assert r.residual.tolist() == [7] and r.edge_flows.tolist() == [3]
""")


@GIL_ONLY
def test_native_phase_allows_heartbeat():
    run("""
assert sys._is_gil_enabled()
sys.setswitchinterval(0.000001)
network = ngc.Graph(2, [0]*3000, [1]*3000, [1000]*3000, [1]*3000, 1.0)
ledger = ngc.FlowLedger(network)
p = policy(ledger)
p.replace(1)
probe = Probe()
operations = [lambda: p.place(1), lambda: p.evaluate(3), lambda: p.replace(3),
              lambda: p.forward(3), lambda: ngc.apply_actions(ledger, [ngc.Replace(p, 3)])]
observations = {}
for name, operation in zip(('place', 'evaluate', 'replace', 'forward', 'ordered_list'), operations, strict=True):
    task = start(lambda: [operation() for _ in range(60)], probe)
    overlaps = 0
    while task[0].is_alive():
        counts = probe.counts()
        overlaps += counts['native_begin'] > counts['native_complete']
    finish(task)
    observations[name] = overlaps

# Each lock-wait entry has a finite native pause, so a retained-GIL plant
# terminates with zero overlap rather than stranding the Python controller.
def active_pause(probe):
    first = probe.counts()
    if first['pause_begin'] <= first['pause_end']:
        return None
    # counts() snapshots while detached. Returning here proves Python ran;
    # the fresh snapshot must still identify the same unfinished pause.
    second = probe.counts()
    if (first['pause_begin'], first['pause_end']) == (second['pause_begin'], second['pause_end']):
        return second
    return None
network = graph()
bare = ngc.FlowLedger(network)
active_ledger = ngc.FlowLedger(network)
owner = policy(active_ledger)
owner.replace(3)
record_id = owner.records()[0].flow_id
cache = ngc.ShortestPathCache(network)
options = ngc.EvaluationOptions(ngc.FlowPolicyConfig(routing=ngc.ShortestPaths()), cache=cache, offer=3)
held = [policy(ngc.FlowLedger(network))]
def destroy():
    held.clear()
cases = [
    ('uid', bare.uid, 'ledger_acquired'),
    ('graph_read', bare.graph, 'ledger_acquired'),
    ('capacities', bare.edge_capacity_units, 'ledger_acquired'),
    ('residual', bare.residual, 'ledger_acquired'),
    ('carried', bare.carried, 'ledger_acquired'),
    ('reserved', bare.reserved, 'ledger_acquired'),
    ('unreserved', bare.unreserved, 'ledger_acquired'),
    ('ledger_records', bare.records, 'ledger_acquired'),
    ('ledger_record', lambda: active_ledger.record(record_id), 'ledger_acquired'),
    ('saturated_boundary', lambda: bare.saturated_boundary_edges(0), 'ledger_acquired'),
    ('reset', bare.reset, 'ledger_acquired'),
    ('copy', bare.copy, 'ledger_acquired'),
    ('program_read', lambda: owner.program, 'ledger_acquired'),
    ('delivered_read', owner.delivered, 'ledger_acquired'),
    ('policy_records', owner.records, 'ledger_acquired'),
    ('policy_record', lambda: owner.record(record_id), 'ledger_acquired'),
    ('construction', lambda: policy(bare), 'ledger_acquired'),
    ('remove', owner.remove, 'ledger_acquired'),
    ('destruction', destroy, 'ledger_acquired'),
    ('freeze', cache.freeze, 'cache_acquired'),
    ('thaw', cache.thaw, 'cache_acquired'),
    ('clear', cache.clear, 'cache_acquired'),
    ('evaluate_demand', lambda: ngc.evaluate_demand(network, 0, 1, options), 'cache_acquired'),
    ('evaluate_pairs', lambda: ngc.evaluate_pairs(network, [(0, 1)], options, workers=1), 'cache_acquired'),
    ('evaluate_failures', lambda: ngc.evaluate_edge_failures(network, 0, 1, options, 1), 'cache_acquired'),
]
for name, operation, stage in cases:
    probe = Probe(pause=stage, auto_release_ms=200)
    task = start(operation, probe)
    overlap = 0
    while task[0].is_alive():
        overlap += active_pause(probe) is not None
    finish(task)
    observations[name] = overlap
assert all(value > 0 for value in observations.values()), observations
# Tunnel-limit preflight grows with the grid and runs before any lease or gate.
# The bounded native pause makes a retained-GIL mistake fail by assertion.
side = 80
us, vs = [], []
for u in range(side*side):
    if u % side + 1 < side: us.append(u); vs.append(u+1)
    if u + side < side*side: us.append(u); vs.append(u+side)
grid = ngc.Graph(side*side, us, vs, [10]*len(us), [1]*len(us), 1.0)
pairs = [(i, side*side-1) for i in range(64)]
for routing in (ngc.DistributedTE(ngc.TunnelCount(1)), ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL, tunnel_count=ngc.TunnelCount(1))):
    options = ngc.EvaluationOptions(ngc.FlowPolicyConfig(routing=routing), design_demand=1, offer=1)
    probe = Probe(pause='tunnel_roster_search', auto_release_ms=200)
    task = start(lambda: ngc.evaluate_pairs(grid, pairs, options, workers=4), probe)
    overlap = 0
    while task[0].is_alive():
        counts = active_pause(probe)
        if counts is not None:
            assert counts['native_begin'] == counts['ledger_wait'] == counts['cache_queued'] == 0
            overlap += 1
    results = finish(task)
    assert overlap > 0, 'tunnel preflight held the GIL'
    assert len(results) == len(pairs) and all(r.delivered == 1 for r in results)
    assert probe.counts()['tunnel_roster_search'] == len(pairs), probe.counts()
    # The constructor's first pass follows the same detached-search rule.
    tunnel_ledger = ngc.FlowLedger(grid)
    constructor = Probe(pause='tunnel_roster_search', auto_release_ms=200)
    task = start(lambda: ngc.FlowPolicy(tunnel_ledger, 0, side*side-1,
                                      ngc.FlowPolicyConfig(routing=routing)), constructor)
    overlap = 0
    while task[0].is_alive():
        counts = active_pause(constructor)
        if counts is not None:
            assert counts['native_begin'] == counts['ledger_wait'] == counts['cache_queued'] == 0
            overlap += 1
    constructed = finish(task)
    assert overlap > 0, 'constructor tunnel preflight held the GIL'
    assert constructor.counts()['tunnel_roster_search'] == 1
    assert constructed.replace(1, design_demand=1).delivered == 1
    small = ngc.Graph(4, [0,1,2], [1,2,3], [10]*3, [1]*3, 1.0)
    for name, call, searches in (
        ('demand', lambda: ngc.evaluate_demand(small, 0, 3, options), 1),
        ('edge_failures', lambda: ngc.evaluate_edge_failures(small, 0, 3, options, 2), 4),
    ):
        preflight = Probe(pause='tunnel_roster_search', auto_release_ms=200)
        task = start(call, preflight)
        overlap = 0
        while task[0].is_alive():
            counts = active_pause(preflight)
            if counts is not None:
                assert counts['native_begin'] == counts['ledger_wait'] == counts['cache_queued'] == 0
                overlap += 1
        result = finish(task)
        assert overlap > 0, (name, 'tunnel preflight held the GIL')
        assert preflight.counts()['tunnel_roster_search'] == searches
        if name == 'demand':
            assert result.delivered == 1
        else:
            assert result.baseline.delivered == 1
            assert [(e.edge, e.delivery_reduction) for e in result.entries] == [(0,1), (1,1), (2,1)]
# Supplied-DAG preparation is the first pass's graph walk, before construction.
from helpers import dag_of
supplied = dag_of(grid, range(grid.num_edges()))
probe = Probe(pause='route_walk', auto_release_ms=200)
task = start(lambda: ngc.ForwardingProgram.from_routes(grid, 0, side*side-1,
                    [ngc.RouteDefinition(0, supplied)]), probe)
overlap = 0
while task[0].is_alive():
    counts = active_pause(probe)
    if counts is not None:
        assert counts['native_begin'] == counts['ledger_wait'] == counts['cache_queued'] == 0
        overlap += 1
program = finish(task)
assert overlap > 0, 'supplied route first-pass walk held the GIL'
assert program.sources == (0,) and program.destinations == (side*side-1,)
assert probe.counts()['route_walk'] == probe.counts()['route_definition_check'] == 1
""")


def test_cache_before_gate_order():
    run("""
cache = ngc.ShortestPathCache(network)
holder = Probe(pause='ledger_acquired')
a = start(lambda: p.place(2), holder)
assert holder.wait('ledger_acquired')
constructor = Probe()
b = start(lambda: policy(ledger, cache=cache), constructor)
assert constructor.wait('cache_acquired')
assert constructor.wait('ledger_wait')
early = constructor.wait('ledger_acquired', seconds=0.1)
try:
    if early:
        q = finish(b)
finally:
    holder.release()
assert finish(a).delivered == 2
if not early:
    q = finish(b)
assert not early, 'construction bypassed the ledger gate'
assert q.place(2).delivered == 2
""")


def test_cached_policies_ignore_busy_original_cache():
    run("""
cache = ngc.ShortestPathCache(network)
q = policy(ledger, cache=cache)
options = ngc.EvaluationOptions(ngc.FlowPolicyConfig(routing=ngc.ShortestPaths()), cache=cache, offer=2)
holder = Probe(pause='cache_acquired')
a = start(lambda: ngc.evaluate_demand(network, 0, 1, options), holder)
assert holder.wait('cache_acquired')
assert q.place(2).delivered == 2
assert ngc.apply_actions(ledger, [ngc.Replace(q, 3)])[0].delivered == 3
holder.release()
assert finish(a).delivered == 2
""")


def test_record_and_workers_do_not_deadlock():
    run("""
p.place(3)
record = p.records()[0]
assert p.record(record.flow_id).delivered == ledger.record(record.flow_id).delivered == 3
options = ngc.EvaluationOptions(ngc.FlowPolicyConfig(routing=ngc.ShortestPaths()), offer=4)
assert [r.delivered for r in ngc.evaluate_pairs(network, [(0, 1)]*8, options, workers=4)] == [4]*8
assert ngc.evaluate_edge_failures(network, 0, 1, options, 3).entries[0].delivery_reduction == 4
""")


def test_same_ledger_calls_serialize():
    run("""
holder, waiter = Probe(pause='ledger_acquired'), Probe()
a = start(lambda: p.place(2), holder)
assert holder.wait('ledger_acquired')
b = start(lambda: p.place(3), waiter)
assert waiter.wait('ledger_wait')
assert waiter.counts()['ledger_acquired'] == 0
holder.release()
assert finish(a).delivered == 2
assert finish(b).delivered == 5
assert ledger.carried().tolist() == [5]
""")


def test_every_ledger_reader_waits_without_gil():
    run("""
p.place(1)
flow_id = p.records()[0].flow_id
readers = [ledger.uid, ledger.graph, ledger.edge_capacity_units, ledger.residual,
           ledger.carried, ledger.reserved, ledger.unreserved, ledger.records,
           lambda: ledger.record(flow_id), lambda: ledger.saturated_boundary_edges([0]),
           p.delivered, p.records, lambda: p.record(flow_id), lambda: p.program,
           ledger.copy, lambda: ledger.reset(), p.remove]
early = []
for index, read in enumerate(readers):
    holder, waiter = Probe(pause='ledger_acquired'), Probe()
    a = start(lambda: p.place(0), holder)
    assert holder.wait('ledger_acquired')
    b = start(read, waiter)
    assert waiter.wait('ledger_wait')
    if waiter.wait('ledger_acquired', seconds=0.1):
        early.append(index)
        print('CHILD: reader', index, 'did not wait', flush=True)
    holder.release()
    finish(a)
    if index in (14, 15):
        with pytest.raises(ValueError, match='bound'):
            finish(b)
    else:
        finish(b)
assert early == [], early
assert ledger.carried().tolist() == [0] and p.program is None
""")


def test_immutable_readers_do_not_wait():
    run("""
cache = ngc.ShortestPathCache(network)
q = policy(ledger, cache=cache)
result = q.place(2)
holder = Probe(pause='ledger_acquired')
a = start(lambda: q.place(0), holder)
assert holder.wait('ledger_acquired')
assert q.config.admission == ngc.Admission.LOSSLESS
assert cache.graph_uid() == network.uid()
assert cache.stats().entries == 1
assert result.delivered == 2
assert result.program.sources == (0,)
assert network.edge_capacity_units().tolist() == [10]
holder.release()
finish(a)
""")


def test_two_ledgers_enter_together():
    run("""
other = policy(ngc.FlowLedger(network))
left, right = Probe(pause='ledger_acquired'), Probe(pause='ledger_acquired')
a = start(lambda: p.place(2), left)
assert left.wait('ledger_acquired')
b = start(lambda: other.place(3), right)
assert right.wait('ledger_acquired')
left.release()
right.release()
assert finish(a).delivered == 2
assert finish(b).delivered == 3
""")


def test_full_calls_and_reads_serialize():
    run("""
p.replace(2)
for operation, expected in zip((lambda: p.evaluate(4), lambda: p.replace(4), lambda: p.forward(4), p.remove), (2, 4, 4, 0)):
    holder, waiter = Probe(pause='ledger_acquired'), Probe()
    a = start(operation, holder)
    assert holder.wait('ledger_acquired')
    b = start(p.delivered, waiter)
    assert waiter.wait('ledger_wait')
    assert waiter.counts()['ledger_acquired'] == 0
    holder.release()
    finish(a)
    assert finish(b) == expected
""")


def test_list_publication_is_atomic():
    run("""
for mode in (ngc.BatchMode.ORDERED, ngc.BatchMode.SIMULTANEOUS):
    ledger = ngc.FlowLedger(network)
    a, b = policy(ledger, lossy=True), policy(ledger, lossy=True)
    holder, waiter = Probe(pause='after_exchange'), Probe()
    writer = start(lambda: ngc.apply_actions(ledger, [ngc.Replace(a, 2), ngc.Replace(b, 3)], mode=mode), holder)
    assert holder.wait('after_exchange')
    reader = start(lambda: (a.program, b.program, ledger.records()), waiter)
    assert waiter.wait('ledger_wait')
    assert waiter.counts()['ledger_acquired'] == 0
    holder.release()
    assert [r.delivered for r in finish(writer)] == [2, 3]
    pa, pb, records = finish(reader)
    assert pa is not None and pb is not None
    assert [r.delivered for r in records] == [2, 3]
""")


def test_reservations_and_program_are_atomic():
    run("""
q = ngc.FlowPolicy(ledger, 0, 1, ngc.FlowPolicyConfig(routing=ngc.DistributedTE(ngc.TunnelCount(1))))
holder, waiter = Probe(pause='after_exchange'), Probe()
a = start(lambda: q.replace(4, design_demand=6), holder)
assert holder.wait('after_exchange')
b = start(lambda: (q.program, q.records(), ledger.reserved()), waiter)
assert waiter.wait('ledger_wait')
assert waiter.counts()['ledger_acquired'] == 0
holder.release()
r = finish(a)
program, records, holds = finish(b)
assert program is not None and r.plan.total_reserved_bandwidth == 6
assert holds.tolist() == [6] and records[0].edge_reservations == [(0, 6)]
""")


CACHE_SETUP = """
cache = ngc.ShortestPathCache(network)
options = ngc.EvaluationOptions(ngc.FlowPolicyConfig(routing=ngc.ShortestPaths()), cache=cache, offer=2)
ngc.evaluate_demand(network, 0, 1, options)
"""


def test_building_cache_serializes_native_calls():
    run(
        CACHE_SETUP
        + """
holder, waiter = Probe(pause='cache_acquired'), Probe()
a = start(lambda: ngc.evaluate_demand(network, 0, 1, options), holder)
assert holder.wait('cache_acquired')
b = start(lambda: ngc.evaluate_demand(network, 0, 1, options), waiter)
assert waiter.wait('cache_queued')
early = waiter.wait('cache_acquired', seconds=0.1)
try:
    if early:
        second = finish(b)
finally:
    holder.release()
first = finish(a)
if not early:
    second = finish(b)
assert first.delivered == second.delivered == 2
assert not early, 'building calls shared a lease'
"""
    )


def test_frozen_cache_allows_overlapping_calls():
    run(
        CACHE_SETUP
        + """
cache.freeze()
left, right = Probe(pause='cache_acquired'), Probe(pause='cache_acquired')
a = start(lambda: ngc.evaluate_demand(network, 0, 1, options), left)
assert left.wait('cache_acquired')
b = start(lambda: ngc.evaluate_demand(network, 0, 1, options), right)
assert right.wait('cache_acquired')
left.release()
right.release()
assert finish(a).delivered == finish(b).delivered == 2
"""
    )


def test_cache_fifo_admission():
    run(
        CACHE_SETUP
        + """
for frozen in (False, True):
    if frozen:
        cache.freeze()
    holder, writer, later = Probe(pause='cache_acquired'), Probe(pause='cache_acquired'), Probe()
    a = start(lambda: ngc.evaluate_demand(network, 0, 1, options), holder)
    assert holder.wait('cache_acquired')
    b = start(cache.thaw, writer)
    assert writer.wait('cache_queued')
    c = start(lambda: ngc.evaluate_demand(network, 0, 1, options), later)
    assert later.wait('cache_queued')
    holder.release()
    assert writer.wait('cache_acquired')
    assert later.counts()['cache_acquired'] == 0
    writer.release()
    finish(a)
    finish(b)
    assert finish(c).delivered == 2
"""
    )


def test_workers_borrow_outer_lease():
    run(
        CACHE_SETUP
        + """
cache.freeze()
holder, lifecycle = Probe(pause='worker_begin'), Probe()
a = start(lambda: ngc.evaluate_pairs(network, [(0, 1)]*6, options, workers=3), holder)
assert holder.wait('worker_begin', count=3)
b = start(cache.thaw, lifecycle)
assert lifecycle.wait('cache_queued')
holder.release()
assert [r.delivered for r in finish(a)] == [2]*6
finish(b)
assert holder.counts()['cache_acquired'] == 1
# With no supplied cache, each worker owns a building cache. Clearing it for
# the next failure state must not enroll a public lifecycle lease either.
network = ngc.Graph(3, [0,0,1,1], [1,1,2,2], [3]*4, [1]*4, 1.0)
options = ngc.EvaluationOptions(ngc.FlowPolicyConfig(routing=ngc.ShortestPaths()), offer=4)
probe = Probe()
result = probe.run(lambda: ngc.evaluate_edge_failures(network, 0, 2, options, workers=2))
assert result.baseline.call_delivered == 4
assert [entry.delivery_reduction for entry in result.entries] == [1]*4
assert probe.counts()['cache_queued'] == 0 and probe.counts()['cache_acquired'] == 0
"""
    )


def test_failed_waiter_unblocks_lifecycle():
    run(
        CACHE_SETUP
        + """
for point in ('cache_queued', 'cache_acquired', 'ledger_wait', 'ledger_acquired'):
    failure = Probe(fail=point)
    before = cache.stats().hits
    with pytest.raises(RuntimeError, match='injected native failure'):
        failure.run(lambda: policy(ledger, cache=cache))
    cache.freeze()
    cache.thaw()
    assert ngc.evaluate_demand(network, 0, 1, options).delivered == 2
"""
    )


def test_worker_validation_uses_leased_mode():
    run(
        CACHE_SETUP
        + """
holder, waiter = Probe(pause='cache_acquired'), Probe()
a = start(cache.freeze, holder)
assert holder.wait('cache_acquired')
b = start(lambda: ngc.evaluate_pairs(network, [(0, 1)]*4, options, workers=2), waiter)
assert waiter.wait('cache_queued')
holder.release()
finish(a)
assert [r.delivered for r in finish(b)] == [2]*4
cache.thaw()
with pytest.raises(ValueError, match='supplied building'):
    ngc.evaluate_pairs(network, [(0, 1)], options, workers=2)
with pytest.raises(ValueError, match='sources'):
    ngc.evaluate_pairs(network, [([], 1)], options, workers=2)
"""
    )


def test_frozen_tally_finishes_before_lifecycle():
    run(
        CACHE_SETUP
        + """
cache.freeze()
before = cache.stats().hits
holder, lifecycle = Probe(pause='cache_tally_complete'), Probe()
a = start(lambda: ngc.evaluate_pairs(network, [(0, 1)]*3, options, workers=3), holder)
assert holder.wait('cache_tally_complete')
assert cache.stats().hits == before + 3
b = start(cache.thaw, lifecycle)
assert lifecycle.wait('cache_queued')
assert lifecycle.counts()['cache_acquired'] == 0
holder.release()
assert len(finish(a)) == 3
finish(b)
assert cache.stats().hits == before + 3
"""
    )


def test_cache_lifetime_is_pinned_during_wait():
    run(
        CACHE_SETUP
        + """
holder, waiter = Probe(pause='cache_acquired'), Probe()
a = start(lambda: ngc.evaluate_demand(network, 0, 1, options), holder)
assert holder.wait('cache_acquired')
keep = [cache]
b = start(lambda: keep[0].freeze(), waiter)
assert waiter.wait('cache_queued')
del cache
gc.collect()
holder.release()
assert finish(a).delivered == 2
finish(b)
assert keep[0].stats().entries == 1
"""
    )


def test_policy_destructor_cannot_deadlock():
    run("""
holder, destructor = Probe(pause='ledger_acquired'), Probe()
ready, drop = threading.Event(), threading.Event()
def lifetime():
    # Create and delete on the same OS thread. Free-threaded CPython may defer
    # cross-thread final decrefs, which would not witness the native destructor.
    owned = policy(ledger)
    owned.place(2)
    ready.set()
    assert drop.wait(5)
    del owned
b = start(lifetime, destructor)
assert ready.wait(5)
before = destructor.counts()
a = start(lambda: p.place(3), holder)
assert holder.wait('ledger_acquired')
drop.set()
assert destructor.wait('ledger_wait', count=before['ledger_wait'] + 1)
early = destructor.wait('ledger_acquired', count=before['ledger_acquired'] + 1, seconds=0.1)
try:
    if early:
        finish(b)
finally:
    holder.release()
assert finish(a).delivered == 3
if not early:
    finish(b)
assert not early, 'destruction bypassed the ledger gate'
assert ledger._state()['orphans']
""")


def test_waiters_pin_every_dependency():
    run("""
holder, waiter = Probe(pause='ledger_acquired'), Probe()
other = [policy(ledger)]
a = start(lambda: p.place(2), holder)
assert holder.wait('ledger_acquired')
method = other[0].place
b = start(lambda: method(3), waiter)
assert waiter.wait('ledger_wait')
other.clear()
gc.collect()
holder.release()
finish(a)
assert finish(b).delivered == 3
assert ledger.carried().tolist() == [5]
""")


def test_identity_order_is_a_serial_order():
    run("""
policies = []
with ThreadPoolExecutor(4) as pool:
    policies = list(pool.map(lambda _: policy(ledger), range(16)))
assert sorted(q._owner() for q in policies) == list(range(1, 17))
with ThreadPoolExecutor(4) as pool:
    results = list(pool.map(lambda q: q.place(0), policies))
assert all(r.delivered == 0 for r in results)
assert ledger._state()['next_flow_id'] == 0
""")


def test_bind_copy_reset_races():
    run("""
del p
holder, waiter = Probe(pause='ledger_acquired'), Probe()
a = start(lambda: policy(ledger), holder)
assert holder.wait('ledger_acquired')
b = start(ledger.copy, waiter)
assert waiter.wait('ledger_wait')
holder.release()
q = finish(a)
with pytest.raises(ValueError, match='bound'):
    finish(b)
with pytest.raises(ValueError, match='bound'):
    ledger.reset()
del q, a
gc.collect()
ledger.reset()
assert ledger.copy().residual().tolist() == [10]
""")


def test_failed_acquisition_releases_waiters():
    run(
        CACHE_SETUP
        + """
for stage in ('cache_queued', 'cache_acquired', 'ledger_wait', 'ledger_acquired'):
    probe = Probe(fail=stage)
    with pytest.raises(RuntimeError, match='injected native failure'):
        probe.run(lambda: policy(ledger, cache=cache))
    assert ledger._state()['bound'] == [p._owner()]
    cache.freeze()
    cache.thaw()
    assert p.place(0).delivered == 0
"""
    )


def test_lowest_worker_exception_survives_overlap():
    run(
        CACHE_SETUP
        + """
cache.freeze()
for entry in ('pairs', 'failures'):
    if entry == 'pairs':
        operation = lambda: ngc.evaluate_pairs(network, [(0, 1)]*4, options, workers=3)
    else:
        big = ngc.Graph(2, [0]*4, [1]*4, [10]*4, [1]*4, 1.0)
        operation = lambda: ngc.evaluate_edge_failures(big, 0, 1, ngc.EvaluationOptions(ngc.FlowPolicyConfig(routing=ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL)), offer=4), 3)
    probe = Probe(workers=[0, 2], reverse=True)
    with pytest.raises(RuntimeError, match='injected worker failure 0$'):
        probe.run(operation)
"""
    )


def test_failure_does_not_strand_waiters():
    run("""
for operation in (lambda: p.place(3), lambda: p.replace(3), lambda: p.forward(3),
                  lambda: ngc.apply_actions(ledger, [ngc.Replace(p, 3)])):
    p.replace(1)
    before = (ledger.carried().tolist(), p.place(0), [r.flow_id for r in p.records()])
    failure = Probe(fail='before_publication')
    with pytest.raises(RuntimeError, match='injected native failure'):
        failure.run(operation)
    assert ledger.carried().tolist() == before[0]
    after = p.place(0)
    assert (after.requested, after.admitted, after.delivered) == (before[1].requested, before[1].admitted, before[1].delivered)
    assert [r.flow_id for r in p.records()] == before[2]
    assert p.replace(2).delivered == 2
""")


def test_zero_empty_and_masked_calls():
    run("""
assert p.place(0).program is None
assert ngc.apply_actions(ledger, []) == []
masked = ngc.FlowPolicy(ledger, 0, 1, ngc.FlowPolicyConfig(routing=ngc.ShortestPaths()), eligibility=network.eligibility(excluded_nodes=[0]))
assert masked.place(5).delivered == 0
before = p.place(3)
for mode in (ngc.BatchMode.ORDERED, ngc.BatchMode.SIMULTANEOUS):
    # A zero LOSSY participant changes no other participant's quantities.
    l = ngc.FlowLedger(network)
    a, z = policy(l, lossy=True), policy(l, lossy=True)
    r = ngc.apply_actions(l, [ngc.Replace(a, 3), ngc.Replace(z, 0)], mode=mode)
    assert [x.delivered for x in r] == [3, 0]
    assert l.carried().tolist() == [3]
assert p.place(0).program is before.program
""")


@GIL_ONLY
def test_immutable_native_operations_overlap_heartbeat():
    run("""
assert sys._is_gil_enabled()
sys.setswitchinterval(0.000001)
big = ngc.Graph(2, [0]*30000, [1]*30000, [10]*30000, [1]*30000, 1.0)
probe = Probe()
edges = list(range(30000))
from helpers import dag_of
mask = big.eligibility(excluded_edges=edges[::2])
dag = dag_of(big, edges)
chain = ngc.Graph(10001, list(range(10000)), list(range(1, 10001)), [10]*10000, [1]*10000, 1.0)
hops = [[i] for i in range(10000)]
path = ngc.Path(chain, hops)
operations = [lambda: ngc.Graph(2, [0]*30000, [1]*30000, [10]*30000, [1]*30000, 1.0),
              lambda: big.eligibility(excluded_edges=edges[::2]),
              lambda: ngc.Eligibility.intersect(mask, mask),
              lambda: ngc.PredecessorDAG.from_csr(big, dag.parent_offsets(), dag.parents(), dag.via_edges()),
              lambda: ngc.Path(chain, hops),
              path.to_dag,
              lambda: ngc.ForwardingProgram.from_dag(big, 0, 1, dag),
              lambda: ngc.ForwardingProgram.from_routes(chain, 0, 10000, [ngc.RouteDefinition(0, path)]),
              lambda: ngc.spf(big, 0, ngc.SpfOptions(ngc.Direction.FORWARD)),
              lambda: ngc.ksp(big, 0, 1, ngc.KspOptions(10)),
              lambda: ngc.resolve_to_paths(big, dag, 0, 1, split_parallel_edges=True, max_paths=10000),
              lambda: ngc.to_units(np.ones(300000), 1.0),
              lambda: ngc.as_float(np.ones(300000, dtype=np.int64), 1.0)]
observations = {}
names = ('graph', 'eligibility', 'intersect', 'csr', 'path', 'path_to_dag', 'program_dag', 'program_routes', 'spf', 'ksp', 'resolve', 'to_units', 'as_float')
for name, operation in zip(names, operations, strict=True):
    task = start(lambda: [operation() for _ in range(8)], probe)
    overlap = 0
    while task[0].is_alive():
        counts = probe.counts()
        overlap += counts['native_begin'] > counts['native_complete']
    finish(task)
    observations[name] = overlap
assert all(value > 0 for value in observations.values()), observations
""")


def test_trace_conversion_is_separate_read():
    run("""
result, trace = p._place_traced(3)
assert result.delivered == 3
before = (ledger.carried().tolist(), p.program, p._state())
first = trace.as_dict()
assert set(first) == {'prepare', 'T'}
assert trace.as_dict() == first
from helpers import plain
assert type(plain(trace)) is dict and plain(trace) == first
# One unconstrained source: the equal-allocation search checks 0, then offer 3.
assert first['T'] == [0, 3]
poisoned = dict(first)
poisoned['T'] = None
assert plain(trace) != poisoned
poisoned['T'] = [0, 4]
assert plain(trace) != poisoned
assert before == (ledger.carried().tolist(), p.program, p._state())
probe = Probe(pause='ledger_acquired')
task = start(lambda: p.place(0), probe)
assert probe.wait('ledger_acquired')
assert trace.as_dict() == first
probe.release()
finish(task)
""")


def test_daemon_exit_releases_native_resources():
    # Py_AtExit runs after Python teardown. Its native waiter owns a different
    # policy on the busy ledger, or queues thaw against the daemon's outer lease.
    for resource in ("ledger", "cache"):
        output = run_python(
            SETUP
            + f"\nresource = {resource!r}\n"
            + """
if resource == 'ledger':
    probe = Probe(pause='ledger_acquired', auto_release_ms=1000)
    ext._at_native_exit(ledger, None, probe)
    operation = lambda: p.replace(3)
    stage = 'ledger_acquired'
else:
    cache = ngc.ShortestPathCache(network)
    options = ngc.EvaluationOptions(ngc.FlowPolicyConfig(routing=ngc.ShortestPaths()), cache=cache, offer=3)
    assert ngc.evaluate_demand(network, 0, 1, options).delivered == 3
    cache.freeze()
    probe = Probe(pause='cache_acquired', auto_release_ms=1000)
    ext._at_native_exit(None, cache, probe)
    operation = lambda: ngc.evaluate_demand(network, 0, 1, options)
    stage = 'cache_acquired'
task = start(operation, probe)
assert probe.wait(stage)
assert probe.counts()['pause_begin'] == 1 and probe.counts()['pause_end'] == 0
print('CHILD: finalizing with daemon holding', resource, flush=True)
""",
            cleanup=False,
        )
        assert "CHILD: native finalizer waiting after Python teardown" in output
        assert "CHILD: native finalizer completed" in output


def test_every_program_representation_uses_the_transaction_gate():
    run("""
from helpers import dag_of, plain
from test_apply_actions import World
from test_contract_examples import mismatch
# Two heads, a branch and a merge, and one shared final link. Each program
# representation describes these same walks; a competing offer makes the final
# link bind. Compare whole serial states, not only aggregate delivery.
network = ngc.Graph(7, [0,1,2,2,3,4,5], [2,2,3,4,5,5,6], [10]*7, [1]*7, 1.0)
dag = dag_of(network, range(7))
paths = [ngc.Path(network, [[first],[branch],[last],[6]]) for first in (0,1) for branch,last in ((2,4),(3,5))]
programs = {
    'pooled_dag': ngc.ForwardingProgram.from_dag(network, [0,1], 6, dag),
    'tunnels': ngc.ForwardingProgram.from_routes(network, [0,1], 6, [ngc.RouteDefinition(i//2, path) for i,path in enumerate(paths)]),
    'split_routes': ngc.ForwardingProgram.from_routes(network, [0,1], 6, [ngc.RouteDefinition(head, dag) for head in (0,1)]),
    'mixed': ngc.ForwardingProgram.from_routes(network, [0,1], 6, [ngc.RouteDefinition(0, dag), ngc.RouteDefinition(0, paths[0]), ngc.RouteDefinition(1, dag), ngc.RouteDefinition(1, paths[2])]),
}
outputs = ngc.OutputOptions(edge_flows=True, residual=True, saturated_boundary_edges=True, edge_reservations=True)
def make(representation):
    world = World(network)
    routing = ngc.ShortestPaths() if representation == 'shortest_paths' else ngc.Pinned(programs[representation])
    config = ngc.FlowPolicyConfig(routing=routing, admission=ngc.Admission.LOSSY)
    for _ in range(2):
        world.add([0,1], [6], config).replace(2)
    return world

def operation(world, entry):
    a,b = world.p
    if entry == 'place':
        return lambda: a.place(7, outputs)
    if entry == 'evaluate':
        return lambda: a.evaluate(9, outputs=outputs)
    if entry == 'replace':
        return lambda: a.replace(9, outputs=outputs)
    if entry == 'forward':
        return lambda: a.forward(9, outputs=outputs)
    mode = ngc.BatchMode.ORDERED if entry == 'ordered' else ngc.BatchMode.SIMULTANEOUS
    return lambda: ngc.apply_actions(world.ledger, [ngc.Replace(a, 9), ngc.Replace(b, 7)], mode=mode, outputs=outputs)

def capture(task):
    try:
        return finish(task)
    except Exception as error:
        # A reader of a partial publication can fail on a retired record as well
        # as return a mixed value. Keep that failure and exercise every other path.
        return ('unexpected task error', type(error).__name__, str(error))

checked, failures = [], []
for representation in ('shortest_paths', *programs):
    entries = ('place', 'evaluate', 'replace', 'ordered', 'simultaneous')
    if representation == 'shortest_paths':
        entries += ('forward',)
    for entry in entries:
        serial, racing = make(representation), make(representation)
        want = plain(operation(serial, entry)())
        want_state = serial.state()
        call = operation(racing, entry)
        if entry != 'evaluate':
            before = racing.state()
            with pytest.raises(RuntimeError, match='injected native failure'):
                Probe(fail='before_publication').run(call)
            assert mismatch(racing.state(), before) is None, (representation, entry, 'failed call published')
        holder = Probe(pause='ledger_acquired' if entry == 'evaluate' else 'after_exchange')
        reader = Probe()
        writer = start(call, holder)
        assert holder.wait('ledger_acquired' if entry == 'evaluate' else 'after_exchange')
        observer = start(racing.state, reader)
        assert reader.wait('ledger_wait')
        early = False
        got_state = None
        try:
            early = reader.wait('ledger_acquired', seconds=0.1)
            if early:
                # A gate-removal plant still makes an owned snapshot while the
                # writer is stopped. Join the reader before resuming writes.
                got_state = capture(observer)
        finally:
            holder.release()
        got = plain(capture(writer))
        if not early:
            got_state = capture(observer)
        differences = (early, mismatch(got, want), mismatch(got_state, want_state))
        if any(differences):
            failures.append((representation, entry, differences))
            print('CHILD: transaction', representation, entry, 'failed', flush=True)
        checked.append((representation, entry))
assert len(checked) == 26
assert failures == [], failures
""")


def test_batch_first_pass_uses_one_detach_and_no_computation_event():
    run("""
chain = ngc.Graph(4, [0,1,2], [1,2,3], [10]*3, [1]*3, 1.0)
for routing in (ngc.DistributedTE(ngc.TunnelCount(1)),
                ngc.CentralizedTE(source_allocation=ngc.SourceAllocation.EQUAL, tunnel_count=ngc.TunnelCount(1))):
    options = ngc.EvaluationOptions(ngc.FlowPolicyConfig(routing=routing), design_demand=1, offer=1)
    for count in (1, 8):
        probe = Probe()
        result = probe.run(ngc.evaluate_pairs, chain, [(0,3)]*count, options, None, 3)
        counts = probe.counts()
        assert len(result) == count and all(r.delivered == 1 for r in result)
        assert counts['argument_begin'] == 1, counts
        assert counts['native_begin'] == counts['native_complete'] == 1, counts
        assert counts['tunnel_roster_search'] == count, counts
    for bad in ((2,2), ([],3), None):
        probe = Probe()
        with pytest.raises((ValueError, TypeError)):
            probe.run(ngc.evaluate_pairs, chain, [(0,3), bad], options, None, 2)
        counts = probe.counts()
        assert counts['argument_begin'] == counts['tunnel_roster_search'] == 1, counts
        assert counts['native_begin'] == counts['native_complete'] == counts['ledger_wait'] == counts['cache_queued'] == 0, counts
    failed = Probe(workers=[0,1], reverse=True)
    with pytest.raises(RuntimeError, match='^injected worker failure 0$'):
        failed.run(ngc.evaluate_pairs, chain, [(0,3), (0,3), ([],3)], options, None, 2)
    assert failed.counts()['native_begin'] == failed.counts()['ledger_wait'] == failed.counts()['cache_queued'] == 0
""")
