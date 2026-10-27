"""Public concurrency promises on ordinary and free-threaded, hooks ON/OFF builds."""

from concurrency_support import run_python


def test_shared_policy_increments():
    run_python("""
        ledger = ngc.FlowLedger(graph(100))
        a, b = policy(ledger), policy(ledger)
        def place(i):
            return (a if i % 2 else b).place(1).call_delivered
        with ThreadPoolExecutor(4) as pool:
            results = list(pool.map(place, range(240)))
        assert sum(results) == a.delivered() + b.delivered() == 100
        assert ledger.carried().tolist() == [100]
        assert ledger.residual().tolist() == [0]
        assert sum(row.delivered for row in ledger.records()) == 100
    """)


def test_replacements_replay_by_flow_id():
    run_python("""
        for mode in (ngc.BatchMode.ORDERED, ngc.BatchMode.SIMULTANEOUS):
            ledger = ngc.FlowLedger(graph())
            a, b = policy(ledger, lossy=True), policy(ledger, lossy=True)
            ngc.apply_actions(ledger, [ngc.Replace(a, 2), ngc.Replace(b, 3)], mode=mode)
            owner_a, owner_b = [row.policy_id for row in ledger.records()]
            def writer():
                for i in range(80):
                    offers = (2, 3) if i % 2 else (6, 4)
                    results = ngc.apply_actions(ledger, [ngc.Replace(a, offers[0]), ngc.Replace(b, offers[1])], mode=mode)
                    assert [row.delivered for row in results] == list(offers)
            def reader():
                for _ in range(200):
                    records = ledger.records()
                    amounts = {row.policy_id: row.delivered for row in records}
                    assert (amounts[owner_a], amounts[owner_b]) in ((2, 3), (6, 4))
                    ids = [row.flow_id for row in records]
                    assert ids == sorted(ids) and len(set(ids)) == 2
            with ThreadPoolExecutor(3) as pool:
                jobs = [pool.submit(writer), pool.submit(reader), pool.submit(reader)]
                for job in jobs:
                    job.result()
    """)


def test_independent_ledgers():
    run_python("""
        network = graph(1000)
        ledgers = [ngc.FlowLedger(network), ngc.FlowLedger(network)]
        policies = [policy(ledger) for ledger in ledgers]
        barrier = threading.Barrier(2)
        def worker(p):
            barrier.wait(5)
            return sum(p.place(1).call_delivered for _ in range(400))
        with ThreadPoolExecutor(2) as pool:
            assert list(pool.map(worker, policies)) == [400, 400]
        assert [ledger.residual().tolist() for ledger in ledgers] == [[600], [600]]
        assert ledgers[0].uid() != ledgers[1].uid()
    """)


def test_immutable_configuration_and_snapshots():
    run_python("""
        network = graph()
        ledger = ngc.FlowLedger(network)
        p = policy(ledger)
        old_result = p.place(3, ngc.OutputOptions(residual=True, edge_flows=True))
        old_record = p.records()[0]
        old_program = p.program
        config = p.config
        def reader(_):
            for _ in range(100):
                assert old_result.delivered == old_record.delivered == 3
                assert old_result.residual.tolist() == [7]
                assert old_result.edge_flows.tolist() == [3]
                assert old_record.route.graph_uid() == network.uid()
                assert old_program.sources == (0,) and old_program.destinations == (1,)
                assert config.admission == ngc.Admission.LOSSLESS
        def writer(_):
            for i in range(50):
                p.replace(i % 10)
        with ThreadPoolExecutor(4) as pool:
            jobs = [pool.submit(reader, i) for i in range(3)] + [pool.submit(writer, 0)]
            for job in jobs:
                job.result()
        assert old_record.delivered == old_result.delivered == 3
    """)


def test_cache_lifecycle_and_retained_values():
    run_python("""
        network = graph()
        cache = ngc.ShortestPathCache(network)
        options = ngc.EvaluationOptions(ngc.FlowPolicyConfig(routing=ngc.ShortestPaths()), cache=cache, offer=3)
        retained = ngc.evaluate_demand(network, 0, 1, options)
        def reader():
            for _ in range(100):
                assert ngc.evaluate_pairs(network, [(0, 1), (0, 1)], options, workers=1)[0].delivered == 3
        def lifecycle():
            for _ in range(100):
                cache.freeze()
                cache.thaw()
        with ThreadPoolExecutor(3) as pool:
            jobs = [pool.submit(reader), pool.submit(reader), pool.submit(lifecycle)]
            for job in jobs:
                job.result()
        assert cache.stats().hits == 400
        cache.clear()
        assert cache.stats().entries == 0
        assert retained.program.sources == (0,) and retained.program.destinations == (1,)
        assert retained.delivered == 3
    """)


def test_public_policy_teardown():
    # The keeper's __del__ runs in interpreter teardown, on the finalizing thread after
    # Py_IsInitialized() is already 0: its native read must retake the GIL there.
    output = run_python(
        """
        class Keeper:
            def __init__(self, ledger):
                self.ledger = ledger

            def __del__(self):
                print("CHILD: teardown read", self.ledger.carried().tolist(), flush=True)

        network = graph()
        ledger = ngc.FlowLedger(network)
        keeper = Keeper(ledger)
        p = policy(ledger)
        p.place(3)
        record = ledger.records()[0]
        del p
        gc.collect()
        assert ledger.records()[0].flow_id == record.flow_id
        assert ledger.carried().tolist() == [3]
        copied = ledger.copy()
        assert copied.carried().tolist() == [3]
        assert copied.uid() != ledger.uid()
        ledger.reset()
        assert ledger.records() == []
        assert ledger.residual().tolist() == [10]
    """,
        cleanup=False,
    )
    assert "CHILD: teardown read [0]" in output.splitlines()


def test_final_policy_rejects_swallowed_constructor_error():
    run_python("""
        with pytest.raises(TypeError):
            class Swallow(ngc.FlowPolicy):
                def __init__(self):
                    try:
                        super().__init__(None, 0, 1, ngc.FlowPolicyConfig(routing=ngc.ShortestPaths()))
                    except Exception:
                        pass
        network = graph()
        ledger = ngc.FlowLedger(network)
        cache = ngc.ShortestPathCache(network)
        cache.freeze()
        shell = ngc.FlowPolicy.__new__(ngc.FlowPolicy)
        with pytest.raises(ngc.ShortestPathCacheMiss):
            ngc.FlowPolicy.__init__(shell, ledger, 0, 1, ngc.FlowPolicyConfig(routing=ngc.ShortestPaths()), cache=cache)
        with pytest.raises(TypeError):
            shell.place(1)
        del shell
        gc.collect()
        ledger.reset()
        class Pretender:
            reads = 0
            @property
            def __class__(self):
                self.reads += 1
                return ngc.FlowPolicy
        fake = Pretender()
        # Python's isinstance accepts this object; its actual native layout is unrelated.
        assert type(fake) is Pretender and isinstance(fake, ngc.FlowPolicy)
        assert fake.reads == 1
        for action in (ngc.Place, ngc.Replace, ngc.Forward):
            with pytest.raises(TypeError, match='policy must be a FlowPolicy'):
                action(fake, 1)
        entries = [
            lambda: ngc.FlowPolicy.config.fget(fake),
            lambda: ngc.FlowPolicy.program.fget(fake),
            lambda: ngc.FlowPolicy.place(fake, 1),
            lambda: ngc.FlowPolicy.evaluate(fake, 1),
            lambda: ngc.FlowPolicy.replace(fake, 1),
            lambda: ngc.FlowPolicy.forward(fake, 1),
            lambda: ngc.FlowPolicy.remove(fake),
            lambda: ngc.FlowPolicy.delivered(fake),
            lambda: ngc.FlowPolicy.records(fake),
            lambda: ngc.FlowPolicy.record(fake, (ledger.uid(), 0)),
        ]
        for entry in entries:
            with pytest.raises(TypeError, match='policy must be a FlowPolicy'):
                entry()
        assert fake.reads == 1
        assert ledger.carried().tolist() == [0]
    """)
