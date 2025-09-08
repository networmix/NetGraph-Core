def test_imports():
    import netgraph_core as ngc

    assert hasattr(ngc, "StrictMultiDiGraph")
    assert hasattr(ngc, "Backend")
    assert hasattr(ngc, "Algorithms")
    assert hasattr(ngc, "EdgeSelection")
    assert hasattr(ngc, "EdgeTieBreak")
