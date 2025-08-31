def test_imports():
    import netgraph_core as ngc

    assert hasattr(ngc, "StrictMultiDiGraph")
    assert hasattr(ngc, "spf")
    assert hasattr(ngc, "ksp")
    assert hasattr(ngc, "max_flow")
