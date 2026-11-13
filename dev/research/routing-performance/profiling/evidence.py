"""Strict complete export comparison; only metadata.workers may differ."""


def normalize(value, parent=None):
    if isinstance(value, dict):
        return {
            k: None if k == "workers" and parent == "metadata" else normalize(v, k)
            for k, v in value.items()
        }
    if isinstance(value, list):
        return [normalize(v) for v in value]
    return value


def equal(a, b, path="root"):
    assert type(a) is type(b), (path, type(a), type(b))
    if isinstance(a, dict):
        assert list(a) == list(b), (path, "keys/order")
        for key in a:
            equal(a[key], b[key], f"{path}.{key}")
    elif isinstance(a, list):
        assert len(a) == len(b), (path, "length")
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            equal(x, y, f"{path}[{i}]")
    else:
        assert a == b, (path, a, b)
