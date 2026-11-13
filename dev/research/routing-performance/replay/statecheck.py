"""Rebuild exact owner state independently from stored observation deltas."""

import base64
import struct
from fractions import Fraction

import numpy as np

from values import dumps


def unpack(value, blobs, core, field=""):
    if not isinstance(value, (dict, list)):
        return value
    if isinstance(value, list):
        return [unpack(v, blobs, core) for v in value]
    if "$array" in value:
        return blobs.get(value)
    if "$float" in value:
        return struct.unpack("!d", bytes.fromhex(value["$float"]))[0]
    if "$fraction" in value:
        return Fraction(*value["$fraction"])
    if "$scalar" in value:
        dtype, raw = value["$scalar"]
        return np.frombuffer(base64.b64decode(raw), dtype=dtype)[0]
    if "$enum" in value:
        cls, name, number = value["$enum"]
        result = getattr(getattr(core, cls), name)
        assert result.value == number
        return result
    if "$tuple" in value or "$list" in value:
        key = "$tuple" if "$tuple" in value else "$list"
        items = [unpack(v, blobs, core) for v in value[key]]
        return tuple(items) if key == "$tuple" else items
    if "$range" in value:
        return range(*value["$range"])
    if "$dict" in value:
        return {
            unpack(k, blobs, core): unpack(v, blobs, core) for k, v in value["$dict"]
        }
    return {k: unpack(v, blobs, core, k) for k, v in value.items()}


def record_key(record):
    uid, number = record["flow_id"]
    return dumps({"$tuple": [{"$uid": {"$list": uid["$uid"]}}, number]})


def apply_delta(previous, delta, blobs, core):
    if delta is None:
        return previous
    if "full" in delta:
        return unpack(delta["full"], blobs, core)
    if "replace" in delta:
        return unpack(delta["replace"], blobs, core)
    if "records" in delta:
        change = delta["records"]
        records = {record_key(v): v for v in previous}
        assert len(records) == len(previous)
        for key in change["removed"]:
            del records[key]
        for key, value in change["added"].items():
            assert key not in records
            records[key] = unpack(value, blobs, core)
            assert record_key(records[key]) == key
        for key, change_value in change["changed"].items():
            records[key] = apply_delta(records[key], change_value, blobs, core)
        assert len(set(change["order"])) == len(change["order"]) and set(
            records
        ) == set(change["order"])
        return [records[key] for key in change["order"]]
    if "indices" in delta or "pair_indices" in delta:
        key = "$ndarray" if "indices" in delta else "$int_pairs"
        indices = blobs.get(delta.get("indices", delta.get("pair_indices")))
        result = {**previous, key: previous[key].copy()}
        result[key].reshape(-1)[indices] = blobs.get(delta["values"])
        return result
    if "fields" in delta:
        result = dict(previous)
        for key, change in delta["fields"].items():
            if change == {"removed": True}:
                del result[key]
            else:
                result[key] = apply_delta(result.get(key), change, blobs, core)
        return result
    if "items" in delta:
        result = list(previous)
        for index, change in delta["items"].items():
            result[int(index)] = apply_delta(result[int(index)], change, blobs, core)
        return tuple(result) if type(previous) is tuple else result
    raise ValueError(f"unknown delta {list(delta)}")


def exact(left, right, path="$"):
    if type(left) is not type(right):
        raise AssertionError((path, type(left).__name__, type(right).__name__))
    if isinstance(left, np.ndarray):
        assert (
            left.dtype == right.dtype
            and left.shape == right.shape
            and left.tobytes() == right.tobytes()
        ), path
    elif type(left) is dict:
        assert left.keys() == right.keys(), (path, "keys")
        for key in left:
            exact(left[key], right[key], path + "." + key)
    elif type(left) in (list, tuple):
        assert len(left) == len(right), (path, "length")
        for i, (a, b) in enumerate(zip(left, right, strict=True)):
            exact(a, b, path + f"[{i}]")
    elif type(left) is float:
        assert struct.pack("!d", left) == struct.pack("!d", right), path
    elif isinstance(left, np.generic):
        assert left.dtype == right.dtype and left.tobytes() == right.tobytes(), path
    else:
        assert left == right, path
