"""Exact values, payload storage, public snapshots and sparse state changes."""

import base64
import hashlib
import json
import operator
from itertools import chain, repeat
import struct
import threading
from fractions import Fraction
from pathlib import Path

import numpy as np

FIELDS = {
    "Graph": "uid resolution num_nodes num_edges max_edge_cost edge_cost_gcd edge_sources edge_destinations edge_capacity_units edge_costs input_positions outgoing_offsets incoming_offsets incoming_edges",
    "Eligibility": "graph_uid uid node_allowed edge_allowed",
    "PredecessorDAG": "graph_uid parent_offsets parents via_edges",
    "Path": "graph_uid hops cost source destination",
    "FlowRecord": "flow_id policy_id sources destinations admission entry_id route min_route_cost offered delivered edge_flows edge_drops edge_reservations",
    "ForwardingProgram": "sources destinations source_weights tunnels",
    "Tunnel": "entry_id head local_index destination edges path_cost status requested_bandwidth reserved_bandwidth planned_bandwidth weight",
    "RouteDefinition": "head route weight",
    "AbsoluteCost": "cost",
    "CostStretch": "factor",
    "TunnelCount": "per_head up_to",
    "ShortestPaths": "allocation source_allocation weight_basis edge_weights tie_break multipath parallel_edges",
    "DistributedTE": "tunnel_count cost_ceiling min_unreserved_bandwidth tie_break parallel_edges",
    "CentralizedTE": "source_allocation source_weights support tunnel_count prefer_headroom",
    "Pinned": "program black_holes",
    "FlowPolicyConfig": "routing admission",
    "OutputOptions": "edge_flows residual saturated_boundary_edges edge_reservations min_cut program by_source",
    "EvaluationOptions": "config eligibility cache outputs design_demand offer",
    "SourceFlow": "source offered admitted delivered source_drops total_edge_drops",
    "TunnelFlow": "entry_id head status assigned delivered total_edge_drops",
    "TunnelPlan": "tunnels_per_head design_demand design_by_head total_reserved_bandwidth heads_with_capacity insufficient_resolution max_fractional_bandwidth_per_tunnel max_bandwidth_per_tunnel planned_bandwidth_per_tunnel total_planned_bandwidth total_tunnel_cost headroom_penalty max_lossless_offer",
    "PlacementResult": "call_delivered call_undelivered termination routes_changed weights_changed reservations_changed delivery_delta call_edge_drops call_source_drops program plan requested offered admitted delivered unadmitted source_drops edge_drops by_source by_tunnel total_cost marginal_cost_curve min_cut edge_flows residual saturated_boundary_edges edge_reservations",
    "ShortestPathCacheStats": "hits misses entries bytes live_values",
    "ShortestPathCache": "graph_uid stats",
    "FlowPolicy": "config program delivered records",
    "FlowLedger": "uid graph edge_capacity_units carried residual reserved unreserved records",
    "Place": "policy offer",
    "Replace": "policy offer design_demand",
    "Forward": "policy offer black_holes",
    "SpfOptions": "direction targets tie_break multipath parallel_edges residual min_residual eligibility",
    "KspOptions": "k tie_break parallel_edges cost_ceiling residual eligibility",
    "SpfResult": "costs dag target_cost",
    "SensitivityEntry": "edge delivery_reduction",
    "SensitivityResult": "baseline entries",
}
FIELDS = {key: value.split() for key, value in FIELDS.items()}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dumps(value):
    return json.dumps(value, separators=(",", ":"), allow_nan=False)


class Blobs:
    def __init__(self, root, write=False, charge=None):
        self.root = Path(root)
        self.charge = charge
        if write:
            self.root.mkdir()
        self.known = set()
        self.lock = threading.RLock()

    def put(self, array):
        array = np.asarray(array)
        if array.dtype.hasobject:
            raise TypeError("object arrays cannot be recorded")
        raw = array.tobytes(order="C")
        key = hashlib.sha256(
            array.dtype.str.encode() + str(array.shape).encode() + raw
        ).hexdigest()
        path = self.root / (key + ".npz")
        with self.lock:
            if key not in self.known:
                if not path.exists():
                    with path.open("xb") as stream:
                        np.savez_compressed(stream, values=array)
                if self.charge is not None:
                    self.charge(path.stat().st_size)
                self.known.add(key)
        return {"$array": key, "dtype": array.dtype.str, "shape": list(array.shape)}

    def get(self, spec):
        with np.load(
            self.root / (spec["$array"] + ".npz"), allow_pickle=False
        ) as archive:
            value = archive["values"]
        if value.dtype.str != spec["dtype"] or list(value.shape) != spec["shape"]:
            raise ValueError("array payload metadata mismatch")
        return value


class Values:
    """Value tags distinguish Python containers, exact floats, enums and arrays."""

    def __init__(self, core, blobs):
        self.core, self.blobs = core, blobs

    def atom(self, value, descend):
        if value is None or type(value) in (bool, int, str):
            return value
        if type(value) is float:
            if np.isnan(value):
                raise ValueError("unexpected NaN")
            return {"$float": struct.pack("!d", value).hex()}
        if type(value) is range:
            return {"$range": [value.start, value.stop, value.step]}
        if isinstance(value, Fraction):
            return {"$fraction": [value.numerator, value.denominator]}
        if isinstance(value, np.generic):
            return {
                "$scalar": [value.dtype.str, base64.b64encode(value.tobytes()).decode()]
            }
        if type(value).__name__ in vars(self.core) and hasattr(
            type(value), "__members__"
        ):
            return {"$enum": [type(value).__name__, value.name, value.value]}
        if type(value) in (tuple, list):
            return {"$" + type(value).__name__: [descend(v) for v in value]}
        if type(value) is dict:
            return {"$dict": [[descend(k), descend(v)] for k, v in value.items()]}
        raise TypeError(f"unhandled value type {type(value).__name__}")

    def literal(self, spec, descend):
        if spec is None or type(spec) in (bool, int, str):
            return spec
        if "$float" in spec:
            return struct.unpack("!d", bytes.fromhex(spec["$float"]))[0]
        if "$range" in spec:
            return range(*spec["$range"])
        if "$fraction" in spec:
            return Fraction(*spec["$fraction"])
        if "$scalar" in spec:
            dtype, raw = spec["$scalar"]
            return np.frombuffer(base64.b64decode(raw), dtype=dtype)[0]
        if "$enum" in spec:
            cls, name, number = spec["$enum"]
            result = getattr(getattr(self.core, cls), name)
            if result.value != number:
                raise ValueError(f"enum value differs: {cls}.{name}")
            return result
        if "$list" in spec:
            return [descend(x) for x in spec["$list"]]
        if "$tuple" in spec:
            return tuple(descend(x) for x in spec["$tuple"])
        if "$dict" in spec:
            return {descend(k): descend(v) for k, v in spec["$dict"]}
        raise ValueError(f"unknown value tag: {list(spec)}")


class Snapshots(Values):
    def __init__(self, core, blobs):
        super().__init__(core, blobs)
        self.uids = {"Graph": {}, "FlowLedger": {}, "Eligibility": {0: "all"}}
        self.previous = {}
        self.owners = {}

    def register(self, obj, oid):
        kind = type(obj).__name__
        if kind in ("FlowPolicy", "FlowLedger", "ShortestPathCache"):
            self.owners[id(obj)] = oid
        if kind in self.uids:
            number = obj.uid()
            self.uids[kind].setdefault(number, oid)

    def uid(self, kind, number):
        if (
            type(number) is not int
            or not 0 <= number < 2**63
            or number == 0
            and kind != "Eligibility"
        ):
            raise ValueError(f"invalid {kind} uid {number!r}")
        if number not in self.uids[kind]:
            raise ValueError(f"unregistered {kind} uid {number}")
        return {"$uid": [kind, self.uids[kind][number]]}

    def snapshot(self, obj, field="", nested=False):
        kind = type(obj).__name__
        if field == "flow_id":
            return (self.uid("FlowLedger", obj[0]), obj[1])
        if field == "graph_uid":
            return self.uid("Graph", obj)
        if isinstance(obj, np.ndarray):
            owner = type(obj.base).__name__ if obj.base is not None else None
            return {
                "$ndarray": obj.copy(),
                "readonly": not obj.flags.writeable,
                "owner_type": owner,
            }
        if field == "edge_flows" and type(obj) is list:
            if not all(map(operator.is_, map(type, obj), repeat(tuple))):
                raise TypeError("edge_flows rows must be tuples")
            if not all(map(operator.eq, map(len, obj), repeat(2))):
                raise TypeError("edge_flows rows must have two entries")
            if not all(
                map(operator.is_, map(type, chain.from_iterable(obj)), repeat(int))
            ):
                raise TypeError("edge_flows values must be Python ints")
            pairs = np.fromiter(
                chain.from_iterable(obj), dtype=np.int64, count=2 * len(obj)
            ).reshape(-1, 2)
            if len(pairs) and (
                np.any(pairs < 0)
                or np.any(pairs[:, 0] >= 2**31)
                or np.any(pairs[:, 1] >= 2**62)
            ):
                raise ValueError("edge_flows exceeds its integer domain")
            return {"$int_pairs": pairs}
        if kind in FIELDS:
            if kind == "Graph" and nested:
                return {"$type": kind, "uid": self.uid(kind, obj.uid())}
            if kind in ("FlowPolicy", "FlowLedger") and nested:
                return {"$owner": [kind, self.owners[id(obj)]]}
            result = {"$type": kind}
            for name in FIELDS[kind]:
                value = getattr(obj, name)
                value = value() if callable(value) else value
                if name == "uid":
                    value = self.uid(kind, value)
                else:
                    value = self.snapshot(value, name, nested=True)
                result[name] = value
            return result
        if type(obj) in (list, tuple):
            return type(obj)(self.snapshot(x, nested=True) for x in obj)
        if type(obj) is dict:
            return {k: self.snapshot(v, nested=True) for k, v in obj.items()}
        # Validate scalar types now; keep them in memory until the delta is made.
        self.atom(obj, lambda x: self.snapshot(x, nested=True))
        return obj

    def pack(self, value):
        if type(value) is dict:
            if "$int_pairs" in value:
                return {"$int_pairs": self.blobs.put(value["$int_pairs"])}
            if "$ndarray" in value:
                return {
                    "$ndarray": self.blobs.put(value["$ndarray"]),
                    "readonly": value["readonly"],
                    "owner_type": value["owner_type"],
                }
            return {k: self.pack(v) for k, v in value.items()}
        return self.atom(value, self.pack)

    def delta(self, old, new, field=""):
        if type(old) is not type(new):
            return {"replace": self.pack(new)}
        if field == "records" and type(old) is list and type(new) is list:

            def keyed(records):
                result = {dumps(self.pack(r["flow_id"])): r for r in records}
                if len(result) != len(records):
                    raise ValueError("duplicate FlowId in records")
                return result

            before, after = keyed(old), keyed(new)
            added = {k: self.pack(v) for k, v in after.items() if k not in before}
            removed = [k for k in before if k not in after]
            changed = {
                k: d
                for k, v in after.items()
                if k in before and (d := self.delta(before[k], v)) is not None
            }
            if added or removed or changed or list(before) != list(after):
                return {
                    "records": {
                        "order": list(after),
                        "added": added,
                        "removed": removed,
                        "changed": changed,
                    }
                }
            return None
        if type(new) is dict:
            if "$int_pairs" in new:
                a, b = old["$int_pairs"], new["$int_pairs"]
                if a.shape != b.shape:
                    return {"replace": self.pack(new)}
                indices = np.flatnonzero(a.reshape(-1) != b.reshape(-1))
                return (
                    {
                        "pair_indices": self.blobs.put(indices),
                        "values": self.blobs.put(b.reshape(-1)[indices]),
                    }
                    if len(indices)
                    else None
                )
            if "$ndarray" in new:
                a, b = old["$ndarray"], new["$ndarray"]
                if (
                    a.dtype != b.dtype
                    or a.shape != b.shape
                    or old["readonly"] != new["readonly"]
                    or old["owner_type"] != new["owner_type"]
                ):
                    return {"replace": self.pack(new)}
                indices = np.flatnonzero(a.reshape(-1) != b.reshape(-1))
                if len(indices):
                    return {
                        "indices": self.blobs.put(indices),
                        "values": self.blobs.put(b.reshape(-1)[indices]),
                    }
                return None
            changes = {}
            for key in old.keys() | new.keys():
                if key not in new:
                    changes[key] = {"removed": True}
                elif key not in old:
                    changes[key] = {"replace": self.pack(new[key])}
                else:
                    item = self.delta(old[key], new[key], key)
                    if item is not None:
                        changes[key] = item
            return {"fields": changes} if changes else None
        if type(new) in (list, tuple):
            if len(old) != len(new):
                return {"replace": self.pack(new)}
            changes = {
                str(i): d
                for i, (a, b) in enumerate(zip(old, new, strict=True))
                if (d := self.delta(a, b)) is not None
            }
            return {"items": changes} if changes else None
        return (
            None
            if self.atom(old, self.pack) == self.atom(new, self.pack)
            else {"replace": self.pack(new)}
        )

    def state(self, oid, obj):
        current = self.snapshot(obj)
        if oid not in self.previous:
            result = {"full": self.pack(current)}
        else:
            result = self.delta(self.previous[oid], current)
        self.previous[oid] = current
        return result


def first_difference(a, b, left_blobs, right_blobs, path="$"):
    if type(a) is not type(b):
        return path, a, b
    if isinstance(a, dict):
        if a.keys() != b.keys():
            return path + ".keys", sorted(a), sorted(b)
        if "$array" in a:
            for key in ("dtype", "shape"):
                if a[key] != b[key]:
                    return path + "." + key, a[key], b[key]
            x, y = left_blobs.get(a), right_blobs.get(b)
            indices = np.flatnonzero(
                x.reshape(-1).view(np.uint8) != y.reshape(-1).view(np.uint8)
            )
            if len(indices):
                index = int(indices[0]) // x.dtype.itemsize
                return (
                    path + f".flat[{index}]",
                    repr(x.reshape(-1)[index]),
                    repr(y.reshape(-1)[index]),
                )
            return None
        for key in sorted(a):
            if result := first_difference(
                a[key], b[key], left_blobs, right_blobs, path + "." + key
            ):
                return result
    elif isinstance(a, list):
        if len(a) != len(b):
            return path + ".length", len(a), len(b)
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            if result := first_difference(
                x, y, left_blobs, right_blobs, path + f"[{i}]"
            ):
                return result
    elif a != b:
        return path, a, b
    return None
