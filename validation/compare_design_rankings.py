"""Compare independently supplied fixed-lift drag changes, with absolute bounds.

No aerodynamic calculation, curve digitization, interpolation, or fitting is
performed here. A successful comparison is not a solver validation certificate.
See INDEPENDENT_RANKING.md for the data contract and current source audit.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path


class RankingInputError(ValueError):
    pass


def require(ok, message):
    if not ok:
        raise RankingInputError(message)


def finite(value, name, positive=False):
    require(type(value) in (int, float) and math.isfinite(value),
            f"{name} must be a finite number")
    require(not positive or value > 0, f"{name} must be positive")
    return float(value)


def direction(interval):
    if interval[1] < 0:
        return "lower_drag"
    if interval[0] > 0:
        return "higher_drag"
    return "unresolved"


def compare(data):
    require(data.get("schema_version") == 1, "unsupported schema_version")
    require(data.get("reference_kind") in ("experiment", "mses", "rans"),
            "reference_kind must distinguish experiment, mses, and rans")
    for key in ("geometry_audited", "conditions_audited", "bounds_audited"):
        require(data.get(key) is True, f"benchmark blocked: {key} is not true")
    for key in ("source", "geometry_audit", "conditions_audit", "bounds_audit"):
        require(isinstance(data.get(key), str) and data[key].strip(), f"missing {key}")
    target = finite(data.get("target_cl"), "target_cl")
    tolerance = finite(data.get("cl_tolerance"), "cl_tolerance", True)
    rows = data.get("geometries")
    require(isinstance(rows, list) and len(rows) >= 2, "need at least two geometries")
    seen = set()
    for row in rows:
        name = row.get("id")
        require(isinstance(name, str) and name and name not in seen,
                "geometry ids must be nonempty and unique")
        seen.add(name)
        digest = row.get("geometry_sha256", "")
        require(isinstance(digest, str) and len(digest) == 64 and
                all(c in "0123456789abcdef" for c in digest),
                f"{name}: missing geometry-file hash")
        require(isinstance(row.get("geometry_file"), str) and row["geometry_file"],
                f"{name}: missing geometry_file")
        for kind in ("reference", "model"):
            point = row.get(kind, {})
            cl = finite(point.get("cl"), f"{name}/{kind}/cl")
            require(abs(cl-target) <= tolerance,
                    f"{name}/{kind}: not at target CL; interpolation is not implicit")
            finite(point.get("cd"), f"{name}/{kind}/cd", True)
            finite(point.get("cd_abs_bound"), f"{name}/{kind}/cd_abs_bound", True)
            require(isinstance(point.get("provenance"), str) and point["provenance"].strip(),
                    f"{name}/{kind}: missing provenance")
        require(row["model"].get("converged") is True, f"{name}: failed model solve")
        require(row["model"].get("local_supersonic") is False,
                f"{name}: local-supersonic diagnostic missing or invalid")
    pairs = []
    for a, b in itertools.combinations(rows, 2):
        result = {"a": a["id"], "b": b["id"], "difference": "B minus A"}
        for kind in ("reference", "model"):
            delta = b[kind]["cd"] - a[kind]["cd"]
            bound = a[kind]["cd_abs_bound"] + b[kind]["cd_abs_bound"]
            interval = [delta-bound, delta+bound]
            result[kind] = {"delta_cd": delta, "delta_drag_counts": delta*1e4,
                            "absolute_bound": bound, "interval": interval,
                            "direction": direction(interval)}
        rd, md = result["reference"]["direction"], result["model"]["direction"]
        result["ranking"] = ("unresolved" if "unresolved" in (rd, md) else
                             "agrees" if rd == md else "disagrees")
        pairs.append(result)
    return {"schema_version": 1, "reference_kind": data["reference_kind"],
            "source": data["source"], "target_cl": target,
            "scope": "Pairwise agreement conditional on supplied audits and bounds; "
                     "not independent validation of the entire solver or optimum.",
            "bound_rule": "Sum of absolute per-geometry bounds; no assumed cancellation "
                          "or statistical independence. Intervals are not confidence intervals.",
            "pairs": pairs}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    try:
        raw = args.input.read_bytes()
        data = json.loads(raw)
        result = compare(data)
        for row in data["geometries"]:
            path = args.input.parent / row["geometry_file"]
            require(hashlib.sha256(path.read_bytes()).hexdigest() == row["geometry_sha256"],
                    f"{row['id']}: geometry file hash mismatch")
        result["input_sha256"] = hashlib.sha256(raw).hexdigest()
    except (RankingInputError, OSError, ValueError, TypeError, AttributeError) as exc:
        result = {"status": "blocked", "reason": str(exc), "pairs": []}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps(result, indent=2, allow_nan=False))
    return 2 if result.get("status") == "blocked" else 0


if __name__ == "__main__":
    raise SystemExit(main())
