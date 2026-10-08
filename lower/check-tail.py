#!/usr/bin/env python3
"""Check the rational rotation-weight dominance used in lower-bound.tex.

Uses only the Python standard library. Run beside grothendieck_1777/,
or pass the path to data/lower_rows.json as the sole argument.
"""

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path


def check(data):
    rows = data["rows"]
    if len(rows) != 35:
        raise ValueError("Expected 35 certificate rows")
    results = []
    for index, row in enumerate(rows):
        angular = row["angular"]
        if row["index"] != index or [a["N"] for a in angular] != [7, 9]:
            raise ValueError(f"Unexpected row or rotation indices in row {index}")
        for item in angular:
            n = item["N"]
            weights = list(map(Fraction, item["weights"]))
            if len(weights) != (n - 1) // 2:
                raise ValueError(f"Wrong number of weights in row {index}, N={n}")
            gap = weights[0] - sum(max(Fraction(0), -w) for w in weights[1:])
            if gap <= 0:
                raise ValueError(f"Dominance failed in row {index}, N={n}: {gap}")
            results.append({"row": index, "N": n, "gap": str(gap)})
    return {"accepted": True, "checks": len(results),
            "minimum_gap": str(min(Fraction(r["gap"]) for r in results)),
            "results": results}


if __name__ == "__main__":
    if len(sys.argv) > 2:
        raise SystemExit("Usage: python check-tail.py [path/to/lower_rows.json]")
    path = Path(sys.argv[1] if len(sys.argv) == 2 else
                "grothendieck_1777/data/lower_rows.json")
    raw = path.read_bytes()
    report = check(json.loads(raw))
    report["input_sha256"] = hashlib.sha256(raw).hexdigest()
    print(json.dumps(report, indent=2))
