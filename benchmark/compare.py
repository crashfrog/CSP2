#!/usr/bin/env python3
"""
Compare two CSP2 equivalence manifests (from normalize.py) and print a verdict.

EQUIVALENT  -> every canonicalized scientific output matches (same files, same hashes).
DIFFERENT   -> lists added/removed files and per-file hash mismatches.

This is the gate for the Nextflow/Python refactor PRs: a refactor that only restructures
orchestration/code MUST leave this reporting EQUIVALENT against the frozen baseline. When
the Phraya rearchitecture begins, mismatches are expected and this report becomes the
diff to reason about (aided by --show to dump the differing canonical content).

Usage:
    python compare.py baseline/manifest.json candidate/manifest.json
    python compare.py base.json cand.json --show snp_distance_matrix   # dump canon diff
Exit code 0 == equivalent, 1 == different, 2 == usage error.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path


def load(p: Path) -> dict:
    return json.loads(Path(p).read_text())


def entries(m: dict) -> dict:
    return m.get("entries", {})


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("baseline", type=Path)
    ap.add_argument("candidate", type=Path)
    ap.add_argument("--show", metavar="SUBSTR", default=None,
                    help="for mismatched files whose path contains SUBSTR, show which "
                         "side differs (hashes + line counts)")
    args = ap.parse_args()

    base, cand = entries(load(args.baseline)), entries(load(args.candidate))
    bkeys, ckeys = set(base), set(cand)

    removed = sorted(bkeys - ckeys)
    added = sorted(ckeys - bkeys)
    common = sorted(bkeys & ckeys)
    mismatched = [k for k in common
                  if base[k].get("sha256") != cand[k].get("sha256")]

    ok = not (removed or added or mismatched)
    print(f"baseline : {args.baseline}  ({len(base)} files)")
    print(f"candidate: {args.candidate}  ({len(cand)} files)")
    print("-" * 60)
    if ok:
        print(f"EQUIVALENT — {len(common)} canonicalized outputs identical.")
        sys.exit(0)

    print("DIFFERENT")
    if removed:
        print(f"\n  MISSING in candidate ({len(removed)}):")
        for k in removed:
            print(f"    - {k}")
    if added:
        print(f"\n  EXTRA in candidate ({len(added)}):")
        for k in added:
            print(f"    + {k}")
    if mismatched:
        print(f"\n  CONTENT DIFFERS ({len(mismatched)}):")
        for k in mismatched:
            b, c = base[k], cand[k]
            print(f"    ~ {k}  [{b.get('kind')}]  "
                  f"{b.get('sha256','?')[:12]}({b.get('lines','?')}ln) -> "
                  f"{c.get('sha256','?')[:12]}({c.get('lines','?')}ln)")
    if args.show:
        hits = [k for k in mismatched if args.show in k]
        print(f"\n  --show '{args.show}': {len(hits)} matching mismatched file(s)")
        for k in hits:
            print(f"    {k}: baseline={base[k]} candidate={cand[k]}")
    sys.exit(1)


if __name__ == "__main__":
    main()
