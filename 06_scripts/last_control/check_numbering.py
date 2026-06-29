#!/usr/bin/env python3
"""
check_numbering.py  --  RUN THIS FIRST.

The single biggest risk in the whole toolkit is residue-numbering offset: if
CHARMM-GUI/OpenMM renumbered residues (e.g. starting from 1, or dropping
unresolved loops), then "resid 347" in your topology is NOT N347/K347 and
every downstream result would be silently wrong.

This script loads each system topology and:
  1. prints the resname found at each key resid (347, 356, 1419, 1423, 1714),
  2. checks them against the expected resnames (LYS/ASN at 347, ASP at the
     acidics, LYS at the DEKA K1419),
  3. if they DON'T match, automatically searches for a constant integer offset
     that makes all of them line up (the typical renumbering case) and tells
     you exactly what to do.

No external binaries needed.
"""

import os, sys
import MDAnalysis as mda

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C


def resname_map(u):
    """resid -> resname (first residue with that resid, protein only)."""
    m = {}
    prot = u.select_atoms("protein")
    for res in prot.residues:
        m.setdefault(int(res.resid), res.resname.strip().upper())
    return m


def ok(resid, resname, is_mut_pos347):
    if resname is None:
        return False
    if resid == 347:
        return resname in (("LYS",) if is_mut_pos347 else ("ASN",)) or resname in ("LYS", "ASN")
    return resname in C.EXPECTED.get(resid, ())


def find_offset(rmap, is_mut):
    """Search for a constant offset o such that resid+o gives expected resnames."""
    keys = [356, 1419, 1423, 1714]                    # offset-invariant identities
    want = {356: {"ASP"}, 1419: {"LYS"}, 1423: {"ASP"}, 1714: {"ASP"}}
    want347 = {"LYS"} if is_mut else {"ASN"}
    present = set(rmap.keys())
    if not present:
        return None
    lo, hi = min(present) - max(keys), max(present) - min(keys)
    hits = []
    for o in range(lo, hi + 1):
        if all(rmap.get(k + o) in want[k] for k in keys) and rmap.get(347 + o) in want347:
            hits.append(o)
    if len(hits) == 1:
        return hits[0]
    return hits  # 0, or >1 ambiguous


def check_system(key):
    top, traj, is_mut = C.SYSTEMS[key]
    if not os.path.exists(top):
        print(f"[skip] {key}: topology not found -> {top}")
        return None
    u = mda.Universe(top)
    rmap = resname_map(u)
    print(f"\n=== {key}   (is_mutant={is_mut};  topology={os.path.basename(top)}) ===")
    allgood = True
    for r in (347, 356, 1419, 1423, 1714):
        nm = rmap.get(r)
        good = ok(r, nm, is_mut)
        allgood &= good
        print(f"   resid {r:>5} -> {str(nm):>6}   {'OK' if good else 'MISMATCH'}")
    if allgood:
        print("   VERDICT: numbering OK -- run the a*/b* scripts as-is.")
        return ("ok", key, 0)

    print("   VERDICT: numbering MISMATCH -- searching for a constant offset...")
    off = find_offset(rmap, is_mut)
    if isinstance(off, int):
        print(f"   >>> Found offset = {off:+d}.")
        print(f"   >>> i.e. your file's resid (347{off:+d})={347+off} is the real residue 347.")
        print(f"   >>> FIX: in config.py set RES_347={347+off}, RES_D356={356+off}, "
              f"RES_K1419={1419+off}, RES_D1423={1423+off}, RES_D1714={1714+off},")
        print(f"            and SITE2_RESIDS=[{1423+off}, {1714+off}, {1400+off}].")
        return ("offset", key, off)
    elif off == [] or off is None:
        print("   >>> No single offset reproduces all expected residues.")
        print("   >>> Your topology may have non-contiguous numbering or missing")
        print("   >>> residues. Inspect manually: open analysis_result.pdb and find")
        print("   >>> the Lys of the DEKA filter (DIII) and the two binding Asp's.")
        return ("manual", key, None)
    else:
        print(f"   >>> Multiple candidate offsets {off}; ambiguous, inspect manually.")
        return ("ambiguous", key, off)


def main():
    print("=" * 70)
    print("Katman 0  --  RESIDUE NUMBERING CHECK (run before everything else)")
    print(f"REPO_ROOT = {C.REPO_ROOT}")
    print("=" * 70)
    results = [check_system(k) for k in C.SYSTEMS]
    results = [r for r in results if r]
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    statuses = {}
    for status, key, off in results:
        statuses.setdefault(status, []).append((key, off))
    if statuses.get("ok") and len(statuses) == 1:
        print("All systems: numbering OK. Proceed with b1 -> a1 -> a3 -> b2 -> b3.")
    else:
        for status, items in statuses.items():
            for key, off in items:
                extra = f" (offset {off:+d})" if isinstance(off, int) and off else ""
                print(f"  {key:<22} {status}{extra}")
        if statuses.get("offset"):
            offs = {o for _, o in statuses["offset"]}
            if len(offs) == 1:
                print(f"\nConsistent offset {offs.pop():+d} across systems -> apply the "
                      "config.py change printed above, then re-run this checker to confirm.")
            else:
                print("\nWARNING: different offsets per system -> handle each system's "
                      "residue ids separately; do not assume a global offset.")


if __name__ == "__main__":
    main()
