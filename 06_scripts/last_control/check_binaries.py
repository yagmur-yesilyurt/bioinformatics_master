#!/usr/bin/env python3
"""
check_binaries.py  --  diagnose why a3/b2/b3 fail.

a3 (APBS), b2 (HOLE) and b3 (fpocket) need external programs that are NOT
part of MDAnalysis. If a3 says "incomplete" immediately, one of these is
missing. This script tells you exactly which, and how to install it.
"""

import shutil, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

CHECKS = [
    ("pdb2pqr30", C.BIN_PDB2PQR, "a3", "conda install -c conda-forge pdb2pqr"),
    ("apbs",      C.BIN_APBS,     "a3", "conda install -c conda-forge apbs"),
    ("fpocket",   C.BIN_FPOCKET,  "b3", "conda install -c conda-forge fpocket"),
    ("hole",      C.BIN_HOLE,     "b2", "see HOLE2 (http://www.holeprogram.org) or `conda install -c conda-forge hole2`"),
]

def main():
    print("=" * 64)
    print("External binary / dependency check")
    print("=" * 64)
    missing = []
    for label, binname, which, how in CHECKS:
        path = shutil.which(binname)
        if path:
            print(f"  [OK ] {label:<10} ({which}) -> {path}")
        else:
            print(f"  [MISSING] {label:<10} ({which}) -> not on PATH")
            print(f"           install: {how}")
            missing.append((label, which))

    # python deps
    try:
        import gridData  # noqa
        print(f"  [OK ] GridDataFormats (a3) -> {gridData.__version__}")
    except Exception:
        print("  [MISSING] GridDataFormats (a3) -> pip install GridDataFormats")
        missing.append(("GridDataFormats", "a3"))

    print("\n" + "=" * 64)
    if not missing:
        print("All good. a3/b2/b3 should run.")
    else:
        affected = sorted({w for _, w in missing})
        print(f"Missing pieces block: {', '.join(affected)}")
        print("Install them (one line, covers everything for a3/b2/b3):")
        print("  conda install -c conda-forge apbs pdb2pqr fpocket hole2")
        print("  pip install GridDataFormats")
        print("\nNote: a1 and b1 need NONE of these and have already run.")


if __name__ == "__main__":
    main()
