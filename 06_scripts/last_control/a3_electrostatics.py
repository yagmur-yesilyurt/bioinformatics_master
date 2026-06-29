#!/usr/bin/env python3
"""
a3_electrostatics.py  --  Katman 0, Task A (electrostatics)

Addresses report sections 3.3 and 3.8 (the "K347+ reshapes the local
electrostatics" claim).

Pipeline, per representative frame (protein only, ligand stripped):
    PDB --pdb2pqr(pH 7.4, PROPKA)--> PQR --apbs--> potential grid (.dx)
then sample the electrostatic potential at the Site-2 centroid (mean CA of
ASP1423 + ASP1714 + VAL1400) and compare WT vs mutant.

Requires external binaries (set paths in config): pdb2pqr30, apbs.
Python dep for grid sampling:  pip install GridDataFormats

DECISION GATE:
  If the mutant potential at the Site-2 centroid is MORE NEGATIVE than WT
  (a deeper cation-attracting well), the K347+ charge measurably reshapes
  the binding-site electrostatics -> supports 3.3 + 3.8.
  If WT and mutant are indistinguishable at the site, the indirect-
  electrostatic mechanism is NOT supported -> rewrite 3.3 honestly.

Compares two pairs by default:
  apoWT vs apoMutant            (cleanest, no ligand)
  holoWT_prot vs holoMutant_prot
"""

import os, sys, subprocess, tempfile, shutil
import numpy as np
import MDAnalysis as mda

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

PAIRS = [("apoWT", "apoMutant"), ("holoWT_prot", "holoMutant_prot")]


def preflight():
    """Return True only if every external piece a3 needs is present.
    Prints a precise diagnosis (this is what 'incomplete' was hiding)."""
    ok = True
    for name, binname in (("pdb2pqr30", C.BIN_PDB2PQR), ("apbs", C.BIN_APBS)):
        if shutil.which(binname) is None:
            print(f"  [MISSING BINARY] {name} not on PATH "
                  f"-> conda install -c conda-forge {'pdb2pqr' if 'pdb2pqr' in name else 'apbs'}")
            ok = False
    try:
        import gridData  # noqa
    except Exception:
        print("  [MISSING PYTHON DEP] GridDataFormats -> pip install GridDataFormats")
        ok = False
    return ok

APBS_TEMPLATE = """read
    mol pqr {pqr}
end
elec
    mg-auto
    dime 97 97 97
    cglen {cglen:.1f} {cglen:.1f} {cglen:.1f}
    fglen {fglen:.1f} {fglen:.1f} {fglen:.1f}
    cgcent mol 1
    fgcent mol 1
    mol 1
    lpbe
    bcfl sdh
    ion charge 1 conc 0.150 radius 2.0
    ion charge -1 conc 0.150 radius 2.0
    pdie 2.0
    sdie 78.54
    srfm smol
    chgm spl2
    sdens 10.0
    srad 1.4
    swin 0.3
    temp 310.0
    calcenergy no
    calcforce no
    write pot dx {dxout}
end
quit
"""


def strip_protein(in_pdb, out_pdb):
    u = mda.Universe(in_pdb)
    prot = u.select_atoms("protein")
    prot.write(out_pdb)
    return u


def run(cmd, cwd=None, **kw):
    print("   $", " ".join(cmd), (f"   (cwd={cwd})" if cwd else ""))
    try:
        return subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=cwd, **kw)
    except subprocess.CalledProcessError as e:
        tail = (e.stderr or e.stdout or "").strip().splitlines()[-10:]
        raise RuntimeError(f"{cmd[0]} failed (exit {e.returncode}):\n      "
                           + "\n      ".join(tail))


def find_dx(workdir, stem):
    """APBS may write the .dx relative to CWD or with/without the abs path.
    Search the likely locations and return the first existing file."""
    import glob
    cands = [os.path.join(workdir, stem + ".dx"),
             os.path.join(os.getcwd(), stem + ".dx")]
    for c in cands:
        if os.path.exists(c):
            return c
    hits = glob.glob(os.path.join(workdir, stem + "*.dx")) + glob.glob(os.path.join(os.getcwd(), stem + "*.dx"))
    if hits:
        return max(hits, key=os.path.getmtime)
    return None


def site2_centroid(pdb):
    u = mda.Universe(pdb)
    sel = u.select_atoms("name CA and (resid " + " ".join(map(str, C.SITE2_RESIDS)) + ")")
    if sel.n_atoms == 0:
        raise RuntimeError("Site-2 residues not found for centroid")
    return sel.center_of_geometry()


def sample_dx(dx_path, point):
    from gridData import Grid
    g = Grid(dx_path)
    return float(g.interpolated(point[0], point[1], point[2]))


def process(key, workdir):
    src = C.REP_FRAMES.get(key)
    if not src or not os.path.exists(src):
        print(f"[skip] {key}: representative frame missing ({src})")
        return None
    base = os.path.join(workdir, key)
    prot_pdb = base + "_prot.pdb"
    pqr = base + ".pqr"
    apbs_in = base + ".in"
    stem = f"{key}_pot"                      # RELATIVE write stem (apbs runs in workdir)

    strip_protein(src, prot_pdb)
    run([C.BIN_PDB2PQR, "--ff=AMBER", "--with-ph=7.4", "--titration-state-method=propka",
         prot_pdb, pqr])

    u = mda.Universe(pqr)
    ext = u.atoms.positions.max(0) - u.atoms.positions.min(0)
    cglen = float(max(ext) * 1.7 + 20)
    fglen = float(max(ext) * 1.2 + 10)
    with open(apbs_in, "w") as fh:
        # reference pqr and output by BASENAME; apbs is launched with cwd=workdir
        fh.write(APBS_TEMPLATE.format(pqr=os.path.basename(pqr), dxout=stem,
                                      cglen=cglen, fglen=fglen))

    res = run([C.BIN_APBS, os.path.basename(apbs_in)], cwd=workdir)
    dx = find_dx(workdir, stem)
    if dx is None:
        tail = (res.stdout or "").strip().splitlines()[-12:]
        raise RuntimeError("APBS produced no .dx (searched workdir and CWD). "
                           "APBS said:\n      " + "\n      ".join(tail))

    pt = site2_centroid(prot_pdb)
    val = sample_dx(dx, pt)   # units: kT/e (APBS default)
    print(f"   {key}: Site-2 centroid potential = {val:+.3f} kT/e   (dx={os.path.basename(dx)})")
    return val


def main():
    workdir = C.ensure_out("electrostatics")
    print("=" * 70)
    print("Katman 0  Task A3: APBS electrostatic potential at Site 2")
    print("=" * 70)
    if not preflight():
        print("\n>>> a3 cannot run until the items above are installed.")
        print(">>> One line:  conda install -c conda-forge apbs pdb2pqr ;  pip install GridDataFormats")
        print(">>> NOTE: given the a1 result (residue 347 ~26 A from Site 2),")
        print(">>>       a3 is now CONFIRMATORY, not load-bearing. You are not blocked.")
        return
    results = {}
    for key in set(sum([list(p) for p in PAIRS], [])):
        try:
            results[key] = process(key, workdir)
        except Exception as e:
            print(f"[error] {key}: {e}")
            results[key] = None

    print("\n" + "=" * 70)
    print("DECISION GATE  (mutant should be MORE NEGATIVE if mechanism holds)")
    print("=" * 70)
    for wt, mut in PAIRS:
        vw, vm = results.get(wt), results.get(mut)
        if vw is None or vm is None:
            print(f"  {wt} vs {mut}: incomplete")
            continue
        d = vm - vw
        verdict = ("SUPPORTS mechanism (mutant deeper negative well)"
                   if d < -1.0 else
                   "NO clear difference -> rewrite 3.3 as negative result"
                   if abs(d) <= 1.0 else
                   "OPPOSITE of expected (mutant less negative) -> investigate")
        print(f"  {wt}={vw:+.2f}  {mut}={vm:+.2f}  delta(mut-WT)={d:+.2f} kT/e  -> {verdict}")


if __name__ == "__main__":
    main()
