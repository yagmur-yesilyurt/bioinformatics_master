#!/usr/bin/env python3
"""
a3b_charge_isolation.py  --  the CONTROLLED version of a3.

Why this exists: a3 compared two DIFFERENT single frames (WT medoid vs mutant
medoid), so any potential difference at Site 2 reflects global structural +
PROPKA-protonation differences, NOT the K347 charge. Worse, a positive charge
makes a distant point MORE POSITIVE, so "mutant more negative" cannot be the
direct K347 field at all.

This script isolates the K347 charge cleanly. On ONE mutant frame it builds two
PQR files with IDENTICAL coordinates and IDENTICAL charges everywhere EXCEPT
the K347 side chain:
    PQR_charged   : K347 side chain keeps its +1 (real Lys).
    PQR_neutral   : K347 side-chain charges zeroed (charge removed, geometry
                    untouched).
APBS on both, sampled at the Site-2 centroid. The difference is the PURE
electrostatic contribution of the K347 side-chain charge at Site 2, with
everything else held fixed.

PREDICTION from a1 (K347 ~26 A from Site 2; Debye length ~8 A at 0.15 M):
  the isolated K347 contribution at Site 2 should be SMALL and, if anything,
  POSITIVE (a +1 charge raises the potential / makes it less negative). A small
  or positive value CONFIRMS that the local-electrostatic mechanism does not
  operate; it cannot deepen a cation-attracting well at Site 2.

Frame used: the mutant frame (has K347). Default holoMutant_prot medoid;
override with --system.
"""

import os, sys, argparse, shutil
import numpy as np
import MDAnalysis as mda

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
import a3_electrostatics as a3   # reuse APBS_TEMPLATE, run, find_dx, sample_dx, site2_centroid, preflight

# Lys side-chain atoms whose charges are zeroed to neutralise K347 (+1).
LYS_SIDECHAIN = {"CB", "HB1", "HB2", "HB3", "CG", "HG1", "HG2", "HG3",
                 "CD", "HD1", "HD2", "HD3", "CE", "HE1", "HE2", "HE3",
                 "NZ", "HZ1", "HZ2", "HZ3"}


def neutralise_k347(pqr_in, pqr_out, resid=347):
    """Copy a PQR, setting the charge of K347 side-chain atoms to 0.000.
    PQR (pdb2pqr) is whitespace-delimited:
      ATOM serial name resname chainID resid x y z charge radius
    """
    n_zeroed = 0
    with open(pqr_in) as fi, open(pqr_out, "w") as fo:
        for line in fi:
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                # locate resid + atom name robustly
                try:
                    name = parts[2]
                    # resid may be parts[4] (with chain) or parts[3] (no chain)
                    rid = None
                    for idx in (4, 3, 5):
                        if idx < len(parts):
                            try:
                                rid = int(parts[idx]); break
                            except ValueError:
                                continue
                except Exception:
                    fo.write(line); continue
                if rid == resid and name in LYS_SIDECHAIN:
                    # charge is the 2nd-to-last field; rewrite it to 0.0000
                    charge = parts[-2]
                    # replace only the last occurrence of that token before radius
                    head, _, tail = line.rpartition(charge)
                    line = head + f"{0.0:.4f}" + tail
                    n_zeroed += 1
            fo.write(line)
    return n_zeroed


def apbs_potential(tag, pqr, workdir):
    u = mda.Universe(pqr)
    ext = u.atoms.positions.max(0) - u.atoms.positions.min(0)
    cglen = float(max(ext) * 1.7 + 20); fglen = float(max(ext) * 1.2 + 10)
    apbs_in = os.path.join(workdir, f"{tag}.in")
    stem = f"{tag}_pot"
    with open(apbs_in, "w") as fh:
        fh.write(a3.APBS_TEMPLATE.format(pqr=os.path.basename(pqr), dxout=stem,
                                         cglen=cglen, fglen=fglen))
    res = a3.run([C.BIN_APBS, os.path.basename(apbs_in)], cwd=workdir)
    dx = a3.find_dx(workdir, stem)
    if dx is None:
        tail = (res.stdout or "").strip().splitlines()[-12:]
        raise RuntimeError(f"APBS produced no .dx for {tag}:\n      " + "\n      ".join(tail))
    return dx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", default="holoMutant_prot",
                    help="mutant system whose medoid frame to use (must have K347)")
    args = ap.parse_args()

    print("=" * 70)
    print("Katman 0  Task A3b: ISOGEOMETRIC isolation of the K347 charge at Site 2")
    print("=" * 70)
    if not a3.preflight():
        print("\n>>> install the missing piece(s) above, then rerun.")
        return

    key = args.system
    src = C.REP_FRAMES.get(key)
    if not src or not os.path.exists(src):
        print(f"[error] {key}: representative frame missing ({src})"); return

    workdir = C.ensure_out("electrostatics_isolation")
    prot_pdb = os.path.join(workdir, f"{key}_prot.pdb")
    a3.strip_protein(src, prot_pdb)

    # sanity: confirm residue 347 is LYS in this frame
    u = mda.Universe(prot_pdb)
    r347 = u.select_atoms("resid 347")
    rn = r347.residues[0].resname if r347.n_atoms else "ABSENT"
    print(f"frame = {key} medoid;  resid 347 -> {rn}")
    if rn != "LYS":
        print("[error] this isolation requires a mutant frame (K347). "
              "Pass --system apoMutant or holoMutant_prot."); return

    pqr_charged = os.path.join(workdir, f"{key}_K347charged.pqr")
    a3.run([C.BIN_PDB2PQR, "--ff=AMBER", "--with-ph=7.4",
            "--titration-state-method=propka", prot_pdb, pqr_charged])

    pqr_neutral = os.path.join(workdir, f"{key}_K347neutral.pqr")
    n = neutralise_k347(pqr_charged, pqr_neutral)
    print(f"neutralised {n} K347 side-chain atom charges (geometry untouched)")
    if n == 0:
        print("[error] no K347 side-chain atoms found to neutralise; check PQR naming.")
        return

    dx_c = apbs_potential(f"{key}_charged", pqr_charged, workdir)
    dx_n = apbs_potential(f"{key}_neutral", pqr_neutral, workdir)

    pt = a3.site2_centroid(prot_pdb)
    phi_c = a3.sample_dx(dx_c, pt)
    phi_n = a3.sample_dx(dx_n, pt)
    contrib = phi_c - phi_n      # pure K347 side-chain contribution at Site 2

    print("\n" + "=" * 70)
    print("RESULT  (pure K347 charge contribution at the Site-2 centroid)")
    print("=" * 70)
    print(f"  phi(Site2) with K347 charge ON  = {phi_c:+.3f} kT/e")
    print(f"  phi(Site2) with K347 charge OFF = {phi_n:+.3f} kT/e")
    print(f"  isolated K347 contribution      = {contrib:+.3f} kT/e")
    print("\nInterpretation:")
    if abs(contrib) < 0.5:
        print("  |contribution| < 0.5 kT/e -> K347's charge has NEGLIGIBLE effect at")
        print("  Site 2. Confirms a1: the local-electrostatic mechanism does not operate.")
    elif contrib > 0:
        print("  contribution POSITIVE -> K347 makes Site 2 LESS negative (as expected")
        print("  for a +1 charge). It cannot deepen a cation-attracting well; the")
        print("  report's 'K347 creates a negative well at Site 2' mechanism is refuted.")
    else:
        print("  contribution NEGATIVE and sizeable -> unexpected for a +1 charge at")
        print("  distance; inspect (PROPKA reassignment near Site 2 in this frame?).")
    print("\nNote: any WT-vs-mutant Site-2 difference seen in a3 that is NOT reproduced")
    print("here comes from global structural/protonation differences between frames,")
    print("not from the K347 charge -> candidate allosteric or single-frame artefact.")


if __name__ == "__main__":
    main()
