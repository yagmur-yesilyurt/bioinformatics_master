#!/usr/bin/env python3
"""
b3_druggability_distribution.py  --  Katman 0, Task B3

Addresses report section 3.4.

The report reports a single fpocket druggability per system (a value that
reached 1.000), computed from ONE representative frame. This script runs
fpocket on an evenly spaced subsample of frames and reports the
druggability of the PHE1760-containing pocket as a DISTRIBUTION, so values
near 1.000 are seen in the context of their frame-to-frame spread rather
than presented as a single high score.

Requires the fpocket executable (config.BIN_FPOCKET).

DECISION GATE:
  If the per-frame druggability of the PHE1760 pocket has a wide spread
  (e.g. std > ~0.15), the single-frame 1.000 is a frame artefact and 3.4
  must report the distribution. If it is consistently high across frames,
  the high druggability is robust. Either way 3.4 becomes defensible. -> 4
"""

import os, sys, glob, subprocess, tempfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PHE1760 = 1760


def write_frame(u, frame_idx, out_pdb):
    u.trajectory[frame_idx]
    u.select_atoms("protein").write(out_pdb)


def run_fpocket(pdb):
    subprocess.run([C.BIN_FPOCKET, "-f", pdb], check=True,
                   capture_output=True, text=True)
    return pdb[:-4] + "_out"


def druggability_of_phe1760(out_dir):
    """Find the pocket whose atoms include PHE1760, return its druggability."""
    info = glob.glob(os.path.join(out_dir, "*_info.txt"))
    pock_dir = os.path.join(out_dir, "pockets")
    if not info or not os.path.isdir(pock_dir):
        return None
    # which pocket contains resid 1760?
    target = None
    for atm in sorted(glob.glob(os.path.join(pock_dir, "pocket*_atm.pdb"))):
        with open(atm) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        resseq = int(line[22:26])
                    except ValueError:
                        continue
                    if resseq == PHE1760:
                        # pocketN_atm.pdb -> N
                        bn = os.path.basename(atm)
                        target = int(bn.replace("pocket", "").split("_")[0])
                        break
        if target is not None:
            break
    if target is None:
        return None
    # parse druggability for "Pocket <target>" from info file
    drug = None
    with open(info[0]) as fh:
        cur = None
        for line in fh:
            s = line.strip()
            if s.lower().startswith("pocket"):
                try:
                    cur = int(s.split()[1])
                except (IndexError, ValueError):
                    cur = None
            if cur == target and "druggability" in s.lower():
                drug = float(s.split(":")[-1].strip())
                break
    return drug


def main():
    import MDAnalysis as mda
    outdir = C.ensure_out("druggability")
    print("=" * 70)
    print("Katman 0  Task B3: multi-frame PHE1760-pocket druggability")
    print("=" * 70)
    for key, (top, traj, _) in C.SYSTEMS.items():
        if not (os.path.exists(top) and os.path.exists(traj)):
            print(f"[skip] {key}: missing files"); continue
        u = mda.Universe(top, traj)
        n = u.trajectory.n_frames
        idxs = np.linspace(0, n - 1, min(C.FPOCKET_N_FRAMES, n)).astype(int)
        vals = []
        with tempfile.TemporaryDirectory() as tmp:
            for fi in idxs:
                pdb = os.path.join(tmp, f"{key}_{fi}.pdb")
                try:
                    write_frame(u, int(fi), pdb)
                    od = run_fpocket(pdb)
                    d = druggability_of_phe1760(od)
                    if d is not None:
                        vals.append(d)
                except Exception as e:
                    print(f"   [warn] frame {fi}: {e}")
        if not vals:
            print(f"  {key}: no PHE1760 pocket found in sampled frames"); continue
        vals = np.array(vals)
        np.savetxt(os.path.join(outdir, f"{key}_druggability.csv"), vals,
                   header="phe1760_pocket_druggability", comments="")
        spread = "ROBUST" if vals.std() <= 0.15 else "FRAME-DEPENDENT"
        print(f"  {key:<22} drugg = {vals.mean():.3f} +/- {vals.std():.3f}  "
              f"(n={len(vals)}, min={vals.min():.3f}, max={vals.max():.3f})  -> {spread}")

        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.hist(vals, bins=20, range=(0, 1), alpha=0.8)
        ax.axvline(vals.mean(), color="k", lw=1)
        ax.set_xlabel("PHE1760 pocket druggability"); ax.set_ylabel("count")
        ax.set_title(f"{key}: druggability over {len(vals)} frames")
        fig.tight_layout(); fig.savefig(os.path.join(outdir, f"{key}_druggability.png"), dpi=130)
        plt.close(fig)

    print("\nDECISION GATE:")
    print("  std <= 0.15  -> high druggability is robust across frames.")
    print("  std  > 0.15  -> the single-frame 1.000 was a frame artefact;")
    print("                 report the distribution in 3.4.")


if __name__ == "__main__":
    main()
