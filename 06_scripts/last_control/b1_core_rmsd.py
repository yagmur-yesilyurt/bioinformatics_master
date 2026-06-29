#!/usr/bin/env python3
"""
b1_core_rmsd.py  --  Katman 0, Task B1  (memory-safe version)

Addresses report section 3.1.

Separates "is the folded core unstable?" from "are the peripheries wagging?":
  1. average structure (streaming, low memory),
  2. align trajectory to it ON DISK (NOT in memory -> no OOM),
  3. per-residue CA-RMSF from the aligned trajectory,
  4. CORE = residues with RMSF below CORE_RMSF_PERCENTILE,
  5. RMSD superposed on the CORE, vs the original global RMSD.

Memory notes:
  - Earlier version used AlignTraj(in_memory=True), which loads all 10000
    frames into RAM and gets OOM-killed. This version writes the aligned
    trajectory to a temp DCD and streams everything. Universes are freed
    between systems.
  - Set ANALYSIS_STEP > 1 to subsample frames if disk/time is tight
    (RMSF/RMSD are robust to modest subsampling).
"""

import os, sys, csv, gc, tempfile
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ANALYSIS_STEP = 1          # use every Nth frame; 1 = all 10000
CA = "protein and name CA"


def core_selection_string(resids):
    if len(resids) == 0:
        return CA
    return CA + " and (resid " + " ".join(str(int(r)) for r in resids) + ")"


def analyze_system(key):
    top, traj, is_mut = C.SYSTEMS[key]
    if not (os.path.exists(top) and os.path.exists(traj)):
        print(f"[skip] {key}: missing files")
        return None
    outdir = C.ensure_out(f"core_rmsd/{key}")
    u = mda.Universe(top, traj)
    n = u.trajectory.n_frames
    print(f"\n=== {key}  (frames={n}, CA={u.select_atoms(CA).n_atoms}, step={ANALYSIS_STEP}) ===", flush=True)

    # 1. average structure (streaming)
    avg = align.AverageStructure(u, u, select=CA, ref_frame=0).run(step=ANALYSIS_STEP)
    ref = avg.results.universe

    # 2. align to average, written to a TEMP DCD on disk (low memory)
    tmp_dcd = os.path.join(outdir, "_aligned_tmp.dcd")
    align.AlignTraj(u, ref, select=CA, filename=tmp_dcd, in_memory=False).run(step=ANALYSIS_STEP)

    # 3. RMSF from the on-disk aligned trajectory (streaming)
    ua = mda.Universe(top, tmp_dcd)
    ca = ua.select_atoms(CA)
    rmsf = rms.RMSF(ca).run().results.rmsf
    resids = ca.resids
    with open(os.path.join(outdir, "rmsf.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["resid", "rmsf_A"])
        for r, f in zip(resids, rmsf): w.writerow([int(r), f"{f:.3f}"])

    # 4. define core
    thr_pct = float(np.percentile(rmsf, C.CORE_RMSF_PERCENTILE))
    core_resids = list(resids[rmsf <= thr_pct])
    core_abs    = list(resids[rmsf <= C.CORE_RMSF_ABS_CUTOFF])
    with open(os.path.join(outdir, "core_residues.txt"), "w") as fh:
        fh.write(f"# percentile core (RMSF <= {thr_pct:.2f} A, p{C.CORE_RMSF_PERCENTILE})\n")
        fh.write(" ".join(map(str, map(int, core_resids))) + "\n")
        fh.write(f"# fixed-cutoff core (RMSF <= {C.CORE_RMSF_ABS_CUTOFF} A)\n")
        fh.write(" ".join(map(str, map(int, core_abs))) + "\n")

    # 5. RMSD (stream over original trajectory)
    Rg = rms.RMSD(u, u, select=CA, ref_frame=0).run(step=ANALYSIS_STEP).results.rmsd
    global_rmsd = Rg[:, 2]
    core_sel = core_selection_string(core_resids)
    Rc = rms.RMSD(u, u, select=core_sel, groupselections=[CA],
                  ref_frame=0).run(step=ANALYSIS_STEP).results.rmsd
    core_rmsd = Rc[:, 2]

    with open(os.path.join(outdir, "rmsd_global_vs_core.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["frame", "global_rmsd", "core_rmsd"])
        for i in range(len(global_rmsd)):
            w.writerow([i, f"{global_rmsd[i]:.3f}", f"{core_rmsd[i]:.3f}"])

    # plots
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(global_rmsd, lw=0.6, label=f"global  {global_rmsd.mean():.2f}\u00b1{global_rmsd.std():.2f}")
    ax.plot(core_rmsd,   lw=0.8, label=f"CORE    {core_rmsd.mean():.2f}\u00b1{core_rmsd.std():.2f}")
    ax.set_xlabel("frame (step %d)" % ANALYSIS_STEP); ax.set_ylabel("RMSD (A)")
    ax.legend(fontsize=8); ax.set_title(f"{key}: global vs core backbone RMSD")
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "rmsd.png"), dpi=130); plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 3.2))
    ax.plot(resids, rmsf, lw=0.6)
    ax.axhline(thr_pct, color="r", ls="--", lw=0.8, label=f"core cutoff p{C.CORE_RMSF_PERCENTILE}={thr_pct:.2f} A")
    ax.set_xlabel("resid"); ax.set_ylabel("CA RMSF (A)"); ax.legend(fontsize=8)
    ax.set_title(f"{key}: per-residue RMSF")
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "rmsf.png"), dpi=130); plt.close(fig)

    res = dict(system=key, n_core=len(core_resids), n_total=ca.n_atoms,
               global_mean=float(global_rmsd.mean()), global_std=float(global_rmsd.std()),
               core_mean=float(core_rmsd.mean()),     core_std=float(core_rmsd.std()))

    # free memory + remove temp aligned trajectory
    del ua, ca, u, Rg, Rc, global_rmsd, core_rmsd, rmsf
    gc.collect()
    try: os.remove(tmp_dcd)
    except OSError: pass
    return res


def main():
    print("=" * 70)
    print("Katman 0  Task B1: core vs global RMSD  (memory-safe)")
    print("=" * 70)
    rows = []
    for key in C.SYSTEMS:
        try:
            r = analyze_system(key)
            if r: rows.append(r)
        except MemoryError:
            print(f"[MemoryError] {key}: raise ANALYSIS_STEP (e.g. 2 or 5) and rerun.")
        except Exception as e:
            print(f"[error] {key}: {e}")
        gc.collect()

    print("\n" + "=" * 70)
    print(f"{'system':<22}{'global RMSD':>16}{'CORE RMSD':>16}{'core/total':>12}")
    print("-" * 66)
    for r in rows:
        print(f"{r['system']:<22}{r['global_mean']:>7.2f}\u00b1{r['global_std']:<6.2f}"
              f"{r['core_mean']:>7.2f}\u00b1{r['core_std']:<6.2f}"
              f"{r['n_core']:>6}/{r['n_total']:<5}")
    print("\nDECISION GATE:")
    print("  CORE RMSD low+flat while global high -> 3.1: high global = periphery,")
    print("    core stable. Rewrite 'stable plateau' as a CORE statement. (-> 4)")
    print("  CORE RMSD also high for holoMutant_prot -> disclose genuine core")
    print("    instability of the central trajectory. (honest -> 4)")


if __name__ == "__main__":
    main()
