#!/usr/bin/env python3
"""
a1_distances_filter.py  --  Katman 0, Task A (distances)

Addresses report sections 3.3 and 3.8.

For each system it measures, per frame:
  (1) K347(NZ) <-> ASP1423 carboxyl   minimum distance   [salt-bridge style]
  (2) K347(NZ) <-> ASP1714 carboxyl   minimum distance
  (3) K347(NZ) <-> D356   carboxyl    minimum distance   [known salt bridge]
  (4) CA-CA geometry of {347, D356, K1419, D1423, D1714}  -> a static
      distance matrix on the representative frame, to map Site 2 relative
      to the DEKA filter (K1419).

In WT (residue 347 = ASN, uncharged) the script automatically falls back to
the side-chain amide N (ND2) and clearly flags that any distance is a
PROXIMITY, not a salt bridge.

DECISION GATE (printed at the end):
  - If min(K347-ASP1423, K347-ASP1714) < GATE_NEAR (~12 A) for a large
    fraction of frames in the protonated mutant, the DI charge is spatially
    close enough to the DIII/DIV acidic residues to support the indirect
    electrostatic mechanism (3.3) and the selectivity-filter-vestibule
    identity (3.8).
  - If it is NOT close, rewrite 3.3 as "indirect mechanism not supported by
    through-space proximity" -- a correct negative result.

Outputs (under OUT_ROOT/distances/<system>/):
  *_timeseries.csv, *_distances.png, and a printed summary table.
"""

import os, sys, csv
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances as mdadist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def cationic_atoms(u, resid):
    """Return the AtomGroup of the side-chain charged/polar tip for resid,
    plus a label describing whether it is a real cation."""
    res = u.select_atoms(f"resid {resid}")
    if res.n_atoms == 0:
        return None, f"resid {resid} ABSENT"
    resname = res.residues[0].resname.strip().upper()
    if resname in ("LYS", "LSN"):
        ag = u.select_atoms(f"protein and resid {resid} and name NZ")
        return ag, f"{resname}{resid} (cation, NZ)"
    if resname == "ARG":
        ag = u.select_atoms(f"protein and resid {resid} and name NH1 NH2 NE")
        return ag, f"{resname}{resid} (cation, guanidinium)"
    if resname in ("ASN", "ASN1"):
        ag = u.select_atoms(f"protein and resid {resid} and name ND2 OD1")
        return ag, f"{resname}{resid} (NEUTRAL amide -> proximity only)"
    # generic fallback: side-chain N/O
    ag = u.select_atoms(f"protein and resid {resid} and not backbone and (name N* or name O*)")
    return ag, f"{resname}{resid} (fallback polar tip)"


def carboxyl_atoms(u, resid):
    res = u.select_atoms(f"resid {resid}")
    if res.n_atoms == 0:
        return None, f"resid {resid} ABSENT"
    resname = res.residues[0].resname.strip().upper()
    ag = u.select_atoms(f"protein and resid {resid} and name OD1 OD2 OE1 OE2")
    if ag.n_atoms == 0:  # not Asp/Glu; use any side-chain O
        ag = u.select_atoms(f"protein and resid {resid} and not backbone and name O*")
    return ag, f"{resname}{resid} (carboxyl)"


def min_dist_series_multi(u, pairs):
    """Single trajectory pass; pairs = {name: (agA, agB)}.
    Returns {name: np.array}. Skips pairs whose groups are empty."""
    valid = {k: v for k, v in pairs.items()
             if v[0] is not None and v[0].n_atoms and v[1] is not None and v[1].n_atoms}
    out = {k: np.empty(u.trajectory.n_frames) for k in valid}
    for i, _ in enumerate(u.trajectory):
        for k, (a, b) in valid.items():
            out[k][i] = mdadist.distance_array(a.positions, b.positions).min()
    return out


def ca_matrix(u, resids, ref_frame=0):
    """Static CA-CA distance matrix on a single frame."""
    u.trajectory[ref_frame]
    cas, labels = [], []
    for r in resids:
        ag = u.select_atoms(f"protein and resid {r} and name CA")
        if ag.n_atoms:
            cas.append(ag.positions[0]); labels.append(str(r))
        else:
            cas.append([np.nan]*3); labels.append(f"{r}?")
    cas = np.array(cas)
    M = mdadist.distance_array(cas, cas)
    return labels, M


def analyze_system(key):
    top, traj, is_mut = C.SYSTEMS[key]
    if not (os.path.exists(top) and os.path.exists(traj)):
        print(f"[skip] {key}: missing {top} or {traj}")
        return None
    u = mda.Universe(top, traj)
    outdir = C.ensure_out(f"distances/{key}")

    # ---- SANITY: print the resname found at each key resid ----
    print(f"\n=== {key}  (frames={u.trajectory.n_frames}, is_mutant={is_mut}) ===")
    print("SANITY (resid -> resname found in YOUR topology):")
    for r in (C.RES_347, C.RES_D356, C.RES_K1419, C.RES_D1423, C.RES_D1714):
        res = u.select_atoms(f"resid {r}")
        nm = res.residues[0].resname if res.n_atoms else "ABSENT"
        print(f"   resid {r:>5} -> {nm}")

    nzA, labA = cationic_atoms(u, C.RES_347)
    d1423, _  = carboxyl_atoms(u, C.RES_D1423)
    d1714, _  = carboxyl_atoms(u, C.RES_D1714)
    d356,  _  = carboxyl_atoms(u, C.RES_D356)
    if nzA is None or nzA.n_atoms == 0:
        print(f"[warn] {key}: residue {C.RES_347} tip not found; skipping distance series")
        return None

    series_all = min_dist_series_multi(u, {
        "K347_ASP1423": (nzA, d1423),
        "K347_ASP1714": (nzA, d1714),
        "K347_D356":    (nzA, d356),
    })
    series = {k: series_all.get(k) for k in ("K347_ASP1423", "K347_ASP1714", "K347_D356")}

    # write CSV
    csv_path = os.path.join(outdir, f"{key}_timeseries.csv")
    cols = [k for k, v in series.items() if v is not None]
    n = u.trajectory.n_frames
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["frame"] + cols)
        for i in range(n):
            w.writerow([i] + [f"{series[c][i]:.3f}" for c in cols])

    # plot
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for c in cols:
        ax[0].plot(series[c], lw=0.6, label=c)
        ax[1].hist(series[c], bins=40, alpha=0.5, label=c)
    ax[0].axhline(C.GATE_NEAR, color="k", ls="--", lw=0.8)
    ax[0].set_xlabel("frame"); ax[0].set_ylabel("min distance (A)"); ax[0].legend(fontsize=8)
    ax[1].axvline(C.GATE_NEAR, color="k", ls="--", lw=0.8)
    ax[1].set_xlabel("min distance (A)"); ax[1].set_ylabel("count"); ax[1].legend(fontsize=8)
    fig.suptitle(f"{key}: residue-347 tip [{labA}] to DIII/DIV acidics")
    fig.tight_layout(); fig.savefig(os.path.join(outdir, f"{key}_distances.png"), dpi=130)
    plt.close(fig)

    # DEKA-filter geometry matrix (static, on frame 0)
    labels, M = ca_matrix(u, [C.RES_347, C.RES_D356, C.RES_K1419, C.RES_D1423, C.RES_D1714])
    np.savetxt(os.path.join(outdir, f"{key}_CAmatrix.csv"),
               M, delimiter=",", header=",".join(labels), comments="")

    # summary + decision gate
    summary = {"system": key, "tip": labA}
    for c in cols:
        v = series[c]
        summary[c] = dict(mean=float(np.mean(v)), std=float(np.std(v)),
                          frac_near=float(np.mean(v < C.GATE_NEAR)),
                          frac_close=float(np.mean(v < C.GATE_CLOSE)))
    return summary, labels, M


def main():
    print("=" * 70)
    print("Katman 0  Task A: distances + DEKA-filter geometry")
    print("=" * 70)
    results = []
    for key in C.SYSTEMS:
        r = analyze_system(key)
        if r: results.append(r)

    print("\n" + "=" * 70)
    print("DECISION GATE  (focus on holoMutant_prot)")
    print("=" * 70)
    for item in results:
        summ = item[0]
        print(f"\n[{summ['system']}]  residue-347 tip = {summ['tip']}")
        for c in ("K347_ASP1423", "K347_ASP1714", "K347_D356"):
            if c in summ:
                s = summ[c]
                verdict = "COUPLED" if s["frac_near"] > 0.5 else "not coupled"
                print(f"   {c:<14} mean={s['mean']:6.2f} A  std={s['std']:5.2f}  "
                      f"frac<{C.GATE_NEAR:.0f}A={s['frac_near']:.2f}  "
                      f"frac<{C.GATE_CLOSE:.0f}A={s['frac_close']:.2f}  -> {verdict}")
        # filter geometry
        labels, M = item[1], item[2]
        try:
            i347, iK1419, iD1423, iD1714 = labels.index("347"), labels.index("1419"), labels.index("1423"), labels.index("1714")
            print(f"   CA-CA: 347-K1419={M[i347,iK1419]:.1f}  347-D1423={M[i347,iD1423]:.1f}  "
                  f"347-D1714={M[i347,iD1714]:.1f}  K1419-D1423={M[iK1419,iD1423]:.1f} A")
        except ValueError:
            pass
    print("\nInterpretation:")
    print("  frac<12A high in holoMutant_prot  -> DI charge spatially coupled to")
    print("  DIII/DIV acidics: supports 3.3 mechanism + 3.8 vestibule identity.")
    print("  frac<12A low                       -> rewrite 3.3 as a negative result.")


if __name__ == "__main__":
    main()
