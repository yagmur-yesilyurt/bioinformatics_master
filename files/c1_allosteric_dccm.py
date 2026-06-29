#!/usr/bin/env python3
"""
c1_allosteric_dccm.py  --  test whether residue 347 is ALLOSTERICALLY coupled
to Site 2 (1423/1714), as a candidate replacement for the refuted local-
electrostatic mechanism.

a1 + a3b showed K347 has no electrostatic effect at Site 2 (26 A, +0.032 kT/e).
But a3 saw the mutant Site 2 more electronegative -- from global structure, not
the charge. If the mutation drives that via long-range correlated motion, we
should see dynamical coupling between 347 and Site 2 in the mutant.

This computes the dynamic cross-correlation matrix (DCCM) of CA fluctuations
from each trajectory (streaming, memory-safe), then for residue 347 reports:
  - correlation to Site 2 (1423, 1714) and to K1419, and to D356 (positive
    control: adjacent, must be strongly coupled);
  - where 347<->Site2 sits in the distribution of ALL residues in the same
    ~26 A distance shell (is the coupling anomalously strong for that distance?);
  - (optional, needs networkx) the strongest correlation PATH 347 -> Site2
    through the contact network, i.e. a concrete allosteric route.

DECISION GATE:
  - mutant 347<->Site2 |correlation| notably higher than WT AND/OR anomalously
    high vs the same-distance shell -> allosteric coupling SUPPORTED. The new
    mechanism for 3.3 becomes "long-range dynamical (allosteric), not local
    electrostatic."
  - no excess coupling, mutant ~ WT ~ background -> NOT allosteric. The a3
    Site-2 difference is then a single-frame/sampling artefact -> requires
    replication; mechanism stays "unestablished / possibly non-specific."

STEP subsamples frames (DCCM is robust to it); raise it if slow.
"""

import os, sys, csv
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STEP = 5
CA = "protein and name CA"
SITE2 = [C.RES_D1423, C.RES_D1714]
SHELL_LO, SHELL_HI = 22.0, 32.0   # distance shell (A) around residue 347
CONTACT = 9.0                     # CA-CA contact cutoff for the network (A)


def dccm_streaming(u):
    ca = u.select_atoms(CA)
    N = ca.n_atoms
    # pass 1: mean positions
    s = np.zeros((N, 3))
    nf = 0
    for _ in u.trajectory[::STEP]:
        s += ca.positions; nf += 1
    mean = s / nf
    # pass 2: covariance of displacements
    Cdot = np.zeros((N, N))
    var = np.zeros(N)
    for _ in u.trajectory[::STEP]:
        d = ca.positions - mean
        Cdot += d @ d.T
        var += np.einsum("ij,ij->i", d, d)
    Cdot /= nf; var /= nf
    denom = np.sqrt(np.outer(var, var))
    denom[denom == 0] = 1e-9
    Cmat = Cdot / denom
    return ca.resids, mean, Cmat


def shortest_path(resids, mean, Cmat, src, dst):
    try:
        import networkx as nx
    except Exception:
        return None
    idx = {int(r): i for i, r in enumerate(resids)}
    if src not in idx or dst not in idx:
        return None
    from MDAnalysis.analysis import distances as D
    dmat = D.distance_array(mean, mean)
    G = nx.Graph()
    N = len(resids)
    for i in range(N):
        for j in range(i + 1, N):
            if dmat[i, j] <= CONTACT:
                c = abs(Cmat[i, j])
                if c > 1e-3:
                    G.add_edge(i, j, weight=-np.log(c))
    try:
        path = nx.shortest_path(G, idx[src], idx[dst], weight="weight")
        length = nx.shortest_path_length(G, idx[src], idx[dst], weight="weight")
        return [int(resids[p]) for p in path], float(length)
    except Exception:
        return None


def analyze_system(key):
    top, traj, is_mut = C.SYSTEMS[key]
    if not (os.path.exists(top) and os.path.exists(traj)):
        print(f"[skip] {key}: missing files"); return None
    u = mda.Universe(top, traj)
    outdir = C.ensure_out(f"allosteric/{key}")
    resids, mean, Cmat = dccm_streaming(u)
    idx = {int(r): i for i, r in enumerate(resids)}

    np.save(os.path.join(outdir, "dccm.npy"), Cmat)

    i347 = idx.get(C.RES_347)
    out = {"system": key, "is_mut": is_mut}
    if i347 is None:
        print(f"[warn] {key}: residue 347 not in CA selection"); return None

    # couplings of interest
    for label, r in [("D356", C.RES_D356), ("K1419", C.RES_K1419),
                     ("ASP1423", C.RES_D1423), ("ASP1714", C.RES_D1714)]:
        j = idx.get(r)
        out[label] = float(Cmat[i347, j]) if j is not None else None

    # same-distance shell background
    from MDAnalysis.analysis import distances as D
    d347 = D.distance_array(mean[i347:i347+1], mean)[0]
    shell = np.where((d347 >= SHELL_LO) & (d347 <= SHELL_HI))[0]
    shell_absC = np.abs(Cmat[i347, shell])
    site2_absC = np.mean([abs(out["ASP1423"] or 0), abs(out["ASP1714"] or 0)])
    pct = float((shell_absC < site2_absC).mean() * 100) if len(shell) else float("nan")
    out["shell_mean_absC"] = float(shell_absC.mean()) if len(shell) else float("nan")
    out["site2_absC"] = float(site2_absC)
    out["site2_percentile_in_shell"] = pct

    # optional allosteric path
    p1 = shortest_path(resids, mean, Cmat, C.RES_347, C.RES_D1423)
    p2 = shortest_path(resids, mean, Cmat, C.RES_347, C.RES_D1714)
    out["path_to_1423"] = p1
    out["path_to_1714"] = p2
    return out


def main():
    print("=" * 70)
    print("Katman 0+  Allosteric coupling test (DCCM 347 <-> Site 2)")
    print("=" * 70)
    rows = [r for r in (analyze_system(k) for k in C.SYSTEMS) if r]

    print(f"\n{'system':<20}{'347-D356':>10}{'347-1419':>10}{'347-1423':>10}{'347-1714':>10}"
          f"{'shell%':>8}")
    print("-" * 68)
    for r in rows:
        def f(x): return f"{x:+.2f}" if isinstance(x, float) else "  NA"
        print(f"{r['system']:<20}{f(r['D356']):>10}{f(r['K1419']):>10}"
              f"{f(r['ASP1423']):>10}{f(r['ASP1714']):>10}"
              f"{r['site2_percentile_in_shell']:>7.0f}%")

    # mutant vs WT comparison for 347<->Site2
    print("\n" + "=" * 70)
    print("DECISION GATE")
    print("=" * 70)
    def site2(r): return r["site2_absC"]
    muts = [r for r in rows if r["is_mut"]]
    wts  = [r for r in rows if not r["is_mut"]]
    if muts and wts:
        mm = np.mean([site2(r) for r in muts]); ww = np.mean([site2(r) for r in wts])
        print(f"  mean |347<->Site2| correlation:  mutant={mm:.3f}  WT={ww:.3f}")
        excess = mm - ww
        shellpct = np.mean([r["site2_percentile_in_shell"] for r in muts])
        verdict = ("ALLOSTERIC COUPLING SUPPORTED" if (excess > 0.10 or shellpct > 80)
                   else "NO clear allosteric coupling -> a3 Site-2 diff is likely a "
                        "single-frame artefact; replication needed")
        print(f"  mutant excess over WT = {excess:+.3f}; mutant Site2 shell percentile "
              f"= {shellpct:.0f}%")
        print(f"  -> {verdict}")
    # report any communication path found
    for r in rows:
        if r["is_mut"] and r.get("path_to_1423"):
            path, length = r["path_to_1423"]
            print(f"\n  [{r['system']}] strongest 347->ASP1423 correlation path "
                  f"(len {length:.2f}, {len(path)} nodes):")
            print("     " + " - ".join(map(str, path)))
    print("\n(Positive control: 347<->D356 should be strongly coupled in mutants;")
    print(" if it is, the DCCM is meaningful. If even D356 coupling is weak, the")
    print(" trajectory alignment/length is the limiting factor.)")
    print("\nNote: install networkx for the allosteric-path output (pip install networkx).")


if __name__ == "__main__":
    main()
