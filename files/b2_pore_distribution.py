#!/usr/bin/env python3
"""
b2_pore_distribution.py  --  Katman 0, Task B2

Addresses report section 3.7.

The report gives a single 500 ns-average pore bottleneck per system. This
script runs HOLE per frame (over an evenly spaced subsample) and reports the
bottleneck radius as a DISTRIBUTION (mean +/- std, full histogram), so the
small inter-system differences (e.g. 1.40 vs 1.56 A) can be tested against
within-system fluctuation instead of asserted from point values.

Requires the HOLE executable (config.BIN_HOLE) and MDAnalysis hole2 wrapper.

DECISION GATE:
  Compute, for each pair of systems, whether the difference in mean
  bottleneck exceeds the pooled std. If the 1.40-vs-1.56 type differences
  are within ~1 std, the categorical "above/below the 1.5 A Na+ radius"
  claim is not supported and 3.7 must be softened.  -> 3.7 toward 4.
"""

import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STRIDE = 25   # analyse every STRIDE-th frame (10000/25 = 400 frames)
NA_RADIUS = 1.5


def bottleneck_series(top, traj, outdir, key):
    import MDAnalysis as mda
    from MDAnalysis.analysis import hole2
    u = mda.Universe(top, traj)
    tmp = os.path.join(outdir, f"hole_tmp_{key}")
    os.makedirs(tmp, exist_ok=True)
    ha = hole2.HoleAnalysis(u, executable=C.BIN_HOLE, prefix=os.path.join(tmp, "h"))
    ha.run(start=0, stop=None, step=STRIDE)
    # min radius per analysed frame
    mins = []
    for profile in ha.results.profiles.values():
        mins.append(float(np.min(profile.radius)))
    ha.delete_temporary_files()
    return np.array(mins)


def main():
    outdir = C.ensure_out("pore")
    print("=" * 70)
    print("Katman 0  Task B2: per-frame pore bottleneck distribution")
    print("=" * 70)
    stats = {}
    for key, (top, traj, _) in C.SYSTEMS.items():
        if not (os.path.exists(top) and os.path.exists(traj)):
            print(f"[skip] {key}: missing files"); continue
        try:
            mins = bottleneck_series(top, traj, outdir, key)
        except Exception as e:
            print(f"[error] {key}: {e}"); continue
        np.savetxt(os.path.join(outdir, f"{key}_bottleneck.csv"), mins,
                   header="bottleneck_radius_A", comments="")
        stats[key] = (float(mins.mean()), float(mins.std()), len(mins))
        print(f"  {key:<22} bottleneck = {mins.mean():.2f} +/- {mins.std():.2f} A  (n={len(mins)})")

        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.hist(mins, bins=30, alpha=0.8)
        ax.axvline(NA_RADIUS, color="r", ls="--", label="dehydrated Na+ 1.5 A")
        ax.axvline(mins.mean(), color="k", ls="-", lw=1)
        ax.set_xlabel("bottleneck radius (A)"); ax.set_ylabel("count")
        ax.set_title(f"{key}: pore bottleneck distribution"); ax.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(outdir, f"{key}_bottleneck.png"), dpi=130)
        plt.close(fig)

    # pairwise significance vs noise
    print("\n" + "=" * 70)
    print("DECISION GATE: is a mean difference larger than the pooled std?")
    print("=" * 70)
    keys = list(stats)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = keys[i], keys[j]
            ma, sa, _ = stats[a]; mb, sb, _ = stats[b]
            diff = abs(ma - mb); pooled = np.hypot(sa, sb)
            tag = "RESOLVABLE" if diff > pooled else "within noise"
            print(f"  {a:<20} vs {b:<20} d={diff:.2f} A  pooled_std={pooled:.2f}  -> {tag}")
    print("\nIf the 1.40-vs-1.56 type contrasts are 'within noise', soften 3.7:")
    print("the categorical above/below-Na+ framing is not supported by point values.")


if __name__ == "__main__":
    main()
