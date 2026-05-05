"""
HOLE ANALYSIS: Nav1.5 Pore Geometry
500 ns trajectory average pore profile

OUTPUT FILE:
------------
06_HOLE_500ns_Average.png - Mean pore radius ± SD across full simulation
"""

import MDAnalysis as mda
from mdahole2.analysis import HoleAnalysis
import matplotlib.pyplot as plt
import numpy as np
import os
import warnings

warnings.filterwarnings('ignore')


# ============================================
# CONFIGURATION
# ============================================
class Config:
    TOPOLOGY        = "step5_input.psf"
    N_FILES         = 10
    HOLE_EXECUTABLE = "hole"
    CVECT           = [0, 0, 1]
    STRIDE          = 10    # every 10th frame → ~1000 analyzed
    OUTPUT_FILENAME = "06_HOLE_500ns_Average.png"


# ============================================
# UTILITY
# ============================================
def generate_dcd_list(n_files):
    return [f"step7_{i}.dcd" for i in range(1, n_files + 1)]


def print_section(title):
    print(f"\n{'='*60}")
    print(title)
    print('='*60)


# ============================================
# MAIN
# ============================================
def main():
    print_section("HOLE ANALYSIS: Nav1.5 — 500 ns Average")

    config = Config()

    # Load trajectory
    dcd_files = generate_dcd_list(config.N_FILES)
    dcd_files = [d for d in dcd_files if os.path.exists(d)]

    if len(dcd_files) == 0:
        print("✗ ERROR: No DCD files found!")
        exit()

    print(f"\n✓ Found {len(dcd_files)} DCD files")

    u = mda.Universe(config.TOPOLOGY, dcd_files)
    total_frames    = len(u.trajectory)
    analyzed_frames = total_frames // config.STRIDE
    print(f"✓ Total frames: {total_frames:,}  |  Analyzed: {analyzed_frames:,}  (stride={config.STRIDE})")

    # Run HOLE
    print(f"\n>>> Running HOLE analysis...")
    H = HoleAnalysis(u, executable=config.HOLE_EXECUTABLE, cvect=config.CVECT)
    H.run(step=config.STRIDE)

    # Average profiles
    all_profiles = list(H.results.profiles.values())

    if len(all_profiles) == 0:
        print("✗ ERROR: No profiles found — check HOLE executable")
        exit()

    print(f"✓ {len(all_profiles)} profiles analyzed")

    # Interpolate onto common z-axis
    z_min    = max(p['rxn_coord'].min() for p in all_profiles)
    z_max    = min(p['rxn_coord'].max() for p in all_profiles)
    z_common = np.linspace(z_min, z_max, 500)

    interpolated = np.array([
        np.interp(z_common, p['rxn_coord'], p['radius'])
        for p in all_profiles
    ])

    mean_radius = np.mean(interpolated, axis=0)
    std_radius  = np.std(interpolated,  axis=0)

    min_r = np.min(mean_radius)
    min_z = z_common[np.argmin(mean_radius)]
    print(f"\n  Bottleneck: {min_r:.3f} Å  at z = {min_z:.1f} Å")

    # Plot
    color      = "steelblue"
    fill_color = "lightsteelblue"

    plt.figure(figsize=(8, 11))

    plt.plot( mean_radius, z_common, color=color, lw=2.5, label="Mean radius")
    plt.plot(-mean_radius, z_common, color=color, lw=2.5)

    plt.fill_betweenx(z_common,
                      -(mean_radius + std_radius),
                        mean_radius + std_radius,
                      color=fill_color, alpha=0.30, label="± 1 SD")
    plt.fill_betweenx(z_common,
                      -(mean_radius - std_radius),
                        mean_radius - std_radius,
                      color='white', alpha=0.40)

    plt.annotate(f'Bottleneck\n{min_r:.2f} Å',
                 xy=(min_r, min_z),
                 xytext=(min_r + 4, min_z + 5),
                 arrowprops=dict(facecolor=color, edgecolor=color, arrowstyle='->'),
                 fontsize=10, fontweight='bold', color=color)

    plt.axvline(x= 1.5, color='black', ls=':', alpha=0.6, label="Dehydrated Na⁺ (1.5 Å)")
    plt.axvline(x=-1.5, color='black', ls=':', alpha=0.6)

    plt.title("Nav1.5 Pore Geometry\n(500 ns Trajectory Average)",
              fontsize=14, fontweight='bold')
    plt.xlabel("Radius (Å)", fontsize=12)
    plt.ylabel("Z Coordinate (Channel Depth) (Å)", fontsize=12)
    plt.xlim(-15, 15)
    plt.legend(loc='upper right')
    plt.grid(True, ls='--', alpha=0.5)
    plt.tight_layout()

    plt.savefig(config.OUTPUT_FILENAME, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n✓ Saved: {config.OUTPUT_FILENAME}")
    print_section("ANALYSIS COMPLETE")


if __name__ == "__main__":
    main()
