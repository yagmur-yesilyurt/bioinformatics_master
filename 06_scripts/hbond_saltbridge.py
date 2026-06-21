"""
HYDROGEN BOND & SALT BRIDGE ANALYSIS: Mexiletine - Nav1.5
Binding interaction characterization at the N347K mutation site

ANALYSES:
---------
1. H-BOND ANALYSIS (Mexiletine ↔ Protein)
   - H-bond occupancy over time
   - Per-residue H-bond frequency
   - Most persistent H-bond donors/acceptors

2. SALT BRIDGE ANALYSIS (N347K mutation effect)
   - K347 (Lys) salt bridges with nearby acidic residues (Asp, Glu)
   - Comparison: WT (N347, neutral) vs Mutant (K347, +1 charge)
   - Salt bridge distance over time
   - Occupancy percentage

OUTPUT FILES:
-------------
1. hbond_timeseries.png        - H-bond count over time
2. hbond_occupancy.png         - Per-residue H-bond frequency bar plot
3. saltbridge_distance.png     - K347 salt bridge distances over time
4. saltbridge_occupancy.png    - Salt bridge occupancy summary
5. hbond_saltbridge_summary.txt - Numerical results for reporting
"""

import MDAnalysis as mda
from MDAnalysis.analysis import hbonds
from MDAnalysis.analysis import distances
import numpy as np
import matplotlib.pyplot as plt
import os
import warnings

warnings.filterwarnings('ignore')
plt.style.use('seaborn-v0_8-darkgrid')


# ============================================
# CONFIGURATION
# ============================================
class Config:
    PSF_FILE       = 'step5_input.psf'
    N_FILES        = 10
    DT             = 0.05    # ns per frame
    STRIDE         = 5       # analyze every 5th frame

    LIGAND_NAME    = 'UNL'
    MUTATION_RESID = 347

    # H-bond geometric criteria (standard)
    HBOND_CUTOFF   = 3.5     # Å donor-acceptor distance
    HBOND_ANGLE    = 150.0   # degrees (D-H...A angle)

    # Salt bridge criteria
    SALTBRIDGE_CUTOFF = 4.0  # Å between charged groups

    # Residues to search for salt bridge partners with K347
    # Acidic residues near P-loop (340-370)
    ACIDIC_RESIDUES = list(range(330, 400))  # search window

    # Binding pocket residues for H-bond analysis
    BINDING_RESIDUES = list(range(340, 361))


# ============================================
# UTILITY FUNCTIONS
# ============================================
def print_section(title):
    print(f"\n{'='*60}")
    print(title)
    print('='*60)


def generate_trajectory_list(n_files):
    return [f"step7_{i}.dcd" for i in range(1, n_files + 1)]


def load_universe(config):
    dcd_files = generate_trajectory_list(config.N_FILES)
    missing   = [f for f in [config.PSF_FILE] + dcd_files if not os.path.exists(f)]
    if missing:
        print(f"✗ ERROR: Missing files: {missing}")
        exit()

    print(f"\n>>> Loading {len(dcd_files)} DCD files (stride={config.STRIDE})...")
    u = mda.Universe(config.PSF_FILE, dcd_files)
    n_analyzed = len(u.trajectory[::config.STRIDE])
    print(f"✓ Loaded: {len(u.trajectory):,} frames total")
    print(f"  Frames to analyze: {n_analyzed:,}")
    return u


# ============================================
# H-BOND ANALYSIS
# ============================================
def analyze_hbonds(universe, config):
    """
    Calculate H-bonds between mexiletine and protein binding pocket.
    Uses MDAnalysis HydrogenBondAnalysis.
    """
    print_section("HYDROGEN BOND ANALYSIS")

    ligand  = universe.select_atoms(f"resname {config.LIGAND_NAME}")
    pocket  = universe.select_atoms(
        f"resid {' '.join(map(str, config.BINDING_RESIDUES))}"
    )

    if len(ligand) == 0:
        print(f"✗ Ligand {config.LIGAND_NAME} not found!")
        return None, None

    print(f"  Ligand atoms:  {len(ligand)}")
    print(f"  Pocket atoms:  {len(pocket)}")

    # Run H-bond analysis
    print("\n>>> Running H-bond analysis...")
    try:
        hbond_analysis = hbonds.HydrogenBondAnalysis(
            universe,
            donors_sel=None,    # auto-detect
            acceptors_sel=None, # auto-detect
            d_a_cutoff=config.HBOND_CUTOFF,
            d_h_a_angle_cutoff=config.HBOND_ANGLE,
            update_selections=False
        )

        # Select between ligand and pocket
        hbond_analysis.donors_sel    = f"(resname {config.LIGAND_NAME} or resid {' '.join(map(str, config.BINDING_RESIDUES))}) and name N* O* S*"
        hbond_analysis.acceptors_sel = f"(resname {config.LIGAND_NAME} or resid {' '.join(map(str, config.BINDING_RESIDUES))}) and name N* O* S*"

        hbond_analysis.run(
            start=0,
            stop=None,
            step=config.STRIDE,
            verbose=True
        )

        results = hbond_analysis.results.hbonds
        print(f"  ✓ H-bond analysis complete")
        return hbond_analysis, results

    except Exception as e:
        print(f"  ⚠ HydrogenBondAnalysis failed: {e}")
        print("  → Falling back to manual distance-based H-bond detection")
        return None, manual_hbond_analysis(universe, config)


def manual_hbond_analysis(universe, config):
    """
    Fallback: simple distance-based H-bond counting between
    mexiletine N/O atoms and protein N/O atoms.
    """
    print("\n>>> Manual distance-based H-bond analysis...")

    ligand_polar = universe.select_atoms(
        f"resname {config.LIGAND_NAME} and name N* O*"
    )
    protein_polar = universe.select_atoms(
        f"resid {' '.join(map(str, config.BINDING_RESIDUES))} and name N* O*"
    )

    print(f"  Ligand polar atoms:  {len(ligand_polar)}")
    print(f"  Protein polar atoms: {len(protein_polar)}")

    hbond_counts  = []
    time_list     = []
    hbond_pairs   = {}   # (lig_atom, prot_resid) → count

    n_frames = len(universe.trajectory[::config.STRIDE])

    for idx, ts in enumerate(universe.trajectory[::config.STRIDE]):
        if idx % 200 == 0:
            print(f"  Frame {idx:5d}/{n_frames} ({idx/n_frames*100:.0f}%)", end='\r')

        if len(ligand_polar) == 0 or len(protein_polar) == 0:
            hbond_counts.append(0)
            time_list.append(ts.frame * config.DT)
            continue

        dist_mat = distances.distance_array(
            ligand_polar.positions,
            protein_polar.positions
        )

        # Count contacts within H-bond distance
        contacts = np.argwhere(dist_mat < config.HBOND_CUTOFF)
        hbond_counts.append(len(contacts))
        time_list.append(ts.frame * config.DT)

        # Track which pairs form contacts
        for li, pi in contacts:
            lig_name  = ligand_polar[li].name
            prot_res  = protein_polar[pi].resid
            prot_name = protein_polar[pi].name
            key = (lig_name, prot_res, prot_name)
            hbond_pairs[key] = hbond_pairs.get(key, 0) + 1

    print(f"\n  ✓ Manual H-bond analysis complete")
    return np.array(time_list), np.array(hbond_counts), hbond_pairs, n_frames


# ============================================
# SALT BRIDGE ANALYSIS
# ============================================
def analyze_salt_bridges(universe, config):
    """
    Analyze salt bridges formed by K347 (N347K mutation).
    In WT: N347 is neutral, no salt bridge possible.
    In Mutant: K347 is +1, can form salt bridge with Asp/Glu.
    """
    print_section("SALT BRIDGE ANALYSIS")

    # K347 NZ atom (lysine side chain amine, +1 charge)
    try:
        k347_nz = universe.select_atoms(f"resid {config.MUTATION_RESID} and name NZ")
        if len(k347_nz) == 0:
            # Try N347 in WT
            k347_nz = universe.select_atoms(
                f"resid {config.MUTATION_RESID} and name ND2"
            )
            is_mutant = False
            print(f"  Residue {config.MUTATION_RESID}: N347 (WT) — no salt bridge expected")
        else:
            is_mutant = True
            print(f"  Residue {config.MUTATION_RESID}: K347 (Mutant) — salt bridge analysis")
    except:
        print("  ✗ Could not find residue 347")
        return None

    # Find acidic residues (Asp/Glu) in search window
    acidic = universe.select_atoms(
        f"resid {' '.join(map(str, config.ACIDIC_RESIDUES))} "
        f"and resname ASP GLU and name OD1 OD2 OE1 OE2"
    )

    if len(acidic) == 0:
        print("  ✗ No acidic residues found in search window")
        return None

    # Get unique acidic residues
    acidic_resids = list(set(acidic.resids))
    print(f"  Acidic residues in window: {acidic_resids}")
    print(f"  K347 NZ atoms: {len(k347_nz)}")

    # Calculate distances over trajectory
    sb_distances = {resid: [] for resid in acidic_resids}
    time_list = []

    n_frames = len(universe.trajectory[::config.STRIDE])
    print(f"\n>>> Calculating salt bridge distances ({n_frames} frames)...")

    for idx, ts in enumerate(universe.trajectory[::config.STRIDE]):
        if idx % 200 == 0:
            print(f"  Frame {idx:5d}/{n_frames} ({idx/n_frames*100:.0f}%)", end='\r')

        time_list.append(ts.frame * config.DT)

        if len(k347_nz) == 0:
            for resid in acidic_resids:
                sb_distances[resid].append(np.nan)
            continue

        for resid in acidic_resids:
            res_atoms = acidic.select_atoms(f"resid {resid}")
            if len(res_atoms) == 0:
                sb_distances[resid].append(np.nan)
                continue

            # Minimum distance between K347 NZ and acidic oxygen atoms
            dist_mat = distances.distance_array(
                k347_nz.positions,
                res_atoms.positions
            )
            sb_distances[resid].append(np.min(dist_mat))

    print(f"\n  ✓ Salt bridge analysis complete")

    # Convert to arrays
    for resid in acidic_resids:
        sb_distances[resid] = np.array(sb_distances[resid])

    return np.array(time_list), sb_distances, acidic_resids, is_mutant if 'is_mutant' in locals() else False


# ============================================
# VISUALIZATION
# ============================================
def plot_hbond_timeseries(time_ns, hbond_counts):
    """H-bond count over time"""
    print("\n>>> Plotting H-bond time series...")

    fig, ax = plt.subplots(figsize=(12, 5))

    ax.plot(time_ns, hbond_counts, color='#2196F3', linewidth=1.0,
            alpha=0.7, label='H-bond contacts')

    # Running average
    window = min(50, len(hbond_counts)//10)
    if window > 1:
        running_avg = np.convolve(hbond_counts,
                                  np.ones(window)/window, mode='valid')
        ax.plot(time_ns[window-1:], running_avg,
                color='#e53935', linewidth=2.0,
                label=f'Running avg ({window} frames)')

    ax.axhline(y=np.mean(hbond_counts), color='black', linestyle='--',
               alpha=0.5, label=f'Mean: {np.mean(hbond_counts):.2f}')

    ax.set_xlabel('Time (ns)', fontsize=13)
    ax.set_ylabel('H-bond contacts (N/O within 3.5 Å)', fontsize=12)
    ax.set_title('Mexiletine–Protein H-bond Contacts Over Time\n'
                 '(Binding pocket residues 340–360)',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig('13_hbond_timeseries.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved: hbond_timeseries.png")


def plot_hbond_occupancy(hbond_pairs, n_frames):
    """Per-residue H-bond frequency"""
    print(">>> Plotting H-bond occupancy...")

    if not hbond_pairs:
        print("  ⚠ No H-bond pairs found")
        return

    # Aggregate by residue
    resid_counts = {}
    for (lig_name, prot_res, prot_name), count in hbond_pairs.items():
        resid_counts[prot_res] = resid_counts.get(prot_res, 0) + count

    # Convert to occupancy %
    resids     = sorted(resid_counts.keys())
    occupancy  = [resid_counts[r] / n_frames * 100 for r in resids]

    # Only show residues with >5% occupancy
    mask      = [o > 5 for o in occupancy]
    resids    = [r for r, m in zip(resids, mask) if m]
    occupancy = [o for o, m in zip(occupancy, mask) if m]

    if not resids:
        print("  ⚠ No residues with >5% H-bond occupancy")
        return

    colors = ['#e53935' if r == 347 else '#1976D2' for r in resids]

    fig, ax = plt.subplots(figsize=(12, 5))
    bars = ax.bar(range(len(resids)), occupancy, color=colors,
                  edgecolor='white', linewidth=0.5)

    ax.set_xticks(range(len(resids)))
    ax.set_xticklabels([str(r) for r in resids], rotation=45, fontsize=10)
    ax.set_xlabel('Residue ID', fontsize=12)
    ax.set_ylabel('H-bond Occupancy (%)', fontsize=12)
    ax.set_title('Per-Residue H-bond Occupancy with Mexiletine\n'
                 '(Red = K347 mutation site)',
                 fontsize=13, fontweight='bold')
    ax.axhline(y=20, color='gray', linestyle='--', alpha=0.5,
               label='20% threshold')
    ax.legend(fontsize=10)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig('hbond_occupancy.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved: hbond_occupancy.png")


def plot_salt_bridges(time_ns, sb_distances, acidic_resids, is_mutant):
    """Salt bridge distances over time"""
    print(">>> Plotting salt bridge distances...")

    # Only plot residues that come within cutoff at least once
    relevant = []
    for resid in acidic_resids:
        dists = sb_distances[resid]
        valid = dists[~np.isnan(dists)]
        if len(valid) > 0 and np.min(valid) < 6.0:
            relevant.append(resid)

    if not relevant:
        print("  ⚠ No salt bridge partners found within 6 Å")
        relevant = acidic_resids[:5]  # show top 5 anyway

    fig, ax = plt.subplots(figsize=(12, 6))

    colors = plt.cm.tab10(np.linspace(0, 1, len(relevant)))
    for resid, color in zip(relevant, colors):
        dists = sb_distances[resid]
        label_resid = f"Resid {resid}"
        ax.plot(time_ns, dists, linewidth=1.0, alpha=0.7,
                color=color, label=label_resid)

    ax.axhline(y=4.0, color='red', linestyle='--', linewidth=1.5,
               label='Salt bridge cutoff (4.0 Å)')

    residue_label = "K347 (N347K)" if is_mutant else "N347 (WT)"
    ax.set_xlabel('Time (ns)', fontsize=13)
    ax.set_ylabel('Distance (Å)', fontsize=12)
    ax.set_title(f'Salt Bridge Analysis: {residue_label} — Acidic Residues\n'
                 f'(Distance between NZ and carboxylate oxygens)',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig('14_saltbridge_distance.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved: saltbridge_distance.png")


def plot_salt_bridge_occupancy(sb_distances, acidic_resids, is_mutant, cutoff=4.0):
    """Salt bridge occupancy bar chart"""
    print(">>> Plotting salt bridge occupancy...")

    occupancies = {}
    for resid in acidic_resids:
        dists = sb_distances[resid]
        valid = dists[~np.isnan(dists)]
        if len(valid) > 0:
            occupancies[resid] = np.sum(valid < cutoff) / len(valid) * 100
        else:
            occupancies[resid] = 0.0

    resids = sorted(occupancies.keys())
    occ    = [occupancies[r] for r in resids]

    fig, ax = plt.subplots(figsize=(12, 5))
    colors = ['#e53935' if o > 20 else '#90A4AE' for o in occ]
    ax.bar(range(len(resids)), occ, color=colors,
           edgecolor='white', linewidth=0.5)

    ax.set_xticks(range(len(resids)))
    ax.set_xticklabels([str(r) for r in resids], rotation=45, fontsize=10)
    ax.axhline(y=20, color='black', linestyle='--', alpha=0.5,
               label='20% threshold')
    ax.set_xlabel('Acidic Residue ID', fontsize=12)
    ax.set_ylabel('Salt Bridge Occupancy (%)', fontsize=12)

    residue_label = "K347 (N347K mutant)" if is_mutant else "N347 (WT)"
    ax.set_title(f'Salt Bridge Occupancy: {residue_label}\n'
                 f'(% frames with distance < {cutoff} Å)',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig('15_saltbridge_occupancy.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved: saltbridge_occupancy.png")


# ============================================
# SUMMARY
# ============================================
def save_summary(time_ns, hbond_counts, hbond_pairs, n_frames,
                 sb_distances, acidic_resids, is_mutant):
    print(">>> Saving summary...")

    lines = [
        "=" * 60,
        "H-BOND & SALT BRIDGE SUMMARY: Nav1.5 + Mexiletine",
        "=" * 60,
        "",
        "--- HYDROGEN BOND ANALYSIS ---",
        f"Mean H-bond contacts/frame: {np.mean(hbond_counts):.2f}",
        f"Max H-bond contacts:        {np.max(hbond_counts):.0f}",
        f"% frames with ≥1 H-bond:    "
        f"{np.sum(hbond_counts >= 1)/len(hbond_counts)*100:.1f}%",
        "",
        "Top H-bond forming residues (>10% occupancy):"
    ]

    resid_counts = {}
    for (lig_name, prot_res, prot_name), count in hbond_pairs.items():
        resid_counts[prot_res] = resid_counts.get(prot_res, 0) + count

    for resid, count in sorted(resid_counts.items(),
                                key=lambda x: -x[1])[:10]:
        occ = count / n_frames * 100
        if occ > 10:
            marker = " ← MUTATION SITE" if resid == 347 else ""
            lines.append(f"  Resid {resid:4d}: {occ:.1f}%{marker}")

    lines += [
        "",
        "--- SALT BRIDGE ANALYSIS ---",
        f"Residue 347: {'K347 (Mutant, +1)' if is_mutant else 'N347 (WT, neutral)'}",
        "",
        "Salt bridge occupancy (cutoff 4.0 Å):"
    ]

    for resid in sorted(acidic_resids):
        dists = sb_distances[resid]
        valid = dists[~np.isnan(dists)]
        if len(valid) > 0:
            occ = np.sum(valid < 4.0) / len(valid) * 100
            mean_d = np.mean(valid)
            lines.append(
                f"  K347–Resid {resid:4d}: "
                f"{occ:.1f}% occupancy, "
                f"mean dist = {mean_d:.2f} Å"
            )

    lines += ["", "=" * 60]

    with open('12_hbond_saltbridge_summary.txt', 'w') as f:
        f.write('\n'.join(lines))

    print("  ✓ Saved: hbond_saltbridge_summary.txt")
    print("\n" + '\n'.join(lines[6:20]))


# ============================================
# MAIN
# ============================================
def main():
    print_section("H-BOND & SALT BRIDGE ANALYSIS: Nav1.5")

    config   = Config()
    universe = load_universe(config)

    # H-bond analysis (holo only — skipped if no ligand)
    hbond_ok = False
    time_ns = hbond_counts = hbond_pairs = n_frames = None

    _, hbond_result = analyze_hbonds(universe, config)

    if isinstance(hbond_result, tuple):
        time_ns, hbond_counts, hbond_pairs, n_frames = hbond_result
        hbond_ok = True
        print("  ✓ H-bond analysis complete")
    else:
        print("  ⊘ H-bond analysis skipped (apo system — no ligand)")

    # Salt bridge analysis (always runs — holo and apo)
    sb_result = analyze_salt_bridges(universe, config)
    if sb_result is None:
        print("  ✗ Salt bridge analysis failed")
        return
    time_sb, sb_distances, acidic_resids, is_mutant = sb_result

    # Plots — hbond only if ligand present
    if hbond_ok:
        plot_hbond_timeseries(time_ns, hbond_counts)
        plot_hbond_occupancy(hbond_pairs, n_frames)
    plot_salt_bridges(time_sb, sb_distances, acidic_resids, is_mutant)
    plot_salt_bridge_occupancy(sb_distances, acidic_resids, is_mutant)

    # Summary
    save_summary(time_ns, hbond_counts, hbond_pairs, n_frames,
                 sb_distances, acidic_resids, is_mutant)

    print_section("ANALYSIS COMPLETE")
    print("Output files:")
    if hbond_ok:
        print("  • 13_hbond_timeseries.png      — H-bond contacts over time")
        print("  • hbond_occupancy.png          — Per-residue H-bond frequency")
    print("  • 14_saltbridge_distance.png   — K347 salt bridge distances")
    print("  • 15_saltbridge_occupancy.png  — Salt bridge occupancy summary")
    print("  • 12_hbond_saltbridge_summary.txt — Numerical results")


if __name__ == "__main__":
    main()
