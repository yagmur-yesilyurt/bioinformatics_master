"""
MD ANALYSIS PIPELINE: Nav1.5
Comprehensive trajectory analysis — works for both HOLO and APO systems

OUTPUT FILES (holo):
--------------------
01_protein_flexibility_rmsf.png
02_ligand_pore_distance.png
03_residue_contact_map.png
04_ligand_conformational_rmsd.png
05_multisite_distances.png
06_ligand_density_map.dx
07_protein_rmsd_global.png
08_protein_rmsd_domains.png
09_protein_rmsd_ploop.png

OUTPUT FILES (apo):
-------------------
01_protein_flexibility_rmsf.png
07_protein_rmsd_global.png
08_protein_rmsd_domains.png
09_protein_rmsd_ploop.png
"""

import MDAnalysis as mda
from MDAnalysis.analysis import distances, rms
from MDAnalysis.analysis.density import DensityAnalysis
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import warnings

warnings.filterwarnings('ignore')
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")


# ============================================
# CONFIGURATION
# ============================================
class Config:
    DCD_PATTERN    = 'step7_*.dcd'
    PSF_FILE       = 'step5_input.psf'
    N_FILES        = 10
    DT             = 0.05   # ns per frame
    STRIDE         = 1
    PORE_RESIDS    = [372, 901, 1422, 1714]
    MUTATION_RESID = 347
    LIGAND_NAME    = 'UNL'

    DOMAINS = {
        'Domain I':   (1,    400),
        'Domain II':  (401,  800),
        'Domain III': (801,  1200),
        'Domain IV':  (1201, 1800)
    }
    DOMAIN_COLORS = {
        'Domain I':   '#e41a1c',
        'Domain II':  '#377eb8',
        'Domain III': '#4daf4a',
        'Domain IV':  '#ff7f00'
    }
    PLOOP_RANGE = (330, 370)


# ============================================
# UTILITY
# ============================================
def print_section(title):
    print(f"\n{'='*60}")
    print(title)
    print('='*60)


def generate_dcd_list(n_files):
    return [f"step7_{i}.dcd" for i in range(1, n_files + 1)]


def load_universe(config):
    dcd_files = sorted(glob.glob(config.DCD_PATTERN))

    if not dcd_files:
        print(f"✗ ERROR: No DCD files found with pattern '{config.DCD_PATTERN}'")
        exit()

    print(f"\n✓ Found {len(dcd_files)} DCD files:")
    for i, f in enumerate(dcd_files, 1):
        size_mb = os.path.getsize(f) / (1024**2)
        print(f"  [{i:2d}] {f:40s} ({size_mb:6.1f} MB)")

    if not os.path.exists(config.PSF_FILE):
        print(f"\n✗ ERROR: {config.PSF_FILE} not found!")
        exit()

    print(f"\nLoading trajectories...")
    try:
        u = mda.Universe(config.PSF_FILE, *dcd_files)
        print(f"✓ Loaded: {len(u.trajectory):,} frames | "
              f"{len(u.trajectory)*config.DT:.1f} ns | "
              f"{len(u.atoms):,} atoms")
        return u
    except Exception as e:
        print(f"✗ ERROR: {e}")
        exit()


def load_reference(config):
    dcd_first = generate_dcd_list(config.N_FILES)[0]
    return mda.Universe(config.PSF_FILE, dcd_first)


def setup_selections(universe, config):
    selections = {}

    lig = universe.select_atoms(f"resname {config.LIGAND_NAME}")
    if len(lig) == 0:
        print(f"⚠ Ligand '{config.LIGAND_NAME}' not found → APO system")
        print("  Ligand-dependent analyses will be skipped.")
        selections['ligand'] = None
    else:
        print(f"✓ Ligand found: {len(lig)} atoms")
        selections['ligand'] = lig

    selections['mutation']   = universe.select_atoms(
        f"resid {config.MUTATION_RESID} and name CA"
    )
    selections['protein']    = universe.select_atoms('protein')
    selections['protein_ca'] = universe.select_atoms('protein and name CA')

    try:
        pore_sel = f"resid {' '.join(map(str, config.PORE_RESIDS))} and name CA"
        pore = universe.select_atoms(pore_sel)
        selections['pore'] = pore if len(pore) > 0 else None
        if selections['pore']:
            print(f"✓ Pore residues found: {len(pore)} atoms")
        else:
            print("⚠ Pore residues not found")
    except:
        selections['pore'] = None

    return selections


# ============================================
# RMSF
# ============================================
def calculate_rmsf(universe, selections, config):
    print("\n[1/9] Calculating RMSF...")

    rmsf_analysis = rms.RMSF(selections['protein_ca']).run()

    fig, ax = plt.subplots(figsize=(14, 6))
    ax.plot(selections['protein_ca'].resids,
            rmsf_analysis.results.rmsf,
            color='#1f77b4', linewidth=1.5)
    ax.axvline(x=config.MUTATION_RESID, color='red', linestyle='--',
               label=f'Mutation Site ({config.MUTATION_RESID})')
    ax.set_xlabel('Residue ID', fontsize=14)
    ax.set_ylabel('RMSF (Å)', fontsize=14)
    ax.set_title('Local Flexibility (RMSF)', fontsize=16, fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('01_protein_flexibility_rmsf.png', dpi=300)
    plt.close()
    print("  ✓ Saved: 01_protein_flexibility_rmsf.png")


# ============================================
# LIGAND ANALYSES (holo only)
# ============================================
def calculate_pore_distance(universe, selections, config, time_array):
    print("\n[2/9] Ligand-Pore Distance...")

    if selections['ligand'] is None:
        print("  ⊘ Skipped (apo)")
        return None
    if selections['pore'] is None:
        print("  ⊘ Skipped (no pore residues)")
        return None

    pore_dists = []
    for ts in universe.trajectory[::config.STRIDE]:
        ligand_com = selections['ligand'].center_of_mass().reshape(1, 3)
        d = distances.distance_array(ligand_com, selections['pore'].positions)
        pore_dists.append(d.min())
    pore_dists = np.array(pore_dists)

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(time_array[::config.STRIDE], pore_dists, color='darkblue', linewidth=1.5)
    ax.axhline(y=15, color='red', linestyle='--', label='Proximity Threshold (15 Å)')
    ax.set_ylabel('Distance to Pore (Å)', fontsize=12)
    ax.set_xlabel('Time (ns)', fontsize=12)
    ax.set_title('Ligand - Pore Cavity Distance', fontweight='bold')
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('02_ligand_pore_distance.png', dpi=300)
    plt.close()
    print("  ✓ Saved: 02_ligand_pore_distance.png")
    return pore_dists


def calculate_contact_map(universe, selections, config):
    print("\n[3/9] Contact Map...")

    if selections['ligand'] is None:
        print("  ⊘ Skipped (apo)")
        return

    res_range = range(340, 361)
    n_frames  = len(universe.trajectory)
    contacts_data = np.zeros((len(res_range), n_frames))

    for i, r in enumerate(res_range):
        ag = universe.select_atoms(f"resid {r} and not name H*")
        if len(ag) > 0:
            for fidx, ts in enumerate(universe.trajectory):
                d = distances.distance_array(ag.positions,
                                             selections['ligand'].positions)
                if (d < 4.5).any():
                    contacts_data[i, fidx] = 1

    fig, ax = plt.subplots(figsize=(12, 6))
    stride_plot = max(1, n_frames // 1000)
    im = ax.imshow(contacts_data[:, ::stride_plot], aspect='auto',
                   cmap='Greys', interpolation='nearest')
    ax.set_yticks(range(len(res_range)))
    ax.set_yticklabels(list(res_range))
    ax.set_xlabel('Time (subsampled)', fontsize=12)
    ax.set_ylabel('Residue ID', fontsize=12)
    ax.set_title('Contact Map (< 4.5 Å)', fontsize=14, fontweight='bold')
    plt.colorbar(im, ax=ax, label='Contact')
    plt.tight_layout()
    plt.savefig('03_residue_contact_map.png', dpi=300)
    plt.close()
    print("  ✓ Saved: 03_residue_contact_map.png")


def calculate_ligand_rmsd(universe, selections):
    print("\n[4/9] Ligand RMSD...")

    if selections['ligand'] is None:
        print("  ⊘ Skipped (apo)")
        return

    lig_rms = rms.RMSD(selections['ligand'], selections['ligand'],
                       select='not name H*', ref_frame=0).run()
    rmsd_data = lig_rms.results.rmsd
    time_ns   = rmsd_data[:, 1] / 1000
    rmsd_vals = rmsd_data[:, 2]

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(time_ns, rmsd_vals, color='darkgreen', linewidth=1.5)
    ax.set_ylabel('Ligand RMSD (Å)', fontsize=12)
    ax.set_xlabel('Time (ns)', fontsize=12)
    ax.set_title('Mexiletine Internal Conformational Changes',
                 fontsize=14, fontweight='bold')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('04_ligand_conformational_rmsd.png', dpi=300)
    plt.close()
    print("  ✓ Saved: 04_ligand_conformational_rmsd.png")


def calculate_multisite_distances(universe, selections, config, time_array):
    print("\n[5/9] Multi-site Distances...")

    if selections['ligand'] is None:
        print("  ⊘ Skipped (apo)")
        return {}

    key_sites = {'K347': config.MUTATION_RESID, 'S344': 344, 'L348': 348}
    if selections['pore'] is not None:
        key_sites['Pore_Center'] = config.PORE_RESIDS[0]

    site_data = {name: [] for name in key_sites}
    for ts in universe.trajectory[::config.STRIDE]:
        ligand_com = selections['ligand'].center_of_mass().reshape(1, 3)
        for name, resid in key_sites.items():
            atom = universe.select_atoms(f"resid {resid} and name CA")
            if len(atom) > 0:
                d = distances.distance_array(ligand_com, atom.positions)[0, 0]
                site_data[name].append(d)
            else:
                site_data[name].append(np.nan)

    fig, ax = plt.subplots(figsize=(14, 6))
    for name, data in site_data.items():
        if len(data) > 0 and not np.all(np.isnan(data)):
            ax.plot(time_array[::config.STRIDE], data,
                    label=name, alpha=0.8, linewidth=1.5)
    ax.set_ylabel('Distance (Å)', fontsize=12)
    ax.set_xlabel('Time (ns)', fontsize=12)
    ax.set_title('Mexiletine Distance to Key Sites', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('05_multisite_distances.png', dpi=300)
    plt.close()
    print("  ✓ Saved: 05_multisite_distances.png")
    return site_data


def calculate_density(universe, selections):
    print("\n[6/9] 3D Density Map...")

    if selections['ligand'] is None:
        print("  ⊘ Skipped (apo)")
        return

    try:
        dens = DensityAnalysis(selections['ligand'], delta=2.0).run()
        dens.density.export('06_ligand_density_map.dx', type='double')
        print("  ✓ Saved: 06_ligand_density_map.dx")
    except Exception as e:
        print(f"  ⚠ Density Error: {e}")


# ============================================
# PROTEIN RMSD
# ============================================
def calculate_global_rmsd(universe, ref, config):
    print("\n[7/9] Global Backbone RMSD...")

    rmsd_analysis = rms.RMSD(
        universe.select_atoms('backbone'),
        ref.select_atoms('backbone'),
        select='backbone', ref_frame=0
    ).run(step=config.STRIDE)

    rmsd_data = rmsd_analysis.results.rmsd
    time_ns   = np.arange(len(rmsd_data)) * config.DT * config.STRIDE
    rmsd_vals = rmsd_data[:, 2]
    mean_r, std_r, max_r = np.mean(rmsd_vals), np.std(rmsd_vals), np.max(rmsd_vals)

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(time_ns, rmsd_vals, color='#1f77b4', linewidth=1.5, label='Backbone RMSD')
    ax.axhline(y=mean_r, color='red', linestyle='--', alpha=0.8,
               label=f'Mean: {mean_r:.2f} Å')
    ax.fill_between(time_ns, mean_r-std_r, mean_r+std_r,
                    alpha=0.15, color='red', label=f'±1 SD ({std_r:.2f} Å)')
    textstr = f'Mean:  {mean_r:.2f} Å\nSD:    {std_r:.2f} Å\nMax:   {max_r:.2f} Å'
    props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5)
    ax.text(0.02, 0.95, textstr, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', bbox=props)
    ax.set_xlabel('Time (ns)', fontsize=13)
    ax.set_ylabel('RMSD (Å)', fontsize=13)
    ax.set_title('Global Protein Backbone RMSD', fontsize=15, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('10_protein_rmsd_global.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved: 07_protein_rmsd_global.png")
    print(f"  Mean: {mean_r:.2f} ± {std_r:.2f} Å  |  Max: {max_r:.2f} Å")
    return time_ns, rmsd_vals


def calculate_domain_rmsd(universe, ref, config):
    print("\n[8/9] Per-Domain Cα RMSD...")

    fig, ax = plt.subplots(figsize=(12, 5))
    domain_results = {}

    for domain_name, (res_start, res_end) in config.DOMAINS.items():
        sel_str = f'name CA and resid {res_start}:{res_end}'
        protein_sel = universe.select_atoms(sel_str)
        ref_sel     = ref.select_atoms(sel_str)

        if len(protein_sel) == 0:
            print(f"  ⚠ No atoms for {domain_name}, skipping")
            continue

        try:
            rmsd_analysis = rms.RMSD(
                protein_sel, ref_sel, select=sel_str, ref_frame=0
            ).run(step=config.STRIDE)

            rmsd_data = rmsd_analysis.results.rmsd
            time_ns   = np.arange(len(rmsd_data)) * config.DT * config.STRIDE
            rmsd_vals = rmsd_data[:, 2]
            domain_results[domain_name] = (time_ns, rmsd_vals)

            label_suffix = ' ← N347K' if domain_name == 'Domain I' else ''
            ax.plot(time_ns, rmsd_vals,
                    color=config.DOMAIN_COLORS.get(domain_name, 'gray'),
                    linewidth=1.5, alpha=0.85,
                    label=f'{domain_name}{label_suffix}  (mean: {np.mean(rmsd_vals):.2f} Å)')
            print(f"  {domain_name}: {np.mean(rmsd_vals):.2f} ± {np.std(rmsd_vals):.2f} Å")

        except Exception as e:
            print(f"  ⚠ {domain_name} failed: {e}")

    ax.set_xlabel('Time (ns)', fontsize=13)
    ax.set_ylabel('Cα RMSD (Å)', fontsize=13)
    ax.set_title('Per-Domain Cα RMSD\n(N347K mutation is in Domain I)',
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('09_protein_rmsd_domains.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved: 08_protein_rmsd_domains.png")
    return domain_results


def calculate_ploop_rmsd(universe, ref, config):
    print("\n[9/9] P-loop RMSD...")

    res_start, res_end = config.PLOOP_RANGE
    sel_str = f'name CA and resid {res_start}:{res_end}'
    protein_sel = universe.select_atoms(sel_str)
    ref_sel     = ref.select_atoms(sel_str)

    if len(protein_sel) == 0:
        print(f"  ⚠ No atoms found for resid {res_start}-{res_end}")
        return None

    rmsd_analysis = rms.RMSD(
        protein_sel, ref_sel, select=sel_str, ref_frame=0
    ).run(step=config.STRIDE)

    rmsd_data = rmsd_analysis.results.rmsd
    time_ns   = np.arange(len(rmsd_data)) * config.DT * config.STRIDE
    rmsd_vals = rmsd_data[:, 2]
    mean_r, std_r, max_r = np.mean(rmsd_vals), np.std(rmsd_vals), np.max(rmsd_vals)

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(time_ns, rmsd_vals, color='#9467bd', linewidth=1.5,
            label=f'P-loop Cα RMSD (resid {res_start}–{res_end})')
    ax.axhline(y=mean_r, color='darkred', linestyle='--', alpha=0.8,
               label=f'Mean: {mean_r:.2f} Å')
    ax.fill_between(time_ns, mean_r-std_r, mean_r+std_r,
                    alpha=0.15, color='darkred', label=f'±1 SD ({std_r:.2f} Å)')
    textstr = f'Mean:  {mean_r:.2f} Å\nSD:    {std_r:.2f} Å\nMax:   {max_r:.2f} Å'
    props = dict(boxstyle='round', facecolor='lavender', alpha=0.5)
    ax.text(0.02, 0.95, textstr, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', bbox=props)
    ax.set_xlabel('Time (ns)', fontsize=13)
    ax.set_ylabel('Cα RMSD (Å)', fontsize=13)
    ax.set_title(
        f'P-loop Region Cα RMSD (resid {res_start}–{res_end})\n'
        f'Mutation Site: N347K (resid {config.MUTATION_RESID})',
        fontsize=14, fontweight='bold'
    )
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('11_protein_rmsd_ploop.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved: 09_protein_rmsd_ploop.png")
    print(f"  P-loop Mean: {mean_r:.2f} ± {std_r:.2f} Å  |  Max: {max_r:.2f} Å")
    return time_ns, rmsd_vals


# ============================================
# SUMMARY
# ============================================
def print_summary(time_array, site_data, pore_dists, global_rmsd, ploop_rmsd):
    print_section("ANALYSIS COMPLETE - SUMMARY")
    print(f"Total simulation time: {time_array[-1]:.1f} ns")

    if global_rmsd is not None:
        _, rmsd_vals = global_rmsd
        mean_r = np.mean(rmsd_vals)
        print(f"Global RMSD: {mean_r:.2f} ± {np.std(rmsd_vals):.2f} Å", end="  ")
        if mean_r < 3.0:
            print("→ STABLE ✓")
        elif mean_r < 5.0:
            print("→ Moderate drift (acceptable for membrane protein)")
        else:
            print("→ ⚠ High drift — check trajectory!")

    if ploop_rmsd is not None:
        _, ploop_vals = ploop_rmsd
        print(f"P-loop RMSD: {np.mean(ploop_vals):.2f} ± {np.std(ploop_vals):.2f} Å")

    if site_data and 'K347' in site_data:
        print(f"Mex-K347 mean distance: {np.nanmean(site_data['K347']):.2f} Å")

    if pore_dists is not None:
        min_pore = np.min(pore_dists)
        print(f"Closest approach to pore: {min_pore:.2f} Å", end="  ")
        if min_pore < 15:
            print("→ MIGRATION RISK — check 02_ligand_pore_distance.png")
        else:
            print("→ No migration detected ✓")

    print("\nDone!")


# ============================================
# MAIN
# ============================================
def main():
    print_section("MD ANALYSIS PIPELINE: Nav1.5")

    config     = Config()
    universe   = load_universe(config)
    ref        = load_reference(config)
    selections = setup_selections(universe, config)

    n_frames   = len(universe.trajectory)
    time_array = np.arange(n_frames) * config.DT

    # Ligand analyses
    calculate_rmsf(universe, selections, config)
    pore_dists = calculate_pore_distance(universe, selections, config, time_array)
    calculate_contact_map(universe, selections, config)
    calculate_ligand_rmsd(universe, selections)
    site_data  = calculate_multisite_distances(universe, selections, config, time_array)
    calculate_density(universe, selections)

    # Protein RMSD
    global_rmsd = calculate_global_rmsd(universe, ref, config)
    calculate_domain_rmsd(universe, ref, config)
    ploop_rmsd  = calculate_ploop_rmsd(universe, ref, config)

    print_summary(time_array, site_data, pore_dists, global_rmsd, ploop_rmsd)


if __name__ == "__main__":
    main()
