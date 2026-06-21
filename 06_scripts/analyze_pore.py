"""
PORE DIMENSIONS ANALYSIS: Nav1.5 Channel
Multi-level pore diameter analysis across complete simulation

OUTPUT FILE:
------------
pore_dimensions_multilevel.png - Three-panel plot showing:
  • Level 1: Selectivity Filter (DEKA Motif)
  • Level 2: Drug Binding Pocket (Mid-Pore)
  • Level 3: Intracellular Gate (Activation)
"""

import MDAnalysis as mda
from MDAnalysis.analysis import distances
import matplotlib.pyplot as plt
import numpy as np


# ============================================
# CONFIGURATION
# ============================================
class Config:
    """Simulation parameters and file paths"""
    TOPOLOGY = "step5_input.pdb"
    N_FILES = 10
    DT_PER_FRAME = 0.1  # ns per frame
    OUTPUT_FILENAME = "07_pore_dimensions_multilevel.png"


# ============================================
# PORE MEASUREMENT DEFINITIONS
# ============================================
PORE_MEASUREMENTS = {
    "Activation Gate": [
        (1466, 1765, "I1466 - V1765 (Main)"),
        (1469, 1768, "A1469 - I1768 (Secondary)")
    ],
    "Mid-Pore (Drug Binding Pocket)": [
        (1462, 1760, "L1462 - F1760 (Mexiletine Site)"),
        (1409, 934,  "I1409 - L934 (Central Cavity)")
    ],
    "Upper Pore (Selectivity Filter)": [
        (372, 1419, "D372 - K1419 (Outer Mouth)"),
        (898, 1712, "E898 - A1712 (Filter Level)")
    ]
}

PLOT_CONFIG = {
    "Activation Gate": {
        "title": "Level 3: Intracellular Gate (Activation)",
        "ylim": (5, 20),
        "position": 2
    },
    "Mid-Pore (Drug Binding Pocket)": {
        "title": "Level 2: Drug Binding Pocket (Mid-Pore)",
        "ylim": None,
        "position": 1
    },
    "Upper Pore (Selectivity Filter)": {
        "title": "Level 1: Selectivity Filter (DEKA Motif)",
        "ylim": None,
        "position": 0
    }
}

COLORS = ['#ff7f0e', '#1f77b4']  # Orange and Blue


# ============================================
# UTILITY FUNCTIONS
# ============================================
def generate_trajectory_list(n_files):
    """Generate list of trajectory file names"""
    return [f"step7_{i}.dcd" for i in range(1, n_files + 1)]


def load_universe(topology, trajectories):
    """Load trajectory files into MDAnalysis Universe"""
    print(">>> Loading simulation files (1-10)...")
    try:
        u = mda.Universe(topology, trajectories)
        return u
    except Exception as e:
        print(f"✗ ERROR: Failed to load files!")
        print(f"  Details: {e}")
        exit()


def create_time_array(universe, dt_per_frame):
    """Create time array for entire simulation"""
    return np.arange(universe.trajectory.n_frames) * dt_per_frame


# ============================================
# DISTANCE CALCULATIONS
# ============================================
def calculate_ca_distance(universe, res1, res2):
    """Calculate CA-CA distance between two residues"""
    atom1 = universe.select_atoms(f"resid {res1} and name CA")
    atom2 = universe.select_atoms(f"resid {res2} and name CA")
    
    if len(atom1) > 0 and len(atom2) > 0:
        d = distances.distance_array(atom1.positions, atom2.positions)[0][0]
        return d
    else:
        return np.nan


def analyze_pore_dimensions(universe, measurements):
    """Calculate pore dimensions across trajectory"""
    print(f"\n>>> Starting analysis: Total {universe.trajectory.n_frames} frames")
    print(f"    Duration: {universe.trajectory.n_frames * 0.1:.1f} ns\n")
    
    results = {level: [] for level in measurements}
    
    for ts in universe.trajectory:
        # Progress indicator
        if ts.frame % 1000 == 0:
            print(f"Processing: {ts.time/1000:.0f} / 500 ns...", end="\r")
        
        for level_name, pairs in measurements.items():
            level_data = []
            for res1, res2, label in pairs:
                d = calculate_ca_distance(universe, res1, res2)
                level_data.append(d)
            results[level_name].append(level_data)
    
    print("\n\n✓ Data collection complete")
    
    # Convert to numpy arrays
    for level in results:
        results[level] = np.array(results[level])
    
    return results


# ============================================
# VISUALIZATION
# ============================================
def create_subplot(ax, time, data, pairs, config, colors):
    """Create individual subplot for pore level"""
    for i in range(len(pairs)):
        ax.plot(time, data[:, i],
                label=pairs[i][2],
                linewidth=1.2,
                alpha=0.9,
                color=colors[i])
    
    ax.set_ylabel("Distance (Å)", fontsize=12)
    ax.set_title(config['title'], fontsize=13, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    if config['ylim'] is not None:
        ax.set_ylim(config['ylim'])


def plot_results(time, results, measurements, output_filename):
    """Create multi-panel figure with all pore measurements"""
    print(">>> Generating figure...")
    
    fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    
    # Plot each level
    for level_name, config in PLOT_CONFIG.items():
        ax = axs[config['position']]
        data = results[level_name]
        pairs = measurements[level_name]
        create_subplot(ax, time, data, pairs, config, COLORS)
    
    # Set x-label only on bottom subplot
    axs[2].set_xlabel("Time (ns)", fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ COMPLETE! Figure saved to: {output_filename}")


# ============================================
# MAIN EXECUTION
# ============================================
def main():
    """Main analysis pipeline"""
    print("="*60)
    print("PORE DIMENSIONS ANALYSIS: Nav1.5 Channel")
    print("="*60)
    
    # Initialize configuration
    config = Config()
    
    # Generate trajectory list
    trajectories = generate_trajectory_list(config.N_FILES)
    
    # Load universe
    universe = load_universe(config.TOPOLOGY, trajectories)
    
    # Create time array
    time = create_time_array(universe, config.DT_PER_FRAME)
    
    # Analyze pore dimensions
    results = analyze_pore_dimensions(universe, PORE_MEASUREMENTS)
    
    # Create visualization
    plot_results(time, results, PORE_MEASUREMENTS, config.OUTPUT_FILENAME)
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)


if __name__ == "__main__":
    main()
