"""
Site 2 Loop (1430-1440) Frame-by-Frame RMSD Distribution & Boxplot
Author: Yagmur Yesilyurt
Description: Computes the structural deviation distribution of the Site 2 loop 
             over time to provide rigid statistical proof of drug clamping.
"""

import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def get_loop_rmsd_distribution(topology, trajectory, system_name):
    print(f"-> Extracting loop RMSD distribution for: {system_name}")
    if not os.path.exists(topology) or not os.path.exists(trajectory):
        print(f"[ERROR] Files missing for {system_name}. Skipping...")
        return None
        
    u = mda.Universe(topology, trajectory)
    
    # 1. Align the entire trajectory to backbone to remove global translation/rotation
    ref = mda.Universe(topology, trajectory)
    alignment = align.AlignTraj(u, ref, select="backbone", in_memory=True)
    alignment.run()
    
    # 2. Select the specific loop region under the microscope (Residues 1430-1440)
    loop_selection = "protein and name CA and resid 1430:1440"
    
    # 3. Compute Frame-by-Frame RMSD against the average structure of this trajectory
    # First, we calculate the average structure coordinates for the selection
    loop_atoms = u.select_atoms(loop_selection)
    n_frames = len(u.trajectory)
    
    # Calculate mean positions
    loop_coords = np.zeros((n_frames, len(loop_atoms), 3))
    for ts in u.trajectory:
        loop_coords[ts.frame] = loop_atoms.positions
    mean_coords = np.mean(loop_coords, axis=0)
    
    # Compute RMSD for each frame relative to this mean conformation
    frame_rmsd = []
    for ts in u.trajectory:
        current_coords = loop_coords[ts.frame]
        # Mass-weighted or unweighted standard RMSD calculation
        diff = current_coords - mean_coords
        rmsd_val = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        frame_rmsd.append(rmsd_val)
        
    return frame_rmsd

if __name__ == "__main__":
    # Exact paths from your verified directory tree
    systems = [
        {
            "topology": "03_md_simulations/mutant/apo/N347K_500NS/analysis_result.pdb",
            "trajectory": "03_md_simulations/mutant/apo/N347K_500NS/analysis_result.dcd",
            "name": "Apo Mutant"
        },
        {
            "topology": "03_md_simulations/mutant/holo/N347K_MEX_PROTONATED_500NS/analysis_result.pdb",
            "trajectory": "03_md_simulations/mutant/holo/N347K_MEX_PROTONATED_500NS/analysis_result.dcd",
            "name": "Holo Mutant"
        },
        {
            "topology": "03_md_simulations/wt/holo/Nav15_MEX_PROTONATED_500NS/analysis_result.pdb",
            "trajectory": "03_md_simulations/wt/holo/Nav15_MEX_PROTONATED_500NS/analysis_result.dcd",
            "name": "Holo WT"
        }
    ]
    
    all_data = []
    
    # Gather data into a Pandas DataFrame for Seaborn plotting
    for sys in systems:
        rmsd_dist = get_loop_rmsd_distribution(sys["topology"], sys["trajectory"], sys["name"])
        if rmsd_dist is not None:
            df_sys = pd.DataFrame({
                "RMSD": rmsd_dist,
                "System": sys["name"]
            })
            all_data.append(df_sys)
            
            # Print descriptive statistics to terminal
            print(f"--- Statistics for {sys['name']} ---")
            print(f"  Mean RMSD: {np.mean(rmsd_dist):.3f} Å")
            print(f"  Std Dev:   {np.std(rmsd_dist):.3f} Å")
            print(f"  Median:    {np.median(rmsd_dist):.3f} Å\n")

    master_df = pd.concat(all_data, ignore_index=True)
    
    # Plotting the Boxplot + Violin Plot overlay for complete distribution visibility
    plt.figure(figsize=(8, 6))
    
    palette = {"Apo Mutant": "gray", "Holo Mutant": "red", "Holo WT": "blue"}
    
    # Violin plot shows the density/shape, Boxplot shows the exact quartiles
    sns.violinplot(x="System", y="RMSD", data=master_df, palette=palette, inner=None, alpha=0.3)
    sns.boxplot(x="System", y="RMSD", data=master_df, palette=palette, width=0.3, linewidth=2, showfliers=False)
    
    plt.title("Structural Deviation (RMSD) Distribution of Site 2 Loop (Resid 1430-1440)", fontsize=12, pad=15)
    plt.xlabel("Simulation System", fontsize=11)
    plt.ylabel("Loop Conformational Deviation from Mean ($\AA$)", fontsize=11)
    plt.grid(True, linestyle=":", alpha=0.5)
    
    plt.tight_layout()
    output_png = "04_postMD_analysis_results/22_site2_loop_boxplot.png"
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    plt.savefig(output_png, dpi=300)
    print(f"[SUCCESS] Statistical boxplot saved to: {output_png}")