"""
Dynamic Cross-Correlation Matrix (DCCM) Analysis Script
Author: Yagmur Yesilyurt
Date: June 2026
Description: Analyzes allosteric communication between Domain I P-loop (resid 330-370)
             and the Domain III/IV Junction Site 2 (resid 1400, 1423, 1424, 1714)
             across Apo Mutant, Holo Mutant, and Holo Wild-Type systems.
"""

import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
import matplotlib.pyplot as plt
import seaborn as sns

def run_dccm_analysis(topology, trajectory, output_png, system_name):
    print(f"\n==================================================")
    print(f"Starting DCCM Analysis for: {system_name}")
    print(f"==================================================")
    
    # 1. Check if files exist
    if not os.path.exists(topology) or not os.path.exists(trajectory):
        print(f"[ERROR] Input files missing for {system_name}. Skipping...")
        print(f"Expected Topology: {topology}")
        print(f"Expected Trajectory: {trajectory}")
        return

    # 2. Load the Universe
    print("-> Loading topology and trajectory...")
    u = mda.Universe(topology, trajectory)
    ref = mda.Universe(topology, trajectory) # Reference structure for alignment
    
    # 3. Trajectory Alignment (RMSD fitting to remove global rotation/translation)
    print("-> Aligning trajectory to the initial frame (backbone based)...")
    alignment = align.AlignTraj(u, ref, select="backbone", in_memory=True)
    alignment.run()
    
    # 4. Atom Selection (C-alpha atoms for accurate backbone dynamics)
    selection_str = (
        "protein and name CA and ("
        "resid 330:370 or "
        "resid 1400 or resid 1423 or resid 1424 or "
        "resid 1714"
    )
    # Handle the fact that residue numbering might slightly vary between WT and Mutant if needed, 
    # but based on your report, identical human Nav1.5 numbering (8VYK) is preserved.
    selection_str += ")"
    
    selection = u.select_atoms(selection_str)
    n_atoms = len(selection)
    resids = selection.resids
    print(f"-> Selected {n_atoms} C-alpha atoms for allosteric network analysis.")
    
    # 5. Extract Coordinates across Trajectory
    stride = 5
    print(f"-> Extracting coordinates (using a stride of {stride} frames)...")
    frames = u.trajectory[::stride]
    n_frames = len(frames)
    
    positions = np.zeros((n_frames, n_atoms, 3))
    for idx, ts in enumerate(frames):
        positions[idx] = selection.positions
        
    # 6. Calculate Atomic Displacements (Delta r = r - r_mean)
    mean_positions = np.mean(positions, axis=0)
    displacements = positions - mean_positions
    
    # 7. Compute Covariance and Normalization to get Correlation Coefficients (C_ij)
    print("-> Computing atomic cross-correlation matrix...")
    covariance = np.zeros((n_atoms, n_atoms))
    for i in range(3): # Loop over X, Y, Z Cartesian coordinates
        covariance += np.dot(displacements[:, :, i].T, displacements[:, :, i]) / n_frames
        
    dccm = np.zeros((n_atoms, n_atoms))
    for i in range(n_atoms):
        for j in range(n_atoms):
            dccm[i, j] = covariance[i, j] / np.sqrt(covariance[i, i] * covariance[j, j])
            
    # 8. Visualizing and Plotting the Heatmap
    print("-> Generating high-resolution publication-quality heatmap...")
    plt.figure(figsize=(12, 10))
    
    sns.heatmap(
        dccm, 
        xticklabels=resids, 
        yticklabels=resids, 
        cmap="RdBu_r", 
        vmin=-1.0, 
        vmax=1.0, 
        center=0.0,
        cbar_kws={'label': 'Correlation Coefficient ($C_{ij}$)'}
    )
    
    plt.title(f"Dynamic Cross-Correlation Matrix ({system_name})\nAllosteric Coupling: DI P-Loop vs. DIII/DIV Junction", fontsize=14, pad=15)
    plt.xlabel("Residue Number", fontsize=12)
    plt.ylabel("Residue Number", fontsize=12)
    
    plt.xticks(rotation=90, fontsize=9)
    plt.yticks(rotation=0, fontsize=9)
    plt.tight_layout()
    
    # Create target output folder if it doesn't exist
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    plt.savefig(output_png, dpi=300)
    plt.close()
    print(f"[SUCCESS] DCCM matrix saved to: {output_png}")

if __name__ == "__main__":
    # Define directory paths relative to the project root folder (bioinformatics_master)
    systems_to_analyze = [
        {
            "topology": "03_md_simulations/mutant/apo/N347K_500NS/analysis_result.pdb",
            "trajectory": "03_md_simulations/mutant/apo/N347K_500NS/analysis_result.dcd",
            "output_png": "04_postMD_analysis_results/apoMutant_postMD/19_allosteric_dccm_apoMutant.png",
            "system_name": "Apo N347K Mutant"
        },
        {
            "topology": "03_md_simulations/mutant/holo/N347K_MEX_PROTONATED_500NS/analysis_result.pdb",
            "trajectory": "03_md_simulations/mutant/holo/N347K_MEX_PROTONATED_500NS/analysis_result.dcd",
            "output_png": "04_postMD_analysis_results/holoMutant_postMD/protonated_mexiletine/19_allosteric_dccm_holoMutant_protonated.png",
            "system_name": "Protonated Mexiletine Holo N347K Mutant"
        },
        {
            "topology": "03_md_simulations/wt/holo/Nav15_MEX_PROTONATED_500NS/analysis_result.pdb",
            "trajectory": "03_md_simulations/wt/holo/Nav15_MEX_PROTONATED_500NS/analysis_result.dcd",
            "output_png": "04_postMD_analysis_results/holoWT_postMD/protonated_mexiletine/19_allosteric_dccm_holoWT_protonated.png",
            "system_name": "Protonated Mexiletine Holo Wild-Type"
        }
    ]
    
    for system in systems_to_analyze:
        run_dccm_analysis(
            topology=system["topology"],
            trajectory=system["trajectory"],
            output_png=system["output_png"],
            system_name=system["system_name"]
        )
    print("\n[FINISH] All DCCM matrix plots successfully generated.")