"""
Per-Residue Root Mean Square Fluctuation (RMSF) Analysis Script
Author: Yagmur Yesilyurt
Date: June 2026
Description: Measures structural flexibility changes across Apo Mutant, 
             Holo Mutant, and Holo WT to prove the allosteric stabilization.
"""

import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
import matplotlib.pyplot as plt

def calculate_rmsf(topology, trajectory, system_name):
    print(f"-> Processing RMSF for: {system_name}")
    if not os.path.exists(topology) or not os.path.exists(trajectory):
        print(f"[ERROR] Files missing for {system_name}. Skipping...")
        return None, None
        
    u = mda.Universe(topology, trajectory)
    ref = mda.Universe(topology, trajectory)
    
    # 1. Align trajectory to remove global translations/rotations
    alignment = align.AlignTraj(u, ref, select="backbone", in_memory=True)
    alignment.run()
    
    # 2. Select C-alpha atoms of the entire protein for a complete profile
    calphas = u.select_atoms("protein and name CA")
    
    # 3. Compute RMSF
    rmsf_analysis = rms.RMSF(calphas).run()
    
    return calphas.resids, rmsf_analysis.results.rmsf

if __name__ == "__main__":
    # Same exact paths from our successful DCCM run
    systems = [
        {
            "topology": "03_md_simulations/mutant/apo/N347K_500NS/analysis_result.pdb",
            "trajectory": "03_md_simulations/mutant/apo/N347K_500NS/analysis_result.dcd",
            "name": "Apo N347K Mutant",
            "color": "black"
        },
        {
            "topology": "03_md_simulations/mutant/holo/N347K_MEX_PROTONATED_500NS/analysis_result.pdb",
            "trajectory": "03_md_simulations/mutant/holo/N347K_MEX_PROTONATED_500NS/analysis_result.dcd",
            "name": "Holo N347K Mutant (Prot. Mexiletine)",
            "color": "red"
        },
        {
            "topology": "03_md_simulations/wt/holo/Nav15_MEX_PROTONATED_500NS/analysis_result.pdb",
            "trajectory": "03_md_simulations/wt/holo/Nav15_MEX_PROTONATED_500NS/analysis_result.dcd",
            "name": "Holo WT (Prot. Mexiletine)",
            "color": "blue"
        }
    ]
    
    plt.figure(figsize=(14, 6))
    
    for sys in systems:
        resids, rmsf_values = calculate_rmsf(sys["topology"], sys["trajectory"], sys["name"])
        if rmsf_values is not None:
            plt.plot(resids, rmsf_values, label=sys["name"], color=sys["color"], alpha=0.7, linewidth=1.5)
            
    plt.title("Per-Residue C-alpha Root Mean Square Fluctuation (RMSF) Comparison", fontsize=14, pad=15)
    plt.xlabel("Residue Number", fontsize=12)
    plt.ylabel("RMSF ($\AA$)", fontsize=12)
    
    # Highlight our functional areas on the plot
    plt.axvspan(330, 370, color='yellow', alpha=0.2, label='Domain I P-Loop')
    plt.axvline(x=1423, color='purple', linestyle='--', alpha=0.5, label='Site 2 Residues (1423, 1714)')
    plt.axvline(x=1714, color='purple', linestyle='--', alpha=0.5)
    
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend(loc="upper right", fontsize=10)
    plt.tight_layout()
    
    output_png = "04_postMD_analysis_results/20_global_rmsf_comparison.png"
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    plt.savefig(output_png, dpi=300)
    print(f"\n[SUCCESS] Comparative RMSF plot saved to: {output_png}")