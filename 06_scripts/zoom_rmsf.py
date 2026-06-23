"""
Zoomed RMSF Verification Script
Author: Yagmur Yesilyurt
Description: Zooms into specific functional regions to eliminate the macro-compression artifact
             and prints precise numerical averages.
"""

import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
import matplotlib.pyplot as plt

def get_rmsf_data(topology, trajectory):
    if not os.path.exists(topology) or not os.path.exists(trajectory):
        return None, None
    u = mda.Universe(topology, trajectory)
    ref = mda.Universe(topology, trajectory)
    alignment = align.AlignTraj(u, ref, select="backbone", in_memory=True)
    alignment.run()
    calphas = u.select_atoms("protein and name CA")
    rmsf_analysis = rms.RMSF(calphas).run()
    return calphas.resids, rmsf_analysis.results.rmsf

if __name__ == "__main__":
    systems = [
        {"topology": "03_md_simulations/mutant/apo/N347K_500NS/analysis_result.pdb",
         "trajectory": "03_md_simulations/mutant/apo/N347K_500NS/analysis_result.dcd",
         "name": "Apo Mutant", "color": "black"},
        {"topology": "03_md_simulations/mutant/holo/N347K_MEX_PROTONATED_500NS/analysis_result.pdb",
         "trajectory": "03_md_simulations/mutant/holo/N347K_MEX_PROTONATED_500NS/analysis_result.dcd",
         "name": "Holo Mutant", "color": "red"},
        {"topology": "03_md_simulations/wt/holo/Nav15_MEX_PROTONATED_500NS/analysis_result.pdb",
         "trajectory": "03_md_simulations/wt/holo/Nav15_MEX_PROTONATED_500NS/analysis_result.dcd",
         "name": "Holo WT", "color": "blue"}
    ]
    
    data = {}
    for sys in systems:
        resids, rmsf = get_rmsf_data(sys["topology"], sys["trajectory"])
        if rmsf is not None:
            data[sys["name"]] = {"resids": resids, "rmsf": rmsf, "color": sys["color"]}
            
    # 1. PRINT EXACT NUMERICAL AVERAGES
    print("\n" + "="*50)
    print("      CRITICAL REGIONS NUMERICAL RMSF (Å)")
    print("="*50)
    for name, d in data.items():
        resids = d["resids"]
        rmsf = d["rmsf"]
        
        # DI P-loop average (330-370)
        ploop_mask = (resids >= 330) & (resids <= 370)
        ploop_avg = np.mean(rmsf[ploop_mask])
        
        # Site 2 critical residues (1423 and 1714)
        site2_mask = np.isin(resids, [1423, 1714])
        site2_avg = np.mean(rmsf[site2_mask])
        
        print(f"{name:12s} -> DI P-loop (330-370) Mean: {ploop_avg:.3f} Å | Site 2 (1423/1714) Mean: {site2_avg:.3f} Å")
    print("="*50)
        
    # 2. PLOT ZOOMED SUBPLOTS
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
    
    # Left Panel: Zoom on Domain I P-loop
    for name, d in data.items():
        mask = (d["resids"] >= 325) & (d["resids"] <= 375)
        ax1.plot(d["resids"][mask], d["rmsf"][mask], label=name, color=d["color"], linewidth=2.5)
    ax1.set_title("Zoomed View: Domain I P-Loop (325-375)", fontsize=12)
    ax1.set_xlabel("Residue Number")
    ax1.set_ylabel("RMSF (Å)")
    ax1.grid(True, linestyle=":")
    ax1.legend()
    
    # Right Panel: Zoom on Site 2 Region (around ASP1423)
    for name, d in data.items():
        mask = (d["resids"] >= 1410) & (d["resids"] <= 1440)
        ax2.plot(d["resids"][mask], d["rmsf"][mask], label=name, color=d["color"], linewidth=2.5)
    ax2.set_title("Zoomed View: Site 2 / DIII-S6 Loop (1410-1440)", fontsize=12)
    ax2.set_xlabel("Residue Number")
    ax2.set_ylabel("RMSF (Å)")
    ax2.grid(True, linestyle=":")
    ax2.legend()
    
    plt.tight_layout()
    output_png = "04_postMD_analysis_results/21_zoomed_rmsf_verification.png"
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    plt.savefig(output_png, dpi=300)
    print(f"\n[SUCCESS] Zoomed visualization saved to: {output_png}\n")