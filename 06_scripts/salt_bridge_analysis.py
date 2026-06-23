"""
K347 - D356 Salt Bridge Distance Time-Series
Author: Yagmur Yesilyurt
Description: Measures the distance between Lys347 (NZ) and Asp356 (OD1/OD2) 
             to track salt bridge stability over time.
"""

import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import matplotlib.pyplot as plt

systems = [
    {
        "topology": "03_md_simulations/mutant/apo/N347K_500NS/analysis_result.pdb",
        "trajectory": "03_md_simulations/mutant/apo/N347K_500NS/analysis_result.dcd",
        "name": "Apo Mutant", "color": "black"
    },
    {
        "topology": "03_md_simulations/mutant/holo/N347K_MEX_PROTONATED_500NS/analysis_result.pdb",
        "trajectory": "03_md_simulations/mutant/holo/N347K_MEX_PROTONATED_500NS/analysis_result.dcd",
        "name": "Holo Mutant (Protonated)", "color": "red"
    }
]

plt.figure(figsize=(12, 5))

for sys in systems:
    if not os.path.exists(sys["topology"]): continue
    u = mda.Universe(sys["topology"], sys["trajectory"])
    
    # Lysine pozitif ucu (NZ) ve Aspartat negatif uçları (OD1, OD2)
    lys347 = u.select_atoms("resid 347 and name NZ")
    asp356 = u.select_atoms("resid 356 and (name OD1 or name OD2)")
    
    bridge_distances = []
    frames = []
    
    for ts in u.trajectory:
        d = distances.distance_array(lys347.positions, asp356.positions).min()
        bridge_distances.append(d)
        frames.append(ts.frame)
        
    plt.plot(frames, bridge_distances, label=sys["name"], color=sys["color"], alpha=0.6)

# Tuz köprüsü üst sınırı (Biyofizikte genellikle 4.0 Angstrom kabul edilir)
plt.axhline(y=4.0, color="green", linestyle="--", label="Salt Bridge Threshold (4.0 Å)")

plt.title("K347 - D356 Salt Bridge Dynamics Over Time", fontsize=13)
plt.xlabel("Frame Index", fontsize=11)
plt.ylabel("Minimum Distance ($\AA$)", fontsize=11)
plt.ylim(2, 8)
plt.grid(True, linestyle=":", alpha=0.5)
plt.legend(loc="upper right")
plt.tight_layout()

output_png = "04_postMD_analysis_results/24_salt_bridge_dynamics.png"
os.makedirs(os.path.dirname(output_png), exist_ok=True)
plt.savefig(output_png, dpi=300)
print(f"[SUCCESS] Tuz köprüsü grafiği kaydedildi: {output_png}")