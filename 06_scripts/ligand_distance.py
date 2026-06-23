"""
Tailored Ligand-Residue Distance Time-Series Script
Author: Yagmur Yesilyurt
Description: Measures exact distances between Mexiletine (UNL) and Site 2 residues.
"""

import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import matplotlib.pyplot as plt

topology = "03_md_simulations/mutant/holo/N347K_MEX_PROTONATED_500NS/analysis_result.pdb"
trajectory = "03_md_simulations/mutant/holo/N347K_MEX_PROTONATED_500NS/analysis_result.dcd"

u = mda.Universe(topology, trajectory)

# Net seçimler: Ligandın adı UNL olarak belirlendi
ligand = u.select_atoms("resname UNL")
asp1423 = u.select_atoms("resid 1423")
asp1714 = u.select_atoms("resid 1714")

print(f"-> Atom Sayıları: UNL: {len(ligand)} | ASP1423: {len(asp1423)} | ASP1714: {len(asp1714)}")

time_series_1423 = []
time_series_1714 = []
frames = []

print("-> Mesafeler hesaplanıyor...")
for ts in u.trajectory:
    d_1423 = distances.distance_array(ligand.positions, asp1423.positions).min()
    d_1714 = distances.distance_array(ligand.positions, asp1714.positions).min()
    
    time_series_1423.append(d_1423)
    time_series_1714.append(d_1714)
    frames.append(ts.frame)

# Grafik Çizimi
plt.figure(figsize=(12, 5))
plt.plot(frames, time_series_1423, label="Distance to ASP1423", color="crimson", alpha=0.8)
plt.plot(frames, time_series_1714, label="Distance to ASP1714", color="darkblue", alpha=0.8)
plt.axhline(y=4.0, color="black", linestyle="--", label="Interaction Threshold (4.0 Å)")

plt.title("Mexiletine (UNL) Bound State Stability at Site 2", fontsize=13)
plt.xlabel("Frame Index", fontsize=11)
plt.ylabel("Minimum Distance ($\AA$)", fontsize=11)
plt.ylim(0, 15)
plt.grid(True, linestyle=":", alpha=0.6)
plt.legend(loc="upper right")
plt.tight_layout()

output_png = "04_postMD_analysis_results/23_mexiletine_site2_distance.png"
os.makedirs(os.path.dirname(output_png), exist_ok=True)
plt.savefig(output_png, dpi=300)
print(f"[SUCCESS] Grafik başarıyla kaydedildi: {output_png}")