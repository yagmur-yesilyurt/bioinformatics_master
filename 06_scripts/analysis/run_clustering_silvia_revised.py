from clustering_revised import get_most_probable_conformation

# --- AYARLAR ---
trajectory_file = "analysis_result.dcd"
topology_file   = "analysis_result.pdb"

print("--- Revize Edilmiş Clustering: Ligand Filtreli ---")

get_most_probable_conformation(
    trajectory_file=trajectory_file,
    topology_file=topology_file,
    output_path=".",
    selection='protein and backbone',
    eps=0.53,
    min_samples=5,
    subsample_rate=1,
    # Yeni parametreler:
    ligand_resname='UNL',
    ploop_resids=[348, 349, 353, 355],  # PRO348, ASP349, THR353, PHE355
    ligand_cutoff_angstrom=5.0         # Bu frame'lerde mexiletine P-loop'a <10 Å
)
