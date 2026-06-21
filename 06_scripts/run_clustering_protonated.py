from clustering_silvia_revised import get_most_probable_conformation
import sys
sys.path.insert(0, '/home/yyesilyurt/PHD_MD/Scripts')

trajectory_file = "analysis_result.dcd"
topology_file   = "analysis_result.pdb"

print("--- Revize Edilmiş Clustering: Ligand Filtreli (Domain IV) ---")

get_most_probable_conformation(
    trajectory_file=trajectory_file,
    topology_file=topology_file,
    output_path=".",
    selection='protein and backbone',
    eps=0.53,
    min_samples=5,
    subsample_rate=1,
    ligand_resname='UNL',
    ploop_resids=[1714, 1423, 1424, 1400, 1420],
    ligand_cutoff_angstrom=5.0
)
