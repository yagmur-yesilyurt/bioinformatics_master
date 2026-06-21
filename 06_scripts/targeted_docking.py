import numpy as np
from pathlib import Path
import subprocess

# ============================================================
# CONFIG — sadece burası değişir
# ============================================================

POCKETS_DIR = "/home/yagmur/Desktop/bioinformatics_master/07_binding_site_analysis/fpocket/apoMutant_MD_rep_frame/pockets"
RECEPTOR_PDBQT = "/home/yagmur/Desktop/bioinformatics_master/02_docking/mexiletine/blind/Mutant/N347K_apoMutant_MD_rep_frame/rep_structure_full_10k_model0.pdbqt"
LIGAND = "/home/yagmur/Desktop/bioinformatics_master/01_structures/ligands/mexiletine/ligand_mexiletine_protonated.pdbqt"
BASE_OUT = "/home/yagmur/Desktop/bioinformatics_master/02_docking/mexiletine/targeted/apoMutant_MD_rep_frame_fpocket"

BOX_SIZE = 25
EXHAUSTIVENESS = 128
NUM_MODES = 20
ENERGY_RANGE = 5

pockets = {
    "pocket2_DomainIV_S6": "pocket2",
    "pocket13_DomainIII_IV_junction": "pocket13",
}

# ============================================================

for name, pocket_id in pockets.items():
    print(f"\n=== Targeted Docking: {name} ===")

    xs, ys, zs = [], [], []
    with open(f"{POCKETS_DIR}/{pocket_id}_vert.pqr") as f:
        for line in f:
            if line.startswith("ATOM"):
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))

    cx, cy, cz = np.mean(xs), np.mean(ys), np.mean(zs)
    print(f"Center: ({cx:.2f}, {cy:.2f}, {cz:.2f})")

    out_dir = Path(f"{BASE_OUT}/{name}")
    out_dir.mkdir(parents=True, exist_ok=True)

    config = out_dir / "config.txt"
    with open(config, "w") as f:
        f.write(f"receptor = {RECEPTOR_PDBQT}\n")
        f.write(f"ligand = {LIGAND}\n")
        f.write(f"center_x = {cx:.3f}\n")
        f.write(f"center_y = {cy:.3f}\n")
        f.write(f"center_z = {cz:.3f}\n")
        f.write(f"size_x = {BOX_SIZE}\n")
        f.write(f"size_y = {BOX_SIZE}\n")
        f.write(f"size_z = {BOX_SIZE}\n")
        f.write(f"exhaustiveness = {EXHAUSTIVENESS}\n")
        f.write(f"num_modes = {NUM_MODES}\n")
        f.write(f"energy_range = {ENERGY_RANGE}\n")
        f.write(f"out = {out_dir}/docked.pdbqt\n")
        f.write(f"log = {out_dir}/log.txt\n")

    print(f"Config saved: {config}")
    subprocess.run(["vina", "--config", str(config)], check=True)
    print(f"✅ {name} done!")