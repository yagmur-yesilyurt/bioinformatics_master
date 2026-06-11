from pathlib import Path
import subprocess

LIGAND = '/home/yagmur/Desktop/bioinformatics_master/01_structures/ligands/mexiletine/ligand_mexiletine_protonated.pdbqt'
RECEPTOR_PDBQT = '/home/yagmur/Desktop/bioinformatics_master/04_postMD_analysis_results/holoMutant_postMD/notr_mexiletine/rep_holoMutant_nolig_model0.pdbqt'
RECEPTOR_PDB = '/home/yagmur/Desktop/bioinformatics_master/04_postMD_analysis_results/holoMutant_postMD/notr_mexiletine/rep_holoMutant_nolig.pdb'
OUT_DIR = Path('/home/yagmur/Desktop/bioinformatics_master/02_docking/mexiletine/blind/Mutant/N347K_holoMutant_notrMD_rep_frame')
EXHAUSTIVENESS = 128
NUM_MODES = 20
ENERGY_RANGE = 5

def compute_box(pdb_file):
    xs, ys, zs = [], [], []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
    center = ((min(xs)+max(xs))/2, (min(ys)+max(ys))/2, (min(zs)+max(zs))/2)
    size = (max(xs)-min(xs)+10, max(ys)-min(ys)+10, max(zs)-min(zs)+10)
    return center, size

OUT_DIR.mkdir(parents=True, exist_ok=True)
center, size = compute_box(RECEPTOR_PDB)
print('Center:', center)
print('Size:', size)
config = OUT_DIR / 'config.txt'
with open(config, 'w') as f:
    f.write(f'receptor = {RECEPTOR_PDBQT}\n')
    f.write(f'ligand = {LIGAND}\n')
    f.write(f'center_x = {center[0]:.3f}\n')
    f.write(f'center_y = {center[1]:.3f}\n')
    f.write(f'center_z = {center[2]:.3f}\n')
    f.write(f'size_x = {size[0]:.3f}\n')
    f.write(f'size_y = {size[1]:.3f}\n')
    f.write(f'size_z = {size[2]:.3f}\n')
    f.write(f'exhaustiveness = {EXHAUSTIVENESS}\n')
    f.write(f'num_modes = {NUM_MODES}\n')
    f.write(f'energy_range = {ENERGY_RANGE}\n')
    f.write(f'out = {OUT_DIR}/docked.pdbqt\n')
    f.write(f'log = {OUT_DIR}/log.txt\n')
subprocess.run(['vina', '--config', str(config)], check=True)
print('Done!')