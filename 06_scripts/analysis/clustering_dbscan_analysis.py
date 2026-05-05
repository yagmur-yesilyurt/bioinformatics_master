import subprocess
from clustering import get_most_probable_conformation
def run_analysis(base, out_base, pos, i):
    command = ['python', f'{base}/functions/analyse.py', '-p', f'{base}/output/complexes/seq_{pos}_{i}/{out_base}_{pos}_{i}_minimised.pdb', '-t', f'{base}/output/complexes/seq_{pos}_{i}/{out_base}_{pos}_{i}_traj.dcd', '-o', f'{base}/output/complexes/seq_{pos}_{i}/{out_base}_{pos}_{i}_reimaged', '-r']
    subprocess.call(command)
    
def reimage_trajectory(base, out_base, pos, i):
    command = ['mdconvert', f'{base}/output/complexes/seq_{pos}_{i}/{out_base}_{pos}_{i}_reimaged.dcd', '-o', f'{base}/output/complexes/seq_{pos}_{i}/{out_base}_{pos}_{i}_traj_reimaged.pdb', '-t', f'{base}/output/complexes/seq_{pos}_{i}/{out_base}_{pos}_{i}_reimaged.pdb']
    subprocess.call(command)
    
    
"""main"""    
run_analysis(args.base_dir, args.out_base, pos, i)
reimage_trajectory(args.base_dir, args.out_base, pos, i)

dock_dir = f'{DOCKING_IN_DIR}/seq_{pos}_{i}'
if not os.path.exists(dock_dir):
    os.makedirs(dock_dir)

# --- Get representative frame & split chains ---
representative_frame, _, _ = get_most_probable_conformation(
    f'{seq_dir}/{args.out_base}_{pos}_{i}_traj_reimaged.pdb',
    pos, i, output_path=seq_dir, subsample_rate=20
    )