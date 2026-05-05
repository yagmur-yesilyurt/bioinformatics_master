import sys, argparse
import mdtraj as md
import plotly.graph_objects as go

parser = argparse.ArgumentParser(description="Memory Efficient Analyse")
parser.add_argument("-p", "--protein", required=True, help="Topology (PDB/PSF)")
parser.add_argument("-t", "--trajectory", nargs='+', required=True, help="DCD files")
parser.add_argument("-o", "--output", required=True, help="Output base name")
parser.add_argument("-r", "--remove-waters", action='store_true', help="Remove solvent/lipids")

args = parser.parse_args()

def process_memory_efficiently():
    processed_chunks = []
    
    # Her bir DCD dosyasını sırayla işle
    for dcd in args.trajectory:
        print(f'Processing {dcd}...')
        # Dosyayı 100 frame'lik parçalar halinde oku (RAM'i korur)
        for chunk in md.iterload(dcd, top=args.protein, chunk=100):
            if args.remove_waters:
                # Sadece protein ve ligandı tut, gerisini hemen sil
                selection = chunk.top.select('not (resname WAT or resname HOH or resname TIP3 or resname POPC or resname SOD or resname CLA)')
                chunk = chunk.atom_slice(selection)
            
            processed_chunks.append(chunk)

    print("Concatenating processed frames...")
    t = md.join(processed_chunks)
    
    print('Imaging and Realigning...')
    t.image_molecules(inplace=True)
    prot_backbone = t.top.select('protein and backbone')
    t.superpose(t[0], atom_indices=prot_backbone)

    print('Saving outputs...')
    t[0].save(args.output + '.pdb')
    t.save(args.output + '.dcd')

    # RMSD Analizi (Hafifletilmiş veri üzerinden)
    lig_atoms = t.top.select("resname UNL")
    rmsds_lig = md.rmsd(t, t, frame=0, atom_indices=lig_atoms) if len(lig_atoms) > 0 else [0]*t.n_frames
    rmsds_bck = md.rmsd(t, t, frame=0, atom_indices=prot_backbone)

    fig = go.Figure()
    fig.add_trace(go.Scatter(y=rmsds_lig, mode='lines', name='Ligand (UNL)'))
    fig.add_trace(go.Scatter(y=rmsds_bck, mode='lines', name='Protein Backbone'))
    fig.write_image(args.output + '.svg')
    print("Done!")

if __name__ == "__main__":
    process_memory_efficiently()
