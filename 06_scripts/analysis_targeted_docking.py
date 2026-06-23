import numpy as np

def get_mode1_coords(pdbqt_file):
    coords = []
    in_model1 = False
    with open(pdbqt_file) as f:
        for line in f:
            if line.startswith("MODEL 1"):
                in_model1 = True
            if line.startswith("ENDMDL") and in_model1:
                break
            if in_model1 and line.startswith(("ATOM", "HETATM")):
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords)

def get_protein_residues(pdb_file):
    residues = {}
    ext = pdb_file.split('.')[-1]
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    resid = int(line[22:26])
                    resn = line[17:20].strip()
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    if resid not in residues:
                        residues[resid] = {"resn": resn, "coords": []}
                    residues[resid]["coords"].append([x, y, z])
                except:
                    pass
    return residues

BASE_D = "/home/yagmur/Desktop/bioinformatics_master/02_docking/mexiletine/targeted"
BASE_B = "/home/yagmur/Desktop/bioinformatics_master/02_docking/mexiletine/blind"

systems = {
    "apoMutant_pocket2_PHE1760": {
        "docked":   f"{BASE_D}/apoMutant_MD_rep_frame_fpocket/pocket2_DomainIV_S6/docked.pdbqt",
        "protein":  f"{BASE_B}/Mutant/N347K_apoMutant_MD_rep_frame/rep_structure_full_10k_model0.pdbqt",
        "score": -6.3, "drugg": 0.970, "site": "PHE1760 / Domain IV S6",
    },
    "apoMutant_pocket13_DIII_IV": {
        "docked":   f"{BASE_D}/apoMutant_MD_rep_frame_fpocket/pocket13_DomainIII_IV_junction/docked.pdbqt",
        "protein":  f"{BASE_B}/Mutant/N347K_apoMutant_MD_rep_frame/rep_structure_full_10k_model0.pdbqt",
        "score": -5.4, "drugg": 0.873, "site": "Domain III/IV junction",
    },
    "apoWT_pocket83_PHE1760": {
        "docked":   f"{BASE_D}/apoWT_MD_rep_frame_fpocket/pocket83/docked.pdbqt",
        "protein":  f"{BASE_B}/WT/8VYK_apoWT_MD_rep_frame/rep_structure_full_10k_model0.pdbqt",
        "score": -999, "drugg": 0.814, "site": "PHE1760 / Domain IV S6",
    },
    "holoMutant_notr_pocket21_PHE1760": {
        "docked":   f"{BASE_D}/holoMutant_notr_MD_rep_frame_fpocket/pocket21/docked.pdbqt",
        "protein":  f"{BASE_B}/Mutant/N347K_holoMutant_notrMex_MD_rep_frame/rep_holoMutant_nolig_model0.pdbqt",
        "score": -6.0, "drugg": 1.000, "site": "PHE1760 / Domain IV S6 (induced fit)",
    },
    "holoMutant_prot_pocket7_PHE1760": {
        "docked":   f"{BASE_D}/holoMutant_protonated_MD_rep_frame_fpocket/pocket7/docked.pdbqt",
        "protein":  f"{BASE_B}/Mutant/N347K_holoMutant_protonatedMex_MD_rep_frame/holoMutant_nolig.pdb",
        "score": -5.6, "drugg": 0.991, "site": "PHE1760 / Domain IV S6",
    },
    "holoMutant_prot_pocket41_DIII_IV": {
        "docked":   f"{BASE_D}/holoMutant_protonated_MD_rep_frame_fpocket/pocket41/docked.pdbqt",
        "protein":  f"{BASE_B}/Mutant/N347K_holoMutant_protonatedMex_MD_rep_frame/holoMutant_nolig.pdb",
        "score": -999, "drugg": 0.010, "site": "Domain III/IV junction",
    },
    "holoWT_prot_pocket1_PHE1760": {
        "docked":   f"{BASE_D}/holoWT_protonated_MD_rep_frame_fpocket/pocket1/docked.pdbqt",
        "protein":  f"{BASE_B}/WT/N347K_holoWT_protonatedMex_MD_rep_frame/holoWT_nolig.pdb",
        "score": -6.0, "drugg": 0.998, "site": "PHE1760 / Domain IV S6",
    },
    "holoWT_prot_pocket91_DIII_IV": {
        "docked":   f"{BASE_D}/holoWT_protonated_MD_rep_frame_fpocket/pocket91/docked.pdbqt",
        "protein":  f"{BASE_B}/WT/N347K_holoWT_protonatedMex_MD_rep_frame/holoWT_nolig.pdb",
        "score": -5.5, "drugg": 0.012, "site": "Domain III/IV junction",
    },
}

output_lines = []

for name, s in systems.items():
    header = f"\n=== {name} | Score: {s['score']} | Drugg: {s['drugg']} | Site: {s['site']} ==="
    print(header)
    output_lines.append(header)
    try:
        lig_coords = get_mode1_coords(s["docked"])
        protein = get_protein_residues(s["protein"])
        contacts = []
        for resid, data in protein.items():
            prot_coords = np.array(data["coords"])
            dists = np.linalg.norm(lig_coords[:, None] - prot_coords[None, :], axis=2)
            if dists.min() < 4.5:
                contacts.append((dists.min(), resid, data["resn"]))
        contacts.sort()
        for d, rid, rn in contacts[:10]:
            line = f"  {rn}{rid}: {d:.2f} Å"
            print(line)
            output_lines.append(line)
    except Exception as e:
        err = f"  ERROR: {e}"
        print(err)
        output_lines.append(err)

outfile = "/home/yagmur/Desktop/bioinformatics_master/04_postMD_analysis_results/targeted_docking_analysis_all7systems.txt"
with open(outfile, "w") as f:
    f.write("\n".join(output_lines))
print(f"\nKaydedildi: {outfile}")