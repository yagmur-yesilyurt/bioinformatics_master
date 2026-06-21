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
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                resid = int(line[22:26])
                resn = line[17:20].strip()
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                if resid not in residues:
                    residues[resid] = {"resn": resn, "coords": []}
                residues[resid]["coords"].append([x, y, z])
    return residues

BASE = "/home/yagmur/Desktop/bioinformatics_master/02_docking/mexiletine/blind"

systems = {
    "N347K_static_PDB": {
        "docked": f"{BASE}/Mutant/N347K_static_PDB/docked.pdbqt",
        "protein": f"{BASE}/Mutant/N347K_static_PDB/receptor_n347k.pdbqt",
    },
    "apoMutant_MD_rep_frame": {
        "docked": f"{BASE}/Mutant/N347K_apoMutant_MD_rep_frame/docked.pdbqt",
        "protein": f"{BASE}/Mutant/N347K_apoMutant_MD_rep_frame/rep_structure_full_10k_model0.pdbqt",
    },
    "holoMutant_notr_MD_rep_frame": {
        "docked": f"{BASE}/Mutant/N347K_holoMutant_notrMex_MD_rep_frame/docked.pdbqt",
        "protein": f"{BASE}/Mutant/N347K_holoMutant_notrMex_MD_rep_frame/rep_holoMutant_nolig_model0.pdbqt",
    },
    "holoMutant_protonated_MD_rep_frame": {
        "docked": f"{BASE}/Mutant/N347K_holoMutant_protonatedMex_MD_rep_frame/docked.pdbqt",
        "protein": f"{BASE}/Mutant/N347K_holoMutant_protonatedMex_MD_rep_frame/holoMutant_nolig.pdb",
    },
    "WT_static_PDB": {
        "docked": f"{BASE}/WT/WT_8VYK_static_PDB/docked.pdbqt",
        "protein": f"{BASE}/WT/WT_8VYK_static_PDB/receptor_wt.pdbqt",
    },
    "apoWT_MD_rep_frame": {
        "docked": f"{BASE}/WT/8VYK_apoWT_MD_rep_frame/docked.pdbqt",
        "protein": f"{BASE}/WT/8VYK_apoWT_MD_rep_frame/rep_structure_full_10k_model0.pdbqt",
    },
    "holoWT_notr_MD_rep_frame": {
        "docked": f"{BASE}/WT/N347K_holoWT_notrMex_MD_rep_frame/docked.pdbqt",
        "protein": f"{BASE}/WT/N347K_holoWT_notrMex_MD_rep_frame/holoWT_nolig.pdb",
    },
    "holoWT_protonated_MD_rep_frame": {
        "docked": f"{BASE}/WT/N347K_holoWT_protonatedMex_MD_rep_frame/docked.pdbqt",
        "protein": f"{BASE}/WT/N347K_holoWT_protonatedMex_MD_rep_frame/holoWT_nolig.pdb",
    },
}

# Sonuçları hem ekrana yaz hem dosyaya kaydet
output_lines = []

for name, s in systems.items():
    header = f"\n=== {name} ==="
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

# Kaydet
outfile = "/home/yagmur/Desktop/bioinformatics_master/02_docking/mexiletine/blind/blind_docking_analysis_all8systems.txt"
with open(outfile, "w") as f:
    f.write("\n".join(output_lines))

print(f"\nKaydedildi: {outfile}")
