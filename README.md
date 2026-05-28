# Bioinformatics Master - Nav1.5 Drug Discovery Project

## Folder Guide

### 📁 `01_structures/`
**Protein and ligand coordinate files**
- `receptors/wt/` — Wild-type Nav1.5 protein structures (PDB format)
- `receptors/mutants/` — Mutant variants (N347K, etc.)
- `ligands/mexiletine/` — Drug molecule structures in multiple formats
- `prepared/` — Protein structures prepared for docking (AutoDock PDBQT format)

### 📁 `02_docking/`
**Molecular docking results and analysis**
- `mexiletine/` — Docking results for mexiletine drug
  - `wt/` — Wild-type binding poses
  - `mutants/` — Mutant protein binding (A344S, D349N, N347K)
  - `blind/` — Blind docking (unbiased search)
- `other_drugs/` — Docking for multidrug compounds (flecainide, lidocaine, propafenone, quinidine)
  - `A344S/` — Docking results for A344S mutant
  - `D349N/` — Docking results for D349N mutant
  - `N347K/` — Docking results for N347K mutant
- `scoring/` — Docking logs, configurations, and scoring reports

### 📁 `03_md_simulations/`
**Molecular dynamics simulation setup and scripts**
- `md.py` — MD simulation analysis and processing scripts

### 📁 `04_postMD_analysis_results/`
**Molecular dynamics simulation analysis and results**
- `holoWT_postMD/` — WT protein + ligand MD analysis
- `apoWT_postMD/` — WT protein alone MD analysis
- `holoMutant_postMD/` — Mutant protein + ligand MD analysis
- `apoMutant_postMD/` — Mutant protein alone MD analysis

Each contains: RMSD plots, RMSF profiles, contact maps, distance plots, clustering analysis

### 📁 `05_literature/`
**Reference documents and protocols**
- `papers/Nav1.5_structure/` — Structural studies of Nav1.5
- `papers/Nav1.5_activators/` — Drug activator research
- `papers/general_docking_md/` — General computational methods
- `books/intelligent_systems/` — Reference textbooks
- `my_protocols/` — Experimental guides and methods (Mutagenesis, Membrane embedding, MD execution, Post-MD analysis, MMGBSA, Visualization)
- `confidential/` — Restricted access documents

### 📁 `06_scripts/`
**Analysis and visualization Python scripts**
- `analysis/` — Data processing and clustering scripts
  - Clustering analysis (DBSCAN)
  - RMSD/RMSF calculations
  - Hydrogen bond and salt bridge analysis
  - Pore dimension analysis
  - Miscellaneous analysis tools
- `reorganize_docking.py` — Docking directory organization automation tool

### 📁 `07_archive/`
**Working notes and test directories**
- `working_notes/` — Backup of experimental runs and test dockings

---

## Quick Navigation

**Start here if you want to...**

- **View protein structures** → Go to `01_structures/receptors/`
- **Check docking results** → Go to `02_docking/mexiletine/mutants/N347K/` or `02_docking/other_drugs/N347K/`
- **Analyze MD simulations** → Go to `04_postMD_analysis_results/holoMutant_postMD/`
- **Read protocols and guides** → Go to `05_literature/my_protocols/`
- **Read papers and references** → Go to `05_literature/papers/`
- **Run analysis scripts** → Go to `06_scripts/analysis/`

---

## File Types Reference

| Extension | Content | Location |
|-----------|---------|----------|
| `.pdb` | Protein structures | `01_structures/receptors/` |
| `.mol2` | Ligand coordinates | `01_structures/ligands/` |
| `.pdbqt` | Prepared structures for docking | `01_structures/prepared/` |
| `.txt` | Logs and analysis results | `02_docking/scoring/` or `04_postMD_analysis_results/` |
| `.py` | Python analysis scripts | `06_scripts/analysis/` |
| `.pdf` | Papers and guides | `05_literature/papers/` or `05_literature/my_protocols/` |

---

## Project Status

✅ **Structures** — 12 protein PDB files (WT + mutants)  
✅ **Docking** — 16 docking setups (Mexiletine × 4 receptors + Multidrug × 3 receptors × 4 ligands)  
✅ **MD Analysis** — 4 simulation sets (holo/apo × WT/Mutant) with plots and statistics  
✅ **Literature** — 25+ reference papers, 6+ experimental protocols, and guides  
✅ **Scripts** — 13+ Python analysis and automation tools  
✅ **Organization** — Restructured docking directories for improved project management  

---

**Questions?** Check the folder-specific README files or examine example files in each directory.
