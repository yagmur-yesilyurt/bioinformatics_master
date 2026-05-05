# Bioinformatics Master - Nav1.5 Drug Discovery Project

## Overview

This project contains computational biology research on **Nav1.5 (voltage-gated sodium channel)** focusing on drug binding, molecular docking, and molecular dynamics simulations. The research explores how different drugs (mexiletine, flecainide, lidocaine, etc.) interact with wild-type and mutant variants of the Nav1.5 protein.

---

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
- `multidrug/` — Docking for other drugs (flecainide, lidocaine, propafenone, quinidine)
- `scoring/` — Docking logs, configurations, and scoring reports

### 📁 `03_md_simulations/`
**Molecular dynamics simulation analysis and results**
- `setup/` — GROMACS configuration files and parameters
- `holoWT_postMD/` — WT protein + ligand MD analysis
- `apoWT_postMD/` — WT protein alone MD analysis
- `holoMutant_postMD/` — Mutant protein + ligand MD analysis
- `apoMutant_postMD/` — Mutant protein alone MD analysis

Each contains: RMSD plots, RMSF profiles, contact maps, distance plots, clustering analysis

### 📁 `04_results/`
**Publication-ready outputs**
- `figures/` — Interaction diagrams and visualization images (PNG)
- `tables/` — Summary data and results tables
- `reports/` — Analysis summaries and reports

### 📁 `05_literature/`
**Reference documents and protocols**
- `papers/Nav1.5_structure/` — Structural studies of Nav1.5
- `papers/Nav1.5_activators/` — Drug activator research
- `papers/general_docking_md/` — General computational methods
- `books/intelligent_systems/` — Reference textbooks
- `protocols/` — Experimental guides and methods
- `confidential/` — Restricted access documents

### 📁 `06_scripts/`
**Analysis and visualization Python scripts**
- `analysis/` — Data processing and clustering scripts
  - Clustering analysis (DBSCAN)
  - RMSD/RMSF calculations
  - Hydrogen bond and salt bridge analysis
  - Pore dimension analysis
- `preprocessing/` — Data preparation tools
- `docking/` — Docking setup and analysis
- `visualization/` — Figure generation scripts

### 📁 `07_archive/`
**Working notes and test directories**
- `working_notes/` — Backup of experimental runs and test dockings

---

## Quick Navigation

**Start here if you want to...**

- **View protein structures** → Go to `01_structures/receptors/`
- **Check docking results** → Go to `02_docking/mexiletine/mutants/N347K/` (or your mutant of interest)
- **Analyze MD simulations** → Go to `03_md_simulations/holoMutant_postMD/`
- **Read papers and guides** → Go to `05_literature/papers/` or `05_literature/protocols/`
- **Run analysis scripts** → Go to `06_scripts/analysis/`

---

## File Types Reference

| Extension | Content | Location |
|-----------|---------|----------|
| `.pdb` | Protein structures | `01_structures/receptors/` |
| `.mol2` | Ligand coordinates | `01_structures/ligands/` |
| `.pdbqt` | Prepared structures for docking | `01_structures/prepared/` |
| `.png` | Figures and plots | `04_results/figures/` |
| `.txt` | Logs and analysis results | `02_docking/scoring/` or `03_md_simulations/` |
| `.py` | Python analysis scripts | `06_scripts/analysis/` |
| `.pdf` | Papers and guides | `05_literature/papers/` or `05_literature/protocols/` |

---

## Project Status

✅ **Structures** — 12 protein PDB files (WT + mutants)  
✅ **Docking** — Results for 3 drug variants, 3+ mutants, blind docking included  
✅ **MD Analysis** — 4 simulation sets (holo/apo × WT/Mutant) with plots and statistics  
✅ **Literature** — 25+ reference papers and experimental protocols  
✅ **Scripts** — 14+ Python analysis tools  

---

**Questions?** Check the folder-specific README files or examine example files in each directory.
