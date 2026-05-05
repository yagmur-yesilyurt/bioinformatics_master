================================================================
  Nav1.5 N347K — MD Analysis Scripts
  Yağmur Yeşilyurt | SIS2E PhD Fellow | UniBicocca
  April 2026
================================================================

ENVIRONMENT
-----------
All scripts run in the md_run conda environment:
  conda activate md_run

For trajectory preparation (analyse_silvia.py):
  conda activate md_run  (same environment)

For MM-GBSA:
  conda activate mmpbsa  (see GUIDE 5)


================================================================
SCRIPTS OVERVIEW
================================================================

1. run_analysis.py  [MASTER RUNNER]
   ─────────────────────────────────
   Run this first for any new simulation. It calls the four
   sub-scripts below automatically after patching the Config.

   Usage:
     cd ~/PHD_MD/Scripts
     python run_analysis.py

   Before running, edit the CONFIG section at the top:
     SYS_NAME   → label for output files (e.g. "holoMutant")
     IS_HOLO    → True if ligand present, False for apo
     SIM_PATH   → full path to the openmm/ directory
     OUTPUT_DIR → where output files will be saved
     MUTATION_RESID → mutation site residue number (347 for N347K)

   Sub-scripts called automatically:
     - all-graphics.py
     - clustering-dbscan.py
     - analyze_pore.py
     - hbond_saltbridge.py

   Note: MM-GBSA is NOT included here. Run it separately
   following GUIDE 5.


----------------------------------------------------------------

2. all-graphics.py  [CALLED BY run_analysis.py]
   ──────────────────────────────────────────────
   Produces all RMSF, RMSD, ligand distance and contact plots.
   Works for both holo (ligand present) and apo systems.

   Outputs:
     01_protein_flexibility_rmsf.png   — RMSF per residue
     02_ligand_pore_distance.png       — mexiletine distance to pore
     03_residue_contact_map.png        — contact map (4.5 Å cutoff)
     04_ligand_conformational_rmsd.png — ligand internal RMSD
     05_multisite_distances.png        — distance to K347, PRO348, pore
     06_ligand_density_map.dx          — 3D density map (VMD-compatible)
     09_protein_rmsd_domains.png       — per-domain RMSD
     10_protein_rmsd_global.png        — global backbone RMSD
     11_protein_rmsd_ploop.png         — P-loop RMSD (resid 330-370)

   Note: ligand outputs (02-06) are skipped automatically for apo.


----------------------------------------------------------------

3. hbond_saltbridge.py  [CALLED BY run_analysis.py]
   ──────────────────────────────────────────────────
   Calculates H-bond contacts between mexiletine and protein,
   and K347-D356 salt bridge occupancy over time.

   Outputs:
     12_hbond_saltbridge_summary.txt  — numerical summary
     13_hbond_timeseries.png          — H-bond contacts over time
     14_saltbridge_distance.png       — K347-D356 distance over time
     15_saltbridge_occupancy.png      — salt bridge occupancy bar chart


----------------------------------------------------------------

4. clustering-dbscan.py  [CALLED BY run_analysis.py]
   ────────────────────────────────────────────────────
   PCA + DBSCAN conformational clustering on Cα atoms.
   Identifies distinct conformational states in the trajectory.

   Outputs:
     08_conformational_clustering_dbscan.png — PCA scatter + cluster labels

   Key parameters (in Config):
     EPS = 10          — DBSCAN neighborhood distance (PCA space)
     MIN_SAMPLES = 50  — minimum points to form a cluster


----------------------------------------------------------------

5. analyze_pore.py  [CALLED BY run_analysis.py]
   ──────────────────────────────────────────────
   Measures pore diameter at three levels along the channel:
     Level 1: Selectivity filter (DEKA motif)
     Level 2: Drug binding pocket (mid-pore)
     Level 3: Intracellular gate (activation gate)

   Outputs:
     07_pore_dimensions_multilevel.png — three-panel time series


----------------------------------------------------------------

6. hole.py  [RUN SEPARATELY]
   ───────────────────────────
   HOLE pore radius profile averaged across the full trajectory.
   Requires the HOLE executable to be installed.

   Usage:
     cd /path/to/openmm
     python ~/PHD_MD/Scripts/hole.py

   Outputs:
     06_HOLE_500ns_Average.png — mean radius ± SD with bottleneck

   Note: Run from the simulation directory (not Scripts/).
   The script reads step7_1.dcd ... step7_10.dcd from the
   current working directory.


----------------------------------------------------------------

7. analyse_silvia.py  [SILVIA MULTARI — RUN BEFORE CLUSTERING]
   ─────────────────────────────────────────────────────────────
   Merges all 10 DCD files, removes water/lipids/ions,
   re-images the system, aligns to frame 0, and saves:
     analysis_result.pdb  — topology reference (first frame)
     analysis_result.dcd  — clean 10,000-frame trajectory

   These files are required by clustering_silvia_revised.py.

   Usage:
     cd /path/to/openmm
     python ~/PHD_MD/Scripts/silvia/analyse_silvia.py \
       -p step5_input.psf \
       -t step7_1.dcd step7_2.dcd step7_3.dcd step7_4.dcd step7_5.dcd \
          step7_6.dcd step7_7.dcd step7_8.dcd step7_9.dcd step7_10.dcd \
       -o analysis_result \
       --remove-waters

   Note: image_molecules may fail with PSF files (known mdtraj issue).
   The script handles this automatically by skipping reimaging.


----------------------------------------------------------------

8. clustering_silvia_revised.py  [CORE CLUSTERING SCRIPT]
   run_clustering_silvia_revised.py  [P-LOOP BINDING — neutral mexiletine]
   run_clustering_protonated.py      [DOMAIN III/IV BINDING — protonated mexiletine]
   ──────────────────────────────────────────────────────────────────────────────────
   Ligand-filtered DBSCAN clustering to select a representative
   frame. Must be run AFTER analyse_silvia.py.

   Method:
     Step 1 — Keep only frames where mexiletine is within
              5.0 Å of specified residues (ligand filter)
     Step 2 — DBSCAN clustering on backbone RMSD of filtered frames
              (eps=0.53 Å, min_samples=5)
     Step 3 — Select medoid of the largest cluster

   run_clustering_silvia_revised.py
     → Use for P-loop binding systems (neutral mexiletine)
     → Ligand filter residues: PRO348, ASP349, THR353, PHE355
     → ploop_resids=[348, 349, 353, 355]

     Usage:
       cd ~/PHD_MD/Scripts/silvia
       python run_clustering_silvia_revised.py

   run_clustering_protonated.py
     → Use for Domain III/IV binding systems (protonated mexiletine)
     → Ligand filter residues: ASP1714, ASP1423, VAL1400, ILE1424, GLY1420
     → ploop_resids=[1714, 1423, 1424, 1400, 1420]

     Usage:
       cd ~/PHD_MD/Scripts/silvia
       python run_clustering_protonated.py

   Outputs (both scripts):
     rep_structure_full_10k.pdb — representative frame (all atoms)
     Trajectory index printed to terminal (needed for MM-GBSA)


================================================================
WORKFLOW SUMMARY
================================================================

For each new simulation, run in this order:

  1. run_analysis.py          (main pipeline)
  2. hole.py                  (from openmm/ directory)
  3. analyse_silvia.py        (from openmm/ directory)
  4. run_clustering_silvia_revised.py  (P-loop binding systems)
     run_clustering_protonated.py       (Domain III/IV binding systems)
  5. MM-GBSA                  (see GUIDE 5, separate environment)

All outputs are saved to OUTPUT_DIR defined in run_analysis.py.
================================================================
