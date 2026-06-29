# Katman 0 — Trajectory-only validation toolkit

These scripts close the **"[Validation required]"** items in the revised report
**without running any new MD**. They analyse the six 500 ns trajectories you
already have. Outputs go to `07_binding_site_analysis/katman0/`.

Configured for your repository layout:
```
03_md_simulations/{wt,mutant}/{apo,holo}/<NAME>_500NS/
    analysis_result.pdb        (topology)
    analysis_result.dcd        (trajectory)
    rep_structure_full_10k.pdb (DBSCAN medoid representative frame)
```

| Script | Report sections | Needs |
|---|---|---|
| `check_numbering.py` | (prerequisite for all) | MDAnalysis only |
| `b1_core_rmsd.py` | 3.1 | MDAnalysis only |
| `a1_distances_filter.py` | 3.3, 3.8 | MDAnalysis only |
| `a3_electrostatics.py` | 3.3, 3.8 | pdb2pqr30 + apbs + GridDataFormats |
| `b2_pore_distribution.py` | 3.7 | HOLE executable |
| `b3_druggability_distribution.py` | 3.4 | fpocket executable |

## 0. Setup

```bash
pip install MDAnalysis numpy matplotlib GridDataFormats
# external binaries (only for a3/b2/b3):
#   conda install -c conda-forge apbs pdb2pqr fpocket
#   HOLE2 suite for b2
```

Point the toolkit at your repo (no path editing needed if you set this env var):
```bash
export BIOINF_ROOT=/absolute/path/to/bioinformatics_master
```
(Or edit `REPO_ROOT` at the top of `config.py`.)

## 1. RUN check_numbering.py FIRST  (non-negotiable)

CHARMM-GUI/OpenMM often renumber residues (e.g. from 1) or drop unresolved
loops. If that happened, `resid 347` is **not** N347/K347 and every result is
silently wrong. `check_numbering.py` loads each topology, prints the resname at
each key resid, and — if they don't match — **auto-detects the constant offset**
and prints the exact `config.py` edit to make.

```bash
python check_numbering.py
```
- All `OK` → proceed.
- `MISMATCH` with a detected offset → apply the printed `RES_*` changes, re-run
  the checker until all `OK`.
- `MISMATCH` with no single offset → numbering is non-contiguous; inspect
  `analysis_result.pdb` manually (find the DEKA Lys and the two binding Asp's).

> Note: the two mutant-holo dirs also contain `rep_structure_full_10k_fixed.pdb`.
> If APBS/fpocket complain about atom names/protonation on those reps, switch
> their `REP_FRAMES` entry in `config.py` to the `*_fixed.pdb` version.

## 1. What each script decides

### a1 — through-space distances + DEKA-filter geometry (3.3, 3.8)
Measures K347(NZ)↔ASP1423/ASP1714 minimum distances per frame and the CA–CA
geometry relative to the DEKA filter residue K1419.

- **Gate:** in `holoMutant_prot`, if `frac<12Å` is high, the DI charge is
  spatially coupled to the DIII/DIV acidics → the indirect mechanism (3.3) and
  the selectivity-filter-vestibule identity (3.8) are **supported by
  measurement**, not assumed.
- If `frac<12Å` is low → rewrite 3.3 as *"indirect mechanism not supported by
  through-space proximity"* (a correct negative result, also worth 4).

### a3 — APBS electrostatic potential at Site 2 (3.3, 3.8)
Builds PQR (pH 7.4, PROPKA) and an APBS potential grid for WT vs mutant, then
samples the potential at the Site-2 centroid.

- **Gate:** mutant potential **more negative** than WT (deeper cation well)
  → K347⁺ measurably reshapes binding-site electrostatics → supports 3.3/3.8.
- WT ≈ mutant → mechanism not supported; soften 3.3.

### b1 — core vs global RMSD (3.1)
Defines the folded core by RMSF percentile (model-agnostic), then reports core
RMSD (core superposition) beside the original global RMSD.

- **Gate:** core RMSD low + flat while global high → the 7–13 Å global values
  are **peripheral motion, not core collapse**; rewrite the "stable plateau"
  claim as a core statement.
- If core RMSD is *also* high for `holoMutant_prot` → disclose genuine core
  instability of the central trajectory (still honest, still 4).

### b2 — pore bottleneck distribution (3.7)
Runs HOLE per frame (subsampled) → bottleneck radius distribution.

- **Gate:** for each system pair, if the mean difference is **within the pooled
  std**, the 1.40-vs-1.56 Å contrasts are noise and the categorical
  "above/below the 1.5 Å Na⁺ radius" framing must be softened.

### b3 — druggability distribution (3.4)
Runs fpocket on subsampled frames → PHE1760-pocket druggability distribution.

- **Gate:** wide spread (std > 0.15) → the single-frame 1.000 was a frame
  artefact; report the distribution. Consistently high → robust.

## 2. Suggested order (fastest score gain first)

0. `check_numbering.py` — confirm/repair residue numbering (prerequisite).
1. `b1_core_rmsd.py` — rescues/qualifies 3.1 (no binaries).
2. `a1_distances_filter.py` — tests 3.3 + 3.8 spatial coupling (no binaries).
3. `a3_electrostatics.py` — confirms the electrostatic mechanism (APBS).
4. `b2_pore_distribution.py` — settles 3.7 noise question (HOLE).
5. `b3_druggability_distribution.py` — settles 3.4 single-frame question (fpocket).

Steps 0–2 need no external binaries and should be done first — they move 3.1,
3.3, 3.7(geometry) and 3.8 toward 4 within a day.

## 3. Writing results back into the report

Each gate has two honest outcomes (supports / refutes), and **both** raise the
section score, because a correctly reported negative result is stronger than an
unvalidated assertion. After running, replace each `[Validation required]` flag
in the report with the measured number and the gate verdict.

> Important: Katman 0 raises 3.1, 3.3, 3.4, 3.7, 3.8 toward 4. The remaining
> ceiling on 3.2 and 3.6 (single-trajectory occupancy and MM-GBSA magnitude)
> needs **replication + FEP/TI** (Katman 1) — that is the separate next step.
