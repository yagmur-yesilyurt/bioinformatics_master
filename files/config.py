"""
config.py  --  Katman 0 toolkit, configured for Yagmur's bioinformatics_master.

Set REPO_ROOT to the absolute path of your repo, then run check_numbering.py
FIRST, then the a*/b* scripts. Outputs go to 07_binding_site_analysis/katman0/.
"""

import os

# ----------------------------------------------------------------------
# 0. Repo root  --  EDIT THIS to the absolute path on the machine you run on.
# ----------------------------------------------------------------------
REPO_ROOT = os.environ.get("BIOINF_ROOT", os.path.abspath("."))

def R(*parts):
    return os.path.join(REPO_ROOT, *parts)

MD = "03_md_simulations"

# Each system dir holds: analysis_result.pdb (topology), analysis_result.dcd
# (trajectory), rep_structure_full_10k.pdb (DBSCAN medoid representative frame).
_DIR = {
    "apoWT":              R(MD, "wt", "apo",  "Nav15_500NS"),
    "apoMutant":          R(MD, "mutant", "apo",  "N347K_500NS"),
    "holoWT_neutral":     R(MD, "wt", "holo", "Nav15_MEX_500NS"),
    "holoMutant_neutral": R(MD, "mutant", "holo", "N347K_MEX_500NS"),
    "holoWT_prot":        R(MD, "wt", "holo", "Nav15_MEX_PROTONATED_500NS"),
    "holoMutant_prot":    R(MD, "mutant", "holo", "N347K_MEX_PROTONATED_500NS"),
}
_IS_MUT = {"apoWT": False, "apoMutant": True, "holoWT_neutral": False,
           "holoMutant_neutral": True, "holoWT_prot": False, "holoMutant_prot": True}

TOPO_NAME = "analysis_result.pdb"
TRAJ_NAME = "analysis_result.dcd"
REP_NAME  = "rep_structure_full_10k.pdb"
# NOTE: the two mutant-holo dirs ALSO contain rep_structure_full_10k_fixed.pdb.
# If APBS/fpocket complain about atom names/protonation on those, switch their
# REP_FRAMES entry below to the *_fixed.pdb version.

SYSTEMS = {k: (os.path.join(_DIR[k], TOPO_NAME),
               os.path.join(_DIR[k], TRAJ_NAME),
               _IS_MUT[k]) for k in _DIR}

REP_FRAMES = {k: os.path.join(_DIR[k], REP_NAME) for k in _DIR}

# ----------------------------------------------------------------------
# 2. Key residues (human Nav1.5 numbering).
#    check_numbering.py verifies these map to the right resnames in YOUR pdb.
# ----------------------------------------------------------------------
RES_347   = 347
RES_D356  = 356
RES_K1419 = 1419
RES_D1423 = 1423
RES_D1714 = 1714
SITE2_RESIDS = [RES_D1423, RES_D1714, 1400]   # ASP1423, ASP1714, VAL1400

# Expected resnames for the sanity check (mutant has LYS at 347, WT has ASN)
EXPECTED = {347: ("LYS", "ASN"), 356: ("ASP",), 1419: ("LYS",),
            1423: ("ASP",), 1714: ("ASP",)}

# ----------------------------------------------------------------------
# 3. Output + tuning
# ----------------------------------------------------------------------
OUT_ROOT = R("07_binding_site_analysis", "katman0")

GATE_NEAR  = 12.0
GATE_CLOSE = 4.0
CORE_RMSF_PERCENTILE = 60.0
CORE_RMSF_ABS_CUTOFF = 3.0
FPOCKET_N_FRAMES = 60

BIN_PDB2PQR = "pdb2pqr30"
BIN_APBS    = "apbs"
BIN_FPOCKET = "fpocket"
BIN_HOLE    = "hole"

def ensure_out(subdir=""):
    d = os.path.join(OUT_ROOT, subdir) if subdir else OUT_ROOT
    os.makedirs(d, exist_ok=True)
    return d
