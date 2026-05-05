"""
╔══════════════════════════════════════════════════════════════╗
║         Nav1.5 MD ANALYSIS PIPELINE — MASTER RUNNER          ║
║                                                              ║
║  Usage:                                                      ║
║      python run_analysis.py                                  ║
║                                                              ║
║  Edit only the CONFIG section below before running.          ║
╚══════════════════════════════════════════════════════════════╝
"""

import os
import sys
import shutil
import subprocess
import tempfile
import time

# ==============================================================================
#   C O N F I G   —   EDIT THIS SECTION FOR EACH SIMULATION
# ==============================================================================

class Config:

    # --- SYSTEM IDENTITY ------------------------------------------------------
    SYS_NAME   = "holoMutant_protonated"       # "holoWT" | "holoMutant" | "apoWT" | "apoMutant"
    IS_HOLO    = True            # True = ligand present | False = apo system

    # --- SIMULATION FILES -----------------------------------------------------
    SIM_PATH   = "/home/yyesilyurt/PHD_MD/N347K_MEX_PROTONATED_500NS/charmm-gui-7521967997/openmm"
    PSF_FILE   = "step5_input.psf"
    N_DCD      = 10              # step7_1.dcd ... step7_N.dcd
    DT         = 0.05            # ns per frame
    STRIDE     = 1               # stride for heavy analyses (contact map, RMSD)
    STRIDE_HBOND = 5             # stride for H-bond / salt bridge

    # --- OUTPUT ---------------------------------------------------------------
    OUTPUT_DIR = "/home/yyesilyurt/PHD_MD/N347K_MEX_PROTONATED_500NS/postMD_analysis"

    # --- BIOLOGY --------------------------------------------------------------
    MUTATION_RESID  = 347
    LIGAND_NAME     = "UNL"
    PORE_RESIDS     = [372, 901, 1422, 1714]
    BINDING_RESIDS  = list(range(340, 361))
    ACIDIC_RESIDS   = list(range(330, 400))

# --- MMGBSA ---------------------------------------------------------------
    RUN_MMGBSA      = False
    MMGBSA_DIR = "/home/yyesilyurt/PHD_MD/N347K_MEX_PROTONATED_500NS/postMD_analysis"
    MMGBSA_START    = 100
    MMGBSA_END      = 1000
    MMGBSA_INTERVAL = 5
    MMGBSA_TEMP     = 310
    RECEPTOR_GROUP  = 1
    LIGAND_GROUP    = 13

    # --- WHICH ANALYSES TO RUN ------------------------------------------------
    RUN_GRAPHICS    = True   # RMSF + ligand analyses + protein RMSD (all-graphics.py)
    RUN_CLUSTERING  = True   # DBSCAN clustering
    RUN_PORE_DIM    = True   # Pore dimensions multilevel
    RUN_HBOND_SB    = True   # H-bond + salt bridge

# ==============================================================================
# END OF CONFIG
# ==============================================================================


# Config values that need to be injected into each sub-script
CONFIG_PATCH = """
import os as _os
# ── Injected by run_analysis.py ──────────────────────────────
class Config:
    DCD_PATTERN     = 'step7_*.dcd'
    PSF_FILE        = '{PSF_FILE}'
    N_FILES         = {N_DCD}
    TOPOLOGY        = '{PSF_FILE}'
    DT              = {DT}
    STRIDE          = {STRIDE}
    MUTATION_RESID  = {MUTATION_RESID}
    LIGAND_NAME     = '{LIGAND_NAME}'
    PORE_RESIDS     = {PORE_RESIDS}
    BINDING_RESIDS  = {BINDING_RESIDS}
    ACIDIC_RESIDS   = {ACIDIC_RESIDS}
    BINDING_RESIDUES = {BINDING_RESIDS}  # alias for hbond script
    ACIDIC_RESIDUES  = {ACIDIC_RESIDS}   # alias for saltbridge script
    PLOOP_RANGE     = (330, 370)
    DOMAINS = {{
        'Domain I':   (1,    400),
        'Domain II':  (401,  800),
        'Domain III': (801,  1200),
        'Domain IV':  (1201, 1800)
    }}
    DOMAIN_COLORS = {{
        'Domain I':   '#e41a1c',
        'Domain II':  '#377eb8',
        'Domain III': '#4daf4a',
        'Domain IV':  '#ff7f00'
    }}
    # hbond_saltbridge specific
    HBOND_CUTOFF    = 3.5
    HBOND_ANGLE     = 150.0
    SALTBRIDGE_CUTOFF = 4.0
    EPS             = 10
    MIN_SAMPLES     = 50
    N_COMPONENTS    = 2
    ALIGN           = True
    SELECTION       = 'name CA'
    TARGET_FRAME    = 1
    DT_PER_FRAME    = {DT}
    OUTPUT_FILENAME = '07_pore_dimensions_multilevel.png'
# ────────────────────────────────────────────────────────────
"""


def print_banner():
    print("\n" + "="*65)
    print(f"  Nav1.5 MD Analysis Pipeline")
    print(f"  System : {Config.SYS_NAME}")
    print(f"  Type   : {'HOLO (ligand present)' if Config.IS_HOLO else 'APO (no ligand)'}")
    print(f"  Input  : {Config.SIM_PATH}")
    print(f"  Output : {Config.OUTPUT_DIR}")
    print("="*65 + "\n")


def setup_output_dir():
    os.makedirs(Config.OUTPUT_DIR, exist_ok=True)
    print(f"✓ Output directory ready: {Config.OUTPUT_DIR}")


def get_script_dir():
    return os.path.dirname(os.path.abspath(__file__))


def build_patch():
    """Fill CONFIG_PATCH template with current Config values."""
    return CONFIG_PATCH.format(
        PSF_FILE        = Config.PSF_FILE,
        N_DCD           = Config.N_DCD,
        DT              = Config.DT,
        STRIDE          = Config.STRIDE,
        MUTATION_RESID  = Config.MUTATION_RESID,
        LIGAND_NAME     = Config.LIGAND_NAME,
        PORE_RESIDS     = Config.PORE_RESIDS,
        BINDING_RESIDS  = Config.BINDING_RESIDS,
        ACIDIC_RESIDS   = Config.ACIDIC_RESIDS,
    )


def run_script(script_name):
    """
    Patch Config into script, run it inside SIM_PATH,
    then move outputs to OUTPUT_DIR with SYS_NAME prefix.
    """
    script_path = os.path.join(get_script_dir(), script_name)
    if not os.path.exists(script_path):
        print(f"  ✗ Script not found: {script_path}")
        return False

    print(f"\n>>> Running: {script_name}")
    t0 = time.time()

    # Read original script, remove its Config class, prepend patched Config
    import re
    with open(script_path, 'r') as f:
        original = f.read()

    # Remove original Config class block only
    # Stops at the first blank line followed by a non-indented line
    original = re.sub(
        r'(?m)^class Config:[^\n]*\n([ \t]+[^\n]*\n|\n)*',
        '',
        original
    )

    patched = build_patch() + "\n" + original

    # Write to temp file
    with tempfile.NamedTemporaryFile(
        mode='w', suffix='.py', delete=False, dir=Config.SIM_PATH
    ) as tmp:
        tmp.write(patched)
        tmp_path = tmp.name

    try:
        result = subprocess.run(
            [sys.executable, tmp_path],
            cwd=Config.SIM_PATH
        )
    finally:
        os.unlink(tmp_path)

    elapsed = time.time() - t0

    # Always collect outputs even if script failed partially
    moved = collect_outputs()

    if result.returncode != 0:
        print(f"  ✗ {script_name} failed (exit code {result.returncode})")
        if moved > 0:
            print(f"  ↑ But {moved} partial output(s) were saved")
        return False

    print(f"  ✓ Done in {elapsed/60:.1f} min | {moved} files → {Config.OUTPUT_DIR}")
    return True


def collect_outputs():
    """Move .png / .txt / .dx / .dat files from SIM_PATH to OUTPUT_DIR.
    Naming: {stem}_{SYS_NAME}.{ext}  e.g. 01_protein_flexibility_rmsf_holoWT.png
    """
    extensions = (".png", ".txt", ".dx", ".dat")
    moved = 0
    os.makedirs(Config.OUTPUT_DIR, exist_ok=True)
    for fname in os.listdir(Config.SIM_PATH):
        if fname.endswith(extensions):
            src = os.path.join(Config.SIM_PATH, fname)
            # Insert SYS_NAME before extension: stem_SYS_NAME.ext
            base, ext = os.path.splitext(fname)
            if base.endswith("_" + Config.SYS_NAME):
                dest_name = fname  # already has suffix
            else:
                dest_name = f"{base}_{Config.SYS_NAME}{ext}"
            dest = os.path.join(Config.OUTPUT_DIR, dest_name)
            try:
                shutil.move(src, dest)
                moved += 1
                print(f"    → {fname}  →  {dest_name}")
            except Exception as e:
                print(f"  ⚠ Could not move {fname}: {e}")
    return moved


def run_mmgbsa():
    """Write mmpbsa.in and run gmx_MMPBSA from MMGBSA_DIR."""
    if not Config.IS_HOLO:
        print("\n[MMGBSA] Skipped — apo system has no ligand")
        return False

    print(f"\n>>> Running MMGBSA in {Config.MMGBSA_DIR}")

    mmpbsa_in = os.path.join(Config.MMGBSA_DIR, "mmpbsa.in")
    with open(mmpbsa_in, "w") as f:
        f.write(
            f'&general\n'
            f'sys_name="{Config.SYS_NAME}",\n'
            f'startframe={Config.MMGBSA_START},\n'
            f'endframe={Config.MMGBSA_END},\n'
            f'interval={Config.MMGBSA_INTERVAL},\n'
            f'PBRadii=2,\n'
            f'temperature={Config.MMGBSA_TEMP},\n'
            f'/\n'
            f'&gb\n'
            f'igb=5,\n'
            f'saltcon=0.150,\n'
            f'/\n'
        )

    out_dat = f"FINAL_RESULTS_{Config.SYS_NAME}.dat"
    out_csv = f"FINAL_RESULTS_{Config.SYS_NAME}.csv"

    cmd = [
        "gmx_MMPBSA", "-O",
        "-i", "mmpbsa.in",
        "-cs", "gromacs_amber.pdb",
        "-ct", "traj.xtc",
        "-ci", "index.ndx",
        "-cg", str(Config.RECEPTOR_GROUP), str(Config.LIGAND_GROUP),
        "-cp", "gromacs.top",
        "-o", out_dat,
        "-eo", out_csv
    ]

    t0 = time.time()
    result = subprocess.run(cmd, cwd=Config.MMGBSA_DIR)
    elapsed = time.time() - t0

    if result.returncode != 0:
        print(f"  ✗ MMGBSA failed")
        return False

    for fname in [out_dat, out_csv]:
        src = os.path.join(Config.MMGBSA_DIR, fname)
        if os.path.exists(src):
            dest = os.path.join(Config.OUTPUT_DIR, fname)
            shutil.copy2(src, dest)
            print(f"  ✓ Copied: {fname}")

    print(f"  ✓ MMGBSA done in {elapsed/60:.1f} min")
    return True


def print_summary(results):
    print("\n" + "="*65)
    print(f"  PIPELINE COMPLETE — {Config.SYS_NAME}")
    print("="*65)
    passed = sum(1 for v in results.values() if v)
    total  = len(results)
    print(f"  {passed}/{total} analyses completed\n")
    for name, ok in results.items():
        print(f"  {'✓' if ok else '✗'}  {name}")
    print(f"\n  Results: {Config.OUTPUT_DIR}")
    print("="*65 + "\n")


# ==============================================================================
# MAIN
# ==============================================================================
def main():
    print_banner()
    setup_output_dir()
 
    results = {}
 
    if Config.RUN_GRAPHICS:
        results["RMSF + Ligand + Protein RMSD"] = run_script("all-graphics.py")
 
    if Config.RUN_CLUSTERING:
        results["DBSCAN Clustering"] = run_script("clustering-dbscan.py")
 
    if Config.RUN_PORE_DIM:
        results["Pore Dimensions"] = run_script("analyze_pore.py")
 
    if Config.RUN_HBOND_SB:
        results["H-bond & Salt Bridge"] = run_script("hbond_saltbridge.py")
 
    if Config.RUN_MMGBSA:
        results["MMGBSA"] = run_mmgbsa()
 
    print_summary(results)
 
 
if __name__ == "__main__":
    main()
