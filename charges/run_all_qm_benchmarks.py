
#!/usr/bin/env python3
"""
run_all_qm_benchmarks.py

Prepares all QM levels using the exact five reaction mappings and the
PDB organization from fast_lmrr_evb_qm_workflow(3).ipynb.

It creates inputs only; it does NOT launch ORCA.

Reference:
    B3LYP/def2-SVP optimization + frequencies

Benchmarks:
    PBE0/def2-SVP single point
    PBE0/def2-TZVP single point
    omegaB97M-V/def2-TZVP single point
    omegaB97X-V/def2-TZVP single point
    PBE0-DH/def2-TZVP single point

All benchmark single points use the optimized B3LYP/def2-SVP geometry
when available. If it is not available, the script searches the existing
xTB-prepared geometry as a fallback.
"""

from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent

SCRIPTS = [
    "01_prepare_b3lyp_reference.py",
    "02_prepare_pbe0_def2svp.py",
    "03_prepare_pbe0_def2tzvp.py",
    "04_prepare_wB97M-V_def2tzvp.py",
    "05_prepare_wB97X-V_def2tzvp.py",
    "06_prepare_pbe0_dh_def2tzvp.py",
]

def main():
    print("="*72)
    print("LmrR QM BENCHMARK PREPARATION")
    print("="*72)
    for script in SCRIPTS:
        path = HERE / script
        print(f"\n>>> {script}")
        subprocess.run([sys.executable, str(path)], check=True)
    print("\n" + "="*72)
    print("ALL QM INPUTS PREPARED")
    print("="*72)
    print("No ORCA calculations were launched.")
    print("Run B3LYP/def2-SVP first. Then run the four benchmark families.")
    print("PBE0-DH/def2-TZVP is the most expensive; test one job first.")

if __name__ == "__main__":
    main()
