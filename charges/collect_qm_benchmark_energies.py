
#!/usr/bin/env python3
"""
collect_qm_benchmark_energies.py

Collects final ORCA single-point/optimization energies from the benchmark
directories and writes a CSV suitable for comparing QM methods.

It looks for:
  FINAL SINGLE POINT ENERGY
  FINAL SINGLE POINT ENERGY (in ORCA output)
and falls back to the last occurrence.

The script does not claim these are activation free energies.
"""

from pathlib import Path
import re
import csv

BASE_DIR = Path(r"D:\PhD_Thesis\LmrR_EVB")
ROOT = BASE_DIR / "charges" / "qm_reference_for_evb"

METHODS = [
    "B3LYP_def2SVP",
    "PBE0_def2SVP",
    "PBE0_def2TZVP",
    "wB97M-V_def2TZVP",
    "wB97X-V_def2TZVP",
    "PBE0-DH_def2TZVP",
]

STEPS = [
    "RS1_1_to_TS1_2",
    "TS1_2_to_PS1_2b",
    "RS1_2b_to_PS1_3",
    "RS2_1_to_TS2_1a",
    "TS2_1a_to_PS2_2",
]

STATES = ["RS","TS","PS"]

def energy_from_out(path):
    text = path.read_text(encoding="utf-8", errors="replace")
    matches = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)", text)
    if matches:
        return float(matches[-1])
    return None

def find_output(state_dir):
    outs = sorted(state_dir.glob("*.out"))
    outs += sorted(state_dir.glob("*.log"))
    if not outs:
        return None
    # Prefer the largest/most recent named output if multiple exist.
    for p in reversed(outs):
        if energy_from_out(p) is not None:
            return p
    return None

def main():
    rows=[]
    for method in METHODS:
        for step in STEPS:
            energies={}
            files={}
            for state in STATES:
                sd = ROOT/method/step/state
                out = find_output(sd)
                e = energy_from_out(out) if out else None
                energies[state]=e
                files[state]=str(out) if out else ""
                rows.append({
                    "method": method,
                    "step": step,
                    "state": state,
                    "energy_hartree": "" if e is None else f"{e:.12f}",
                    "output_file": files[state],
                })
            rs,ts,ps = energies["RS"],energies["TS"],energies["PS"]
            if rs is not None and ts is not None:
                rows.append({
                    "method": method,
                    "step": step,
                    "state": "BARRIER_TS_minus_RS",
                    "energy_hartree": f"{ts-rs:.12f}",
                    "output_file": "",
                })
            if rs is not None and ps is not None:
                rows.append({
                    "method": method,
                    "step": step,
                    "state": "REACTION_PS_minus_RS",
                    "energy_hartree": f"{ps-rs:.12f}",
                    "output_file": "",
                })

    out = ROOT/"qm_benchmark_energy_summary.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as f:
        w=csv.DictWriter(f, fieldnames=["method","step","state","energy_hartree","output_file"])
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote: {out}")

if __name__ == "__main__":
    main()
