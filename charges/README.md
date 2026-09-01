
# LmrR QM benchmark code package

These scripts were built from the uploaded
`fast_lmrr_evb_qm_workflow(3).ipynb`.

The notebook uses five reaction mappings and the following PDB files:

- RS1_1_to_TS1_2 -> step_1_1_RS/TS/PS.pdb
- TS1_2_to_PS1_2b -> step_1_2_RS/TS/PS.pdb
- RS1_2b_to_PS1_3 -> step_1_3_RS/TS/PS.pdb
- RS2_1_to_TS2_1a -> step_2_1_RS/TS/PS.pdb
- TS2_1a_to_PS2_2 -> step_2_2a_RS.pdb, step_2_2b_TS.pdb, step_2_2b_PS.pdb

## Important correction

The uploaded notebook actually generates:
    B3LYP D3BJ / def2-SVP
for the final ORCA/Gaussian QM optimization/frequency inputs.

Therefore B3LYP is retained as the reference calculation.

## QM hierarchy

1. B3LYP-D3BJ/def2-SVP
   Main reference geometry optimization + frequencies.

2. PBE0-D3BJ/def2-SVP
   Single-point benchmark on the B3LYP optimized geometry.

3. PBE0-D3BJ/def2-TZVP
   Higher-basis single-point benchmark.

4. omegaB97M-V/def2-TZVP
   Independent modern range-separated hybrid benchmark.
   VV10 is part of the functional; D3BJ is NOT added.

5. omegaB97X-V/def2-TZVP
   Independent modern range-separated hybrid benchmark.
   VV10 is part of the functional; D3BJ is NOT added.

6. PBE0-DH/def2-TZVP
   Expensive double-hybrid benchmark.
   Test one representative structure first.

## Run

Open PowerShell in this package directory:

    python run_all_qm_benchmarks.py

This only prepares ORCA input files.

It does NOT run ORCA.

## Recommended order

First:
    B3LYP/def2-SVP optimization + frequencies

Then:
    PBE0/def2-SVP SP
    PBE0/def2-TZVP SP
    omegaB97M-V/def2-TZVP SP
    omegaB97X-V/def2-TZVP SP

Finally:
    PBE0-DH/def2-TZVP SP

## Energy interpretation

QM:
    Delta E_QM^‡ = E_TS - E_RS

EVB:
    Delta G_EVB^‡ = G_TS - G_RS

These are not the same physical quantity.

Your experimental WT activation free energy:
    20.8 kcal mol^-1

is the EVB calibration target, not a value that should be inserted directly as a QM energy difference.

## Output organization

charges/
└── qm_reference_for_evb/
    ├── B3LYP_def2SVP/
    ├── PBE0_def2SVP/
    ├── PBE0_def2TZVP/
    ├── wB97M-V_def2TZVP/
    ├── wB97X-V_def2TZVP/
    └── PBE0-DH_def2TZVP/

Each method contains:
    step/
        RS/
            input/
            geometry/
            energy/
            charges/
            frequencies/
        TS/
        PS/

The existing fast workflow is not overwritten.
