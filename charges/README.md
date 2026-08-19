# LmrR QM → EVB reference pipeline

This is the modified version of the uploaded `fast_lmrr_evb_qm_workflow.ipynb` logic, reorganized specifically to create an EVB-ready QM reference dataset.

## Output structure

```text
qm_reference_for_evb/
├── step_1_1/
│   ├── RS/
│   │   ├── geometry/
│   │   ├── energy/
│   │   ├── charges/
│   │   └── frequencies/
│   ├── TS/
│   └── PS/
├── step_1_2/
├── step_1_3/
├── step_2_1/
├── step_2_2/
├── summaries/
│   ├── qm_state_summary.csv
│   ├── qm_barriers_for_evb_comparison.csv
│   ├── qm_mulliken_charges_long.csv
│   ├── qm_frequencies_long.csv
│   └── qm_mapping_validation.csv
└── evb_reference/
    └── qm_vs_evb_comparison_template.csv
```

## Recommended use

1. Put the script in `D:\PhD_Thesis\LmrR_EVB\charges`.
2. Run:

```powershell
python qm_reference_pipeline.py --prepare
```

3. Run the generated DFT inputs for RS/TS/PS.
4. Put/copy the completed output files and optimized XYZ geometries into the corresponding state folders.
5. Run:

```powershell
python qm_reference_pipeline.py --collect
python qm_reference_pipeline.py --validate
```

Or run everything after the QM outputs exist:

```powershell
python qm_reference_pipeline.py --all
```

The key table for later EVB comparison is:

`qm_reference_for_evb/summaries/qm_barriers_for_evb_comparison.csv`

Do not equate QM ΔE‡ with EVB ΔG‡. Use the QM values as the electronic/reference layer and the EVB values as the protein/explicit-solvent free-energy layer.
