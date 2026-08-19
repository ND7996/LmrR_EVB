# step_2_2 / RS

## Starting structure

`step_2_2a_RS.pdb`

## QM input

Main DFT:

`input/step_2_2_RS_DFT.inp`

Higher-level benchmark:

`input/step_2_2_RS_DFT_def2-TZVP_SP.inp`

## Expected outputs

After running the main DFT job, the pipeline expects:

- QM energy
- optimized geometry
- Mulliken charges
- vibrational frequencies

The optimized geometry should preferably be stored as:

`geometry/step_2_2_RS_DFT_opt.xyz`

The DFT output file can be placed in:

`input/`

or elsewhere underneath this state directory.

The collector searches recursively for:

- `.out`
- `.log`
- `.txt`

## EVB relevance

This state is one of the QM reference states used before EVB parameterization.

The QM energy is an electronic/reference quantity.

It should not be interpreted directly as the final EVB activation free energy.
