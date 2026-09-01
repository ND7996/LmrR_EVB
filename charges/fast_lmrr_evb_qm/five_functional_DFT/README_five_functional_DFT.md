# Five-functional final-QM workflow

Project:
`D:\PhD_Thesis\LmrR_EVB`

Charges:
`D:\PhD_Thesis\LmrR_EVB\charges`

Final-QM methods:

1. TPSSh-D3(BJ)
2. PBEh-3c
3. B3PW91-D3(BJ)
4. PBE0-D3(BJ) [REFERENCE]
5. M06

For every reaction step and every functional, separate RS, TS and PS
ORCA input files are generated.

The original PDB-preserving XYZ files, PDB-to-XYZ maps, xTB runners,
TS constraint files and bond maps remain in the original step folders.

Functional-specific final-QM files are stored under:

five_functional_DFT/<STEP>/<FUNCTIONAL>/

Files:
- <STEP>_RS.inp
- <STEP>_TS.inp
- <STEP>_PS.inp
- corresponding .out files after ORCA is run
- corresponding .gjf files for methods with a direct Gaussian route

PBEh-3c is run through native ORCA and therefore does not receive a
silently substituted Gaussian method.

Results:
- energies.csv
- activation_barriers.csv
- reaction_energies.csv
- deviation_from_PBE0.csv
- five_functional_DFT_summary.csv
