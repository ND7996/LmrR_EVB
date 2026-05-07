#start stepwise mechanism structure prep

1- start with tripeptide generation

**Preparing a Tripeptide In Maestro**

1. Create the amino acid peptide chain (add fragments) using 3D builder of desired amino acids (add fragments - amino acid)
2. Mutate the Cys to Selenocystiene using the - Set element of selected atom in 3D builder
3. Select other edits - change atom properties - change residue name to SEL

4. All Parameter file -  **gpxqoplsaa_all.prm**
5. One PDB file with the protein and the substrate - **GPX6_h2o2.pdb**
Creating the params for selenol and seleninic acid


2- Parameters for HOX-PAF


set SCHRODINGER=C:\Program Files\Schrodinger2022-1

"%SCHRODINGER%\ligprep" -imae HOX.mae -omae HOX_prep.mae -ph 7.0 -pht 0.0 -ff OPLS4 -WAIT

"%SCHRODINGER%\epik" -imae HOX_prep.mae -omae HOX_epik.mae -ph 7.0 -pht 0.0 -WAIT

"%SCHRODINGER%\bmin" HOX_epik.mae -WAIT

"%SCHRODINGER%\ffld_server" -imae maestro HOX_epik.mae -version 14 -print_parameters -out_file HOX_PARAM.log

"%SCHRODINGER%\jaguar" run -jobname HOX_charges HOX_epik.mae -dft b3lyp -basis 6-31gs -igeopt 1 -WAIT

3- Parameters for Substrates

HEX_ENAL
WATER
imminium
2_methylindole
iminium enamine tautomerization



