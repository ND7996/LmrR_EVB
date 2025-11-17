#!/bin/bash
# Run qprep5 + q_genfeps for WT and mutants
# Paths
PARAMS_DIR="/home/hp/nayanika/github/LmrR_EVB/parameters"
INPUT_DIR="/home/hp/nayanika/github/LmrR_EVB/input"
RESULTS_DIR="/home/hp/results/LMRR"
MUTATIONS_DIR="/home/hp/nayanika/github/LmrR_EVB/structures/mutations"
RS_SCRIPT="/home/hp/nayanika/github/LmrR_EVB/cluster_scripts/run_qdyn_5.sh"

# Mutants list (WT included at end)
MUTANTS="ASN14A GLU103A VAL98A ARG10A LYS100A LEU9A ILE15A"

for NAME in $MUTANTS; do
    echo "🔹 Processing $NAME ..."
    FOLDER="${RESULTS_DIR}/${NAME}"
    TOPFILE=$(ls "${FOLDER}"/*.top 2>/dev/null | head -n 1)
    RXFILE="${FOLDER}/relax_012.re"
    
    # Find the solvated PDB file for this mutant
    # WT2 is special: LMRR_WT2_solvated.pdb (no mutation name)
    # Others: LMRR_WT2_{MUTATION}_solvated.pdb
    if [ "$NAME" = "WT2" ]; then
        PDBFILE="${MUTATIONS_DIR}/LMRR_WT2_solvated.pdb"
    else
        PDBFILE="${MUTATIONS_DIR}/LMRR_WT2_${NAME}_solvated.pdb"
    fi
    
    # Check if required files exist
    if [ ! -f "$TOPFILE" ]; then
        echo "❌ Missing .top file in $FOLDER, skipping"
        continue
    fi
    if [ ! -f "$RXFILE" ]; then
        echo "❌ Missing relax_012.re in $FOLDER, skipping"
        continue
    fi
    if [ ! -f "$PDBFILE" ]; then
        echo "❌ Missing $PDBFILE, skipping"
        continue
    fi
    
    # Change to the mutant directory
    cd "$FOLDER"
    
    # Run qprep5 with direct command line input
    echo "   ➡️ Running qprep5 for $NAME (using $(basename $PDBFILE))"
    qprep5 <<EOF
readlib ${PARAMS_DIR}/qoplsaa.lib
readlib ${PARAMS_DIR}/LMRR.lib
readprm ${PARAMS_DIR}/qoplsaa_all.prm
readtop ${TOPFILE}
rx ${RXFILE}
writepdb ${PDBFILE} y
quit
EOF
    
    # Remove existing replica folders if they exist
    echo "   🧹 Cleaning up existing replica folders..."
    rm -rf replica*
    
    # Run q_genfeps.py using the solvated.pdb
    echo "   ➡️ Running q_genfeps.py for $NAME"
    q_genfeps.py "${INPUT_DIR}/genfeps.proc" \
        --pdb "$PDBFILE" \
        relax_012.inp relax \
        --repeats 5 \
        --frames 51 \
        --fromlambda 1.0 \
        --prefix replica \
        --rs "$RS_SCRIPT"
    
    # Return to original directory
    cd - > /dev/null
done

echo "✅ All jobs completed."