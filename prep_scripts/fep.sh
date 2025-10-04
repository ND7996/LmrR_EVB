#!/bin/bash
# Run qprep5 + q_genfeps for WT and mutants
# Paths
PARAMS_DIR="/home/hp/nayanika/github/LmrR_EVB/parameters"
INPUT_DIR="/home/hp/nayanika/github/LmrR_EVB/input"
RESULTS_DIR="/home/hp/results/LMRR_PAF"
RS_SCRIPT="/home/hp/nayanika/github/LmrR_EVB/cluster_scripts/run_qdyn_5.sh"

# Mutants list (WT included at end)
MUTANTS="ASN18A  ASN87A  ASP99A  GLU7A  LEU17A  LYS21A  MET88A  PHE92A  SER94A  SER96A TRP95A"

for NAME in $MUTANTS; do
    echo "🔹 Processing $NAME ..."
    FOLDER="${RESULTS_DIR}/${NAME}"
    TOPFILE=$(ls "${FOLDER}"/*.top 2>/dev/null | head -n 1)
    RXFILE="${FOLDER}/relax_012.re"
    PDBFILE="${FOLDER}/minim.pdb"
    
    # Check if required files exist
    if [ ! -f "$TOPFILE" ]; then
        echo "❌ Missing .top file in $FOLDER, skipping"
        continue
    fi
    if [ ! -f "$RXFILE" ]; then
        echo "❌ Missing relax_012.re in $FOLDER, skipping"
        continue
    fi
    
    # Change to the mutant directory
    cd "$FOLDER"
    
    # Run qprep5 with direct command line input
    echo "   ➡️ Running qprep5 for $NAME"
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
    
    # Run q_genfeps.py if pdb was generated
    if [ -f "$PDBFILE" ]; then
        echo "   ➡️ Running q_genfeps.py for $NAME"
        q_genfeps.py "${INPUT_DIR}/genfeps.proc" \
            --pdb "$PDBFILE" \
            relax_012.inp relax \
            --repeats 5 \
            --frames 51 \
            --fromlambda 1.0 \
            --prefix replica \
            --rs "$RS_SCRIPT"
    else
        echo "❌ No minim.pdb for $NAME after qprep5, skipping genfeps"
    fi
    
    # Return to original directory
    cd - > /dev/null
done

echo "✅ All jobs completed."