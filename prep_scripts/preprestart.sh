#!/bin/bash
# Run qprep5 to generate minim.pdb for WT and mutants from relax_012.re
# It will use the .top file that is already present in each folder.

PARAMS_DIR="/home/hp/nayanika/github/LmrR_EVB/parameters"
RESULTS_DIR="/home/hp/results/LMRR_PAF"

# Mutants list (each has its own folder in RESULTS_DIR)
MUTANTS="ASN18A ASN87A ASP99A GLU7A LEU17A LYS21A MET88A PHE92A SER94A SER96A WT TRP95A"

for NAME in $MUTANTS; do
    echo "Processing $NAME ..."

    FOLDER="${RESULTS_DIR}/${NAME}"

    TOPFILE=$(ls "$FOLDER"/*.top 2>/dev/null)
    RXFILE="${FOLDER}/relax_012.re"
    OUTPDB="${FOLDER}/minim.pdb"
    INPFILE="${FOLDER}/preprestart.inp"

    if [ ! -f "$TOPFILE" ]; then
        echo "❌ No .top file found in $FOLDER, skipping."
        continue
    fi

    if [ ! -f "$RXFILE" ]; then
        echo "❌ No relax_012.re file in $FOLDER, skipping."
        continue
    fi

    # Create input file for qprep5
    cat > "$INPFILE" <<EOF
readlib ${PARAMS_DIR}/qoplsaa.lib
readlib ${PARAMS_DIR}/LMRR.lib
readprm ${PARAMS_DIR}/qoplsaa_all.prm
readtop ${TOPFILE}
rx ${RXFILE}
writepdb ${OUTPDB} y
quit
EOF

    # Run qprep5
    qprep5 "$INPFILE"

done

echo "✅ All minim.pdb files generated (where inputs existed)."
