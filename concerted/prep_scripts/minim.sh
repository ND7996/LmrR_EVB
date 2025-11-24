#!/bin/bash
# Run genrelax for WT and all mutants, results stored in /home/hp/results/LMRR_PAF

run_genrelax() {
    NAME=$1
    TOP=$2
    PDB=$3

    OUTDIR="/home/hp/results/LMRR/${NAME}"

    # Remove old directory if exists
    if [ -d "$OUTDIR" ]; then
        echo "Removing existing directory: $OUTDIR"
        rm -rf "$OUTDIR"
        # Wait until it's really gone
        while [ -d "$OUTDIR" ]; do
            sleep 1
        done
    fi

    # Run genrelax (will create OUTDIR itself)
    q_genrelax.py /home/hp/nayanika/github/LmrR_EVB/input/genrelax.proc \
        --top "$TOP" \
        --pdb "$PDB" \
        --fep /home/hp/nayanika/github/LmrR_EVB/input/fep/LMRR_WT4.fep \
        --outdir "$OUTDIR" \
        --rest top \
        --rs /home/hp/nayanika/github/LmrR_EVB/cluster_scripts/run_qdyn_5.sh
}

# Wild-type
run_genrelax WT \
    /home/hp/nayanika/github/LmrR_EVB/structures/LMRR_WT2_solvated.top \
    /home/hp/nayanika/github/LmrR_EVB/structures/LMRR_WT2_solvated.pdb

# Original mutants
run_genrelax ASN18A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_ASN18A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_ASN18A_solvated.pdb
run_genrelax ASN87A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_ASN87A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_ASN87A_solvated.pdb
run_genrelax GLU7A  /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_GLU7A_solvated.top  /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_GLU7A_solvated.pdb
run_genrelax LEU17A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_LEU17A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_LEU17A_solvated.pdb
run_genrelax LYS21A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_LYS21A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_LYS21A_solvated.pdb
run_genrelax MET88A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_MET88A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_MET88A_solvated.pdb
run_genrelax PHE92A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_PHE92A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_PHE92A_solvated.pdb
run_genrelax SER94A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_SER94A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_SER94A_solvated.pdb
run_genrelax SER96A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_SER96A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_SER96A_solvated.pdb
run_genrelax ASP99A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_ASP99A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_ASP99A_solvated.pdb
run_genrelax TRP95A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_TRP95A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_TRP95A_solvated.pdb

# New mutants
run_genrelax ASN14A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_ASN14A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_ASN14A_solvated.pdb
run_genrelax LEU9A  /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_LEU9A_solvated.top  /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_LEU9A_solvated.pdb
run_genrelax ARG10A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_ARG10A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_ARG10A_solvated.pdb
run_genrelax LYS100A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_LYS100A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_LYS100A_solvated.pdb
run_genrelax GLU103A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_GLU103A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_GLU103A_solvated.pdb
run_genrelax ILE15A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_ILE15A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_ILE15A_solvated.pdb
run_genrelax VAL98A /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_VAL98A_solvated.top /home/hp/nayanika/github/LmrR_EVB/structures/mutations/LMRR_WT2_VAL98A_solvated.pdb