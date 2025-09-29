#!/bin/bash

# Batch solvation script for mutant PDB files
# Author: Automated solvation pipeline

# Set directories
INPUT_DIR="/home/hp/nayanika/github/vaccum/mutations"
OUTPUT_DIR="/home/hp/nayanika/github/vaccum/mutations"  # Changed to save in mutations folder
PARAM_DIR="/home/hp/nayanika/github/LmrR_EVB/parameters"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if qprep5 exists
if ! command -v qprep5 &> /dev/null; then
    echo "Error: qprep5 command not found!"
    echo "Make sure Q5 is properly installed and in your PATH"
    exit 1
fi

echo "Starting batch solvation process..."
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "================================================"

# Counter for processed files
count=0

# Change to output directory to avoid path issues
cd "$OUTPUT_DIR" || exit 1

# Process each PDB file in the input directory
for pdb_file in "$INPUT_DIR"/*.pdb; do
    # Check if files exist (in case no .pdb files found)
    if [ ! -e "$pdb_file" ]; then
        echo "No PDB files found in $INPUT_DIR"
        exit 1
    fi
    
    # Extract basename without extension
    basename=$(basename "$pdb_file" .pdb)
    
    # Skip problematic A11L and A91E mutations for now
    if [[ "$basename" == *"A11L"* ]] || [[ "$basename" == *"A91E"* ]]; then
        echo "Skipping $basename (incomplete sidechain - needs manual fixing)"
        echo "  --------------------------------"
        continue
    fi
    
    echo "Processing: $basename"
    
    # Create Q input file for this mutant (use full path for readpdb)
    cat > "${basename}_solvation.inp" << EOF
readlib $PARAM_DIR/qoplsaa.lib
readlib $PARAM_DIR/LMRR.lib
readprm $PARAM_DIR/qoplsaa_all.prm
readpdb $pdb_file
set solvent_pack 2.7
boundary sphere 8:CA 25.
solvate 8:CA 25. grid HOH
maketop ${basename}_solvated.top
writetop ${basename}_solvated.top
writepdb ${basename}_solvated.pdb y
quit
EOF
    
    echo "  Created input file: ${basename}_solvation.inp"
    
    # Run qprep5 
    echo "  Running solvation for $basename..."
    qprep5 "${basename}_solvation.inp" > "${basename}_solvation.log" 2>&1
    
    # Check if solvation was successful
    if [ -f "${basename}_solvated.pdb" ] && [ -f "${basename}_solvated.top" ]; then
        echo "  ✓ Successfully created ${basename}_solvated.pdb"
        echo "  ✓ Files saved: ${basename}_solvated.pdb, ${basename}_solvated.top, ${basename}_solvation.log"
        ((count++))
    else
        echo "  ✗ Error: Solvation failed for $basename"
        echo "    Check log file: ${basename}_solvation.log"
        
        # Show last few lines of log for debugging
        echo "    Last few lines of log:"
        tail -10 "${basename}_solvation.log" | sed 's/^/      /'
    fi
    
    # Keep the input file for debugging
    echo "  Input file kept: ${basename}_solvation.inp"
    
    echo "  Completed: $basename"
    echo "  --------------------------------"
done

echo "================================================"
echo "Batch solvation completed!"
echo "Successfully processed: $count files"
echo "All files saved to: $OUTPUT_DIR"
echo ""
echo "Note: A11L and A91E mutations were skipped due to incomplete sidechains."
echo "To fix these, use a structure building tool to add missing LEU/GLU atoms."

# Optional: Create summary file
echo "Solvation Summary - $(date)" > "$OUTPUT_DIR/solvation_summary.txt"
echo "Processed $count PDB files" >> "$OUTPUT_DIR/solvation_summary.txt"
echo "Skipped A11L and A91E mutations (incomplete sidechains)" >> "$OUTPUT_DIR/solvation_summary.txt"
echo "Files created:" >> "$OUTPUT_DIR/solvation_summary.txt"
ls -1 "$OUTPUT_DIR"/*_solvated.pdb >> "$OUTPUT_DIR/solvation_summary.txt" 2>/dev/null

echo "Summary saved to: $OUTPUT_DIR/solvation_summary.txt"