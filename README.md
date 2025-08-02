
# LmrR_pAF Alanine Scanning — Free Energy Perturbation (FEP) Workflow

This repository documents an automated and modular pipeline for **computational alanine scanning** of the **artificial enzyme LmrR_pAF** involved in catalyzing the **Friedel–Crafts alkylation** and **hydrazone formation** reactions. The project combines **Empirical Valence Bond (EVB)** simulations using the **Q software suite** with structural mutagenesis and analysis of catalytic hotspots.

 *This study aims to reproduce and extend findings from experimental alanine scanning (Roelfes Lab) by evaluating the impact of active-site residue mutations on activation free energies and mechanistic performance.*

---

## Project Scope

**Protein**: LmrR (Lactococcal multidrug resistance regulator)  
**Variant**: LmrR_pAF (p-aminophenylalanine at position V15)  
**Target Reactions**:  
-  *Friedel–Crafts Alkylation* (primary EVB model validation)  
-  *Hydrazone Formation* (future application)

---

## Workflow Overview

The pipeline performs alanine scanning by systematically mutating catalytic site residues of LmrR to alanine, generating structures and calculating reaction profiles using the EVB method. It comprises 7 major steps:

---

### STEP 1: Structure Preparation

- Input: WT LmrR structure (PDB: `6I8N`)
- Mutation scripts apply **single-point alanine substitutions** using PyMOL automation.
- Outputs: Mutant PDB files with consistent naming (`LmrR_MUTX.pdb`)

---

### STEP 2: Solvation with Qprep

- Generates solvated systems and topologies using `qprep5`
- Applies explicit solvation and periodic boundary conditions
- Output: `${base_name}_solvated.pdb`, `${base_name}_solvated.top`

---

### STEP 3: FEP Setup

- Script: `makeFEP.py`
- Uses `qmap` files tailored to LmrR mutants
- Outputs lambda maps and perturbation parameters for EVB

---

### STEP 4: Input Generation for Relaxation

- Script: `q_genrelax.py`
- Produces relaxation `.inp` files and execution scripts
- Inputs: `--genrelax.proc`, `--top`, `--pdb`, `--fep`
- Output directory: `minim/`

---

### STEP 5: PDB Minimization for FEP

- Performs relaxation runs (`rx relax_012.re`)
- Reconstructs minimized structure via `qprep5`
- Outputs: `minim.pdb` (used for EVB)

---

### STEP 6: FEP Input Generation

- Script: `q_genfeps.py`
- Generates replicas for FEP sampling (default: 3 replicas, 51 frames)
- Inputs: `--genfeps.proc`, `minim.pdb`, `--prefix replica`, `--rs run_qdyn_5.sh`

---

### STEP 7: FEP Analysis

- Tools: `q_mapper.py`, `q_analysefeps.py`
- Analysis steps:
  - Mapping FEP data across replicas
  - Extracting statistics from `test.out`
  - Exporting structured outputs (`stats_output.txt`, CSV, LaTeX tables)

---

