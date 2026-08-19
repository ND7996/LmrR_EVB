# LmrR QM reference dataset for EVB

This directory contains the QM reference layer used before EVB parameterization.

## Workflow

PDB -> XYZ -> DFT -> energy/geometry/charges/frequencies -> QM reference -> EVB calibration -> EVB/MD

## State organization

Each reaction step contains:

- RS
- TS
- PS

Each state contains:

- input/
- geometry/
- energy/
- charges/
- frequencies/

## QM reference quantities

- electronic energy
- optimized geometry
- Mulliken charges
- frequencies
- QM activation energy
- QM reaction energy

## Important interpretation

QM DeltaE^‡ is an electronic/reference-model quantity. EVB DeltaG^‡ is an ensemble free-energy quantity obtained after sampling the protein/solvent environment. They should not be treated as identical observables.

## Commands

```powershell
python qm_reference_pipeline.py --prepare
python qm_reference_pipeline.py --collect
python qm_reference_pipeline.py --validate
```
