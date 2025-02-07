# Quasi-Harmonic Approximation (QHA) workflow

## Overview
This guide outlines the steps required to estimate temperature-dependent mechanical properties of periodic systems using the Quasi-Harmonic Approximation (QHA) method. The workflow involves structure relaxation, deformation, force constant extraction, and processing of thermal properties.

## Prerequisites
- You have to have VASP license.
- Required packages
  - conda environment - call it phonon.
  - Install phonopy from Atsushi Togo.
  - Pymatgen, ASE and Pandas
- Install and activate the phonon environment:
  ```sh
  conda activate phonon
  ```
- Ensure `VASP` is installed and properly configured with mpirun.
- Set the curresponding paths and libraries in sh files with VASP run.
- Python scripts (`preprocess.py`, `postprocess.py`, `QHA_phonopy_process.py`, `QHA_manual_process.py`, `QHA_elastic_constants.py`) should be in the working directory.

## Workflow

### Step 1: Relax Input Structure
```sh
relax input structure under vasp_input/GO/
cp vasp_input/GO/CONTCAR POSCAR_unitcell
```

### Step 2: Generate Deformed Structures
```sh
python preprocess.py
```
- Creates bi-axial deformation structures.
- Deformation magnitude and direction can be controlled in `preprocess.py`.

### Step 3: Static Calculations
```sh
sh exec_vasp_run_for_static.sh
```
- Creates `deformed_structures/static_run/` and runs static calculations.

### Step 4: DFPT Calculations for Vibrational Energy
```sh
sh exec_vasp_run_for_DFPT.sh
```
- Creates a DFPT supercell and runs DFPT calculations in VASP.

### Step 5: Extract Force Constants
```sh
sh post_process_DFPT.sh
```
- Collects force constants from `vasprun.xml`.
- Generates `thermal_properties.yaml`.

### Step 6: Generate Energy-Volume Data
```sh
python postprocess.py
```
- Creates `QHA_phonpy_process/e-v.dat`.

### Step 7: Process Thermal Properties from phonopy
```sh
python QHA_phonopy_process.py
```
- Generates thermal output files under `QHA_phonopy_results/`.

### Step 8: Compute temperature dependant Lattice Parameter
```sh
python QHA_manual_process.py
```
- Computes lattice constant's variations with temperature.

### Step 9: Compute temperature dependant Elastic constants
```sh
python QHA_elastic_constants.py
```
- Computes lattice constant's variations with temperature.


## Output Files
- `thermal_properties.yaml`: Contains extracted force constants.
- `QHA_phonopy_results/`: Stores processed thermal properties.
- `QHA_phonpy_process/e-v.dat`: Energy-volume relationship data.
- Computed lattice and elastic property variations.

## Notes
- Modify `preprocess.py` to control deformation parameters.
- Ensure VASP is configured correctly for DFPT calculations.

## License
This guide is intended for research purposes. Use at your own discretion.
