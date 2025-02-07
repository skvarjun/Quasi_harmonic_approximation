#!/usr/bin/env python
import os
from ase.io import read
from pymatgen.io.vasp.outputs import Outcar
import shutil

def main():
    os.makedirs('QHA_phonopy_process', exist_ok=True)
    with open("QHA_phonopy_process/e-v.dat", "w") as f:
        i = 1
        qha_dirs = sorted([float(d.replace('qha_', '')) for d in os.listdir(os.getcwd()) if os.path.isdir(d) and d.startswith("qha")], reverse=False)
        for folder in [f"qha_{j:.3f}" for j in qha_dirs]:
            full_path = os.path.join(os.getcwd(), folder)
            if os.path.isdir(full_path) and folder.startswith("qha"):
                print(f"System {folder}")
                strt = read(f"{folder}/static_run/CONTCAR", format="vasp")
                outcar = Outcar(f"{folder}/static_run/OUTCAR")
                vol = strt.get_volume()
                e_tot = outcar.final_energy
                f.write(f"{vol:20.10f} {e_tot:20.10e}\n")
                shutil.copy(f"{folder}/static_run/CONTCAR", f"QHA_phonopy_process/POSCAR-{i}")
                shutil.copy(f"{folder}/thermal_properties.yaml", f"QHA_phonopy_process/thermal_properties.yaml-{i}")
                i += 1


if __name__ == "__main__":
    main()
