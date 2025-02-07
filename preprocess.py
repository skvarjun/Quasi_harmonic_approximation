#!/usr/bin/env python
import shutil
from pathlib import Path
from ase.geometry import cell_to_cellpar
from ase.io import read
from pymatgen.core import Structure
import numpy as np
import matplotlib.pyplot as plt
font = {'family': 'serif',
        'weight': 'normal',
        'size': 14}
plt.rc('font', **font)
plt.rcParams["figure.figsize"] = [8, 6]
plt.rcParams['text.usetex'] = True


def pre():
    unitcell = Structure.from_file('POSCAR_unitcell')
    supercell = unitcell.make_supercell([1, 1, 1])
    str_name = "POSCAR_1x1x1"
    supercell.to(fmt='poscar', filename=f'{str_name}.vasp')

    initial_atoms = read(f"{str_name}.vasp", format="vasp")

    # Initial volume
    print(f"Initial volume: {initial_atoms.get_volume():.3f} \u212B^3")

    for ii, strain in enumerate(scale_factors):
        bulk = initial_atoms.copy()
        cell = bulk.get_cell()

        # Apply biaxial strain: scale only the first two lattice vectors (x and y), leave the z-axis unchanged
        deformation_tensor = np.eye(3)
        deformation_tensor[0, 0] = 1 + strain  # x-axis deformation
        deformation_tensor[1, 1] = 1 + strain  # y-axis deformation

        new_cell = np.dot(deformation_tensor, cell)
        bulk.set_cell(new_cell, scale_atoms=True)  # Set the new deformed lattice

        # Calculate the new volume after deformation
        volume = bulk.get_volume()
        print(f"Strain: {strain * 100:.1f}%, Volume: {volume:.3f} Å^3")

        workdir = Path(f"qha_{volume:.3f}")
        workdir.mkdir(exist_ok=True)

        bulk.write(str(workdir / "POSCAR"), format="vasp")

        shutil.copy("vasp_input/INCAR", workdir)
        shutil.copy("vasp_input/KPOINTS", workdir)
        shutil.copy("vasp_input/POTCAR", workdir)

        print(f"Volume of sample {ii:4d}:    {bulk.get_volume():.3f} \u212B^3")
        print(f"Input files written to:   {workdir}")


def input_vol():
    initial_atoms = read("POSCAR_unitcell", format="vasp")
    a, b, c, alpha, beta, gamma = [], [], [], [], [], []
    for strain in scale_factors:
        bulk = initial_atoms.copy()
        cell = bulk.get_cell()

        deformation_tensor = np.eye(3)
        deformation_tensor[0, 0] = 1 + strain  # x-axis deformation
        deformation_tensor[1, 1] = 1 + strain  # y-axis deformation

        new_cell = np.dot(deformation_tensor, cell)
        bulk.set_cell(new_cell, scale_atoms=True)
        volume = bulk.get_volume()

        contcar = read(f"qha_{volume:.3f}/POSCAR", format="vasp")
        cell = contcar.get_cell()
        cell_parameters = cell_to_cellpar(cell)

        a.append(cell_parameters[0])
        b.append(cell_parameters[1])
        c.append(cell_parameters[2])
        alpha.append(cell_parameters[3])
        beta.append(cell_parameters[4])
        gamma.append(cell_parameters[5])

    fig, axs = plt.subplots(2, 3, figsize=(10, 7))
    axs[0, 0].plot(scale_factors, a, '-o', label='a', color='blue')
    axs[0, 0].set_title("a")
    axs[0, 0].set_xlabel("Scale factors")
    axs[0, 0].set_ylabel("Lattice constants")
    axs[0, 0].legend()

    axs[0, 1].plot(scale_factors, b, '-s', label='b', color='orange')
    axs[0, 1].set_title("b")
    axs[0, 1].set_xlabel("Scale factors")
    axs[0, 1].set_ylabel("Lattice constants")
    axs[0, 1].legend()

    axs[0, 2].plot(scale_factors, c, '-v', label='c', color='green')
    axs[0, 2].set_title("c")
    axs[0, 2].set_xlabel("Scale factors")
    axs[0, 2].set_ylabel("Lattice constants")
    axs[0, 2].legend()

    axs[1, 0].plot(scale_factors, alpha, '-o', label=r"$\alpha$", color='red')
    axs[1, 0].set_title(r"$\alpha$")
    axs[1, 0].set_xlabel("x")
    axs[1, 0].set_ylabel("Lattice angle")
    axs[1, 0].legend()

    axs[1, 1].plot(scale_factors, beta, '-s', label=r"$\beta$", color='purple')
    axs[1, 1].set_title(r"$\beta$")
    axs[1, 1].set_xlabel("x")
    axs[1, 1].set_ylabel("Lattice angle")
    axs[1, 1].legend()

    axs[1, 2].plot(scale_factors, gamma, '-v', label=r"$\gamma$", color='brown')
    axs[1, 2].set_title(r"$\gamma$")
    axs[1, 2].set_xlabel("x")
    axs[1, 2].set_ylabel("Lattice angle")
    axs[1, 2].legend()

    plt.suptitle("Variation of lattice parameters")
    plt.tight_layout()
    plt.legend()
    plt.savefig("lattice_deformation_info.png", dpi=300)
    #plt.show()


if __name__ == "__main__":
    scale_factors = np.arange(-0.050, 0.051, 0.005)
    pre()
    input_vol()