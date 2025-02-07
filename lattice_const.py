# Author: Arjun S Kulathuvayal. Intellectual property. Copyright strictly restricted
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from ase.geometry import cell_to_cellpar
from ase.io import read
font = {'family': 'serif',
        'weight': 'normal',
        'size': 14}
plt.rc('font', **font)
plt.rcParams["figure.figsize"] = [8, 6]
plt.rcParams['text.usetex'] = True


def main():
    data = np.loadtxt("helmholtz-volume_fitted.dat", skiprows=28, max_rows=201)
    T = np.arange(0, 1000, 100)

    for i in range(1, 10):
        plt.plot(data[:, 0], data[:, i], '-o', label=f"{T[i]} K")
    plt.xlim(1000, 1050)
    plt.ylim(-0.90, 0.30)
    plt.xlabel(r"Volume ($A^3$)")
    plt.ylabel("Free energy (eV)")
    plt.legend(loc="lower left")
    plt.savefig("helmholtz-volume_fitted.pdf")

def vol_temp():
    a, b, c, al, be, gm, vol, strain = input_vol()
    data = np.loadtxt("helmholtz-volume_fitted.dat", skiprows=219, max_rows=10)

    plt.plot(strain, vol, '-o', label=r"$\mathcal{V}$")

    eq_volume, strains = [], []
    for each in data[:, 0]:
        plt.axhline(y=each, color='r', linestyle='--', linewidth=2)

        for i in range(len(vol) - 1):
            # Check if there's an intersection between vol[i] and vol[i+1] with the axhline
            if (vol[i] - each) * (vol[i + 1] - each) <= 0:
                # Calculate the x-value (strain) of the intersection using linear interpolation
                x_intersect = strain[i] + (each - vol[i]) * (strain[i + 1] - strain[i]) / (vol[i + 1] - vol[i])
                print(f"Intersection at y = {each} occurs at x = {x_intersect:.3f}")
                eq_volume.append(each)
                strains.append(x_intersect)
    plt.xlabel("Strain")
    plt.ylabel("Volume")
    # plt.xlim(0.01, 0.02)
    # plt.ylim(399, 400.75)
    plt.grid()
    plt.savefig("figures/vol_strain_lin_fit.png")
    plt.show()
    plt.close()



    a_0 = 3.07068
    a = []
    T = np.arange(0, 1000, 100)
    for i, eachEqVol in enumerate(eq_volume):
        a.append((1+strains[i])*a_0)

    plt.plot(T, a, '-o', label=r"$a$")
    plt.xlabel("Temperature (K)")
    plt.ylabel(r"Lattice constant a ($\AA$)")
    plt.grid()
    plt.savefig("figures/lattice_temp.png")
    plt.show()


def input_vol():
    initial_atoms = read(f"POSCAR_1x1x1.vasp", format="vasp")
    scale_factors = np.arange(-0.030, 0.031, 0.005)
    print(f"Initial volume: {initial_atoms.get_volume():.3f} \u212B^3")
    a, b, c, alpha, beta, gamma, vol = [], [], [], [], [], [], []
    for ii, strain in enumerate(scale_factors):
        bulk = initial_atoms.copy()
        cell = bulk.get_cell()

        # Apply biaxial strain: scale only the first two lattice vectors (x and y), leave the z-axis unchanged
        deformation_tensor = np.eye(3)
        deformation_tensor[0, 0] = 1 + strain  # x-axis deformation
        deformation_tensor[1, 1] = 1 + strain  # y-axis deformation
        # z-axis (third lattice vector) is left untouched

        new_cell = np.dot(deformation_tensor, cell)
        bulk.set_cell(new_cell, scale_atoms=True)  # Set the new deformed lattice

        volume = bulk.get_volume()
        print(f"Strain: {strain * 100:.1f}%, Volume: {volume:.3f} Å^3")

        contcar = read(f"qha_{volume:.3f}/CONTCAR", format="vasp")
        cell = contcar.get_cell()

        cell_parameters = cell_to_cellpar(cell)
        a.append(cell_parameters[0])
        b.append(cell_parameters[1])
        c.append(cell_parameters[2])
        alpha.append(cell_parameters[3])
        beta.append(cell_parameters[4])
        gamma.append(cell_parameters[5])
        vol.append(contcar.get_volume())

    # plt.plot(scale_factors, a, '-o', label=r"a")
    # plt.plot(scale_factors, b, '-o', label=r"b")
    # plt.plot(scale_factors, c, '-o', label=r"c")
    # plt.plot(scale_factors, alpha, '-o', label=r"alpha")
    # plt.plot(scale_factors, beta, '-o', label=r"beta")
    # plt.plot(scale_factors, gamma, '-o', label=r"gamma")
    # plt.legend()
    # plt.show()
    return a, b, c, alpha, beta, gamma, vol, scale_factors




if __name__ == '__main__':
    vol_temp()
