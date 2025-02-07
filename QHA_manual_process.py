# Author: Arjun S Kulathuvayal. Intellectual property. Copyright strictly restricted
import sys, yaml, os
from os import getcwd
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from ase.io import read
from pymatgen.core import Structure

font = {'family': 'serif',
        'weight': 'normal',
        'size': 14}
plt.rc('font', **font)
plt.rcParams["figure.figsize"] = [8, 6]
plt.rcParams['text.usetex'] = True

pd.set_option('display.max_rows', 14)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)


def extract_thermal_yaml():
    for folder_name in deformation_strs:
        data = yaml.load(open(f'{folder_name}/thermal_properties.yaml'), Loader=yaml.FullLoader)
        g = open(f'{folder_name}/static_run/OSZICAR', 'r')
        g1 = Structure.from_file(f'{folder_name}/static_run/POSCAR')
        natoms = len(g1.sites)
        l = [line for line in g]
        energy = float(l[-1].split()[2])
        zpe = energy / natoms

        struct = Structure.from_file(f'{folder_name}/POSCAR-unitcell')
        natoms = len(struct.sites)
        df = pd.DataFrame()
        temp_array = [];
        energy_array = [];
        cv_array = []
        data_th = data['thermal_properties']

        for m in range(len(data_th)):
            temp = data_th[m]['temperature']
            en1 = (data_th[m]['free_energy'] * 0.01) / (natoms)
            cv = (data_th[m]['heat_capacity']) / (natoms)
            en1 = en1
            en = en1 + zpe
            temp_array.append(temp)
            cv_array.append(cv)
            energy_array.append(en)

        df['temperature'] = temp_array
        df['helm_energy'] = energy_array
        df['Cv'] = cv_array
        df.to_csv(f'{folder_name}/energy_temp.dat', index=False, sep='\t')


def data_fetch():
    df2 = pd.DataFrame()
    for folder_name in deformation_strs:
        df_en = pd.read_csv(f'{folder_name}/energy_temp.dat', usecols=[0, 1],
                            names=['temperature', f'helm_energy_{folder_name}'], skiprows=1, sep='\t')
        if len(df2) == 0:
            df2 = df_en.copy()
        else:
            df2 = df2.merge(df_en, on='temperature')
    return df2


def graph():
    df = data_fetch()

    energies_temp = df.iloc[-1][1:].values

    coefficients = np.polyfit(strain_list, energies_temp, 2)
    polynomial = np.poly1d(coefficients)
    x_fit = np.linspace(min(strain_list), max(strain_list), 100)
    y_fit = polynomial(x_fit)

    plt.scatter(strain_list, energies_temp, c='r', label='Actual data')
    plt.plot(x_fit, y_fit, label='Quadratic Fit', linestyle=':', color='blue')
    plt.xlabel("Strain")
    plt.ylabel("Energy (eV/atom)")
    plt.axhline(y=min(y_fit), linestyle='--', color='black', label=f'Min free energy = {min(y_fit):.2f}')
    x_min = -coefficients[1] / (2 * coefficients[0])
    plt.axvline(x=x_min, color='black', linestyle='--', label=f'Min strain = {x_min:.4f}')
    plt.title(f"Quadratic Fit a = {coefficients[0]:.2f}, b = {coefficients[1]:.2f}, c = {coefficients[2]:.2f}")
    plt.legend()
    plt.savefig("QHA_manual_results/energy_strain.png")
    plt.close()
    #plt.show()


def min_strain_for_eachT():
    df = data_fetch()
    min_energies, min_strains =[], []
    for index, row in df.iterrows():
        energy_for_T = df.iloc[index, 1:].to_numpy()
        coefficients = np.polyfit(strain_list, energy_for_T, 2)
        polynomial = np.poly1d(coefficients)
        x_fit = np.linspace(min(strain_list), max(strain_list), 100)
        y_fit = polynomial(x_fit)
        min_energy = min(y_fit)
        min_strain = -coefficients[1] / (2 * coefficients[0]) #x axis
        min_energies.append(min_energy)
        min_strains.append(min_strain)
    return min_energies, min_strains


def QHA_lattice_change():
    min_energy_list, min_strain_list = min_strain_for_eachT()
    a_0 = 3.16312
    a = []
    T = np.arange(0, 1010, 10)
    for i, temp in enumerate(T):
        a.append((1+min_strain_list[i])*a_0)

    plt.plot(T, a, '-o', c='r', label=r"$a$")
    plt.xlabel("Temperature (K)")
    plt.ylabel(r"Lattice constant a ($\AA$)")
    plt.grid()
    plt.legend()
    plt.savefig("QHA_manual_results/lattice_change.png")
    plt.close()


if __name__ == '__main__':
    os.makedirs('QHA_manual_results', exist_ok=True)
    strain_list = np.arange(-0.050, 0.051, 0.005)
    qha_dirs = sorted([float(d.replace('qha_', '')) for d in os.listdir(os.getcwd()) if os.path.isdir(d) and d.startswith("qha")], reverse=False)
    deformation_strs = [f"qha_{j:.3f}" for j in qha_dirs]
    extract_thermal_yaml()
    graph()
    QHA_lattice_change()
