# Author: Arjun S Kulathuvayal. Intellectual property. Copyright strictly restricted
import sys, yaml, os
from os import getcwd
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from ase.io import read
from phonopy import PhonopyQHA
from pymatgen.core import Structure
from scipy.interpolate import griddata

font = {'family': 'serif',
        'weight': 'normal',
        'size': 14}
plt.rc('font', **font)
plt.rcParams["figure.figsize"] = [8, 6]
plt.rcParams['text.usetex'] = True
from yaml import CLoader as Loader
pd.set_option('display.max_rows', 14)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)

def get_vol_en():
    volumes = []
    energies = []
    for line in open("QHA_phonopy_process/e-v.dat"):
        v, e = line.split()
        volumes.append(float(v))
        energies.append(float(e))

    entropy = []
    cv = []
    fe = []
    for index in range(1, len(volumes) + 1):
        filename = "QHA_phonopy_process/thermal_properties.yaml-%d" % index
        print("Reading %s" % filename)
        thermal_properties = yaml.load(open(filename), Loader=Loader)["thermal_properties"]
        temperatures = [v["temperature"] for v in thermal_properties]
        cv.append([v["heat_capacity"] for v in thermal_properties])
        entropy.append([v["entropy"] for v in thermal_properties])
        fe.append([v["free_energy"] for v in thermal_properties])
    return volumes, energies, temperatures, cv, entropy, fe


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


def energy_strain():
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
    a_0 = 3.09247
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


def phase_data():
    volumes, energies, temperatures, cv, entropy, fe = get_vol_en()
    P_vals = np.arange(0, 100, 1)
    GE = []
    df = pd.DataFrame(columns=['P']+temperatures[:-1])

    for p in P_vals:
        qha = PhonopyQHA(
            volumes,
            energies,
            pressure=p,
            eos='murnaghan',
            temperatures=temperatures,
            free_energy=np.transpose(fe),
            cv=np.transpose(cv),
            entropy=np.transpose(entropy),
            #t_max=4010,
            verbose=False,
        )
        x = qha.gibbs_temperature
        GE.append(x)
        df.loc[len(df)] = np.insert(x, 0, p)
    df.to_csv("QHA_manual_results/phase_data.csv")
    GE = np.array(GE)
    T_values = np.array(temperatures[:-1])
    T_mesh, P_mesh = np.meshgrid(T_values, P_vals)

    T_flat = T_mesh.flatten()
    P_flat = P_mesh.flatten()
    GFE_flat = GE.flatten()

    T_fine = np.linspace(T_values.min(), T_values.max(), 100)
    P_fine = np.linspace(P_vals.min(), P_vals.max(), 100)
    T_fine_mesh, P_fine_mesh = np.meshgrid(T_fine, P_fine)

    GFE_fine = griddata((T_flat, P_flat), GFE_flat, (T_fine_mesh, P_fine_mesh), method='cubic')
    contour = plt.contourf(T_fine_mesh, P_fine_mesh, GFE_fine, levels=50, cmap="viridis")

    cbar = plt.colorbar(contour)
    cbar.set_label("Gibbs Free Energy (eV)")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Pressure (GPa)")
    plt.title("Gibbs Free Energy in PT space")
    plt.grid(True, linestyle="--", alpha=0.3)
    plt.savefig("QHA_manual_results/phase_data.png", dpi=300)
    #plt.show()

if __name__ == '__main__':
    os.makedirs('QHA_manual_results', exist_ok=True)
    strain_list = np.arange(-0.050, 0.051, 0.005)
    qha_dirs = sorted([float(d.replace('qha_', '')) for d in os.listdir(os.getcwd()) if os.path.isdir(d) and d.startswith("qha")], reverse=False)
    deformation_strs = [f"qha_{j:.3f}" for j in qha_dirs]
    extract_thermal_yaml()
    energy_strain()
    QHA_lattice_change()
    phase_data()
