# Author: Arjun S Kulathuvayal. Intellectual property. Copyright strictly restricted
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
font = {'family': 'serif',
        'weight': 'normal',
        'size': 14}
plt.rc('font', **font)
plt.rcParams["figure.figsize"] = [8, 6]
plt.rcParams['text.usetex'] = True

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



def tem_dep_ela_con():
    df = data_fetch()
    avg_sec_order_strain_deriv = []
    for index, row in df.iterrows():
        helm_E = row[1:].tolist()
        helm_en_per_area = []
        for i, energy in enumerate(helm_E):
            unitcel_str = Poscar.from_file(f"{deformation_strs[i]}/static_run/CONTCAR")
            lattice = unitcel_str.structure.lattice.matrix
            energy_per_atom = energy / unitcel_str.structure.num_sites
            en_per_area = energy_per_atom / np.linalg.norm(np.cross(lattice[0], lattice[1]))
            helm_en_per_area.append(en_per_area)

        coefficients = np.polyfit(strain_list, helm_en_per_area, 2)

        polynomial = np.poly1d(coefficients)
        x_fit = np.linspace(min(strain_list), max(strain_list), 100)
        y_fit = polynomial(x_fit)

        dy_dx = np.gradient(y_fit, x_fit)
        d2y_dx2 = np.gradient(dy_dx, x_fit)

        sec_der_avg = np.average(d2y_dx2)
        avg_sec_order_strain_deriv.append(sec_der_avg*160.217) #convert to GPa


    plt.plot(df['temperature'], avg_sec_order_strain_deriv, '-o', c='b', label='Pls correct')
    plt.xlabel("Temperature (K)")
    plt.ylabel("Elastic constant (GPa)")
    plt.grid(True)
    plt.savefig("QHA_manual_results/elastic_constant.png", dpi=300)
    plt.close()


def another_temp_dep():
    df = data_fetch()
    vol_dat = np.loadtxt("QHA_phonopy_results/volume_temperature.dat")

    curvatures_aka_ela_constants = []
    temps = []
    for index, row in df.iterrows():
        helm_E = row[1:].tolist()
        helm_en_per_area = []
        for energy in helm_E:
            helm_en_per_area.append(energy)

        coefficients = np.polyfit(strain_list, helm_en_per_area, 2)

        curvatures_aka_ela_con = np.round((coefficients[0] * 160.217 / vol_dat[index][1]), 4)
        curvatures_aka_ela_constants.append(curvatures_aka_ela_con)
        temps.append(vol_dat[index][0])

    plt.plot(temps, curvatures_aka_ela_constants, '-s', c='r', label='Pls correct')
    plt.xlabel("Temperature (K)")
    plt.ylabel("Elastic constant (GPa)")
    plt.grid(True)
    plt.savefig("QHA_phonopy_results/elastic_constants_try_2.png", dpi=300)
    plt.close()

if __name__ == '__main__':
    strain_list = np.arange(-0.050, 0.051, 0.005)
    qha_dirs = sorted(
        [float(d.replace('qha_', '')) for d in os.listdir(os.getcwd()) if os.path.isdir(d) and d.startswith("qha")],
        reverse=False)
    deformation_strs = [f"qha_{j:.3f}" for j in qha_dirs]
    tem_dep_ela_con()
    #another_temp_dep()
