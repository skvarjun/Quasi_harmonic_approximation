"""Example of QHA calculation by Al."""
import os
import numpy as np
import pandas as pd
import yaml
from yaml import CLoader as Loader
from phonopy import PhonopyQHA
import matplotlib.pyplot as plt
font = {'family': 'serif',
        'weight': 'normal',
        'size': 14}
plt.rc('font', **font)
plt.rcParams["figure.figsize"] = [8, 6]
plt.rcParams['text.usetex'] = True


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

def main():
    volumes, energies, temperatures, cv, entropy, fe = get_vol_en()
    qha = PhonopyQHA(
        volumes,
        energies,
        pressure=1,
        temperatures=temperatures,
        free_energy=np.transpose(fe),
        cv=np.transpose(cv),
        entropy=np.transpose(entropy),
        t_max=1000,
        verbose=True,
    )

    os.makedirs("QHA_phonopy_results", exist_ok=True)
    qha.plot_helmholtz_volume().savefig("QHA_phonopy_results/plot_helmholtz_volume.png")
    qha.plot_volume_temperature().savefig("QHA_phonopy_results/plot_volume_temperature.png")
    qha.plot_thermal_expansion().savefig("QHA_phonopy_results/plot_thermal_expansion.png")
    # plot = qha.plot_volume_expansion()
    # if plot:
    #     plot.show()
    qha.plot_gibbs_temperature().savefig("QHA_phonopy_results/plot_gibbs_temperature.png")
    qha.plot_bulk_modulus_temperature().savefig("QHA_phonopy_results/plot_bulk_modulus_temperature.png")
    qha.plot_heat_capacity_P_numerical().savefig("QHA_phonopy_results/plot_heat_capacity_P_numerical.png")
    qha.plot_heat_capacity_P_polyfit().savefig("QHA_phonopy_results/plot_heat_capacity_P_polyfit.png")
    qha.plot_gruneisen_temperature().savefig("QHA_phonopy_results/plot_gruneisen_temperature.png")

    qha.write_volume_temperature("QHA_phonopy_results/volume_temperature.dat")


def phase_data():
    volumes, energies, temperatures, cv, entropy, fe = get_vol_en()
    P_vals = np.linspace(0, 150, 100)

    GE = []
    for p in P_vals:
        qha = PhonopyQHA(
            volumes,
            energies,
            pressure=p,
            temperatures=temperatures,
            free_energy=np.transpose(fe),
            cv=np.transpose(cv),
            entropy=np.transpose(entropy),
            t_max=1000,
            verbose=False,
        )

        x = qha.gibbs_temperature
        GE.append(x)


    P_grid, T_grid = np.meshgrid(P_vals, range(0, 1000, 10))
    G = np.array(GE)

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(P_grid, T_grid, G, cmap='viridis',  alpha=0.8)

    ax.set_xlabel('Pressure (GPa)')
    ax.set_ylabel('Temperature (K)')
    ax.set_zlabel('Gibbs Free Energy (eV)')
    ax.set_title('Finding minima of Gibbs Free Energy')
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
    plt.savefig('QHA_phonopy_results/phase_data.png', dpi=300)


if __name__ == "__main__":
    main()
    phase_data()
