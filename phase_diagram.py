# Author: Arjun S Kulathuvayal. Intellectual property. Copyright strictly restricted
from itertools import combinations
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from sympy.physics.quantum.gate import Phase

font = {'family': 'serif',
        'weight': 'normal',
        'size': 14}
plt.rc('font', **font)
plt.rcParams["figure.figsize"] = [8, 6]
plt.rcParams['text.usetex'] = True


def find_boundary(A, B, Temperature):
    for i in range(len(A) - 1):
        if (A[i] - B[i]) * (A[i + 1] - B[i + 1]) < 0:
            t = (B[i] - A[i]) / ((A[i + 1] - A[i]) - (B[i + 1] - B[i]))
            boundary_point = A[i] + t * (A[i + 1] - A[i])
            corresponding_temperature = Temperature[i] + t * (Temperature[i + 1] - Temperature[i])

            return boundary_point, corresponding_temperature
    print("No intersection found")
    return None, None


def expt():
    T = range(0, 4000, 200)
    P = range(0, 20, 1)

    bcc = pd.read_csv('BCC/QHA_manual_results/phase_data.csv')
    fcc = pd.read_csv('FCC/QHA_manual_results/phase_data.csv')
    hcp = pd.read_csv('HCP/QHA_manual_results/phase_data.csv')
    for i, p in enumerate(P):
        bcc_GFE = bcc.iloc[i, 2:].to_numpy()
        fcc_GFE = fcc.iloc[i, 2:].to_numpy()
        #hcp_GFE = hcp.iloc[i, 2:].to_numpy()/2 #2atoms present in HCP eV/mol
        plt.plot(T, bcc_GFE, '-', label='BCC', color='red')
        plt.plot(T, fcc_GFE, '-', label='FCC', color='blue')
        #plt.plot(hcp_GFE, T, '-', label='HCP', color='green')
        plt.xlabel("Temperature (K)")
        plt.ylabel("Gibbs Free Energy (eV/mol)")
    #plt.legend(2)
    plt.show()

def calc_phase(p1, p2, p1_ph, p2_ph):
    p1_avg = np.average(p1)
    p2_avg = np.average(p2)
    if p1_avg < p2_avg:
        return p1_ph
    else:
        return p2_ph


def any_phase_change(x, y, threshold):
    dy_dx = np.gradient(y, x)
    d2y_dx2 = np.gradient(dy_dx, x)
    kink_indices = np.where(np.abs(d2y_dx2) > threshold)[0]
    if len(kink_indices) > 0:
        print("Lines crosses")
        return "Yes", kink_indices[0]
    else:
        min_index = np.argmin(y)
        return "Yes", min_index


def main(phase_1, phase_2, n_atom_1, n_atom_2):
    T = range(0, 4010, 10)
    P = range(0, 100, 1)
    p1 = pd.read_csv(f'{phase_1}/QHA_manual_results/phase_data.csv')
    p2 = pd.read_csv(f'{phase_2}/QHA_manual_results/phase_data.csv')

    boundary, boundary_temp, boundary_pressure = [], [], []
    first_phase, second_phase = [], []
    for i, p in enumerate(P):
        #p = 44
        p1_GFE = p1.iloc[p, 2:].to_numpy()/n_atom_1
        p2_GFE = p2.iloc[p, 2:].to_numpy()/n_atom_2
        phase = np.abs(p1_GFE - p2_GFE)
        phase_change, T_index = any_phase_change(T, phase, threshold=0.000001)
        print(phase_change, T_index, T[T_index])

        fig, axs = plt.subplots(1, 2, figsize=(8, 4))
        axs[0].plot(T, p1_GFE, '-', label=phase_1, color='red')
        axs[0].plot(T, p2_GFE, '-', label=phase_2, color='blue')
        axs[0].set_title(f"Pressure = {p} GPa")
        axs[0].set_xlabel("Temperature (K)")
        axs[0].set_ylabel("G (eV/mol)")
        axs[0].legend()
        axs[0].grid()
        axs[1].plot(T, phase, '-', label=r"$\Delta G$", color='green')
        axs[1].set_xlabel("Temperature (K)")
        if phase_change == "Yes":
            phase_status_in = calc_phase(p1_GFE[0:T_index], p2_GFE[0:T_index], phase_1, phase_2)
            phase_status_fi = calc_phase(p1_GFE[T_index:], p2_GFE[T_index:], phase_1, phase_2)
            axs[1].set_title(f"[0, {T[T_index]:.2f}] : {phase_status_in} \n[{T[T_index]:.2f}, {max(T)}] : {phase_status_fi}")
        else:
            phase_status  = calc_phase(p1_GFE, p2_GFE)
            if phase_status == 'P1':
                axs[1].set_title(f"Phase throughout is {phase_1}")
            elif phase_status == 'P2':
                axs[1].set_title(f"Phase throughout is {phase_2}")

        axs[1].set_ylabel(r"$\|G_{{{phase_1}}} - G_{{{phase_2}}}\|$".format(phase_1=phase_1, phase_2=phase_2))
        axs[1].legend()
        axs[1].grid()
        plt.tight_layout()
        plt.legend()
        plt.show()
        plt.close()

        boundary.append((p1_GFE[T_index]+p2_GFE[T_index])/2)
        boundary_temp.append(T[T_index])
        boundary_pressure.append(p)
        first_phase.append(phase_status_in)
        second_phase.append(phase_status_fi)


    plt.plot(P, boundary_temp, '-o', label=f'{phase_1}/{phase_2} boundary', color='green' )
    plt.ylabel("Temperature (K)")
    plt.xlabel("Pressure (GPa)")
    plt.legend()
    plt.savefig(f"phase_diagram_{phase_1}_{phase_2}.png", dpi=300)
    plt.close()
    df = pd.DataFrame({'P':p1['P'], 'T':boundary_temp, 'crt_GFE': boundary, '1st_phase':first_phase, '2nd_phase':second_phase})
    df.to_csv(f"phase_diagram_{phase_1}_{phase_2}.csv", index=False)


def conti_phase_diagram(phase_1='BCC', phase_2='HCP', n_atom_1=1, n_atom_2=2):
    T = range(0, 4010, 10)
    P = range(0, 60, 1)
    p1 = pd.read_csv(f'{phase_1}/QHA_manual_results/phase_data.csv')
    p2 = pd.read_csv(f'{phase_2}/QHA_manual_results/phase_data.csv')

    boundary_data = []
    for i, p in enumerate(P):
        p1_GFE = p1.iloc[p, 2:].to_numpy() / n_atom_1
        p2_GFE = p2.iloc[p, 2:].to_numpy() / n_atom_2
        boundary_data.append(np.abs(p1_GFE - p2_GFE))
        # plt.plot(T, p1_GFE, '-', label=phase_1, color='red')
        # plt.plot(T, p2_GFE, '-', label=phase_2, color='blue')
        # plt.legend()
        # plt.show()
        # exit()


    plt.pcolormesh(T, P, boundary_data, shading='auto', cmap='coolwarm')
    plt.colorbar(label=r'$\Delta$ G (eV/mol)')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Pressure (GPa)')
    plt.text(x=2200, y=14, s="BCC", color='white', rotation=0)
    plt.text(x=2200, y=2, s="HCP", color='white', rotation=0)
    plt.title(f'Phase Diagram using QHA\nBlueish boundary separating {phase_1} and {phase_2}')
    plt.grid(True, linestyle="--", alpha=0.3)
    plt.savefig(f"phase_diagram_{phase_1}_{phase_2}_GFE_trend.png", dpi=300)
    plt.show()


def final_merge():
    df = pd.read_csv('phase_diagram_BCC_HCP.csv')
    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.scatter(df['T'], df['crt_GFE'], c='r', label='y = sin(x)')
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('Gibbs Free Energy (eV/mol)')
    ax1.tick_params(axis='x', labelcolor='b')

    ax2 = ax1.twiny()
    ax2.scatter(df['P'], df['crt_GFE'], c='b', label='z = exp(x)')
    ax2.set_xlabel('Pressure (GPa)', color='r')
    ax2.tick_params(axis='x', labelcolor='r')
    plt.title('Values of Gibbs Free Energy at Phase Boundary')
    plt.grid(True)

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    # expt()
    # main('BCC', 'HCP', 1, 2)
    # conti_phase_diagram()
    final_merge()
