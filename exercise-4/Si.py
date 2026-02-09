import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def make_data(filename1, filename2):
    """
    Read energies and lattice constants from file then write new file with the 
    volume and energy data for plotting and curve fitting.
    
    Args:
        filename1: Path to the input file
        filename2: Path to the output file
    """
    lattice = np.empty(0)
    energies = np.empty(0)
    
    with open(filename1, 'r') as f:
        for line in f:
            # Split the line by whitespace
            items = line.split()
            
            if len(items) >= 2:
                # Extract the integer number from the filename
                filename_str = items[0]
                # Extract digits from the filename (e.g., "mt-data/100.out" -> 100)
                num_str = ''.join(filter(str.isdigit, filename_str))
                if num_str:
                    lattice = np.append(lattice, int(num_str))
                # Extract the second last item (output)
                energies = np.append(energies, float(items[-2]))

    # sort both arrays
    sorted_indices = np.argsort(lattice)
    lattice = lattice[sorted_indices] * 0.01             # convert to Bohrs
    energies = energies[sorted_indices]

    lattice = lattice * 0.529177210903                   # convert to Angstroms
    volumes = lattice**3                                 # volume of unit cell
                                                         # note: 8 atoms per unit cell
    # write data to file 2
    df = pd.DataFrame({
        "Volume": volumes,
        "Energy": energies
    })

    df.to_csv(filename2, index=False, header=True)

def fit_Murnaghan(filename1):
    """
    Read energies as volumes from file, then fit the Murnaghan equation of state to the data and
    write the parameters to a file for subsequent plotting.

    Arguments:
        filename1: Path to the input file
    """
    data = pd.read_csv(filename1)
    V = data["Volume"].values 
    E = data["Energy"].values

    B0_trial = 100 * 1e9 * 1e-30 / (13.6056931229947 * 1.60218e-19) / 8 # convert 100 GPa to Ry/Angstrom^3, via energy per atom
    p0 = [-19.19270380, 20*8, B0_trial, 4] # initial guess for the parameters (E0, V0, B0, B0_prime)

    # Fit Murnaghan EOS
    popt_murn, _ = curve_fit(murnaghan, V, E, p0=p0, maxfev=100000)
    df = pd.DataFrame([popt_murn], columns=["E0", "V0", "B0", "B0_prime"])
    
    df.to_csv("Murnaghan_params.txt", index=False, header=True) # units: Rydbergs and Angstroms

    plot_EOS(filename1, "Murnaghan_params.txt")


def murnaghan(V, E0, V0, B0, B0_prime):
    eta = (V / V0)**(1 / 3)
    return E0 + (B0 * V0 / B0_prime) * eta**3 * (eta**(- 3 * B0_prime) / (B0_prime - 1) + 1) - (B0 * V0 / (B0_prime - 1))

def fit_Birch(filename1):
    """
    Read energies as volumes from file, then fit the Birch-Murnaghan equation of state to the data and
    write the parameters to a file for subsequent plotting.

    Arguments:
        filename1: Path to the input file
    """
    data = pd.read_csv(filename1)
    V = data["Volume"].values
    E = data["Energy"].values

    B0_trial = 100 * 1e9 * 1e-30 / (13.6056931229947 * 1.60218e-19) / 8 # convert 100 GPa to Ry/Angstrom^3
    p0 = [-19.19270380, 20*8, B0_trial, 4] # initial guess for the parameters (E0, V0, B0, B0_prime)

    # Fit Birch-Murnaghan EOS
    popt_birch, _ = curve_fit(birch_murnaghan, V, E, p0=p0, maxfev=100000)
    df = pd.DataFrame([popt_birch], columns=["E0", "V0", "B0", "B0_prime"])

    df.to_csv("Birch_params.txt", index=False, header=True) # units: Rydbergs and Angstroms

    plot_EOS(filename1, "Birch_params.txt")


def birch_murnaghan(V, E0, V0, B0, B0_prime):
    eta = (V / V0)**(1 / 3)
    return E0 + (9 * B0 * V0 / 16) * ((eta**(- 2) - 1)**3 * B0_prime + (eta**(- 2) - 1)**2 * (6 - 4 * eta**(- 2)))

def fit_Vinet(filename1):
    """
    Read energies as volumes from file, then fit the Vinet equation of state to the data and
    write the parameters to a file for subsequent plotting.

    Arguments:
        filename1: Path to the input file
    """
    data = pd.read_csv(filename1)
    V = data["Volume"].values
    E = data["Energy"].values

    B0_trial = 100 * 1e9 * 1e-30 / (13.6056931229947 * 1.60218e-19) / 8 # convert 100 GPa to Ry/Angstrom^3
    p0 = [-19.19270380, 20*8, B0_trial, 4] # initial guess for the parameters (E0, V0, B0, B0_prime)

    # Fit Vinet EOS
    popt_vinet, _ = curve_fit(vinet, V, E, p0=p0, maxfev=100000)
    df = pd.DataFrame([popt_vinet], columns=["E0", "V0", "B0", "B0_prime"])
    
    df.to_csv("Vinet_params.txt", index=False, header=True) # units: Rydbergs and Angstroms

    plot_EOS(filename1, "Vinet_params.txt")


def vinet(V, E0, V0, B0, B0_prime):
    eta = (V / V0)**(1 / 3)
    return E0 + (4 * B0 * V0 / (B0_prime - 1)**2) + (2 * B0 * V0 / (B0_prime - 1)**2) * np.exp(1.5 * (B0_prime - 1) * (1 - eta)) * (3 * (B0_prime - 1) * (1 - eta) - 2)


def plot_EOS(filename1, filename2):
    """
    Plot one of the fitted equations of state along with the original data points.
    
    Arguments:
        filename1: Path to the file containing the original data points
        filename2: Path to the file containing the fitted parameters
    """

    data1 = pd.read_csv(filename1)
    V     = data1["Volume"].values
    E     = data1["Energy"].values

    data2    = pd.read_csv(filename2)
    E0       = data2["E0"].values[0]
    V0       = data2["V0"].values[0]
    B0       = data2["B0"].values[0]
    B0_prime = data2["B0_prime"].values[0]

    if "Murnaghan" in filename2:
        E_fit = murnaghan(V, E0, V0, B0, B0_prime)
        label = "Murnaghan EOS"
    elif "Birch" in filename2:
        E_fit = birch_murnaghan(V, E0, V0, B0, B0_prime)
        label = "Birch-Murnaghan EOS"
    elif "Vinet" in filename2:
        E_fit = vinet(V, E0, V0, B0, B0_prime)
        label = "Vinet EOS"

    plt.plot(V, E, 'o', label='Data points')
    plt.plot(V, E_fit, label=label)
    plt.xlabel(r'Volume [$\AA^3$]')
    plt.ylabel(r'Energy [Ry]')
    plt.legend()
    plt.show()

def convert(filename1, filename2, filename3):
    """convert units of fitted parameters to Angstroms and Rydbergs for volume and energy, GPa for bulk modulus"""
    data1 = pd.read_csv(filename1)
    E01      = data1["E0"].values[0]
    V01      = data1["V0"].values[0] / 8 # volume per atom
    B01      = data1["B0"].values[0] * 8 * (13.6056931229947 * 1.60218e-19) / 1e-30 / 1e9 # convert from Ry/Angstrom^3 to GPa
    B0_prime1= data1["B0_prime"].values[0]

    data2 = pd.read_csv(filename2)
    E02      = data2["E0"].values[0]
    V02      = data2["V0"].values[0] / 8 # volume per atom
    B02      = data2["B0"].values[0] * 8 * (13.6056931229947 * 1.60218e-19) / 1e-30 / 1e9 # convert from Ry/Angstrom^3 to GPa
    B0_prime2= data2["B0_prime"].values[0]

    data3 = pd.read_csv(filename3)
    E03      = data3["E0"].values[0]
    V03      = data3["V0"].values[0] / 8 # volume per atom
    B03      = data3["B0"].values[0] * 8 * (13.6056931229947 * 1.60218e-19) / 1e-30 / 1e9 # convert from Ry/Angstrom^3 to GPa
    B0_prime3= data3["B0_prime"].values[0]

    df = pd.DataFrame(index = ["Murnaghan", "Birch-Murnaghan", "Vinet"], columns=["E0 [Ry/atom]", "V0 [A^3/atom]", "B0 [GPa]", "B0'"])
    df.loc["Murnaghan"] = [E01, V01, B01, B0_prime1]
    df.loc["Birch-Murnaghan"] = [E02, V02, B02, B0_prime2]
    df.loc["Vinet"] = [E03, V03, B03, B0_prime3]
    df.to_csv("EOS_parameters_converted.txt", index=False, header=True)

def plot_4(filename1, filename2, filename3, filename4):
    """plot the original data and the three fitted EOS on 4 subplots"""
    data1 = pd.read_csv(filename1)
    V     = data1["Volume"].values
    E     = data1["Energy"].values

    data2    = pd.read_csv(filename2)
    E0_1       = data2["E0"].values[0]
    V0_1       = data2["V0"].values[0]
    B0_1       = data2["B0"].values[0]
    B0_prime_1 = data2["B0_prime"].values[0]

    data3    = pd.read_csv(filename3)
    E0_2       = data3["E0"].values[0]
    V0_2       = data3["V0"].values[0]
    B0_2       = data3["B0"].values[0]
    B0_prime_2 = data3["B0_prime"].values[0]

    data4    = pd.read_csv(filename4)
    E0_3       = data4["E0"].values[0]
    V0_3       = data4["V0"].values[0]
    B0_3       = data4["B0"].values[0]
    B0_prime_3 = data4["B0_prime"].values[0]

    E_1 = murnaghan(V, E0_1, V0_1, B0_1, B0_prime_1)
    E_2 = birch_murnaghan(V, E0_2, V0_2, B0_2, B0_prime_2)
    E_3 = vinet(V, E0_3, V0_3, B0_3, B0_prime_3)

    fig, ax = plt.subplots(2, 2, figsize=(12, 10), sharex='col', sharey='row')

    ax[0, 0].plot(V, E, '.', label='QE data')
    ax[0, 0].set_ylabel(r'Energy [Ry]', fontsize=16)
    ax[0, 0].legend(fontsize=12)

    ax[0, 1].plot(V, E, '.', label='QE data')
    ax[0, 1].plot(V, E_1, color='red', label='Murnaghan')
    ax[0, 1].legend(fontsize=12)

    ax[1, 0].plot(V, E, '.', label='QE data')
    ax[1, 0].plot(V, E_2, color='red',label='Birch-Murnaghan')
    ax[1, 0].set_xlabel(r'Volume [$\AA^3$]', fontsize=16)
    ax[1, 0].set_ylabel(r'Energy [Ry]', fontsize=16)
    ax[1, 0].legend(fontsize=12)

    ax[1, 1].plot(V, E, '.', label='QE data')
    ax[1, 1].plot(V, E_3, color='red', label='Vinet')
    ax[1, 1].set_xlabel(r'Volume [$\AA^3$]', fontsize=16)
    ax[1, 1].legend(fontsize=12)

    plt.tight_layout()
    plt.show()

def bulk_properties(filename):
    """Determine the bulk properties of silicon from the fitted parameters 
    of the three equations of state and write to file"""
    data = pd.read_csv(filename)
    V0_M = data["V0 [A^3/atom]"].values[0] * 1e-30 # convert from A^3 to m^3
    B0_M = data["B0 [GPa]"].values[0]      * 1e9 # convert from GPa to Pa
    V0_B = data["V0 [A^3/atom]"].values[1] * 1e-30 # convert from A^3 to m^3
    B0_B = data["B0 [GPa]"].values[1]      * 1e9 # convert from GPa to Pa
    V0_V = data["V0 [A^3/atom]"].values[2] * 1e-30 # convert from A^3 to m^3
    B0_V = data["B0 [GPa]"].values[2]      * 1e9 # convert from GPa to Pa

    h_cross = 1.054571726e-34 # Planck's constant in Js
    m_si = 28.0855 * 1.66053906660e-27 # mass of silicon atom in kg
    kB = 1.3806488e-23 # Boltzmann constant in J/K

    # Calculate the mass density of Si from each EOS
    rho_M = m_si / V0_M
    rho_B = m_si / V0_B
    rho_V = m_si / V0_V

    # Calculate the sound velocity in Si from each EOS
    cL_M = np.sqrt(B0_M / rho_M)
    cL_B = np.sqrt(B0_B / rho_B)
    cL_V = np.sqrt(B0_V / rho_V)

    # Calculate the Debye temperature of Si from each EOS
    kD_M = (6 * np.pi**2 / V0_M)**(1 / 3)
    kD_B = (6 * np.pi**2 / V0_B)**(1 / 3)
    kD_V = (6 * np.pi**2 / V0_V)**(1 / 3)

    theta_M = h_cross * cL_M * kD_M / kB
    theta_B = h_cross * cL_B * kD_B / kB
    theta_V = h_cross * cL_V * kD_V / kB

    # round results to 6 decimal places and write to file
    rho_M = round(rho_M, 6)
    rho_B = round(rho_B, 6)
    rho_V = round(rho_V, 6)
    cL_M = round(cL_M, 6)
    cL_B = round(cL_B, 6)
    cL_V = round(cL_V, 6)
    theta_M = round(theta_M, 6)
    theta_B = round(theta_B, 6)
    theta_V = round(theta_V, 6)

    df = pd.DataFrame(index=["Murnaghan", "Birch-Murnaghan", "Vinet"], columns=["Density [kg/m^3]", "Sound velocity [m/s]", "Debye temperature [K]"])
    df.loc["Murnaghan"] = [rho_M, cL_M, theta_M]
    df.loc["Birch-Murnaghan"] = [rho_B, cL_B, theta_B]
    df.loc["Vinet"] = [rho_V, cL_V, theta_V]
    df.to_csv("Si_bulk.txt", index=True, header=True)

if __name__ == "__main__":
    filename1 = "Si-RyA.txt"
    filename2 = "Murnaghan_params.txt"
    filename3 = "Birch_params.txt"
    filename4 = "Vinet_params.txt"
    filename5 = "EOS_parameters_converted.txt"
    bulk_properties(filename5)