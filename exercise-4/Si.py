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

    B0_trial = 100 * 1e9 * 1e-30 / (13.6056931229947 * 1.60218e-19) # convert 100 GPa to Ry/Angstrom^3
    p0 = [-19.19270234, 20*8, B0_trial, 4] # initial guess for the parameters (E0, V0, B0, B0_prime)

    # Fit Murnaghan EOS
    popt_murn, _ = curve_fit(murnaghan, V, E, p0=p0, maxfev=100000, bounds=([-100, 1, 0, 0], [0, 400, 1, 10]))
    df = pd.DataFrame([popt_murn], columns=["E0", "V0", "B0", "B0_prime"])
    
    df.to_csv("Murnaghan_params.txt", index=False, header=True) # units: Rydbergs and Angstroms

    plot_EOS(filename1, "Murnaghan_params.txt")


def murnaghan(V, E0, V0, B0, B0_prime):
    eta = (V / V0)**(1 / 3)
    return E0 + (B0 * V / B0_prime) * eta**3 * (eta**(- 3 * B0_prime) / (B0_prime - 1) + 1) - (B0 * V0 / (B0_prime - 1))

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

    B0_trial = 100 * 1e9 * 1e-30 / (13.6056931229947 * 1.60218e-19) # convert 100 GPa to Ry/Angstrom^3
    p0 = [-19.19270234, 20*8, B0_trial, 4] # initial guess for the parameters (E0, V0, B0, B0_prime)

    # Fit Birch-Murnaghan EOS
    popt_birch, _ = curve_fit(birch_murnaghan, V, E, p0=p0, maxfev=100000, bounds=([-100, 1, 0, 0], [0, 400, 1, 10]))
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
    E = data["Energy"].valuesbounds=([-100, 1, 0, 0], [0, 400, 1, 10])

    B0_trial = 100 * 1e9 * 1e-30 / (13.6056931229947 * 1.60218e-19) # convert 100 GPa to Ry/Angstrom^3
    p0 = [-19.19270234, 20*8, B0_trial, 4] # initial guess for the parameters (E0, V0, B0, B0_prime)

    # Fit Vinet EOS
    popt_vinet, _ = curve_fit(vinet, V, E, p0=p0, maxfev=100000, bounds=([-100, 1, 0, 0], [0, 400, 1, 10]))
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


if __name__ == "__main__":
    filename1 = "Si-RyA.txt"
    fit_Birch(filename1)