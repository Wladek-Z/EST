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
    volumes = lattice**3

    energies = energies * 13.6056931229947 * 1.60218e-19 # convert to Joules
    
    # write data to file 2
    df = pd.DataFrame({
        "Volume": volumes,
        "Energy": energies
    })

    df.to_csv(filename2, index=False, header=True)

def curve_fit_data(filename1, filename2):
    """
    Read energies as volumes from file, then fit 3 different equations of state to the data and
    write the fitted data to new files for subsequent plotting.

    Arguments:
        filename1: Path to the input file
        filename2: Path to the output file containing each fit's parameters
    """
    data = pd.read_csv(filename1)
    V = data["Volume"].values
    E = data["Energy"].values

    V = V * 1e-30 # convert to m^3

    df = pd.DataFrame(index=["Murnaghan", "Birch-Murnaghan", "Vinet"], columns=["E0", "V0", "B0", "B0_prime"])
    
    p0 = [-19.19270234 * 13.6056931229947 * 1.60218e-19, 20**3 * 1e-30, 100e9, 4] # initial guess for the parameters (E0, V0, B0, B0_prime)

    # Fit Murnaghan EOS
    popt_murn, _ = curve_fit(murnaghan, V, E, p0=p0, maxfev=1000000)
    df.loc["Murnaghan"] = popt_murn

    # Fit Birch-Murnaghan EOS
    popt_birch, _ = curve_fit(birch_murnaghan, V, E, p0=p0, maxfev=1000000)
    df.loc["Birch-Murnaghan"] = popt_birch

    # Fit Vinet EOS
    popt_vinet, _ = curve_fit(vinet, V, E, p0=p0, maxfev=1000000)
    df.loc["Vinet"] = popt_vinet

    df.to_csv(filename2, index=True, header=True) # units: m^3 and J


def murnaghan(V, E0, V0, B0, B0_prime):
    eta = (V / V0)**(1 / 3)
    return E0 + (B0 * V / B0_prime) * eta**3 * (eta**(- 3 * B0_prime) / (B0_prime - 1) + 1) - (B0 * V0 / (B0_prime - 1))

def birch_murnaghan(V, E0, V0, B0, B0_prime):
    eta = (V / V0)**(1 / 3)
    return E0 + (9 * B0 * V0 / 16) * ((eta**(- 2) - 1)**3 * B0_prime + (eta**(- 2) - 1)**2 * (6 - 4 * eta**(- 2)))

def vinet(V, E0, V0, B0, B0_prime):
    eta = (V / V0)**(1 / 3)
    return E0 + (4 * B0 * V0 / (B0_prime - 1)**2) + (2 * B0 * V0 / (B0_prime - 1)**2) * np.exp(1.5 * (B0_prime - 1) * (1 - eta)) * (3 * (B0_prime - 1) * (1 - eta) - 2)


if __name__ == "__main__":
    curve_fit_data("Si_AJ.txt", "parameters2.txt")