import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob

def get_data(filename1):
    """
    Read input and output data from file and return them as arrays.
    
    Args:
        filename: Path to the data file

    Returns:
        x: Array of input values
        y: Array of output values
    """
    x = np.empty(0)
    y = np.empty(0)
    
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
                    x = np.append(x, int(num_str))
                # Extract the second last item (output)
                y = np.append(y, float(items[-2]))

    # sort y array according to x values
    sorted_indices = np.argsort(x)
    x = x[sorted_indices]
    y = y[sorted_indices]

    return x, y

def read_stress(filename):
    """
    Read stress data from file and return it as an array.
    
    Args:
        filename: Path to the stress data file
    Returns:
        e33: array of strain rate values
        s11: array of stress values
        s33: array of stress values
    """
    data = pd.read_csv(filename)
    ca = data['ca10000'].values * 0.0001
    s11 = data['s11'].values 
    s33 = data['s33'].values 
    # calculate the actual strain rate
    e33 = ca - 1
    # sort the data
    sort = np.argsort(e33)
    e33 = e33[sort]
    s11 = s11[sort]
    s33 = s33[sort]
    # units: none, Ry/A^3, Ry/A^3
    return e33, s11, s33

def plot_stress(filename):
    """
    Plot stress data from file.
    
    Args:
        filename: Path to the stress data file
    """
    e33, s11, s33 = read_stress(filename)

    x = np.linspace(min(e33), max(e33), 100)
    line = lambda m, c: m * x + c
    m11, c11, m33, c33 = fit_stress(filename)

    print(f"m11: {m11:.6f}, m33: {m33:.6f}")

    y11 = line(m11, c11)
    y33 = line(m33, c33)
    
    plt.figure(figsize=(10, 6))
    plt.scatter(e33, s11, label=r'$\sigma_{11}$ data', color='cyan', s=20, marker='x')
    plt.plot(x, y11, label=r'$\sigma_{11}$ fit', color='b', linestyle='-', linewidth=2)
    plt.scatter(e33, s33, label=r'$\sigma_{33}$ data', color='lime', s=20, marker='x')
    plt.plot(x, y33, label=r'$\sigma_{33}$ fit', color='g', linestyle='-', linewidth=2)
    plt.xlabel(r'Strain Rate $\epsilon_{33}$', fontsize=12)
    plt.ylabel(r'Stress $\sigma_{ij}$ [Ry/Bohr$^{3}$]', fontsize=12)
    plt.legend()
    plt.show()

def fit_stress(filename):
    """
    Fit stress data from file to a linear model and return the fitted parameters.
    
    Args:
        filename: Path to the stress data file
    Returns:
        m11: Gradient of the s11 fitted line
        c11: Intercept of the s11 fitted line
        m33: Gradient of the s33 fitted line
        c33: Intercept of the s33 fitted line
    """
    e33, s11, s33 = read_stress(filename)
    
    # Fit a linear model to the data
    m11, c11 = np.polyfit(e33, s11, 1)
    m33, c33 = np.polyfit(e33, s33, 1)
    
    return m11, c11, m33, c33

def elsatic_constants(filename):
    """
    Calculate elastic constants from fitted stress data and write to file.
    
    Args:
        filename: Path to the stress data file
    """
    m11, _, m33, _ = fit_stress(filename)

    # convert from Ry/Bohr^3 to GPa
    c12 = -m11 / (0.529177249e-10)**3 * 2.1798741e-18 / 1e9
    c11 = -m33 / (0.529177249e-10)**3 * 2.1798741e-18 / 1e9
    
    with open("task 3.1.2/elastic_constants.txt", 'w') as f:
        f.write(f"c12: {c12:.6f} GPa\n")
        f.write(f"c11: {c11:.6f} GPa")

def dos(folder, fermifile):
    """
    Plot density of states for each k-point grid size data in the folder.
    
    Arguments:
        folder: path to the folder containing dos data files for varying k-point grid sizes
        fermifile: path to the fermi energy data file
    """
    # get files from folder and sort them by the integer number in the filename
    files = glob.glob(f"{folder}/*.dat")
    files.sort(key=lambda x: int(''.join(filter(str.isdigit, x))))
    files = files[::2]

    _, Ef = get_data(fermifile)
    Ef = Ef[::2]

    for i, f in enumerate(files):
        data = np.loadtxt(f, comments='#')
        E = data[:, 0] - Ef[i]  # Energy values
        D = data[:, 1]  # Density of states values
        plt.plot(E, D, label=f"N = {''.join(filter(str.isdigit, f.split('/')[-1].split('.')[0]))}")
    
    plt.xlabel(r"Energy ($E - E_F$) [eV])", fontsize=12)
    plt.ylabel("Density of States", fontsize=12)
    plt.axvline(x=0, color='r', linestyle='--')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    folder="task 3.2/ASn-tetrahedron-data"
    filename = "task 3.2/ASn-fermi.txt"
    dos(folder, filename)