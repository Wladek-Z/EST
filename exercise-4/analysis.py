import numpy as np
import matplotlib.pyplot as plt

def read_energies(filename):
    """
    Read energies from a file where each row contains the cutoff energy (first item)
    and the total energy (second last item).
    
    Args:
        filename: Path to the input file
    
    Returns:
        tuple: (ecutwfc, energies) - two lists containing the extracted values
    """
    ecutwfc = np.empty(0)
    energies = np.empty(0)
    
    with open(filename, 'r') as f:
        for line in f:
            # Split the line by whitespace
            items = line.split()
            
            if len(items) >= 2:
                # Extract the integer number from the filename
                filename_str = items[0]
                # Extract digits from the filename (e.g., "mt-data/100.out" -> 100)
                num_str = ''.join(filter(str.isdigit, filename_str))
                if num_str:
                    ecutwfc = np.append(ecutwfc, int(num_str))
                # Extract the second last item (energy)
                energies = np.append(energies, float(items[-2]))

    # sort both arrays
    sorted_indices = np.argsort(ecutwfc)
    ecutwfc = ecutwfc[sorted_indices]
    energies = energies[sorted_indices]
    
    return ecutwfc, energies

def plot_convergence(filename):
    """
    Plot the change in total energy with respect to the last value as a function of
    cutoff energy to find convergence.

    Args:
        filename: file containing the energies data
    """
    Ec, Et = read_energies(filename)

    # difference in total energy with respect to the last value
    dE = np.abs(Et - Et[-1])

    # find meV tolerance in Ry 
    tol_10meV = 10 / 13605.6931229947  # 10 meV in Ry
    tol_1meV = 1 / 13605.6931229947    # 1 meV in Ry
    tol_5meV = 5 / 13605.6931229947    # 5 meV in Ry

    plt.plot(Ec, dE, marker='o', color='b', label=r'$\alpha$-Sn')
    plt.axhline(y=tol_5meV, color='r', linestyle='--', label='5 meV')
    #plt.axhline(y=tol_1meV, color='g', linestyle='--', label='1 meV')
    plt.yscale('log')
    plt.xlabel(r'$E_c$ [Ry]')
    plt.ylabel(r'$\Delta E$ [Ry]')
    plt.legend()
    plt.xticks(Ec[::2])
    plt.grid()
    plt.show()

def energy_diff(file1, file2):
    """
    Compute the difference in total energy between two datasets.

    Args:
        file1: first energy data file
        file2: second energy data file
    """
    Ec1, Et1 = read_energies(file1)
    Ec2, Et2 = read_energies(file2)

    # Ensure both datasets have the same cutoff energies for comparison
    if not np.array_equal(Ec1, Ec2):
        raise ValueError("Cutoff energies do not match between the two datasets.")

    dE = Et1 - Et2
    tol_5meV = 5 / 13605.6931229947    # 5 meV in Ry
    
    plt.plot(Ec1, dE, marker='o', color='b', label=r'$\Delta E$')
    plt.axhline(y=tol_5meV, color='r', linestyle='--', label='5 meV')
    plt.legend()
    plt.xlabel(r'$E_c$ [Ry]')
    plt.ylabel(r'$E_{\alpha} - E_{\beta}$ [Ry]')
    plt.xticks(Ec1[::2])
    plt.grid()
    plt.show()



if __name__ == "__main__":
    filename1 = "a-energies.txt"
    filename2 = "b-energies.txt"
    plot_convergence(filename1)