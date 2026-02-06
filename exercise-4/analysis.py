import numpy as np
import matplotlib.pyplot as plt

def read_energies(filename, pos=-2):
    """
    Read energies from a file where each row contains some input (first item)
    and an output (second last item).
    
    Args:
        filename: Path to the input file
        pos: Position of the output value in the line (negative index)
    
    Returns:
        tuple: (numbers, outputs) - two lists containing the extracted values
    """
    numbers = np.empty(0)
    outputs = np.empty(0)
    
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
                    numbers = np.append(numbers, int(num_str))
                # Extract the second last item (output)
                outputs = np.append(outputs, float(items[pos]))

    # sort both arrays
    sorted_indices = np.argsort(numbers)
    numbers = numbers[sorted_indices]
    outputs = outputs[sorted_indices]
    
    return numbers, outputs


def plot_convergence(filename, pos=-2, name='alpha'):
    """
    Plot the change in total energy with respect to the last value as a function of
    cutoff energy to find convergence.

    Args:
        filename: file containing the energies data
        pos: Position of the output value in the line of the input file (negative index)
    """
    x, y = read_energies(filename, pos=pos)

    # difference in total energy with respect to the last value
    dy = np.abs(y - y[-1])

    # find meV tolerance in Ry 
    tol_10meV = 10 / 13605.6931229947  # 10 meV in Ry
    tol_1meV = 1 / 13605.6931229947    # 1 meV in Ry
    tol_5meV = 5 / 13605.6931229947    # 5 meV in Ry

    plt.plot(x, dy, marker='o', color='b', label=fr'$\{name}$-Sn')
    plt.axhline(y=tol_5meV, color='r', linestyle='--', label='5 meV')
    #plt.axhline(y=tol_1meV, color='r', linestyle='--', label='1 meV')
    plt.yscale('log')
    plt.xlabel(r'$E_c$ [Ry]')
    plt.ylabel(r'$\Delta E$ [Ry]')
    #plt.ylabel('k-points used')
    plt.legend()
    #plt.xticks(x[::4])
    plt.grid()
    plt.show()

def k_plot(filename1, filename2):
    """
    Plot the number of k-points used against the k-point grid size
    
    Arguments:
        filename1: file containing k-point data for quartz
        filename2: file containing k-point data for alpha-tin
    """
    x1, y1 = read_energies(filename1, pos=-1)
    x2, y2 = read_energies(filename2, pos=-1)
    y_3 = (x1**3)

    plt.loglog(x1, y1, color='b', label=r'SiO$_2$')
    plt.loglog(x2, y2, color='orange', label=r'$\alpha$-Sn')
    plt.loglog(x1, y_3, linestyle='--', color='r', label=r'$N^3$')
    plt.xlabel('k-point grid size')
    plt.ylabel('k-points used')
    #plt.grid()
    plt.legend()
    plt.show()

def energy_diff(file1, file2, pos=-2):
    """
    Compute the difference in total energy between two datasets.

    Args:
        file1: first energy data file
        file2: second energy data file
        pos: Position of the output value in the line of the input file (negative index)
    """
    Ec1, Et1 = read_energies(file1, pos=pos)
    Ec2, Et2 = read_energies(file2, pos=pos)

    # Ensure both datasets have the same cutoff energies for comparison
    if not np.array_equal(Ec1, Ec2):
        raise ValueError("Cutoff energies do not match between the two datasets.")

    dE1 = np.abs(Et1 - Et2)
    dE = np.abs(dE1 - dE1[-1])
    tol_1meV = 1 / 13605.6931229947    # 1 meV in Ry
    
    plt.plot(Ec1, dE, marker='o', color='b', label=r'$\Delta E$')
    plt.axhline(y=dE[-1]+tol_1meV, color='r', linestyle='--', label='1 meV')
    plt.legend()
    plt.yscale('log')
    plt.xlabel(r'$E_c$ [Ry]')
    plt.ylabel(r'$|E_{\alpha} - E_{\beta}|$ [Ry]')
    #plt.xticks(Ec1[::4])
    plt.grid()
    plt.show()


if __name__ == "__main__":
    filename1 = "A-energies.txt"
    filename2 = "B-energies.txt"
    plot_convergence(filename1, name='alpha')
    plot_convergence(filename2, name='beta')
    energy_diff(filename1, filename2)
    