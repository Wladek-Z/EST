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


def plot_convergence(filename, pos=-2):
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
    tol_5meV = 5 / 13605.6931229947    # 5 meV in Ry

    plt.plot(x, dy/12, marker='o', color='b', label='ice-II')
    plt.axhline(y=tol_5meV, color='r', linestyle='--', label='5 meV')
    plt.yscale('log')
    plt.xlabel(r'Energy cutoff $E_c$ [Ry]', fontsize=12)
    plt.ylabel(r'$\Delta E$ [Ry/molecule]', fontsize=12)
    #plt.xlabel(r'k-point grid size $(N_k, N_k, N_k)$', fontsize=12)
    plt.legend()
    #plt.xticks(x[::2])
    plt.grid()
    plt.show()

if __name__ == "__main__":
    filename = 'iceII-energies.txt'
    plot_convergence(filename)