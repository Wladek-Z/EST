import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import fsolve

def get_data(filename1):
    """
    Read x and y data from file and return y
    
    Args:
        filename1: Path to the data file

    Returns:
        y: Array of values corresponding to some inputs
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
    y = y[sorted_indices]

    return y

def collect_data(filename1, filename2, output):
    """
    combine the energy and volume data into a single file for plotting and curve fitting.
    
    Arguments:
        filename1: path to energy data file
        filename2: path to volume data file
        output: name of output file
    """

    y1 = get_data(filename1)
    y2 = get_data(filename2)

    # convert volume data from au^3 (bohr^3) to angstrom^3
    y2 = y2 * (0.529177210903)**3

    df = pd.DataFrame({
        "Volume": y2,
        "Energy": y1
    })

    df.to_csv(output, index=False, header=True)

def vinet(V, E0, V0, B0, B0_prime):
    eta = (V / V0)**(1 / 3)
    return E0 + (4 * B0 * V0 / (B0_prime - 1)**2) + (2 * B0 * V0 / (B0_prime - 1)**2) * np.exp(1.5 * (B0_prime - 1) * (1 - eta)) * (3 * (B0_prime - 1) * (1 - eta) - 2)


def fit_Vinet(filename1, filename2, molecules):
    """
    Read energies as volumes from file, then fit the Vinet equation of state to the data and
    write the parameters to a file for subsequent plotting.

    Arguments:
        filename1: Path to the input file
        filename2: Path to the output file where fitted parameters will be saved
        molecules: number of molecules per unit cell
    """
    data = pd.read_csv(filename1)
    V = data["Volume"].values / molecules
    E = data["Energy"].values / molecules

    B0_trial = 100 * 1e9 * 1e-30 / (13.6056931229947 * 1.60218e-19) / 2 # convert 100 GPa to Ry/Angstrom^3
    p0 = [-19.19270380, 40, B0_trial, 4] # initial guess for the parameters (E0, V0, B0, B0_prime)

    # Fit Vinet EOS
    popt, _ = curve_fit(vinet, V, E, p0=p0, maxfev=100000)
    df = pd.DataFrame([popt], columns=["E0", "V0", "B0", "B0_prime"])

    # Plot the fit against the data
    V_fit = np.linspace(V[0], V[-1], 100)
    E_fit = vinet(V_fit, popt[0], popt[1], popt[2], popt[3])

    fig = plt.figure(figsize=[8, 6])
    plt.plot(V, E, 'o', label='data')
    plt.plot(V_fit, E_fit, color='red', label='Vinet')
    plt.xlabel(r"Volume [$\AA^3$/molecule]", fontsize=12)
    plt.ylabel("Energy [Ry/molecule]", fontsize=12)
    plt.legend()
    plt.tight_layout()
    plt.show()

    
    df.to_csv(filename2, index=False, header=True) # units: Rydbergs and Angstroms

def plot(filename, molecules):
    data = pd.read_csv(filename)
    V = data["Volume"].values / molecules
    E = data["Energy"].values / molecules
    fig = plt.figure(figsize=[8, 6])
    plt.plot(V, E, 'o')
    plt.xlabel(r"Volume [$\AA^3$/molecule]", fontsize=12)
    plt.ylabel("Energy [Ry/molecule]", fontsize=12)
    plt.show()

class Phase:
    """Class for calculating the phase transition pressure"""

    def __init__(self, params1, params2):
        """
        read in the vinet parameters
        """
        data1 = pd.read_csv(params1)
        self.E01, self.V01, self.B01, self.B0_prime1 = data1.iloc[0]

        data2 = pd.read_csv(params2)
        self.E02, self.V02, self.B02, self.B0_prime2 = data2.iloc[0]

    def common_tangent(self, filename1, filename2, molecules1, molecules2):
        """
        Read the Vinet parameters from file, then calculate the common tangent to the two curves and plot.

        Arguments:
            filname1: Path to the raw volume-energy data for the first phase
            filename2: Path to the raw volume-energy data for the second phase
            molecules1: number of molecules in the first phase
            molecules2: number of molecules in the second phase
        """
        data1 = pd.read_csv(filename1)
        V1 = data1["Volume"].values / molecules1
        E1 = data1["Energy"].values / molecules1

        data2 = pd.read_csv(filename2)
        V2 = data2["Volume"].values / molecules2
        E2 = data2["Energy"].values / molecules2

        # Calculate the Vinet fits
        E1_fit = vinet(V1, self.E01, self.V01, self.B01, self.B0_prime1)
        E2_fit = vinet(V2, self.E02, self.V02, self.B02, self.B0_prime2)

        # Find the common tangent of the data between the two fits
        initial_guess = [24, 19, -0.01]
        V1_t, V2_t, P = fsolve(self._common_tangent_equations, initial_guess)

        # Compute energy at transition volume
        E1_t = vinet(V1_t, self.E01, self.V01, self.B01, self.B0_prime1)

        # gradient and intercept
        m = -P
        c = E1_t - m * V1_t

        # create x data for common tangent line
        x = np.linspace(min(V1_t, V2_t)-4, max(V1_t, V2_t)+4, 4)

        # common tangent line
        tangent = m * x + c

        # Convert from Rydbergs/Angstrom^3 to GPa
        P_GPa = P * (13.6056931229947 * 1.60218e-19) / 1e-30 / 1e9
        # print results
        print(f"Pressure = {P_GPa:.6f} GPa")
        print(f"V1 = {V1_t:.6f} Angstrom^3/molecule, V2 = {V2_t:.6f} Angstrom^3/molecule")
        print(f"Volume difference = {((V1_t - V2_t)/V1_t * 100):.6f}%")

        # Plot the data, the fits, and the common tangent
        plt.figure(figsize=(8, 6))
        plt.plot(V1, E1, 'o', label='ice-Ih data', color='cyan')
        plt.plot(V2, E2, 'o', label='ice-VIII data', color='lime')
        plt.plot(V1, E1_fit, label='ice-Ih fit', color='blue')
        plt.plot(V2, E2_fit, label='ice-VIII fit', color='green')
        plt.plot(x, tangent, label='common tangent', linestyle='--', color='black')

        plt.xlabel(r"Volume [$\AA^3$/molecule]", fontsize=12)
        plt.ylabel("Energy [Ry/molecule]", fontsize=12)
        plt.legend()
        plt.show()


    def _common_tangent_equations(self, vars):
        """
        Function to minimise for points through which the common tangent passes and the gradient (-pressure) of this line
        """
        V1, V2, P = vars
        # P = - (E1 - E2) / (V1  - V2)
        eq1 = self.derivative(V1, self.E01, self.V01, self.B01, self.B0_prime1) + P
        eq2 = self.derivative(V2, self.E02, self.V02, self.B02, self.B0_prime2) + P
        # Make eq3 symmetric: slope should be the same regardless of order
        eq3 = (vinet(V1, self.E01, self.V01, self.B01, self.B0_prime1) - vinet(V2, self.E02, self.V02, self.B02, self.B0_prime2)) / (V1 - V2) + P
        return [eq1, eq2, eq3]

    def derivative(self, V, E0, V0, B0, B0_prime):
        """
        Calculate the numerical derivative of the Vinet equation of state with respect to volume.
        """
        h = 1e-5
        r = V + h
        l = V - h
        return (vinet(r, E0, V0, B0, B0_prime) - vinet(l, E0, V0, B0, B0_prime)) / (2 * h)
        

    
if __name__ == "__main__":
    Ih_VIII = Phase("iceIh-params.txt", "iceVIII-params.txt")
    Ih_VIII.common_tangent("iceIh-data.txt", "iceVIII-data.txt", 8, 8)
