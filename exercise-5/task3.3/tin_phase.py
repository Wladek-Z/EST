import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.optimize import fsolve

def get_data(filename1, filename2, filename3):
    """
    Read x and y data from file and return y
    
    Args:
        filename1: Path to the energy data file 
        filename2: Path to the volume data file
        filename3: Path to the output file where processed data will be saved

    Returns:
        y: Array of values corresponding to some inputs
    """
    x1 = np.empty(0)
    y1 = np.empty(0)
    
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
                    x1 = np.append(x1, int(num_str))
                # Extract the second last item (output)
                y1 = np.append(y1, float(items[-2]))

    # sort y array according to x values
    sorted_indices = np.argsort(x1)
    y1 = y1[sorted_indices]

    x2 = np.empty(0)
    y2 = np.empty(0)
    
    with open(filename2, 'r') as f:
        for line in f:
            # Split the line by whitespace
            items = line.split()
            
            if len(items) >= 2:
                # Extract the integer number from the filename
                filename_str = items[0]
                # Extract digits from the filename (e.g., "mt-data/100.out" -> 100)
                num_str = ''.join(filter(str.isdigit, filename_str))
                if num_str:
                    x2 = np.append(x2, int(num_str))
                # Extract the second last item (output)
                y2 = np.append(y2, float(items[-2]))

    
    sorted_indices = np.argsort(x2)
    y2 = y2[sorted_indices]
    # Note: energy in Ry, volume in Angstrom^3

    # write data to file 3
    df = pd.DataFrame({
        "Volume": y2 * (0.529177210903)**3, # convert from Bohr^3 to Angstrom^3
        "Energy": y1
    })

    df.to_csv(filename3, index=False, header=True)


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

    B0_trial = 100 * 1e9 * 1e-30 / (13.6056931229947 * 1.60218e-19) / 2 # convert 100 GPa to Ry/Angstrom^3
    p0 = [-19.19270380, 40, B0_trial, 4] # initial guess for the parameters (E0, V0, B0, B0_prime)

    # Fit Vinet EOS
    popt_vinet, _ = curve_fit(vinet, V, E, p0=p0, maxfev=100000)
    df = pd.DataFrame([popt_vinet], columns=["E0", "V0", "B0", "B0_prime"])
    
    df.to_csv("B-params.txt", index=False, header=True) # units: Rydbergs and Angstroms

def plot(filename):
    data = pd.read_csv(filename)
    V = data["Volume"].values
    E = data["Energy"].values
    plt.plot(V, E, 'o')
    plt.xlabel(r"Volume [$\AA^3$]", fontsize=12)
    plt.ylabel("Energy [Ry]")
    plt.show()


class Phase:

    def __init__(self, filenameA, filenameB):
        """
        read in the vinet parameters
        """
        dataA = pd.read_csv(filenameA)
        self.E0A, self.V0A, self.B0A, self.B0_primeA = dataA.iloc[0]

        dataB = pd.read_csv(filenameB)
        self.E0B, self.V0B, self.B0B, self.B0_primeB = dataB.iloc[0]

    def common_tangent(self, filename1, filename2):
        """
        Read the Vinet parameters from file, then calculate the common tangent to the two curves and plot.

        Arguments:
            filname1: Path to the raw volume-energy data for the first phase
            filename2: Path to the raw volume-energy data for the second phase
        """
        dataA = pd.read_csv(filename1)
        VA = dataA["Volume"].values
        EA = dataA["Energy"].values

        dataB = pd.read_csv(filename2)
        VB = dataB["Volume"].values
        EB = dataB["Energy"].values

        # Calculate the Vinet fits
        EA_fit = vinet(VA, self.E0A, self.V0A, self.B0A, self.B0_primeA)
        EB_fit = vinet(VB, self.E0B, self.V0B, self.B0B, self.B0_primeB)

        # Find the common tangent of the data between the two fits

        initial_guess = [50, 60, -0.01]
        VA_t, VB_t, P = fsolve(self._common_tangent_equations, initial_guess)

        # Compute energies at transition volumes
        EA_t = vinet(VA_t, self.E0A, self.V0A, self.B0A, self.B0_primeA)
        EB_t = vinet(VB_t, self.E0B, self.V0B, self.B0B, self.B0_primeB)

        # gradient
        m = -P

        # intercept
        c = EA_t - m * VA_t

        # create x data for common tangent line
        x = np.linspace(min(VA_t, VB_t)-10, max(VA_t, VB_t)+10, 200)

        # common tangent line
        tangent = m * x + c

        # Convert from Rydbergs/Angstrom^3 to GPa
        P_GPa = P * (13.6056931229947 * 1.60218e-19) / 1e-30 / 1e9
        # print results
        print(f"Pressure = {P_GPa:.6f} GPa")
        print(f"VA = {VA_t:.6f} Angstrom^3, VB = {VB_t:.6f} Angstrom^3")
        print(f"Volume difference = {((VA_t - VB_t)/VA_t * 100):.6f}%")

        # Plot the data, the fits, and the common tangent
        plt.figure(figsize=(8, 6))
        plt.plot(VA, EA, 'o', label=r'$\alpha$-Sn data', color='cyan')
        plt.plot(VB, EB, 'o', label=r'$\beta$-Sn data', color='lime')
        plt.plot(VA, EA_fit, label=r'$\alpha$-Sn fit', color='blue')
        plt.plot(VB, EB_fit, label=r'$\beta$-Sn fit', color='green')
        plt.plot(x, tangent, label='common tangent', linestyle='--', color='black')

        plt.xlabel(r"Volume [$\AA^3$]", fontsize=12)
        plt.ylabel("Energy [Ry]")
        plt.legend()
        plt.grid(True)
        plt.show()


    def _common_tangent_equations(self, vars):
        """
        Function to minimise for points through which the common tangent passes and the gradient (-pressure) of this line
        """
        AV, BV, P = vars
        # P = - (EA - EB) / VA  - VB
        eq1 = self.derivative(AV, self.E0A, self.V0A, self.B0A, self.B0_primeA) + P
        eq2 = self.derivative(BV, self.E0B, self.V0B, self.B0B, self.B0_primeB) + P
        # Make eq3 symmetric: slope should be the same regardless of order
        eq3 = (vinet(AV, self.E0A, self.V0A, self.B0A, self.B0_primeA) - vinet(BV, self.E0B, self.V0B, self.B0B, self.B0_primeB)) / (AV - BV) + P
        return [eq1, eq2, eq3]

    def derivative(self, V, E0, V0, B0, B0_prime):
        """
        Calculate the numerical derivative of the Vinet equation of state with respect to volume.
        """
        h = 1e-5
        r = V + h
        l = V - h
        return (vinet(r, E0, V0, B0, B0_prime) - vinet(l, E0, V0, B0, B0_prime)) / (2 * h)
        

def vinet(V, E0, V0, B0, B0_prime):
    eta = (V / V0)**(1 / 3)
    return E0 + (4 * B0 * V0 / (B0_prime - 1)**2) + (2 * B0 * V0 / (B0_prime - 1)**2) * np.exp(1.5 * (B0_prime - 1) * (1 - eta)) * (3 * (B0_prime - 1) * (1 - eta) - 2)


if __name__ == "__main__":
    phase = Phase("A-params.txt", "B-params.txt")
    phase.common_tangent("A-data.txt", "B-data.txt")