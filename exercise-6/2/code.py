import numpy as np
import matplotlib.pyplot as plt

def phonon_dispersion(filename):
    """Plot the phonon dispersion from gp file. with speeds of each branch."""
    # Calculate velocities from gradients
    mt = 363.045 #373.995
    ml = 527.17 #542.785
    # Perform unit conversion
    vt = mt * 100 * 2.99792458e8 * 0.529177210544 * 10.2 / 10**10
    vl = ml * 100 * 2.99792458e8 * 0.529177210544 * 10.2 / 10**10
    print(rf"Transverse velocities = {vt:.6f} m/s")
    print(rf"Longitudinal velocity = {vl:.6f} m/s")


    # Read in the data
    data = np.loadtxt(filename)

    nbands = data.shape[1] - 1
    for band in range(nbands):
        plt.plot(data[:, 0], data[:, band + 1], linewidth=1, alpha=0.5, color='k')

    # Plot gradients near 0
    x = np.linspace(0, 1, 4)
    yt = mt * x
    yl = ml * x
    plt.plot(x, yt, linestyle="--", color='b', label=f"Transverse $v_s$\ngradient = {mt}")
    plt.plot(x, yl, linestyle="--", color='r', label=f"Longitudinal $v_s$\ngradient = {ml}")
    # High symmetry k-points from disp.in file
    plt.axvline(x=data[0,0], linewidth=0.5, color='k', alpha=0.5)
    plt.axvline(x=data[50,0], linewidth=0.5, color='k', alpha=0.5)
    plt.axvline(x=data[100,0], linewidth=0.5, color='k', alpha=0.5)
    plt.axvline(x=data[150,0], linewidth=0.5, color='k', alpha=0.5)
    plt.axvline(x=data[200,0], linewidth=0.5, color='k', alpha=0.5)
    plt.axvline(x=data[250,0], linewidth=0.5, color='k', alpha=0.5)
    
    plt.xticks(ticks=[0, data[50, 0], data[100, 0], data[150, 0], data[200, 0], data[250, 0]], \
               labels=[r'$\Gamma$', 'X', 'W', 'K', r'$\Gamma$', 'L'])
    plt.ylabel(r'Frequency $\omega$ (cm$^{-1}$)')
    plt.xlim(data[0, 0], data[-1, 0])
    plt.ylim(0, )
    plt.legend(loc="upper right")
    plt.show()
    

if __name__ == "__main__":
    phonon_dispersion("disp7.gp")