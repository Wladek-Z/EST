import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

def phonon_dispersion(filename):
    """Plot the phonon dispersion from gp file."""
    # Read in the data
    data = np.loadtxt(filename)

    nbands = data.shape[1] - 1
    for band in range(nbands):
        plt.plot(data[:, 0], data[:, band + 1], linewidth=1, alpha=0.5, color='k')
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
    plt.show()

def phonon_dos(filename):
    """Plot the phonon DOS from dat file."""
    # Read in the data 
    freq, dos, pdos_1, pdos_2 = np.loadtxt(filename, unpack=True)

    # Plot the DOS
    plt.plot(freq, dos, c='k', lw=0.5, label='Total')
    plt.plot(freq, pdos_1, c='b', lw=0.5, label='Si (1)')
    plt.plot(freq, pdos_2, c='r', lw=0.5, label='si (2)')
    plt.xlabel(r'Frequency $\omega$ (cm$^{-1}$)')
    plt.ylabel(r'Phonon DOS (states/cm$^{-1}$/unit-cell)')
    plt.legend(frameon=False, loc='upper left')
    plt.xlim(freq[0], freq[-1])
    plt.show()

def dispersion_dos(bandfile, dosfile):
    plt.rcParams["figure.dpi"] = 150
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["figure.figsize"] = (8, 6)

    # ==========================
    # ---- BAND STRUCTURE ----
    # ==========================

    band_data = np.loadtxt(bandfile)
    nbands = band_data.shape[1] - 1

    # ==========================
    # -------- DOS ------------
    # ==========================

    freq, dos, pdos_1, pdos_2 = np.loadtxt(dosfile, unpack=True)

    # ==========================
    # -------- PLOTTING -------
    # ==========================

    fig = plt.figure()
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1], wspace=0.05)

    ax_band = plt.subplot(gs[0])
    ax_dos  = plt.subplot(gs[1], sharey=ax_band)
    
    # Remove duplicate y-axis labels
    plt.setp(ax_dos.get_yticklabels(), visible=False)

    # ---- Plot bands ----
    for band in range(nbands):
        ax_band.plot(band_data[:, 0], band_data[:, band + 1], linewidth=1, alpha=0.5, color='k')

    # High symmetry k-points from disp.in file
    ax_band.axvline(x=band_data[0,0], linewidth=0.5, color='k', alpha=0.5)
    ax_band.axvline(x=band_data[50,0], linewidth=0.5, color='k', alpha=0.5)
    ax_band.axvline(x=band_data[100,0], linewidth=0.5, color='k', alpha=0.5)
    ax_band.axvline(x=band_data[150,0], linewidth=0.5, color='k', alpha=0.5)
    ax_band.axvline(x=band_data[200,0], linewidth=0.5, color='k', alpha=0.5)
    ax_band.axvline(x=band_data[250,0], linewidth=0.5, color='k', alpha=0.5)

    ax_band.set_xticks([0, band_data[50, 0], band_data[100, 0], band_data[150, 0], band_data[200, 0], band_data[250, 0]])
    ax_band.set_xticklabels([r'$\Gamma$', 'X', 'W', 'K', r'$\Gamma$', 'L'])
    ax_band.set_ylabel(r'Frequency $\omega$ (cm$^{-1}$)')
    ax_band.set_xlabel(r'Wave vector $k$')
    ax_band.set_xlim(band_data[0, 0], band_data[-1, 0])
    ax_band.set_ylim(0, freq[-1])

    # ---- Plot DOS ----
    ax_dos.plot(dos, freq, c='k', lw=0.5, label='Total')
    ax_dos.plot(pdos_1, freq, c='b', lw=0.5, label='Si (1)')
    ax_dos.plot(pdos_2, freq, c='r', lw=0.5, label='si (2)')
    ax_dos.set_xlabel('Phonon DOS')
    ax_dos.set_xticks([0, 0.05, 0.1])
    ax_dos.set_xticklabels([0, 0.05, 0.1])
    ax_dos.legend(frameon=False, loc='lower right')

    plt.show()



if __name__ == "__main__":
    filename1 = "si.ph-disp.gp"
    filename2 = "si.ph-dos.dat"
    #phonon_dispersion(filename1)
    dispersion_dos(filename1, filename2)