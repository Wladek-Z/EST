import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib import gridspec

def band_dos(bandfile, dosfile, sym_points, labels):
    plt.rcParams["figure.dpi"] = 150
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["figure.figsize"] = (8, 6)

    # ==========================
    # ---- BAND STRUCTURE ----
    # ==========================

    band_data = np.loadtxt(bandfile)

    k = np.unique(band_data[:,0])
    nk = len(k)

    total_points = len(band_data)
    nbands = total_points // nk

    bands = band_data[:,1][:nbands*nk].reshape((nbands, nk))

    # ==========================
    # -------- DOS ------------
    # ==========================

    filename = dosfile

    with open(filename, "r") as f:
        first_line = f.readline()

    EF = float(re.search(r'EFermi\s*=\s*([0-9.]+)', first_line).group(1))

    dos_data = np.loadtxt(filename, comments="#")
    E = dos_data[:,0]
    DOS = dos_data[:,1]

    # ==========================
    # -------- PLOTTING -------
    # ==========================

    fig = plt.figure()
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1], wspace=0.05)

    ax_band = plt.subplot(gs[0])
    ax_dos  = plt.subplot(gs[1], sharey=ax_band)

    # ---- Plot bands ----
    for band in range(nbands):
        ax_band.plot(k, bands[band,:], linewidth=1, alpha=0.75, color='midnightblue')

    # Fermi energy line
    ax_band.axhline(EF, linestyle=(0, (5, 5)), linewidth=0.75, color='k')
    
    # Add EF to y-axis ticks, moving any ticks that are too close
    yticks = list(ax_band.get_yticks())
    threshold = 1.0  # If a tick is within 1 unit of EF, move it
    
    # Check and adjust ticks that are too close to EF
    for i, tick in enumerate(yticks):
        if abs(tick - EF) < threshold:
            # Move the tick away from EF by 2 units
            if tick < EF:
                yticks[i] = tick - 2
            else:
                yticks[i] = tick + 2
    
    yticks.append(EF)
    ax_band.set_yticks(sorted(yticks))
    yticklabels = [f'{tick:.0f}' if tick != EF else r'$E_F$' for tick in sorted(yticks)]
    ax_band.set_yticklabels(yticklabels)

    for x in sym_points:
        ax_band.axvline(x, linewidth=0.75, color='k', alpha=0.5)

    ax_band.set_xticks(sym_points)
    ax_band.set_xticklabels(labels)
    ax_band.set_ylabel("Energy [eV]")
    ax_band.set_xlim(min(k), max(k))
    ax_band.set_xlabel("Wave vector")


    # ---- Plot DOS ----
    ax_dos.plot(DOS, E, color='midnightblue', linewidth=1)
    ax_dos.axhline(EF, linestyle=(0, (5, 5)), linewidth=0.75, color='k')

    ax_dos.set_xlabel("DOS")
    ax_dos.set_xlim(left=0)

    # Remove duplicate y-axis labels
    plt.setp(ax_dos.get_yticklabels(), visible=False)

    #plt.suptitle("Quartz Band Structure and DOS")
    plt.show()

if __name__ == "__main__":
    bandfile_quartz = 'Q-bands/quartz-bands.dat.gnu'
    dosfile_quartz = 'Q8.dat'
    sym_points_quartz = [0, 0.5774, 0.9107, 1.5773, 2.0319, 2.6092, 2.9426, 3.7494]
    labels_quartz = [r'$\Gamma$', 'M', 'K', r'$\Gamma$', 'A', 'L', 'H', r'$\Gamma$']

    bandfile_alpha = 'ASn-bands-nscf/alpha-bands.dat.gnu'
    dosfile_alphs = 'ASn14.dat'
    sym_points_alpha = [0, 1, 1.5, 1.8536, 2.9142, 3.7802]
    labels_alpha = [r'$\Gamma$', 'X', 'W', 'K', r'$\Gamma$', 'L']


    bandfile_beta = 'BSn-bands-nscf/beta-bands.dat.gnu'
    dosfile_beta = 'BSn14.dat'
    sym_points_beta = [0, 1.0439, 2.4896, 4.8046, 5.4534, 6.6826, 7.1826, 8.6862, 10.8726]
    labels_beta = [r'$\Gamma$', 'X', 'M', r'$\Gamma$', 'Z', 'N', 'P', r'$Z_1$', r'$\Gamma$']


    #band_dos(bandfile_quartz, dosfile_quartz, sym_points_quartz, labels_quartz)
    band_dos(bandfile_alpha, dosfile_alphs, sym_points_alpha, labels_alpha)
    #band_dos(bandfile_beta, dosfile_beta, sym_points_beta, labels_beta)