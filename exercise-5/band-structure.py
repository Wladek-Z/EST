import matplotlib.pyplot as plt
from matplotlib import rcParamsDefault
import numpy as np
#matplotlib inline

def alpha_bands():
    plt.rcParams["figure.dpi"]=150
    plt.rcParams["figure.facecolor"]="white"
    plt.rcParams["figure.figsize"]=(8, 6)

    # load data
    data = np.loadtxt('ASn-bands/alpha-bands.dat.gnu')


    k = np.unique(data[:,0])
    nk = len(k)

    total_points = len(data)
    nbands = total_points // nk  # integer division

    print("nk =", nk)
    print("nbands =", nbands)

    bands = data[:,1][:nbands*nk].reshape((nbands, nk))


    for band in range(len(bands)):
        plt.plot(k, bands[band, :], linewidth=1, alpha=0.5, color='k')
    plt.xlim(min(k), max(k))

    # Fermi energy
    plt.axhline(10.812, linestyle=(0, (5, 5)), linewidth=0.75, color='k', alpha=0.5)
    # High symmetry k-points (check bands_pp.out)0.0000
    plt.axvline(0.0000, linewidth=0.75, color='k', alpha=0.5)
    plt.axvline(1.0000, linewidth=0.75, color='k', alpha=0.5)
    plt.axvline(1.5000, linewidth=0.75, color='k', alpha=0.5)
    plt.axvline(1.8536, linewidth=0.75, color='k', alpha=0.5)
    plt.axvline(2.9142, linewidth=0.75, color='k', alpha=0.5)
    plt.axvline(3.7802, linewidth=0.75, color='k', alpha=0.5)
    #plt.axvline(0.8660, linewidth=0.75, color='k', alpha=0.5)
    #plt.axvline(1.8660, linewidth=0.75, color='k', alpha=0.5)
    #plt.axvline(2.2196, linewidth=0.75, color='k', alpha=0.5)
    # text labels
    plt.xticks(ticks= [0, 1.0, 1.5,1.8536, 2.9142, 3.7802], \
            labels=['$\Gamma$', 'X', 'W','K', '$\Gamma$','L'])
    plt.ylabel("Energy (eV)")
    plt.text(2.3, 5.6, 'Fermi energy')
    plt.show()


def quartz_bands():
    ef = 2.556

    plt.rcParams["figure.dpi"]=150
    plt.rcParams["figure.facecolor"]="white"
    plt.rcParams["figure.figsize"]=(8, 6)

    # load data
    data = np.loadtxt('Q-bands/quartz-bands.dat.gnu')


    k = np.unique(data[:,0])
    nk = len(k)

    total_points = len(data)
    nbands = total_points // nk  # integer division

    print("nk =", nk)
    print("nbands =", nbands)

    bands = data[:,1][:nbands*nk].reshape((nbands, nk))


    for band in range(len(bands)):
        plt.plot(k, bands[band, :], linewidth=1, alpha=0.75, color='midnightblue')
    plt.xlim(min(k), max(k))

    # Fermi energy
    plt.axhline(ef, linestyle=(0, (5, 5)), linewidth=0.75, color='k')
    # High symmetry k-points (check bands_pp.out)0.0000
    plt.axvline(0.0000, linewidth=0.75, color='k')
    plt.axvline(0.5774, linewidth=0.75, color='k')
    plt.axvline(0.9107, linewidth=0.75, color='k')
    plt.axvline(1.5773, linewidth=0.75, color='k')
    plt.axvline(2.0319, linewidth=0.75, color='k')
    plt.axvline(2.6092, linewidth=0.75, color='k')
    plt.axvline(2.9426, linewidth=0.75, color='k')
    plt.axvline(3.7494, linewidth=0.75, color='k')
    # text labels
    plt.xticks(ticks= [0, 0.5774, 0.9107,1.5773, 2.0319, 2.6092, 2.9426, 3.7494], \
            labels=[r'$\Gamma$', 'M', 'K', r'$\Gamma$', 'A', 'L', 'H', r'$\Gamma$'])
    plt.ylabel("Energy [eV]")
    #plt.text(2.3, 5.6, 'Fermi energy')
    ax = plt.gca()
    # Add Fermi energy to y-axis ticks
    yticks = list(ax.get_yticks())
    yticks.append(ef)
    ax.set_yticks(sorted(yticks))

    # Add 'EF' label for the Fermi energy tick
    yticklabels = [f'{tick:.0f}' if tick != ef else r'$E_F$' for tick in sorted(yticks)]
    ax.set_yticklabels(yticklabels)
    plt.show()

if __name__ == '__main__':
    quartz_bands()