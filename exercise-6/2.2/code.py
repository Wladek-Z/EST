import numpy as np
import matplotlib.pyplot as plt

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
    #plt.axvline(x=data[200,0], linewidth=0.5, color='k', alpha=0.5)
    
    plt.xticks(ticks=[0, data[50, 0], data[100, 0], data[150, 0]], \
               labels=[r'$\Gamma$', 'X', r'$\Gamma$', 'L'])
    plt.ylabel(r'Frequency $\omega$ (cm$^{-1}$)')
    plt.xlim(data[0, 0], data[-1, 0])
    plt.ylim(0, )
    plt.show()

def neutron(filename):
    """Plot the phonon dispersion relation against the neutron data"""
    # Read in DFT data
    data = np.loadtxt(filename)

    nbands = data.shape[1] - 1
    for band in range(nbands):
        plt.plot(data[:, 0], data[:, band + 1], linewidth=1, alpha=0.5, color='k')
    # High symmetry k-points from disp.in file
    plt.axvline(x=data[0,0], linewidth=0.5, color='k', alpha=0.5)
    plt.axvline(x=data[50,0], linewidth=0.5, color='k', alpha=0.5)
    plt.axvline(x=data[100,0], linewidth=0.5, color='k', alpha=0.5)
    plt.axvline(x=data[150,0], linewidth=0.5, color='k', alpha=0.5)

    # Read in neutron data
    data1 = np.loadtxt("100-1.gp")
    path1 = data1[:, 0] * data[50, 0]
    freq1 = data1[:, 1] * 1e12 / 2.99792458e8 / 100 

    data2 = np.loadtxt("100-2.gp")
    path2 = data2[:, 0] * data[50, 0]
    freq2 = data2[:, 1] * 1e12 / 2.99792458e8 / 100 

    data3 = np.loadtxt("100-3.gp")
    path3 = data3[:, 0] * data[50, 0]
    freq3 = data3[:, 1] * 1e12 / 2.99792458e8 / 100 

    data4 = np.loadtxt("110-1.gp")
    path4 = data[100, 0] - data4[:, 0] * (data[100, 0] - data[50, 0])
    freq4 = data4[:, 1] * 1e12 / 2.99792458e8 / 100 

    data5 = np.loadtxt("110-2.gp")
    path5 = data[100, 0] - data5[:, 0] * (data[100, 0] - data[50, 0])
    freq5 = data5[:, 1] * 1e12 / 2.99792458e8 / 100 

    data6 = np.loadtxt("110-3.gp")
    path6 = data[100, 0] - data6[:, 0] * (data[100, 0] - data[50, 0])
    freq6 = data6[:, 1] * 1e12 / 2.99792458e8 / 100 

    data7 = np.loadtxt("110-4.gp")
    path7 = data[100, 0] - data7[:, 0] * (data[100, 0] - data[50, 0])
    freq7 = data7[:, 1] * 1e12 / 2.99792458e8 / 100 

    data8 = np.loadtxt("111-1.gp")
    path8 = data8[:, 0] * 2 * (data[150, 0] - data[100, 0]) + data[100, 0]
    freq8 = data8[:, 1] * 1e12 / 2.99792458e8 / 100 

    data9 = np.loadtxt("111-2.gp")
    path9 = data9[:, 0] * 2 * (data[150, 0] - data[100, 0]) + data[100, 0]
    freq9 = data9[:, 1] * 1e12 / 2.99792458e8 / 100 

    data10 = np.loadtxt("111-3.gp")
    path10 = data10[:, 0] * 2 * (data[150, 0] - data[100, 0]) + data[100, 0]
    freq10 = data10[:, 1] * 1e12 / 2.99792458e8 / 100 

    plt.scatter(path1, freq1, marker="^", color='k')
    plt.scatter(path2, freq2, marker="^", color='k')
    plt.scatter(path3, freq3, marker="^", color='k')
    plt.scatter(path4, freq4, marker="^", color='k')
    plt.scatter(path5, freq5, marker="^", color='k')
    plt.scatter(path6, freq6, marker="^", color='k')
    plt.scatter(path7, freq7, marker="^", color='k')
    plt.scatter(path8, freq8, marker="^", color='k')
    plt.scatter(path9, freq9, marker="^", color='k')
    plt.scatter(path10, freq10, marker="^", color='k')

    
    plt.xticks(ticks=[0, data[50, 0], data[100, 0], data[150, 0]], \
               labels=[r'$\Gamma$', 'X', r'$\Gamma^{\prime}$', r'L$^{\prime}$'])
    plt.ylabel(r'Frequency $\omega$ (cm$^{-1}$)')
    plt.xlim(data[0, 0], data[-1, 0])
    plt.ylim(0, )
    plt.show()
    
if __name__ == "__main__":
    filename = "disp7.gp"
    neutron(filename)
