import numpy as np
import matplotlib.pyplot as plt
import os

def phonon_dos(folder):
    """Plot all phonon DOS from dat files in a single folder"""
    # Read in the data 
    freq = []
    dos = []
    pdos_1 = []
    pdos_2 = []
    q = []
    for file in os.listdir(folder):
        filename = os.path.join(folder, file)
        freq_, dos_, pdos_1_, pdos_2_ = np.loadtxt(filename, unpack=True)
        freq.append(freq_)
        dos.append(dos_)
        pdos_1.append(pdos_1_)
        pdos_2.append(pdos_2_)
        # Find the number in the filename
        q.append(int(file.split('dfpt')[1].split('.')[0]))

    # Sort the data by q
    q, freq, dos, pdos_1, pdos_2 = zip(*sorted(zip(q, freq, dos, pdos_1, pdos_2)))

    # Plot each DOS
    fig = plt.figure(figsize=(10, 10))
    num = len(q) - 1
    for i in range(len(q)):
        plt.plot(freq[num-i], dos[num-i] + 0.1*num - i*0.1, label=rf'$q$ = {num+2-i}')
        plt.axhline(0.1*num - i*0.1, color='k', linestyle='--', linewidth=0.5, alpha=0.5)

    plt.xlabel(r'Frequency $\omega$ (cm$^{-1}$)', fontsize=12)
    plt.ylabel(r'Phonon DOS (states/cm$^{-1}$/u.c.)', fontsize=12)
    plt.legend(loc='upper left')
    plt.xlim(0, freq[0][-1])
    plt.tight_layout()
    plt.show()

def phonon_dispersion(folder):
    """Plot all phonon dispersions from folder on the same figure"""
    # create custom colour list for 10 phonon dispersions
    colours = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    # Read in the data 
    counter = -1
    data_list = []
    q_list = []
    fig = plt.figure(figsize=(10,8))
    #for file in os.listdir(folder):
    for i in np.arange(10, 0, -1):
        counter += 1
        file = f"disp{i}.gp"
        filename = os.path.join(folder, file)
        data = np.loadtxt(filename)
        nbands = data.shape[1] - 1
        data_list.append(data)
        q_list.append(i)
        plt.plot(data[:, 0], data[:, 1], linewidth=1.5, alpha=0.5, color=colours[counter], label=rf'$q$ = {i}')
        for band in range(1, nbands):
            plt.plot(data[:, 0], data[:, band + 1], linewidth=1.5, alpha=0.5, color=colours[counter])

    # sort data by q
    q_list, data_list = zip(*sorted(zip(q_list, data_list)))
    data_10 = data_list[-1]
        
    # High symmetry k-points from disp.in file
    # Use converged data
    plt.axvline(x=data_10[0,0], linewidth=0.5, color='k', alpha=0.5)
    plt.axvline(x=data_10[50,0], linewidth=0.5, color='k', alpha=0.5)
    plt.axvline(x=data_10[100,0], linewidth=0.5, color='k', alpha=0.5)
    plt.axvline(x=data_10[150,0], linewidth=0.5, color='k', alpha=0.5)
    plt.axvline(x=data_10[200,0], linewidth=0.5, color='k', alpha=0.5)
    plt.axvline(x=data_10[250,0], linewidth=0.5, color='k', alpha=0.5)
    
    plt.xticks(ticks=[0, data_10[50, 0], data_10[100, 0], data_10[150, 0], data_10[200, 0], data_10[250, 0]], \
               labels=[r'$\Gamma$', 'X', 'W', 'K', r'$\Gamma$', 'L'])
    plt.ylabel(r'Frequency $\omega$ (cm$^{-1}$)', fontsize=12)
    plt.xlabel(r'Wave vector $k$', fontsize=12)
    plt.xlim(data_10[0, 0], data_10[-1, 0])
    plt.ylim(0, )
    plt.legend()
    plt.show()

if __name__ == "__main__":
    folder = 'disp-data'
    phonon_dispersion(folder)
