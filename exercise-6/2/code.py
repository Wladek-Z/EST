import numpy as np
import matplotlib.pyplot as plt

def speeds(filename):
    """Calculate the speed of sound from the phonon dispersion data"""
    data = np.loadtxt(filename)
    q = data[:, 0]
    omega = data[:, 1]
    # Calculate the speed of sound as the slope of the dispersion at q=0
    v_s = omega[1] / q[1]
    return v_s