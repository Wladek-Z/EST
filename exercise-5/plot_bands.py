#!/usr/bin/env python3
"""
Plot band structure along a k-point path in the Brillouin zone.
Reads data from quartz-bands.dat.gnu and identifies high-symmetry points.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def read_band_data(filename):
    """
    Read band structure data from .gnu file.
    Returns: k_distance (distance along path), energies (in eV)
    """
    data = np.loadtxt(filename)
    k_distance = data[:, 0]
    energies = data[:, 1]
    return k_distance, energies

def get_high_symmetry_labels():
    """
    Get the labels and relative positions for the known k-point path.
    Path: Gamma - M - K - Gamma - A - L - H - Gamma
    Each segment has 50 k-points.
    Returns: indices and labels of high-symmetry points
    """
    # K-point path with 50 points per segment
    npts_per_segment = 50
    path_labels = [
        (0 * npts_per_segment, 'Γ'),
        (1 * npts_per_segment, 'M'),
        (2 * npts_per_segment, 'K'),
        (3 * npts_per_segment, 'Γ'),
        (4 * npts_per_segment, 'A'),
        (5 * npts_per_segment, 'L'),
        (6 * npts_per_segment, 'H'),
        (7 * npts_per_segment, 'Γ'),
    ]
    
    indices = [idx for idx, _ in path_labels]
    labels = [label for _, label in path_labels]
    
    return indices, labels

def plot_band_structure(data_dir, prefix='Q', band_file_prefix='quartz'):
    """
    Plot the band structure for the specified material.
    
    Args:
        data_dir: Directory containing the band structure files
        prefix: Prefix for directory name (e.g., 'Q' for Q-bands/)
        band_file_prefix: Prefix for actual band data files (e.g., 'quartz')
    """
    gnu_file = Path(data_dir) / f'{prefix}-bands' / f'{band_file_prefix}-bands.dat.gnu'
    
    if not gnu_file.exists():
        print(f"Error: {gnu_file} not found")
        return
    
    # Read band energies
    k_distance, energies = read_band_data(str(gnu_file))
    
    # Get high-symmetry points from the known k-point path
    special_indices, labels = get_high_symmetry_labels()
    
    # Convert indices to distances along the k-point path
    special_distances = []
    for idx in special_indices:
        if idx < len(k_distance):
            special_distances.append(k_distance[idx])
        else:
            # For points beyond the data, use the last distance
            special_distances.append(k_distance[-1])
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Plot all bands - scatter plot to show individual bands clearly
    unique_distances = np.unique(k_distance)
    
    # Group energies by k-point (all bands at each k-point)
    for i, x in enumerate(unique_distances):
        mask = np.isclose(k_distance, x)
        y_vals = energies[mask]
        ax.scatter([x] * len(y_vals), y_vals, c='blue', s=1, alpha=0.7)
    
    # Mark high-symmetry points with vertical lines and labels
    for dist, label in zip(special_distances, labels):
        ax.axvline(x=dist, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    
    # Set x-axis with high-symmetry point labels
    ax.set_xticks(special_distances)
    ax.set_xticklabels(labels, fontsize=12, fontweight='bold')
    
    # Labels and formatting
    ax.set_ylabel('Energy (eV)', fontsize=12)
    ax.set_title(f'{prefix} Band Structure: Γ-M-K-Γ-A-L-H-Γ', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add Fermi level at 0
    ax.axhline(y=0, color='red', linestyle=':', linewidth=1.5, alpha=0.7, label='Fermi level')
    ax.legend(fontsize=10)
    
    plt.tight_layout()
    
    output_path = Path(data_dir) / f'{prefix}_band_structure.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()  # Close the figure instead of showing it
    
    print(f"Band structure plot saved to {output_path}")

if __name__ == '__main__':
    import sys
    
    # Set the data directory (task 3.2 in exercise-5)
    data_dir = Path('/home/nairne/Desktop/EST/exercise-5/task 3.2')
    
    # Plot for quartz (Q)
    print("Plotting Quartz (Q) band structure...")
    plot_band_structure(data_dir, prefix='Q', band_file_prefix='quartz')
    
    # Uncomment to plot other materials
    # plot_band_structure(data_dir, prefix='ASn', band_file_prefix='alpha')
    # plot_band_structure(data_dir, prefix='BSn', band_file_prefix='beta')
