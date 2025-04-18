import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as mcolors

df_mesh1 = pd.read_csv('BigMesh_Li4SiO4_0.0_beryllium_0.0_lead_0.0_2.0_1.csv')
df_mesh2 = pd.read_csv('BigMesh_Li4SiO4_0.0_lead_10.0_lead_0.0_2.0_1.csv')

r_grid = np.linspace(0, 700, num=100) #[NEW] original: (25, 200, num=25), 0, 600
z_grid = np.linspace(-400, 400, num=120) #[NEW] original: (-200, 200, num=50), -700, 700
x_values = r_grid[:-1]  # Take the left edges of bins for x
z_values = z_grid[:-1]  # Take the lower edges of bins for z
vv_points = np.loadtxt("/home/hice1/dcox67/TBR/data/" + 'arc_vv.txt')
#print("vv_points:", vv_points)
x_overlay, z_overlay = vv_points[:, 0], vv_points[:, 1]

reflector_points = np.array([
    [620.0, 0.0],
    [593.80, 107.15],
    [524.64, 202.69],
    [435.08, 276.26],
    [348.42, 319.90],
    [279.26, 328.87],
    [231.67, 302.21],
    [202.76, 242.79],
    [187.29, 157.06],
    [180.76, 54.32],
    [180.76, -54.32],
    [187.29, -157.06],
    [202.76, -242.79],
    [231.67, -302.21],
    [279.26, -328.87],
    [348.42, -319.90],
    [435.08, -276.26],
    [524.64, -202.69],
    [593.80, -107.15]
])
x_reflector, z_reflector = reflector_points[:, 0], reflector_points[:, 1]

channel_inner = np.array([
    [ 5.21338306e+02,  1.33839943e-16],
    [ 5.06909882e+02,  5.90019801e+01],
    [ 4.68939665e+02,  1.11456372e+02],
    [ 4.19847355e+02,  1.51786136e+02],
    [ 3.72249711e+02,  1.75752759e+02],
    [ 3.33880377e+02,  1.80729805e+02],
    [ 3.07192779e+02,  1.65774183e+02],
    [ 2.91166267e+02,  1.32837778e+02],
    [ 2.82681159e+02,  8.58277161e+01],
    [ 2.79113597e+02,  2.96682763e+01],
    [ 2.79113597e+02, -2.96682763e+01],
    [ 2.82681159e+02, -8.58277161e+01],
    [ 2.91166267e+02, -1.32837778e+02],
    [ 3.07192779e+02, -1.65774183e+02],
    [ 3.33880377e+02, -1.80729805e+02],
    [ 3.72249711e+02, -1.75752759e+02],
    [ 4.19847355e+02, -1.51786136e+02],
    [ 4.68939665e+02, -1.11456372e+02],
    [ 5.06909882e+02, -5.90019801e+01]
])
x_inner, z_inner = channel_inner[:, 0], channel_inner[:, 1]

channel_outer = np.array([
[ 5.23397238e+02,  3.39747548e-16],
[ 5.08759609e+02,  5.98574810e+01],
[ 4.70408515e+02,  1.12838039e+02],
[ 4.20945244e+02,  1.53472551e+02],
[ 3.72844863e+02,  1.77692315e+02],
[ 3.33480616e+02,  1.82798416e+02],
[ 3.05668505e+02,  1.67212621e+02],
[ 2.89247100e+02,  1.33464661e+02],
[ 2.80692490e+02,  8.60695398e+01],
[ 2.77113597e+02,  2.97317380e+01],
[ 2.77113597e+02, -2.97317380e+01],
[ 2.80692490e+02, -8.60695398e+01],
[ 2.89247100e+02, -1.33464661e+02],
[ 3.05668505e+02, -1.67212621e+02],
[ 3.33480616e+02, -1.82798416e+02],
[ 3.72844863e+02, -1.77692315e+02],
[ 4.20945244e+02, -1.53472551e+02],
[ 4.70408515e+02, -1.12838039e+02],
[ 5.08759609e+02, -5.98574810e+01]
])
x_outer, z_outer = channel_outer[:, 0], channel_outer[:, 1]

merged = pd.merge(
    df_mesh1, df_mesh2,
    on=['mesh 1 z', 'mesh 1 x', 'energy low [eV]', 'energy high [eV]'],
    suffixes=('_1', '_2')
)
merged['diff'] = merged['mean_1'] - merged['mean_2']

energy_ranges = merged[['energy low [eV]', 'energy high [eV]']].drop_duplicates().values
#print(repr(energy_ranges))

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(25, 8), dpi=600)  # 1 row, 5 columns
axes = axes.flatten()  # Ensure it's a flat iterable list

for i, (energy_low, energy_high) in enumerate(energy_ranges):
    subset = merged[
        (merged['energy low [eV]'] == energy_low) &
        (merged['energy high [eV]'] == energy_high)
    ]
    grouped = subset.groupby(['mesh 1 z', 'mesh 1 x'])['diff'].mean().unstack()
    X, Z = np.meshgrid(x_values, z_values)
    grouped.replace([0, np.inf, -np.inf], np.nan, inplace=True)
    ax = axes[i]
    im = ax.pcolormesh(X, Z, grouped.values, shading='auto', cmap='coolwarm')
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Flux Ratio")

    # Set axis labels and title
    ax.set_xlabel("r")
    ax.set_ylabel("z")
    ax.set_title(f"{energy_low:.2e} eV - {energy_high:.2e} eV")

    # Overlay markers and lines
    ax.scatter(x_overlay, z_overlay, color='cyan', marker='o', s=5, alpha=0.5)
    ax.plot(x_overlay, z_overlay, color='cyan', alpha=0.5, linewidth=1, label='PFC')

    ax.scatter(x_inner, z_inner, color='green', marker='o', s=5, alpha=0.5)
    ax.plot(x_inner, z_inner, color='green', alpha=0.5, linewidth=1, label="Inner Channel")

    ax.scatter(x_outer, z_outer, color='blue', marker='o', s=5, alpha=0.6)
    ax.plot(x_outer, z_outer, color='blue', alpha=0.6, linewidth=1, label="Outer Channel")

    ax.scatter(x_reflector, z_reflector, color='lime', s=5, alpha=0.3)
    ax.plot(x_reflector, z_reflector, color='lime', alpha=0.3, linewidth=1, label="Reflector")

    ax.legend(
        loc='upper center', 
        bbox_to_anchor=(0.5, -0.15),  # Move below subplot
        fancybox=True, shadow=True, ncol=2
    )


# Optimize spacing
plt.subplots_adjust(bottom=0.2)  # Adds space at the bottom for legends
plt.tight_layout()

# Save the figure
fig.savefig("Ratiomap.png", dpi=600)
