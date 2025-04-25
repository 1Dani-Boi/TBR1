import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm

df_mesh1 = pd.read_csv('BigMesh_lead_0.0_lead_0.0_lead_0.0_2.0_1.csv')
df_mesh2 = pd.read_csv('BigMesh_lead_50.0_lead_10.0_lead_0.0_2.0_1.csv')

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

merged['mesh 1 x'] = merged['mesh 1 x']*7
merged['mesh 1 z'] = merged['mesh 1 z']*(800/120)
#print(f'mreged r:{merged['mesh 1 x']}')
merged['diff'] = ((merged['mean_2'] - merged['mean_1']) / merged['mean_1']) #* 100
merged['std'] =  merged['std. dev._2'] / merged['std. dev._1'] #np.abs( ((merged['std. dev._2'] - merged['std. dev._1']) / merged['std. dev._1']) )
# merged.loc[merged['diff'] >= 1, 'diff'] = 0

energy_ranges = merged[['energy low [eV]', 'energy high [eV]']].drop_duplicates().values
#print(repr(energy_ranges))

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(25, 8), dpi=600)  # 1 row, 5 columns
axes = axes.flatten()  # Ensure it's a flat iterable list

unique_z = np.sort(merged['mesh 1 z'].unique())
middle_z = 60#unique_z[len(unique_z) // 2]  # pick middle z

# Step 4: Filter for that z value
middle_slice = merged[merged['mesh 1 z'] == middle_z]
middle_slice = middle_slice.sort_values(by='mesh 1 x')
cleaned = middle_slice.dropna(subset=['diff'])
filtered = cleaned.loc[cleaned.groupby('mesh 1 x')['diff'].idxmax()]
plt.figure(figsize=(10, 6), dpi=300)

plt.errorbar(filtered['mesh 1 x'], filtered['diff'], yerr=filtered['std'] , fmt='x', label='Percent difference between cases', color='black', capsize=3, elinewidth=0.5, ecolor='red')
plt.axvline(523.3, color='red', linestyle='--', label='Channel Inner')
plt.axvline(620, color='blue', linestyle='--', label='Tank Inner')
plt.title(" ")
plt.xlabel("Radial position (cm)")
plt.ylabel("Difference in flux (%)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.xlim([520, 650])
plt.ylim([-10, 10])
plt.savefig("Ratio_vs_x_middle_z.png", dpi=300)
plt.show()

# for i, (energy_low, energy_high) in enumerate(energy_ranges):
#     subset = merged[
#         (merged['energy low [eV]'] == energy_low) &
#         (merged['energy high [eV]'] == energy_high)
#     ]
#     grouped = subset.groupby(['mesh 1 z', 'mesh 1 x'])['diff'].mean().unstack()
#     X, Z = np.meshgrid(x_values, z_values)
#     grouped.replace([0, np.inf, -np.inf], np.nan, inplace=True)
#     ax = axes[i]
#     im = ax.pcolormesh(X, Z, grouped.values, shading='auto', cmap='coolwarm')
#     cbar = fig.colorbar(im, ax=ax)
#     cbar.set_label("Flux Ratio")

#     # Set axis labels and title
#     ax.set_xlabel("r")
#     ax.set_ylabel("z")
#     ax.set_title(f"{energy_low:.2e} eV - {energy_high:.2e} eV")

#     # Overlay markers and lines
#     ax.scatter(x_overlay, z_overlay, color='cyan', marker='o', s=5, alpha=0.5)
#     ax.plot(x_overlay, z_overlay, color='cyan', alpha=0.5, linewidth=1, label='PFC')

#     ax.scatter(x_inner, z_inner, color='green', marker='o', s=5, alpha=0.5)
#     ax.plot(x_inner, z_inner, color='green', alpha=0.5, linewidth=1, label="Inner Channel")

#     ax.scatter(x_outer, z_outer, color='blue', marker='o', s=5, alpha=0.6)
#     ax.plot(x_outer, z_outer, color='blue', alpha=0.6, linewidth=1, label="Outer Channel")

#     ax.scatter(x_reflector, z_reflector, color='lime', s=5, alpha=0.3)
#     ax.plot(x_reflector, z_reflector, color='lime', alpha=0.3, linewidth=1, label="Reflector")

#     ax.legend(
#         loc='upper center', 
#         bbox_to_anchor=(0.5, -0.15),  # Move below subplot
#         fancybox=True, shadow=True, ncol=2
#     )


# # Optimize spacing
# plt.subplots_adjust(bottom=0.2)  # Adds space at the bottom for legends
# plt.tight_layout()

# # Save the figure
# fig.savefig("Ratiomap.png", dpi=600)


# # Individual Plots for verification:

# energy_ranges = df_mesh1[['energy low [eV]', 'energy high [eV]']].drop_duplicates().values

# fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(25, 8), dpi=600)  # 1 row, 5 columns
# axes = axes.flatten()  # Ensure it's a flat iterable list

# for i, (energy_low, energy_high) in enumerate(energy_ranges):
#     subset = df_mesh1[(df_mesh1['energy low [eV]'] == energy_low) & (df_mesh1['energy high [eV]'] == energy_high)]

#     grouped = subset.groupby(['mesh 1 z', 'mesh 1 x'])['mean'].mean().unstack()
#     if energy_high == 2.00e+07:
#         grouped[grouped >= 0.01] = np.nan
#     #X, Z = np.meshgrid(x_values, z_values)
#     grouped.replace([0, np.inf, -np.inf], np.nan, inplace=True)
#     #ax = axes[i]
#     #im = ax.pcolormesh(X, Z, grouped.values, shading='auto', cmap='coolwarm')
#     # im = ax.pcolormesh(
#     #     X, Z, grouped.values,
#     #     shading='auto',
#     #     cmap='coolwarm',
#     #     norm=LogNorm(vmin=np.nanmin(grouped.values[grouped.values > 0]), vmax=np.nanmax(grouped.values))
#     # )
#     # cbar = fig.colorbar(im, ax=ax)
#     # cbar.set_label("Flux")

#     # Set axis labels and title
#     ax.set_xlabel("r")
#     ax.set_ylabel("z")
#     ax.set_title(f"{energy_low:.2e} eV - {energy_high:.2e} eV")

#     # Overlay markers and lines
#     ax.scatter(x_overlay, z_overlay, color='cyan', marker='o', s=5, alpha=0.5)
#     ax.plot(x_overlay, z_overlay, color='cyan', alpha=0.5, linewidth=1, label='PFC')

#     ax.scatter(x_inner, z_inner, color='green', marker='o', s=5, alpha=0.5)
#     ax.plot(x_inner, z_inner, color='green', alpha=0.5, linewidth=1, label="Inner Channel")

#     ax.scatter(x_outer, z_outer, color='blue', marker='o', s=5, alpha=0.6)
#     ax.plot(x_outer, z_outer, color='blue', alpha=0.6, linewidth=1, label="Outer Channel")

#     ax.scatter(x_reflector, z_reflector, color='lime', s=5, alpha=0.3)
#     ax.plot(x_reflector, z_reflector, color='lime', alpha=0.3, linewidth=1, label="Reflector")

#     ax.legend(
#         loc='upper center', 
#         bbox_to_anchor=(0.5, -0.15),  # Move below subplot
#         fancybox=True, shadow=True, ncol=2
#     )


# # Optimize spacing
# plt.subplots_adjust(bottom=0.2)  # Adds space at the bottom for legends
# plt.tight_layout()

# # Save the figure
# fig.savefig("Base_Case.png", dpi=600)





# energy_ranges = df_mesh2[['energy low [eV]', 'energy high [eV]']].drop_duplicates().values

# fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(25, 8), dpi=600)  # 1 row, 5 columns
# axes = axes.flatten()  # Ensure it's a flat iterable list

# for i, (energy_low, energy_high) in enumerate(energy_ranges):
#     subset = df_mesh2[(df_mesh2['energy low [eV]'] == energy_low) & (df_mesh2['energy high [eV]'] == energy_high)]

#     grouped = subset.groupby(['mesh 1 z', 'mesh 1 x'])['mean'].mean().unstack()
#     X, Z = np.meshgrid(x_values, z_values)
#     grouped.replace([0, np.inf, -np.inf], np.nan, inplace=True)
#     ax = axes[i]
#     im = ax.pcolormesh(X, Z, grouped.values, shading='auto', cmap='coolwarm')
#     cbar = fig.colorbar(im, ax=ax)
#     cbar.set_label("Flux")

#     # Set axis labels and title
#     ax.set_xlabel("r")
#     ax.set_ylabel("z")
#     ax.set_title(f"{energy_low:.2e} eV - {energy_high:.2e} eV")

#     # Overlay markers and lines
#     ax.scatter(x_overlay, z_overlay, color='cyan', marker='o', s=5, alpha=0.5)
#     ax.plot(x_overlay, z_overlay, color='cyan', alpha=0.5, linewidth=1, label='PFC')

#     ax.scatter(x_inner, z_inner, color='green', marker='o', s=5, alpha=0.5)
#     ax.plot(x_inner, z_inner, color='green', alpha=0.5, linewidth=1, label="Inner Channel")

#     ax.scatter(x_outer, z_outer, color='blue', marker='o', s=5, alpha=0.6)
#     ax.plot(x_outer, z_outer, color='blue', alpha=0.6, linewidth=1, label="Outer Channel")

#     ax.scatter(x_reflector, z_reflector, color='lime', s=5, alpha=0.3)
#     ax.plot(x_reflector, z_reflector, color='lime', alpha=0.3, linewidth=1, label="Reflector")

#     ax.legend(
#         loc='upper center', 
#         bbox_to_anchor=(0.5, -0.15),  # Move below subplot
#         fancybox=True, shadow=True, ncol=2
#     )


# # Optimize spacing
# plt.subplots_adjust(bottom=0.2)  # Adds space at the bottom for legends
# plt.tight_layout()

# # Save the figure
# fig.savefig("MultDopant_Case.png", dpi=600)




# # Uncertainty plot:

# energy_ranges = df_mesh1[['energy low [eV]', 'energy high [eV]']].drop_duplicates().values

# fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(25, 8), dpi=600)  # 1 row, 5 columns
# axes = axes.flatten()  # Ensure it's a flat iterable list

# for i, (energy_low, energy_high) in enumerate(energy_ranges):
#     subset = df_mesh1[(df_mesh1['energy low [eV]'] == energy_low) & (df_mesh1['energy high [eV]'] == energy_high)]

#     grouped = subset.groupby(['mesh 1 z', 'mesh 1 x'])['std. dev.'].mean().unstack()
#     groupedmean = subset.groupby(['mesh 1 z', 'mesh 1 x'])['mean'].mean().unstack()
#     percent = (grouped / groupedmean) / 100 
#     X, Z = np.meshgrid(x_values, z_values)
#     grouped.replace([0, np.inf, -np.inf], np.nan, inplace=True)
#     groupedmean.replace([0, np.inf, -np.inf], np.nan, inplace=True)
#     ax = axes[i]
#     im = ax.pcolormesh(X, Z, percent.values, shading='auto', cmap='coolwarm')
#     cbar = fig.colorbar(im, ax=ax)
#     cbar.set_label("Uncertainty (Percentage)")

#     # Set axis labels and title
#     ax.set_xlabel("r")
#     ax.set_ylabel("z")
#     ax.set_title(f"{energy_low:.2e} eV - {energy_high:.2e} eV")

#     # Overlay markers and lines
#     ax.scatter(x_overlay, z_overlay, color='cyan', marker='o', s=5, alpha=0.5)
#     ax.plot(x_overlay, z_overlay, color='cyan', alpha=0.5, linewidth=1, label='PFC')

#     ax.scatter(x_inner, z_inner, color='green', marker='o', s=5, alpha=0.5)
#     ax.plot(x_inner, z_inner, color='green', alpha=0.5, linewidth=1, label="Inner Channel")

#     ax.scatter(x_outer, z_outer, color='blue', marker='o', s=5, alpha=0.6)
#     ax.plot(x_outer, z_outer, color='blue', alpha=0.6, linewidth=1, label="Outer Channel")

#     ax.scatter(x_reflector, z_reflector, color='lime', s=5, alpha=0.3)
#     ax.plot(x_reflector, z_reflector, color='lime', alpha=0.3, linewidth=1, label="Reflector")

#     ax.legend(
#         loc='upper center', 
#         bbox_to_anchor=(0.5, -0.15),  # Move below subplot
#         fancybox=True, shadow=True, ncol=2
#     )


# # Optimize spacing
# plt.subplots_adjust(bottom=0.2)  # Adds space at the bottom for legends
# plt.tight_layout()

# # Save the figure
# fig.savefig("Base_Case_Uncertainty_Percentage.png", dpi=600)

