import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm

df_mesh1 = pd.read_csv('BigMesh_lead_0.0_lead_0.0_lead_0.0_2.0_1.csv')
df_mesh2 = pd.read_csv('BigMesh_lead_50.0_lead_10.0_lead_0.0_2.0_1.csv')



flux_df.columns = [' '.join(filter(None, map(str, col))).strip() for col in flux_df.columns]

print(flux_df.columns)
energy_ranges = flux_df[['energy low [eV]', 'energy high [eV]']].drop_duplicates().values
plt.figure()
for i, (energy_low, energy_high) in enumerate(energy_ranges):
    subset = flux_df[(flux_df['energy low [eV]'] == energy_low) & (flux_df['energy high [eV]'] == energy_high)]
    #print('subset:' , subset.columns)
    voxel_volume = 2*20*1  # cm^3

    flux_std = (subset['std. dev.'] / voxel_volume) * neutrons_per_second
    radial_positions = subset['mesh 2 x']*2-20
    plt.errorbar( radial_positions,  (subset['mean'] / voxel_volume) * neutrons_per_second, label=f'{energy_low:.2e} to {energy_high:.2e} eV' , yerr=flux_std, capsize=3)
    #plt.plot(radial_positions, (subset['mean'] / voxel_volume) * neutrons_per_second, label=f'{energy_low:.2e} to {energy_high:.2e} eV')
    #plt.fill_between(subset['mesh 2 x']*2-20, (subset['mean'] / voxel_volume) * neutrons_per_second - flux_std, (subset['mean'] / voxel_volume) * neutrons_per_second + flux_std, alpha=0.3)
    
plt.yscale('log')
plt.xlabel('Radial Distance from Blanket Edge (cm)')
plt.ylabel('Flux [n/cm^2-s]')
plt.title(f'Flux vs X for {energy_low} to {energy_high} eV')
plt.grid(True, which='both', ls='--')
plt.axvline(x=0, color='red', linestyle='--', linewidth=1, label='PFC')
plt.axvline(x=3.3, color='black', linestyle='--', linewidth=1, label='Channel Inner')
plt.axvline(x=23.3, color='orange', linestyle='--', linewidth=1, label='Channel Outer')
plt.axvline(x=100, color='cyan', linestyle='--', linewidth=1, label='Blanket Outer')
plt.legend()
plt.xlim(-24, 150)
plt.ylim(1e-10,0.4e-5)
plt.savefig(f'Flux vs X for {energy_low:.2e} to {energy_high:.2e} eV together 2.png')