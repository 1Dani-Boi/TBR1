import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors as mcolors



spectrum_mesh1 = pd.read_csv('Spectrum_lead_0.0_lead_0.0_lead_0.0_2.0_1.csv')
spectrum_mesh1 = pd.read_csv('Spectrum_lead_50.0_lead_10.0_lead_0.0_2.0_1.csv')



plt.figure()
# trim the last energy filter bin edge to make the number of x and y values the same
plt.step(cell_energy_filter.values[:-1], channel_flux_adjusted, label='Flux in FLiBe Channel')
plt.step(cell_energy_filter.values[:-1], blanket_flux_adjusted, label='Flux in FLiBe Tank')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Flux [n/cm^2-s]')
plt.xlabel('Energy [eV]')
plt.grid(True, which='both', ls='--')
plt.legend()
plt.savefig('Comparison_CHannel_Blanket_Spectra.png')






