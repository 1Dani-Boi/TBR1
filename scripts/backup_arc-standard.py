import arc_2 as anp
import openmc
import openmc.data
import numpy as np
import os
import sys
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import subprocess
import pandas as pd
import shutil
import matplotlib.colors as mcolors
import arc_2.materials_2 as mat

openmc.config["cross_sections"] = '/home/hice1/dcox67/endfb-viii.0-hdf5/cross_sections.xml'

# ==============================================================================
# Geometry
# ==============================================================================
def create_arc(Li6_enrichment, dopant, dopant_mass, multiplier_material, multiplier_thickness, reflector_material, reflector_thickness, channel_thickness, order):
    device = anp.generate_device(Li6_enrichment = Li6_enrichment, dopant = dopant, dopant_mass = dopant_mass,
                                multiplier_material = multiplier_material, multiplier_thickness = multiplier_thickness, reflector_material = reflector_material,
                                reflector_thickness = reflector_thickness, channel_thickness = channel_thickness, order = order)
    
    # Plotting
    # plot = openmc.Plot()
    # plot.filename = 'geometry_plot'
    # plot.basis = 'xz'
    # plot.origin = (350, 0, 0)
    # plot.width = (700, 800)
    # plot.pixels = (plot.width[0]*10, plot.width[1]*10)
    # plot.color_by = 'cell'

    # plot.cell_colors = { device._cells[0].id: 'beige', device._cells[1].id: 'lightcoral', device._cells[2].id: 'yellow', device._cells[3].id: 'orange', device._cells[4].id: 'lime', device._cells[5].id: 'navy', device._cells[6].id: 'lightcyan', device._cells[7].id: 'black'}
    # color = ['beige', 'lightcoral', 'yellow', 'orange', 'lime', 'navy', 'lightcyan', 'black']
    # color_dict = {cell.id: color[i] for i, cell in enumerate(device._cells)}
    # plot.colors = color_dict
    # for count, cell in enumerate(device._cells):
    #     print(f'Device name: {cell.name} with color: {color[count]}')
    
    # plot.highlight_domains(geometry=device.geometry, domains=device._cells)
    
    # plots = openmc.Plots([plot])
    # plots.export_to_xml()
    
    # ==============================================================================
    # Settings
    # ==============================================================================
    
    """ Source Definition """
    source = openmc.Source()
    source.space = openmc.stats.CylindricalIndependent(openmc.stats.Discrete(400, 1), openmc.stats.Uniform(a=-np.pi/18, b=np.pi/18), openmc.stats.Discrete(0, 1)) # original openmc.stats.Discrete(450, 1), openmc.stats.Uniform(a=-np.pi/18, b=np.pi/18)
    source.angles = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([14.1E6], [1.0])
    
    device.settings.source = source
    # energy filter
    # energy_filter = openmc.EnergyFilter.from_group_structure("CCFE-709")
    # ==============================================================================
    # Blanket Material
    # ==============================================================================
    # breeding_material = openmc.Material(material_id = 56)  # Pb84.2Li15.8
    # breeding_material.add_element('Pb', 84.2)
    # breeding_material.add_element('Li', 15.8)
    # breeding_material.set_density('g/cm3', 11.)
    # device.blanket.fill = breeding_material
    # ==============================================================================
    # Tallies
    # ==============================================================================
    # """ Cylindrical Mesh Tally """
    r_grid = np.linspace(0, 700, num=100) #[NEW] original: (25, 200, num=25), 0, 600
    z_grid = np.linspace(-400, 400, num=120) #[NEW] original: (-200, 200, num=50), -700, 700
    mesh = openmc.CylindricalMesh(r_grid=r_grid, z_grid=z_grid) #[NEW]
    mesh.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
    mesh_filter = openmc.MeshFilter(mesh)

    energy_filter = openmc.EnergyFilter(np.array([1,10.0e6,14.0e6,20.0e6]))
    neutron_particle_filter = openmc.ParticleFilter(['neutron'])
    
    # device.add_tally('Mesh Tally', ['flux', '(n,Xt)', 'heating-local', 'absorption'], filters=[mesh_filter])
    device.add_tally('Mesh Tally', ['flux'], filters=[mesh_filter, energy_filter, neutron_particle_filter])
    
    # """ TBR Tally """
    # tbr_filter1 = openmc.MaterialFilter(anp.tungsten)
    # device.add_tally('Tbr Plasma-facing Component Tally ', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[tbr_filter1])
    
    # tbr_filter2 = openmc.MaterialFilter(device.vcrti_VV)
    # device.add_tally('Tbr Vacuum Vessel Tally ', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[tbr_filter2])
    
    # tbr_filter3 = openmc.MaterialFilter(device.doped_flibe_channels)
    tbr_filter3 = openmc.CellFilter(device.get_cell(name = 'channel'))
    device.add_tally('TBR Channel Tally', ['(n,Xt)'], nuclides = ['Li6', 'Li7'], filters = [tbr_filter3])
    
    # tbr_filter4 = openmc.MaterialFilter(device.vcrti_BI)
    # device.add_tally('Tbr Tank Inner Tally ', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[tbr_filter4])

    tbr_filter5 = openmc.CellFilter(device.get_cell(name = 'blanket'))
    # tbr_filter5 = openmc.MaterialFilter(device.doped_flibe_blanket) #device.doped_flibe, initially wanted 'doped_mat' which doesn't exist
    device.add_tally('TBR Blanket Tally', ['(n,Xt)'], nuclides = ['Li6', 'Li7'], filters = [tbr_filter5])
    
    cell_energy_filter = openmc.EnergyFilter.from_group_structure('CCFE-709')
    device.add_tally('Channel Flux Spectrum', ['flux'], filters = [tbr_filter3, cell_energy_filter, neutron_particle_filter])
    device.add_tally('Blanket Flux Spectrum', ['flux'], filters = [tbr_filter5, cell_energy_filter, neutron_particle_filter])
    # tbr_filter6 = openmc.MaterialFilter(device.vcrti_BO)
    # device.add_tally('Tbr Tank Outer Tally ', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[tbr_filter6])
    
    r_grid2 = np.linspace(500, 700, num=100) #[NEW] original: (25, 200, num=25), 0, 600
    z_grid2 = np.array([-10,10]) #[NEW] original: (-200, 200, num=50), -700, 700
    mesh2 = openmc.CylindricalMesh(r_grid=r_grid2, z_grid=z_grid2) #[NEW]
    mesh2.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
    mesh2_filter = openmc.MeshFilter(mesh2)

    #device.add_tally('1D Flux', ['flux'], filters=[mesh2_filter, energy_filter, neutron_particle_filter])
    device.add_tally('1D Flux', ['flux'], filters=[mesh2_filter, neutron_particle_filter])

    # ---------- Dosage: -----------------------------------
    # energies = np.array([1e-5, 1e-3, 1.0, 10.0])  # in MeV
    # dose_factors = np.array([1e-8, 1e-7, 1e-6])   # in Sv cm² (example only)

    # energy_func_filter = openmc.EnergyFunctionFilter(energies, dose_factors)
    x0, y0, z0 = 621.0, 1, 0  # coordinates 
    delta = 10  # cm, small width
    dose_mesh = openmc.RegularMesh()
    dose_mesh.dimension = [1, 1, 1]
    dose_mesh.lower_left = [x0-delta/2, y0-delta/2, z0-delta/2]
    dose_mesh.upper_right = [x0+delta/2, y0+delta/2, z0+delta/2]

    dose_mesh_filter = openmc.MeshFilter(dose_mesh)

    # Neutron ambient dose equivalent
    n_energies, n_factors = openmc.data.dose_coefficients('neutron', 'AP')  # 'h' for ambient dose equivalent

    # Photon ambient dose equivalent
    p_energies, p_factors = openmc.data.dose_coefficients('photon', 'AP')

    # Add the energy function filter to tallies
    n_filter = openmc.EnergyFunctionFilter(n_energies, n_factors)
    p_filter = openmc.EnergyFunctionFilter(p_energies, p_factors)

    device.add_tally('Neutron Dose', ['flux'], filters = [n_filter, dose_mesh_filter])
    device.add_tally('Photon Dose', ['flux'], filters = [p_filter, dose_mesh_filter])

    #energyeet = 14e6
    #sigma_t = mat.tungsten.macroscopic_total_xs(energyeet)
    #mfp = 1.0 / sigma_t
    #print(f"MFP of W at {energyeet:.2e} eV: {mfp:.5f} cm")

    # ==============================================================================
    # Run
    # ==============================================================================
    
    device.settings.photon_transport = True
    device.survival_biasing = True
    device.build()
    device.export_to_xml(remove_surfs=True)
    
    # geometry_plot = "geometry_plot.png"
    # if not os.path.exists(geometry_plot):
    #      openmc.plot_geometry() #path_input = 'plots.xml'
    
    # set run parameters
    # device.settings.threads = 24
    device.settings.particles = int(1e6)
    device.settings.batches = 10  
    device.settings.inactive = 1  

    # ww = openmc.WeightWindows(
    #     mesh=my_mesh,  # openmc.RegularMesh or openmc.CylindricalMesh
    #     lower_ww=lower_bounds,  # ndarray of shape (nx, ny, nz)
    #     upper_ww=upper_bounds,
    #     energy_bounds=energy_bounds  # optional
    # )
    # device.settings.weight_windows = [ww]

    # device.settings.variance_reduction = True

    return device
# ================================================================================
# additional runs
# ================================================================================

def make_materials_geometry_tallies(Li6_enrichment, dopant = str(sys.argv[1]), dopant_mass = float(sys.argv[2]),
                                multiplier_material = sys.argv[3], multiplier_thickness = float(sys.argv[4]), reflector_material = sys.argv[5],
                                reflector_thickness = float(sys.argv[6]), channel_thickness = float(sys.argv[7]), order = sys.argv[8]):
    """Makes a neutronics model of a blanket and simulates the TBR value.

    Arguments:
        enrichment (float): the enrichment percentage of Li6 in the breeder material
    
    Returns:
        resutsl (dict): simulation tally results for TBR along with the standard deviation and enrichment
    """
    # RUN OPENMC
    device = create_arc(Li6_enrichment, dopant, dopant_mass, multiplier_material, multiplier_thickness, reflector_material, reflector_thickness, channel_thickness, order)
    print(device.Li6_enrichment)
    
    #remove old output files
    for file in os.listdir('.'):
        if file.endswith('.h5'):
            os.remove(file)
    print('running with a lot of particles!!!!')
    sp_filename = device.run()  # runs with reduced amount of output printing, output = false

    # OPEN OUPUT FILE
    sp = openmc.StatePoint(sp_filename)

    tbr_tally = sp.get_tally(name='TBR Blanket Tally')
    tbr_tally_2 = sp.get_tally(name='TBR Channel Tally')
    flux_tally = sp.get_tally(name='Mesh Tally')
    channel_tally = sp.get_tally(name='Channel Flux Spectrum')
    blanket_tally = sp.get_tally(name='Blanket Flux Spectrum')
    flux_tally_21 = sp.get_tally(name='1D Flux')
    cell_energy_filter = openmc.EnergyFilter.from_group_structure('CCFE-709')

    df = tbr_tally.get_pandas_dataframe()
    df2 = tbr_tally_2.get_pandas_dataframe()
    df_mesh = flux_tally.get_pandas_dataframe()
    channel_df = channel_tally.get_pandas_dataframe()
    blanket_df = blanket_tally.get_pandas_dataframe()
    flux_df = flux_tally_21.get_pandas_dataframe()


    channel_flux = channel_tally.mean.flatten()
    blanket_flux = blanket_tally.mean.flatten()

    vol_calc_load = openmc.VolumeCalculation.from_hdf5('/home/hice1/dcox67/TBR/data/arc-1_volumes.h5')
    blanket_volume = vol_calc_load.volumes[8].n
    channel_volume = vol_calc_load.volumes[5].n

    reactor_power = 500e6 # in Watt  (from JBall)
    neutrons_per_J = 3.546e11 # 1MJ/17.6 Mev  (from JBall)
    neutrons_per_second = reactor_power / neutrons_per_J  # [Watt / Joule] = [1/s]

    channel_flux_adjusted = (channel_flux / channel_volume) * neutrons_per_second 
    blanket_flux_adjusted = (blanket_flux / blanket_volume) * neutrons_per_second 

    new_df = pd.DataFrame({
        'channel_flux_adjusted': channel_flux_adjusted,
        'blanket_flux_adjusted': blanket_flux_adjusted,
        'energy_filter_values': cell_energy_filter.values[:-1]
    })
    
    new_df.to_csv(f'Spectrum_{str(sys.argv[1])}_{float(sys.argv[2])}_{sys.argv[3]}_{float(sys.argv[4])}_{sys.argv[5]}_{float(sys.argv[6])}_{float(sys.argv[7])}_{sys.argv[8]}_particles:{device.settings.particles}.csv', index=False)


    combined_df = pd.concat([df, df2], ignore_index = True)
    combined_df.to_csv(f'{device.Li6_enrichment}%.csv', index = False)
    tbr_tally_result = df['mean'].sum() + df2['mean'].sum()
    tbr_tally_std_dev = df['std. dev.'].sum() + df2['mean'].sum()


    


    # ========================
    # Plotting               =
    # ========================

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

    df_mesh.columns = [' '.join(col).strip() for col in df_mesh.columns]
    energy_ranges = df_mesh[['energy low [eV]', 'energy high [eV]']].drop_duplicates().values

    print(df_mesh.columns)

    df_mesh.to_csv(f'BigMesh_{str(sys.argv[1])}_{float(sys.argv[2])}_{sys.argv[3]}_{float(sys.argv[4])}_{sys.argv[5]}_{float(sys.argv[6])}_{float(sys.argv[7])}_{sys.argv[8]}_particles:{device.settings.particles}.csv', index=False)

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(25, 8), dpi=600)  # 1 row, 5 columns
    axes = axes.flatten()  # Ensure it's a flat iterable list

    for i, (energy_low, energy_high) in enumerate(energy_ranges):
        # Filter data for the specific energy range
        subset = df_mesh[(df_mesh['energy low [eV]'] == energy_low) & (df_mesh['energy high [eV]'] == energy_high)]

        # Group by x and y, then take the mean in case of duplicates
        grouped = subset.groupby(['mesh 1 z', 'mesh 1 x'])['mean'].mean().unstack()
        X, Z = np.meshgrid(x_values, z_values)

        # norm = mcolors.LogNorm(vmin=grouped.values.min(), vmax=grouped.values.max())

        # Ensure there are no zero or negative values in the data
        grouped.replace(0, np.nan, inplace=True)  # Replace zeros with NaN to avoid log issues
        
        # Select subplot
        ax = axes[i]

        # Plot
        im = ax.pcolormesh(X, Z, grouped.values, shading='auto', cmap='plasma')
        
        # Add a colorbar for the image
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Mean Flux")

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
    fig.savefig(f"flux_all_energy_ranges_dopant:{float(sys.argv[2])}_mult:{float(sys.argv[6])}_particles:{device.settings.particles}.png", dpi=600)

    
    # =================== EDGE PLOT =================================================================================


    flux_df.columns = [' '.join(filter(None, map(str, col))).strip() for col in flux_df.columns] #.columns = [' '.join(col).strip() for col in df_mesh.columns]
    # flux_df.to_csv(f'Edge_CSV_dopant:{float(sys.argv[2])}.csv')
    # print(flux_df.columns)
    #energy_ranges = flux_df[['energy low [eV]', 'energy high [eV]']].drop_duplicates().values
    plt.figure()
    # for i, (energy_low, energy_high) in enumerate(energy_ranges):
    #     subset = flux_df[(flux_df['energy low [eV]'] == energy_low) & (flux_df['energy high [eV]'] == energy_high)]
    #     #print('subset:' , subset.columns)
    #     voxel_volume = 2*20*1  # cm^3

    #     flux_std = (subset['std. dev.'] / voxel_volume) * neutrons_per_second
    #     radial_positions = subset['mesh 2 x']*2-20
    #     plt.errorbar( radial_positions,  (subset['mean'] / voxel_volume) * neutrons_per_second, label=f'{energy_low:.2e} to {energy_high:.2e} eV' , yerr=flux_std, capsize=3)
    #     #plt.plot(radial_positions, (subset['mean'] / voxel_volume) * neutrons_per_second, label=f'{energy_low:.2e} to {energy_high:.2e} eV')
    #     #plt.fill_between(subset['mesh 2 x']*2-20, (subset['mean'] / voxel_volume) * neutrons_per_second - flux_std, (subset['mean'] / voxel_volume) * neutrons_per_second + flux_std, alpha=0.3)
    print(flux_df.columns)
    voxel_volume = 2*20*1  # cm^3
    flux_std = (flux_df['std. dev.'] / voxel_volume) * neutrons_per_second
    radial_positions = flux_df['mesh 2 x']*2-20
    plt.errorbar( radial_positions,  (flux_df['mean'] / voxel_volume) * neutrons_per_second, yerr=flux_std, capsize=3)  
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
    plt.savefig(f'Flux vs X all energies together for mult:{float(sys.argv[6])}.png')

    # DOSE STUFF: ---------------------------------------------------------------------------------------------------------
    tally_n = sp.get_tally(name='Neutron Dose')
    tally_p = sp.get_tally(name='Photon Dose')

    dose_n = tally_n.mean.ravel()[0]
    dose_p = tally_p.mean.ravel()[0]

    print(f"Neutron dose per source particle: {dose_n:.3e} Sv/source")
    print(f"Photon dose per source particle: {dose_p:.3e} Sv/source")

    photons_per_second = 2* neutrons_per_second
    occupancy_seconds = 8 * 5 * 52 * 3600

    dose_rate_n = dose_n * neutrons_per_second
    dose_rate_p = dose_p * photons_per_second
    total_dose_rate = dose_rate_n + dose_rate_p

    # Annual dose
    annual_dose = total_dose_rate * occupancy_seconds
    print(f"Annual dose at point: {annual_dose:.3e} Sv/year")






    return {'enrichment': device.Li6_enrichment,
            'tbr_tally_result': tbr_tally_result,
            'tbr_tally_std_dev': tbr_tally_std_dev}

results = []
for Li6_enrichment in [7.5]:  # percentage enrichment from 0% Li6 to 100% Li6, 0.01, 7.5, 15, 25, 50, 75, 99.99
    results.append(make_materials_geometry_tallies(Li6_enrichment, str(sys.argv[1]), float(sys.argv[2]),
                    sys.argv[3], float(sys.argv[4]), sys.argv[5],
                    float(sys.argv[6]), float(sys.argv[7]), sys.argv[8]))
print(results)
#results.append(make_materials_geometry_tallies(7.5))
# PLOTS RESULTS
x = [entry['enrichment'] for entry in results]
y = [entry['tbr_tally_result'] for entry in results]
#error_y = {'array': [entry['tbr_tally_std_dev'] for entry in results]}
plt.plot(x, y)
plt.title="TBR as a function of Li6 enrichment",
plt.xtitle="Li6 enrichment (%)",
plt.ytitle="TBR"
plt.grid()
plt.show()
# ================================================================================
print(str(sys.argv[9]))

def move_files(source_dir, target_dir):
    # Check if the source directory exists
    if not os.path.exists(source_dir):
        print(f"Error: Source directory '{source_dir}' does not exist.")
        return
    
    # Create the target directory if it doesn't exist
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    
    # Loop through all files in the source directory
    for filename in os.listdir(source_dir):
        source_file = os.path.join(source_dir, filename)
        
        # Skip directories (we only want to move files)
        if os.path.isdir(source_file):
            continue
        
        if filename == 'arc-standard.py':
            continue
        if filename == 'backup_arc-standard.py':
            continue
        if filename == 'slurm_running.sh':
            continue
        # Create the full target file path
        target_file = os.path.join(target_dir, filename)
        
        # Move the file
        try:
            shutil.move(source_file, target_file)
            print(f"Moved '{filename}' to '{target_dir}'")
        except Exception as e:
            print(f"Error moving '{filename}': {e}")

os.mkdir(str(sys.argv[9]))
move_files('/home/hice1/dcox67/skibidi/TBR/scripts', str(sys.argv[9]))
print("OpenMC files moved to new directory:", str(sys.argv[9]))
'''
try:
    if sys.argv[9] is not None:
        os.mkdir(str(sys.argv[9]))
        device.move_files(str(sys.argv[9]))
        print("OpenMC files moved to new directory:", str(sys.argv[9]))

except:
    print("No directory specified, using this one")
'''
# =============================================
# tally plot
# =============================================
# out_file = ""
# for file in os.listdir('.'):
#     if file.endswith('.h5'):
#         if file != "summary.h5":
#             out_file = file
# command = ["openmc-plot-mesh-tally", out_file]
# # Run the command
# subprocess.run(command)