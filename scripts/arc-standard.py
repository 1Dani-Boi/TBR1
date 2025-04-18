import arc_2 as anp
import openmc
import numpy as np
import os
import sys
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import subprocess
import pandas as pd
import shutil
import matplotlib.colors as mcolors
from openmc_plasma_source import TokamakSource
import arc_2.materials_2 as mat

# print('HERE:', dir(mat) )
# print(mat.vcrti)
# print(openmc.Material(mat.vcrti))
openmc.config["cross_sections"] = '/home/hice1/dcox67/endfb-viii.0-hdf5/cross_sections.xml'

# ==============================================================================
# Geometry
# ==============================================================================
def create_arc(Li6_enrichment, dopant, dopant_mass, multiplier_material, multiplier_thickness, reflector_material, reflector_thickness, channel_thickness, order):
    device = anp.generate_device(Li6_enrichment = Li6_enrichment, dopant = dopant, dopant_mass = dopant_mass,
                                multiplier_material = multiplier_material, multiplier_thickness = multiplier_thickness, reflector_material = reflector_material,
                                reflector_thickness = reflector_thickness, channel_thickness = channel_thickness, order = order)
    
    # Plotting
    plot = openmc.Plot()
    plot.filename = 'geometry_plot'
    plot.basis = 'xz'
    plot.origin = (350, 0, 0)
    plot.width = (700, 800)
    plot.pixels = (plot.width[0]*10, plot.width[1]*10)
    plot.color_by = 'cell'

    # plot.cell_colors = { device._cells[0].id: 'beige', device._cells[1].id: 'lightcoral', device._cells[2].id: 'yellow', device._cells[3].id: 'orange', device._cells[4].id: 'black', device._cells[5].id: 'navy', device._cells[6].id: 'red', device._cells[7].id: 'black'}
    color = ['black', 'black', 'yellow', 'orange', 'red', 'navy', 'lightcyan', 'cyan', 'pink', 'blue', 'green']
    color_dict = {cell.id: color[i] for i, cell in enumerate(device._cells)}
    plot.colors = color_dict
    for count, cell in enumerate(device._cells):
        print(f'Device name: {cell.name} with color: {color[count]}')
    
    plot.highlight_domains(geometry=device.geometry, domains=device._cells)
    
    plots = openmc.Plots([plot])
    # plot.plot_geometry(geometry, colors=color_dict)

    plots.export_to_xml()
    #openmc.plot_geometry()
    # ==============================================================================
    # Settings
    # ==============================================================================
    
    """ Source Definition """
    source = openmc.Source()
    source.space = openmc.stats.CylindricalIndependent(openmc.stats.Discrete(400, 1), openmc.stats.Uniform(a=-np.pi/18, b=np.pi/18), openmc.stats.Discrete(0, 1)) # original openmc.stats.Discrete(450, 1), openmc.stats.Uniform(a=-np.pi/18, b=np.pi/18)
    source.angles = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([14.1E6], [1.0])
    
    device.settings.source = source


    # source = TokamakSource(

    # elongation=1.5, # paper says 1.85

    # ion_density_centre=1.8e20,# paper

    # ion_density_peaking_factor=1,# source paper

    # ion_density_pedestal=1.8e20,# paper has these as the same n_p = n_0

    # ion_density_separatrix=4.95e19,# ratiod

    # ion_temperature_centre=27e3,# paper

    # ion_temperature_peaking_factor=4.74e3, #ratiod

    # ion_temperature_pedestal=3.58e3, #ratiod

    # ion_temperature_separatrix=0.1e3, # kept same as demo

    # major_radius=400,#3.3 in paper

    # minor_radius=120, #1.15 in paper

    # pedestal_radius=0.8 * 120, # jball

    # mode="H",

    # shafranov_factor=0.24789,# good enough to divide by 2

    # triangularity=0.5,# jball

    # ion_temperature_beta=2,#bruh moment

    # sample_size=1000, # the number of individual sources to use to make a combined source

    # angles=( 0 , np.pi/18) # angle in radians

    # ).make_openmc_sources()



    # device.settings.source = source



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

    r_grid2 = np.linspace(500, 700, num=100) #[NEW] original: (25, 200, num=25), 0, 600
    z_grid2 = np.array([-10,10]) #[NEW] original: (-200, 200, num=50), -700, 700
    mesh2 = openmc.CylindricalMesh(r_grid=r_grid2, z_grid=z_grid2) #[NEW]
    mesh2.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
    mesh2_filter = openmc.MeshFilter(mesh2)

    detector_points = [(520.0, 0.0)]
    mesh1 = openmc.Mesh()
    mesh1.dimension = [50, 50]
    mesh1.lower_left = [400, -100]
    mesh1.upper_right = [600, 100]
    mesh1_filter = openmc.MeshFilter(mesh1)





    # point_filters = [openmc.CellFilter(openmc.Position(*point)) for point in detector_points]

    # #creates an empty tally object
    # my_tallies = openmc.Tallies()

    # # sets up filters for the tallies
    # neutron_particle_filter = openmc.ParticleFilter(['neutron'])

    # # creates an array of energy bins to use for the tally
    # # these are not linearly spaced as they have extra bins in key energy ranges
    # # A full list of energy structures is available here
    # # https://github.com/openmc-dev/openmc/blob/6254be37582e09acff038f5656332b89e53e4eae/openmc/mgxs/__init__.py#L50-L420
    # energy_filter = openmc.EnergyFilter.from_group_structure('CCFE-709')

    # # setup the filters for the cell tally
    # inner_surface_filter = openmc.SurfaceFilter(inner_surface)
    # mid_surface_filter = openmc.SurfaceFilter(mid_surface) 

    # # in openmc a surface current includes positive contributions from neutrons
    # # passing in one direction and negative contributions from neutrons passing
    # # in the other direction. A openmc.CellFromFilter can be used to get
    # # contributions from a single direction.

    # # create the tally
    # inner_surface_spectra_tally = openmc.Tally(name='inner_surface_spectra_tally')
    # inner_surface_spectra_tally.scores = ['current']
    # inner_surface_spectra_tally.filters = [inner_surface_filter, neutron_particle_filter, energy_filter]
    # my_tallies.append(inner_surface_spectra_tally)

    # mid_surface_spectra_tally = openmc.Tally(name='mid_surface_spectra_tally')
    # mid_surface_spectra_tally.scores = ['current']
    # mid_surface_spectra_tally.filters = [mid_surface_filter, neutron_particle_filter, energy_filter]
    # my_tallies.append(mid_surface_spectra_tally)

    # # open the results file
    # results = openmc.StatePoint(results_filename)

    # #extracts the tally values from the simulation results
    # inner_surface_tally = results.get_tally(name='inner_surface_spectra_tally')
    # mid_surface_tally = results.get_tally(name='mid_surface_spectra_tally')

    # # these are the widths of each energy bin (energy bins vary in size to get detail in specific areas of the spectrum)
    # bin_boundaries = energy_filter.lethargy_bin_width

    # inner_current = inner_surface_tally.mean.flatten()
    # mid_current = mid_surface_tally.mean.flatten()

    # normalised_inner_current = inner_current / bin_boundaries
    # normalised_mid_current = mid_current / bin_boundaries

    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.step(energy_filter.values[:-1], normalised_mid_current, label='mid surface')
    # plt.step(energy_filter.values[:-1], normalised_inner_current, label='inner surface')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.legend()
    # plt.ylabel('Neutron current per unit lethargy')
    # plt.xlabel('Neutron Energy [eV]')
    # plt.show()




    # mesh = openmc.CylindricalMesh.from_domain(
    #     domain=my_geometry, # the corners of the mesh are being set automatically to surround the geometry
    #     dimension=[10, 20, 30] # 100 voxels in each axis direction (r, z, phi)
    # )
    # my_tallies = openmc.Tallies()
    # # Create mesh filter for tally
    # mesh_filter = openmc.MeshFilter(mesh)

    # # Create mesh tally to score tritium production
    # mesh_tally_1 = openmc.Tally(name='tbr_on_mesh')
    # mesh_tally_1.filters = [mesh_filter]
    # mesh_tally_1.scores = ['(n,Xt)']  # where X is a wildcard
    # my_tallies.append(mesh_tally_1)

    # # Create mesh tally to score heating
    # mesh_tally_2 = openmc.Tally(name='heating_on_mesh')
    # mesh_tally_2.filters = [mesh_filter]
    # mesh_tally_2.scores = ['heating']
    # my_tallies.append(mesh_tally_2)

    # # loads up the output file from the simulation
    # statepoint = openmc.StatePoint(sp_filename)

    # # extracts the mesh tally by name
    # my_tbr_tally = statepoint.get_tally(name='tbr_on_mesh')

    # # converts the tally result into a VTK file
    # mesh.write_data_to_vtk(
    #     filename="tbr_tally_on_cy_mesh.vtk",
    #     datasets={"mean": my_tbr_tally.mean}  # the first "mean" is the name of the data set label inside the vtk file
    # )



    
    cell_energy_filter = openmc.EnergyFilter.from_group_structure('CCFE-709')
    #cell_bin_boundaries = energy_filter.lethargy_bin_width



    energy_filter = openmc.EnergyFilter(np.array([1,10.0e6,14.0e6,20.0e6]))   #np.linspace(1e-1,20e6,50)

    energy1_bins = np.array([1,10.0e6,14.0e6,20.0e6])
    energy1_filter = openmc.EnergyFilter(energy1_bins)

    neutron_particle_filter = openmc.ParticleFilter(['neutron'])



    device.add_tally('Mesh Tally', ['flux'], filters=[mesh_filter, energy1_filter, neutron_particle_filter])
    device.add_tally('1D Flux', ['flux'], filters=[mesh2_filter, energy_filter, neutron_particle_filter])

    tbr_filter3 = openmc.CellFilter(device.get_cell(name = 'channel'))
    device.add_tally('TBR Channel Tally', ['(n,Xt)'], nuclides = ['Li6', 'Li7'], filters = [tbr_filter3])
    device.add_tally('Channel Flux Spectrum', ['flux'], filters = [tbr_filter3, cell_energy_filter, neutron_particle_filter])

    tbr_filter5 = openmc.CellFilter(device.get_cell(name = 'blanket'))
    device.add_tally('TBR Blanket Tally', ['(n,Xt)'], nuclides = ['Li6', 'Li7'], filters = [tbr_filter5])
    device.add_tally('Blanket Flux Spectrum', ['flux'], filters = [tbr_filter5, cell_energy_filter, neutron_particle_filter])



    # device.pfc = openmc.Cell(region=pfc, fill=tungsten, name='PFC')
    # device.vv = openmc.Cell(region=vv, fill=vcrti_VV, name='VV')
    # device.channel = openmc.Cell(region=channel, fill=doped_flibe_channels, name='channel')
    # device.multiplier = openmc.Cell(region = multiplier, fill = multiplier_material, name = 'multiplier')
    # device.tank_inner = openmc.Cell(region=tank_inner, fill=vcrti_BI, name='tank inner')

    #vcrti2121 = vcrti.clone()
    #vcrti2122 = openmc.Material(material_id=4)
    #print(openmc.Material(material_id=4))
    #print(vcrti2121)
    #print('Skippledi', vcrti2122)
    #vcrti_skibidi = mat.vcrti.clone()
    print(mat.vcrti)
    vv_abs_filter = openmc.MaterialFilter(mat.vcrti.clone())
    #device.add_tally('VCRTI Absorption', ['absorption'], filters = [vv_abs_filter])

    

    #device.add_tally('Flux Mesh High E Resolution', ['flux'], filters=[mesh_filter, energy_filter, neutron_particle_filter] )
    device.add_tally('Heating Local', ['heating-local'], filters=[mesh_filter])

    # region_dict = device.region_dict


    # # Create a tally for each surface
    # for name, region in region_dict.items():
    #     print("About to Flux Tally Add On Foenem")
    #     # Get all surfaces that define the boundary of the region
    #     surfaces = region.get_surfaces()  # Returns dict of {surface_id: surface}

    #     # Combine all surfaces into a single SurfaceFilter
    #     surface_filter = openmc.SurfaceFilter(list(surfaces.values()))
    #     device.add_tally(f'Flux on {name}', ['flux'], filters=[surface_filter, energy1_filter, neutron_particle_filter] )
    #     print("Flux Tally Added On Foenem")
    


    # ==============================================================================
    # Run
    # ==============================================================================
    
    device.settings.photon_transport = True
    device.survival_biasing = True
    device.build()
    device.export_to_xml(remove_surfs=True)
    
    geometry_plot = "geometry_plot.png"
    if not os.path.exists(geometry_plot):
         openmc.plot_geometry() #path_input = 'plots.xml'
    


    # all_cells_dict = device.geometry.get_all_cells()
    # print("ON FOENEM:", all_cells_dict)

    # regionpfc = device.pfc.region
    # surfacespfc = regionpfc.get_surfaces()
    # print("SURFACCES:", list(surfacespfc.values()))
    # surface_filter_pfc = openmc.SurfaceFilter(list(surfacespfc.values()))
    # device.add_tally('PFC Current', ['current'], filters=[surface_filter_pfc, energy1_filter, neutron_particle_filter] )
    


    # set run parameters
    # device.settings.threads = 24
    device.settings.particles = int(1e3)
    device.settings.batches = 10  
    device.settings.inactive = 1  

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
    sp_filename = device.run()  # runs with reduced amount of output printing, output = false

    # OPEN OUPUT FILE
    sp = openmc.StatePoint(sp_filename)

    # mesh_tally_2 = sp.get_tally(name='Mesh Tally 2')
    tbr_tally = sp.get_tally(name='TBR Blanket Tally')
    tbr_tally_2 = sp.get_tally(name='TBR Channel Tally')
    mesh_tally = sp.get_tally(name='Mesh Tally')
    heat_tally = sp.get_tally(name='Heating Local')
    flux_tally = sp.get_tally(name='1D Flux')
    #abs_tally = sp.get_tally(name='VCRTI Absorption')
    
    channel_tally = sp.get_tally(name='Channel Flux Spectrum')
    blanket_tally = sp.get_tally(name='Blanket Flux Spectrum')
    cell_energy_filter = openmc.EnergyFilter.from_group_structure('CCFE-709')
    # pfc_tally = sp.get_tally(name='PFC Current')
    # E_tally = sp.get_tally(name='Flux Mesh High E Resolution')

    # # mesh_df = mesh_tally_2.get_pandas_dataframe()
    mesh_df = mesh_tally.get_pandas_dataframe()
    df = tbr_tally.get_pandas_dataframe()
    df2 = tbr_tally_2.get_pandas_dataframe()
    heat_df = heat_tally.get_pandas_dataframe()
    flux_df = flux_tally.get_pandas_dataframe()
    #abs_df = abs_tally.get_pandas_dataframe()
    #print(abs_df.columns)

    # channel_df = channel_tally.get_pandas_dataframe()
    # blanket_df = blanket_tally.get_pandas_dataframe()
    # E_df = E_tally.get_pandas_dataframe()
    # pfc_df = pfc_tally.get_pandas_dataframe()

    reactor_power = 500e6 # in Watt  (from JBall)
    neutrons_per_J = 3.546e11 # 1MJ/17.6 Mev  (from JBall)
    neutrons_per_second = reactor_power / neutrons_per_J  # [Watt / Joule] = [1/s]

    #mesh_df.to_csv('BigMesh.csv', index=False)
    #print(mesh_df.columns)
    mesh_df.columns = [' '.join(col).strip() for col in mesh_df.columns]
    energy_ranges = mesh_df[['energy low [eV]', 'energy high [eV]']].drop_duplicates().values

    #openmc.plot_geometry()
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

    # Create plots for each energy range
    for energy_low, energy_high in energy_ranges:
        # Filter data for the specific energy range
        subset = mesh_df[(mesh_df['energy low [eV]'] == energy_low) & (mesh_df['energy high [eV]'] == energy_high)]

        # Group by x and y, then take the mean in case of duplicates
        grouped = subset.groupby(['mesh 1 z', 'mesh 1 x'])['mean'].mean().unstack()
        X, Z = np.meshgrid(x_values, z_values)

        # Plot
        fig, ax = plt.subplots(figsize=(12,12), dpi=300)
        #im = ax.imshow(grouped, origin='lower', cmap='inferno', aspect='auto')
        im = ax.pcolormesh(X,Z,grouped.values, shading='auto',cmap='plasma')
        # Add a colorbar for the image
        cbar = fig.colorbar(im, ax=ax, label="Mean Flux")

        # Set axis labels and title
        ax.set_xlabel("r")
        ax.set_ylabel("z")
        ax.set_title(f"Mean Flux for Energy Range {energy_low:.2e} eV - {energy_high:.2e} eV")

        #ax.axvline(x=520, color='cyan', linestyle='--', linewidth=2)
        ax.scatter(x_overlay, z_overlay, color='cyan', marker='o',  s=5, alpha=0.5)
        ax.plot(x_overlay, z_overlay, color='cyan', alpha=0.5, linewidth=1, label='PFC')

        ax.scatter(x_inner, z_inner, color='green', marker='o', s=5, alpha=0.5)
        ax.plot(x_inner, z_inner, color='green', alpha=0.5, linewidth=1, label="Inner Channel")

        ax.scatter(x_outer, z_outer, color='blue', marker='o', s=5, alpha=0.6)
        ax.plot(x_inner, z_inner, color='blue', alpha=0.6, linewidth=1, label="Outer Channel")

        ax.scatter(x_reflector, z_reflector, color='lime', s=5, alpha=0.3)
        ax.plot(x_reflector, z_reflector, color='lime', alpha=0.3, linewidth=1, label="Reflector")

        ax.legend()

        # Save the figure
        fig.savefig(f"mean_flux_{energy_low:.2e}_to_{energy_high:.2e}_cylsource2.png")


    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(25, 8), dpi=600)  # 1 row, 5 columns
    axes = axes.flatten()  # Ensure it's a flat iterable list

    for i, (energy_low, energy_high) in enumerate(energy_ranges):
        # Filter data for the specific energy range
        subset = mesh_df[(mesh_df['energy low [eV]'] == energy_low) & (mesh_df['energy high [eV]'] == energy_high)]

        # Group by x and y, then take the mean in case of duplicates
        grouped = subset.groupby(['mesh 1 z', 'mesh 1 x'])['mean'].mean().unstack()
        X, Z = np.meshgrid(x_values, z_values)

        norm = mcolors.LogNorm(vmin=grouped.values.min(), vmax=grouped.values.max())

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
    fig.savefig("mean_flux_all_energy_ranges3.png", dpi=600)


    # PLOT FLUX ENERGY SPECTRA (ONE FOR EACH ENERGY BIN): ----------------------------------------------------------------------------------------------------
     
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

   


    channel_flux = channel_tally.mean.flatten()
    blanket_flux = blanket_tally.mean.flatten()

    vol_calc_load = openmc.VolumeCalculation.from_hdf5('/home/hice1/dcox67/TBR/data/arc-1_volumes.h5')
    blanket_volume = vol_calc_load.volumes[8].n
    channel_volume = vol_calc_load.volumes[5].n

    channel_flux_adjusted = (channel_flux / channel_volume) * neutrons_per_second 
    blanket_flux_adjusted = (blanket_flux / blanket_volume) * neutrons_per_second 

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







    


    #heat_df.columns = [' '.join(col).strip() for col in mesh_df.columns]
    #grouped = subset.groupby(['mesh 1 z', 'mesh 1 x'])['mean'].mean().unstack()
    #fig, ax = plt.subplots()
    #im = ax.imshow(grouped, origin='lower', cmap='inferno', aspect='auto')
    #ax.scatter(x_overlay, z_overlay, color='cyan', marker='o', label="Overlay Points", s=5, alpha=0.5)
    #ax.plot(x_overlay, z_overlay, color='cyan', alpha=0.5, linewidth=1)

    #ax.scatter(x_reflector, z_reflector, color='lime', s=5, alpha=0.3, label="Reflector Points")
    #ax.plot(x_reflector, z_reflector, color='lime', alpha=0.3, linewidth=1)

    # Add a colorbar for the image
    #cbar = fig.colorbar(im, ax=ax, label="Mean Local Heat")

    # Set axis labels and title
    #ax.set_xlabel("r")
    #ax.set_ylabel("z")
    #ax.set_title("Mean Heat Local")

    # Save the figure
    #fig.savefig("Heat_SLURMED.png")




    #print(mesh1_df)
    #print(mesh1_df.head(20))

    # print(pfc_df.head(10))
    sum_means = []
    sum_std = []
    energy_low = []
    energy_high = []

    # energy_filter = openmc.EnergyFilter.from_group_structure('CCFE-709')
    # energy_bins = energy_filter.bins

    # # Loop over each starting index (0, 1, 2, 3, 4)
    # for start in range(711):
    #     # Select every 5th element starting from `start`
    #     subset = E_df.iloc[start::711]
        
    #     # Sum the 'mean' and 'std. dev.' columns for this subset
    #     sum_means.append(subset['mean'].sum())
    #     sum_std.append(subset['std. dev.'].sum())

    #     energy_low.append( E_df['energy low [eV]'][start] )
    #     energy_high.append( E_df['energy high [eV]'][start] )
    
    # plt.figure()
    # for i in range(len(sum_means)):
    #     plt.hlines(sum_means[i], energy_low[i], energy_high[i], color='b', linewidth=3)

    # # Add labels, grid, and title
    # plt.xlabel('Energy (eV)', fontsize=12)
    # plt.ylabel('Sum of Means', fontsize=12)
    # plt.title('Step Graph of Sum of Means vs Energy Range', fontsize=14)
    # plt.grid(True)
    # plt.xscale('log')
    # plt.savefig("onjah.png")


    # # Get the tally by ID or name
    # tally = sp.get_tally(name="Mesh Tally")

    # # Extract mesh filter
    # mesh_filter = tally.find_filter(openmc.MeshFilter)
    # mesh = mesh_filter.mesh

    # # Extract energy filter and get index of energy of interest
    # energy_filter = tally.find_filter(openmc.EnergyFilter)
    # energy_bins = energy_filter.bins
    # energy_idx = 4  # Choose appropriate energy index

    # # Get tally data
    # data = tally.get_reshaped_data()
    # print(data)
    # num_mesh_bins = np.prod(mesh.dimension)
    # num_energy_bins = len(energy_bins) - 1  # Energy bins are given as bin edges

    # # Ensure the data shape matches expectation
    # assert data.shape == (num_mesh_bins, num_energy_bins), f"Unexpected data shape: {data.shape}"

    # # Extract the data for the desired energy bin
    # data = data[:, energy_idx]

    # # Reshape to match the mesh (r, z) dimensions
    # r_bins, z_bins = mesh.dimension
    # data = data.reshape((r_bins, z_bins))

    # # Plot
    # fig, ax = plt.subplots()
    # c = ax.imshow(data.T, origin='lower', extent=[mesh.lower_left[0], mesh.upper_right[0],
    #                                             mesh.lower_left[1], mesh.upper_right[1]],
    #             aspect='auto', cmap='inferno')

    # plt.colorbar(c, label="Tally Value")
    # ax.set_xlabel("r [cm]")
    # ax.set_ylabel("z [cm]")
    # ax.set_title("R-Z Mesh Tally at {} MeV".format(energy_bins[energy_idx]))

    # # Save the figure instead of displaying it
    # plt.savefig("mesh_tally.png")

    
    #bin_boundaries = np.logspace(np.log10(1e0), np.log10(20e6), num=6)
    #normalised_current = current / bin_boundaries

    # print("Flux detector: ", mesh1_df)
    
    # mesh_df.to_csv('mesh.csv', index = False)
    # print(mesh_df.columns)
    # skibidi_df = mesh_df
    # skibidi_df.columns = ['_'.join(filter(None, col)).strip() for col in skibidi_df.columns]
    # filtered_df = skibidi_df[(skibidi_df['mesh 1_x'] == 20) & (skibidi_df['mesh 1_z'] == 28)]
    # skibidi = filtered_df
    # x = []
    # y = []

    # for idx, row in skibidi.iterrows():
    #     x.extend([row['energy low [eV]'], row['energy high [eV]']])
    #     y.extend([row['mean'], row['mean']])

    # # Plotting the step plot
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.step(x, y, where='post', label='Mean Values', color='blue')
    # plt.xlabel('Energy')
    # plt.ylabel('Mean')
    # plt.title('Step Plot of Mean Values for Energy Ranges')
    # plt.grid(True)
    # plt.legend()
    # plt.savefig('skippledi.png')
    # print('FIGURE OUTPUT!')

    #print('dataframe result:', df)
    #print('dataframe result 2', df2)
    #print(mesh_df.head(200))

    combined_df = pd.concat([df, df2], ignore_index = True)
    combined_df.to_csv(f'{device.Li6_enrichment}%.csv', index = False)
    tbr_tally_result = df['mean'].sum() + df2['mean'].sum()
    tbr_tally_std_dev = df['std. dev.'].sum() + df2['mean'].sum()
    tbr_from_channel = df2['mean'].sum()
    tbr_from_blanket = df['mean'].sum()

    #absorption = abs_df['mean'].sum()


    print('Total TBR: ', tbr_tally_result)
    print('Channel TBR: ', tbr_from_channel)
    print('Blanket TBR: ', tbr_from_blanket)
    #print('Absorption: ', absorption)

    # these are the widths of each energy bin (energy bins vary in size to get detail in specific areas of the spectrum)
    # energy_filter = openmc.EnergyFilter.from_group_structure('CCFE-709')
    # bin_boundaries = energy_filter.lethargy_bin_width

    # print(mesh_tally_2)
    # current = mesh_tally_2.mean.flatten()
    # print(np.shape(current))
    # normalized_current = current / bin_boundaries




    # mesh_tally_result = mesh_df['mean'].sum() / 

    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.step(energy_filter.values[:-1], normalized_current, label='bruh moment')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.legend()
    # plt.ylabel('Neutron current per unit lethargy')
    # plt.xlabel('Neutron Energy [eV]')
    # plt.savefig('skibidi.png')

    # command = ["openmc-plot-mesh-tally", sp_filename]
    # # Run the command
    # subprocess.run(command)
    
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
# print(str(sys.argv[9]))









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
        # Create the full target file path
        target_file = os.path.join(target_dir, filename)
        
        # Move the file
        try:
            shutil.move(source_file, target_file)
            print(f"Moved '{filename}' to '{target_dir}'")
        except Exception as e:
            print(f"Error moving '{filename}': {e}")




import h5py

# Path to the .h5 file (replace with the actual path to your file)
h5_file_path = "../../../endfb-viii.0-hdf5/neutron/V51.h5"

# Open the HDF5 file
with h5py.File(h5_file_path, 'r') as h5file:
    # List all datasets in the file
    datasets = list(h5file.keys())
    print(f"Datasets in the file: {datasets}")

    # Look for temperature datasets
    if 'kTs' in h5file:
        temperatures = h5file['kTs']
        print("Available temperatures (in MeV):")
        for temp in temperatures:
            print(f"{temp} MeV")

        # Convert temperatures from MeV to Kelvin
        kelvin_temperatures = [float(temp) * 11604.525 for temp in temperatures]
        print("\nAvailable temperatures (in Kelvin):")
        for temp in kelvin_temperatures:
            print(f"{temp:.2f} K")
    else:
        print("No temperature data ('kTs') found in this file.")




# os.mkdir(str(sys.argv[9]))
# move_files('/home/hice1/dcox67/TBR/scripts', str(sys.argv[9]))
# print("OpenMC files moved to new directory:", str(sys.argv[9]))
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