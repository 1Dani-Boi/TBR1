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

openmc.config["cross_sections"] = '/home/hice1/dcox67/endfb-viii.0-hdf5/cross_sections.xml'

# ==============================================================================
# Geometry
# ==============================================================================
def create_arc(Li6_enrichment, dopant, dopant_mass, multiplier_material, multiplier_thickness, reflector_material, reflector_thickness, channel_thickness, order,
                gap_thickness, multiplier_2_material, multiplier_2_thickness):
    device = anp.generate_device(Li6_enrichment, dopant, dopant_mass,
                                multiplier_material, multiplier_thickness, reflector_material,
                                reflector_thickness, channel_thickness, order, gap_thickness, multiplier_2_material, multiplier_2_thickness)
    
    # Plotting
    plot = openmc.Plot()
    plot.filename = 'geometry_plot'
    plot.basis = 'xz'
    plot.origin = (350, 0, 0)
    plot.width = (700, 800)
    plot.pixels = (plot.width[0]*10, plot.width[1]*10)
    plot.color_by = 'cell'

    # plot.cell_colors = { device._cells[0].id: 'beige', device._cells[1].id: 'lightcoral', device._cells[2].id: 'yellow', device._cells[3].id: 'orange', device._cells[4].id: 'lime', device._cells[5].id: 'navy', device._cells[6].id: 'lightcyan', device._cells[7].id: 'black'}
    # color = ['beige', 'lightcoral', 'yellow', 'orange', 'lime', 'navy', 'lightcyan', 'black']
    # color_dict = {cell.id: color[i] for i, cell in enumerate(device._cells)}
    # plot.colors = color_dict
    # for count, cell in enumerate(device._cells):
    #     print(f'Device name: {cell.name} with color: {color[count]}')
    
    # plot.highlight_domains(geometry=device.geometry, domains=device._cells)
    
    plots = openmc.Plots([plot])
    plots.export_to_xml()
    
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
    r_grid = np.linspace(0, 620, num=620) #[NEW] original: (25, 200, num=25), 0, 600
    z_grid = np.array([-10,10]) #[NEW] original: (-200, 200, num=50), -700, 700
    mesh = openmc.CylindricalMesh(r_grid=r_grid, z_grid=z_grid) #[NEW]
    mesh.phi_grid = np.array([0, (2 * np.pi)/(18 * 2)])
    mesh_filter = openmc.MeshFilter(mesh)
    
    device.add_tally('Mesh Tally', ['flux', '(n,Xt)', 'heating-local', 'absorption'], filters=[mesh_filter])
    
    # """ TBR Tally """
    # tbr_filter1 = openmc.MaterialFilter(anp.tungsten)
    # device.add_tally('Tbr Plasma-facing Component Tally ', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[tbr_filter1])
    
    # tbr_filter2 = openmc.MaterialFilter(device.vcrti_VV)
    # device.add_tally('Tbr Vacuum Vessel Tally ', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[tbr_filter2])
    
    # tbr_filter3 = openmc.MaterialFilter(device.doped_flibe_channels)
    if int(order) == 1:
        tbr_filter3 = openmc.CellFilter(device.get_cell(name = 'channel'))
        device.add_tally('TBR Channel Tally', ['(n,Xt)'], nuclides = ['Li6', 'Li7'], filters = [tbr_filter3])

        tbr_filter5 = openmc.CellFilter(device.get_cell(name = 'blanket'))
        device.add_tally('TBR Blanket Tally', ['(n,Xt)'], nuclides = ['Li6', 'Li7'], filters = [tbr_filter5])

    if int(order) == 2:
        tbr_filter3 = openmc.CellFilter(device.get_cell(name = 'channel'))
        device.add_tally('TBR Channel Tally', ['(n,Xt)'], nuclides = ['Li6', 'Li7'], filters = [tbr_filter3])

        tbr_filter4 = openmc.CellFilter(device.get_cell(name = 'blanket 2a'))
        device.add_tally('TBR Blanket 2a Tally', ['(n,Xt)'], nuclides = ['Li6', 'Li7'], filters = [tbr_filter4])

        tbr_filter5 = openmc.CellFilter(device.get_cell(name = 'blanket 2b'))
        device.add_tally('TBR Blanket 2b Tally', ['(n,Xt)'], nuclides = ['Li6', 'Li7'], filters = [tbr_filter5])
    # tbr_filter4 = openmc.MaterialFilter(device.vcrti_BI)
    # device.add_tally('Tbr Tank Inner Tally ', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[tbr_filter4])

    # tbr_filter5 = openmc.MaterialFilter(device.doped_flibe_blanket) #device.doped_flibe, initially wanted 'doped_mat' which doesn't exist
    
    # tbr_filter6 = openmc.MaterialFilter(device.vcrti_BO)
    # device.add_tally('Tbr Tank Outer Tally ', ['(n,Xt)', 'fission', 'kappa-fission', 'fission-q-prompt', 'fission-q-recoverable', 'heating', 'heating-local'], filters=[tbr_filter6])
    
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
    
    # set run parameters
    # device.settings.threads = 24
    device.settings.particles = int(1e5)
    device.settings.batches = 10  
    device.settings.inactive = 1  

    return device
# ================================================================================
# additional runs
# ================================================================================

def make_materials_geometry_tallies(Li6_enrichment, dopant, dopant_mass,
                                multiplier_material, multiplier_thickness, reflector_material,
                                reflector_thickness, channel_thickness, order, gap_thickness, multiplier_2_material, multiplier_2_thickness):
    """Makes a neutronics model of a blanket and simulates the TBR value.

    Arguments:
        enrichment (float): the enrichment percentage of Li6 in the breeder material
    
    Returns:
        resutsl (dict): simulation tally results for TBR along with the standard deviation and enrichment
    """
    # RUN OPENMC
    device = create_arc(Li6_enrichment, dopant, dopant_mass, multiplier_material, multiplier_thickness, reflector_material, reflector_thickness, channel_thickness, order,
                        gap_thickness, multiplier_2_material, multiplier_2_thickness)
    #remove old output files
    for file in os.listdir('.'):
        if file.endswith('.h5'):
            os.remove(file)
    sp_filename = device.run(output = False)  # runs with reduced amount of output printing, output = false
    # OPEN OUPUT FILE
    sp = openmc.StatePoint(sp_filename)
    # export mesh tally to csv
    mesh_tally = sp.get_tally(name='Mesh Tally')
    df_mesh = mesh_tally.get_pandas_dataframe()
    df_mesh.to_csv(f'{dopant_mass}wppm_mesh.csv', index = False)
    #export tbr tallies to csv
    if int(order) == 1:
        tbr_tally_blanket = sp.get_tally(name='TBR Blanket Tally')
        tbr_tally_channel = sp.get_tally(name='TBR Channel Tally')
        df_blanket = tbr_tally_blanket.get_pandas_dataframe()
        df_channel = tbr_tally_channel.get_pandas_dataframe()
        print(df_blanket)
        print(df_channel)
        combined_df = pd.concat([df_blanket, df_channel], ignore_index = True)
        combined_df.to_csv(f'{dopant_mass}wppm.csv', index = False)
        tbr_tally_result = df_blanket['mean'].sum() + df_channel['mean'].sum()
        tbr_tally_std_dev = df_blanket['std. dev.'].sum() + df_channel['mean'].sum()
    elif int(order) == 2:
        tbr_tally_channel = sp.get_tally(name='TBR Channel Tally')
        tbr_tally_blanket_2a = sp.get_tally(name='TBR Blanket 2a Tally')
        tbr_tally_blanket_2b = sp.get_tally(name='TBR Blanket 2b Tally')        
        df_channel = tbr_tally_channel.get_pandas_dataframe()
        df_blanket_2a = tbr_tally_blanket_2a.get_pandas_dataframe()
        df_blanket_2b = tbr_tally_blanket_2b.get_pandas_dataframe()
        print(df_channel)
        print(df_blanket_2a)
        print(df_blanket_2b)
        combined_df = pd.concat([df_channel, df_blanket_2a, df_blanket_2b], ignore_index = True)
        combined_df.to_csv(f'{dopant_mass}wppm.csv', index = False)
        tbr_tally_result = df_channel['mean'].sum() + df_blanket_2a['mean'].sum() + df_blanket_2b['mean'].sum() 
        tbr_tally_std_dev = df_channel['mean'].sum() + df_blanket_2a['std. dev.'].sum() + df_blanket_2b['std. dev.'].sum()
    # command = ["openmc-plot-mesh-tally", sp_filename]
    # # Run the command
    # subprocess.run(command)
    
    return {'dopant mass': dopant_mass,
            'tbr_tally_result': tbr_tally_result,
            'tbr_tally_std_dev': tbr_tally_std_dev}

results = []
for dopant_wppm in np.arange(0, 100, 50):  # dopant wppm; start at 0wppm, end at 900wppm, step of 50wppm
    print(dopant_wppm)
    results.append(make_materials_geometry_tallies(7.5, str(sys.argv[1]), dopant_wppm, 
                    sys.argv[2], float(sys.argv[3]), sys.argv[4],
                    float(sys.argv[5]), float(sys.argv[6]), sys.argv[7], float(sys.argv[8]), sys.argv[9], float(sys.argv[10])))
print(results)
# ================================================================================
print(str(sys.argv[11]))

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

os.mkdir(str(sys.argv[11]))
move_files('/home/hice1/dcox67/skibidi/TBR/scripts', str(sys.argv[11]))
print("OpenMC files moved to new directory:", str(sys.argv[11]))
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