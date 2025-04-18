import numpy as np
import openmc
from arc_2.components_2 import *
from arc_2.materials_2 import *
from arc_2.constants_2 import neutrons_per_MJ
import shutil

"""
Class for simulating an ARC-class reactor blanket doped with either
U-238 or Th-233 for the purpose of breeding fissile material.

Instances of the class used for simulations are generated using the 
generate_device function which allows for a device instance to be generated 
with a given Li-6 enrichment and fertile inventory.

This class also stores as a class variable its fusion power and total neutron 
source rate so that these values can be easily accessed by both depletion scripts
and post processing scripts

"""
class Device(openmc.model.Model):

    fusion_power = 500 # MW
    neutron_source_rate = fusion_power * neutrons_per_MJ

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

        self._cells = []
        self._components = []
        self._tallies = []

        self.materials = []

        self.dopant = None
        self.dopant_mass = 0

        self.name = 'Untitled'

        # Define simulation boundary
        self.boundary = openmc.Sphere(r=1000, boundary_type='vacuum')
        self.domain = openmc.Cell(region=-self.boundary, name='domain')

        self._init_settings()

    def __setattr__(self, name, value):
        # Override __setattr__ to add components to the internal list
        if isinstance(value, Component):
            self._components.append(value)
        elif isinstance(value, openmc.Cell):
            self._cells.append(value)
        super().__setattr__(name, value)

    def _init_settings(self):
        """Initializes a settings object and sets it to fixed source mode"""
        settings = openmc.Settings()
        settings.run_mode = 'fixed source'
        self.settings = settings
        settings.temperature = {
            'default': 900.0,
            'method': 'interpolation'  # Ensure interpolation is enabled if needed
        }

    def build(self):
        """
        Builds the model geometry from specifically named components.
        Must be called before model can be run.
        """
        comp_cells = [cell for comp in self._components for cell in comp.cells]
        comp_regs = [~comp for comp in self._components if comp.exclude]

        #Exclude areas which have been defined as being components from total sim domain
        new_domain_reg = self.domain.region & openmc.Intersection(comp_regs)
        self.domain.region = new_domain_reg

        all_cells = comp_cells + self._cells
        self.univ = openmc.Universe(cells=all_cells)

        self.geometry = openmc.Geometry(self.univ)
        self.tallies = openmc.Tallies(self._tallies)

        if self.materials is []:
            self.materials = openmc.Materials(self.geometry.get_all_materials().values())
        
        super().export_to_xml()

    def add_tally(self, name, scores, nuclides=None, filters=None):
        """
        Creates a tally from given kwargs and adds it to the tallies object

        Parameters
        ----------
        name : str
            Name to be used for the tally
        scores : str
            tally scores, e. g. flux, absorption, (n,Xt)
        nuclides : list of str
            nuclides to individually track
        filters: list of openmc.Filter
                openmc filters like openmc.CellFilter to include in the tally

        Returns
        -------
        openmc.Tally, the generated tally object
        """
        tally = openmc.Tally(name=name)
        
        if nuclides is not None:
            tally.nuclides = nuclides

        tally.filters = filters
        tally.scores = scores

        self._tallies.append(tally)

        return tally

    def get_cell(self, name):
        comp_cells = [cell for comp in self._components for cell in comp.cells]
        all_cells = comp_cells + self._cells

        for cell in all_cells:
            if cell.name == name:
                return cell
        
        print("WARNING: Cell with name:", name, "not found. returning None.")
        return None

def generate_device(Li6_enrichment = 7.5, dopant = "Li4SiO4", dopant_mass = 0, multiplier_material = beryllium, multiplier_thickness = 0, reflector_material = lead, reflector_thickness = 0, channel_thickness = 2, order = 1, vv_file = 'arc_vv.txt', blanket_file = "arc_blanket.txt", dopant_mass_units = "wppm"):
    """
    Generates a device object with specified fertile inventory and Li-6 enrichment.

    Parameters
    ----------
    dopant : str
        "U" for U-238 -> Pu-239, "Th" for Th-232 -> U233
    dopant_mass : float
        mass of fertile material to dope into the blanket in units of kg
    Li6_enrichment : float
        The percent of the lithium which is Li-6 instead of Li-7.
    vv_file : str
        name of file in the API data folder which contains the set of points which define the vv boundary
    blanket_file : str
        name of file in the API data folder which contains the set of points which define the blanket tank boundary

    Returns
    -------
    arc_nonproliferation.Device, the generated device object.
    """
    #inputs
    if multiplier_material == 'beryllium':
        multiplier_material = beryllium
    else:
        multiplier_material = lead
    reflector_material = lead

    device = Device()
    device.dopant = dopant
    device.Li6_enrichment = Li6_enrichment

    # print(Li6_enrichment)
    # print(dopant)
    # print(dopant_mass)
    # print(multiplier_material)
    # print(type(multiplier_material))
    # print(multiplier_thickness)
    # print(reflector_material)
    # print(reflector_thickness)
    # print(channel_thickness)
    # print(order)
    # ==============================================================================
    # Geometry
    # ==============================================================================
    if int(order) == 1:    
        """ PFCs and Vacuum Vessel """

        vv_points = np.loadtxt("/home/hice1/dcox67/TBR/data/" + vv_file)

        r_mid = np.mean(vv_points[:, 0])  # Average of all r values
        z_mid = np.mean(vv_points[:, 1])  # Average of all z values
        print(f"Midpoint of Vacuum Vessel (r, z): ({r_mid:.2f}, {z_mid:.2f})")

        pfc_polygon = openmc.model.Polygon(vv_points, basis='rz')
        vv_inner_edge = pfc_polygon.offset(0.3) #other lit says 0.1cm [PFC]
        vv_channel_inner = vv_inner_edge.offset(1.0) #VV [STR1]
        channel_outer = vv_channel_inner.offset(channel_thickness) #FLiBe channels [FLiBe1]
        multiplier_outer = channel_outer.offset(multiplier_thickness) #[*]Be multiplier [Be]
        vv_channel_outer = multiplier_outer.offset(3.0) #Channel shell [STR2]

        vertices_in = vv_channel_inner.points
        print("channel inner vertices:", vertices_in)

        vertices_out = channel_outer.points
        print("channel outer vertices:", vertices_out)

        # Make surfaces for surface tally:
        # pfc_surface = pfc_polygon.boundary  # Original VV surface
        # vv_inner_surface = vv_inner_edge.boundary
        # vv_channel_inner_surface = vv_channel_inner.boundary
        # channel_outer_surface = channel_outer.boundary
        # multiplier_outer_surface = multiplier_outer.boundary
        # vv_channel_outer_surface = vv_channel_outer.boundary
        
        # surfaces = [
        #     (pfc_surface, "PFC Surface"),
        #     (vv_inner_surface, "VV Inner Surface"),
        #     (vv_channel_inner_surface, "VV Channel Inner Surface"),
        #     (channel_outer_surface, "Channel Outer Surface"),
        #     (multiplier_outer_surface, "Multiplier Outer Surface"),
        #     (vv_channel_outer_surface, "VV Channel Outer Surface")
        # ]   
        # print( list(surfaces.values()) )
        # device.surfaces = surfaces


        """ Blanket and Outer Blanket Tank """

        blanket_points = np.loadtxt("/home/hice1/dcox67/TBR/data/" + blanket_file)

        blanket_inner = openmc.model.Polygon(blanket_points, basis='rz')
        reflector_outer = blanket_inner.offset(reflector_thickness)
        vertices = reflector_outer.points
        print("vertices:", vertices)
        # gap = blanket_inner.offset(1.0)
        blanket_outer = reflector_outer.offset(3.0) #Blanket tank outer [FLiBe2] originially offset of 2 from gap

        regions = openmc.model.subdivide([pfc_polygon,
                                        vv_inner_edge, vv_channel_inner,
                                        channel_outer, multiplier_outer, vv_channel_outer,
                                        blanket_inner, reflector_outer, blanket_outer]) #[*]

        plasma, pfc, vv, channel, multiplier, tank_inner, salt, reflector, tank_outer, outside = regions #[*]
    
        # region_dict = {
        #     "Plasma": plasma,
        #     "PFC": pfc,
        #     "VV": vv,
        #     "Channel": channel,
        #     "Multiplier": multiplier,
        #     "Tank Inner": tank_inner,
        #     "Salt": salt,
        #     "Reflector": reflector,
        #     "Tank Outer": tank_outer
        # }
        # device.region_dict = region_dict
    
    
    
    # Read volume calc file
    vol_calc_load = openmc.VolumeCalculation.from_hdf5('/home/hice1/dcox67/TBR/data/arc-1_volumes.h5')
    flibe_volume = vol_calc_load.volumes[8].n
    channels_volume = vol_calc_load.volumes[5].n


    if dopant_mass_units == "wppm": #units used to be kg (switched it)
        doped_flibe = make_doped_flibe(dopant, 
                                        dopant_mass, 
                                        volume=flibe_volume + channels_volume, 
                                        Li6_enrichment=Li6_enrichment, 
                                        name="doped flibe blanket")
        
    elif dopant_mass_units == 'kg': #units used to be wppm (switched it)
        doped_flibe = make_impure_flibe(dopant_mass, 
                                        name="doped flibe blanket")

    else:
        raise ValueError("Invalid units specified for dopant mass")
        
    device.doped_flibe = doped_flibe

    doped_flibe_channels = doped_flibe.clone()
    doped_flibe_channels.volume = channels_volume
    doped_flibe_channels.name = "doped flibe channels"
    device.doped_flibe_channels = doped_flibe_channels

    doped_flibe_blanket = doped_flibe.clone()
    doped_flibe_blanket.volume = flibe_volume
    doped_flibe_blanket.name = "doped flibe blanket"
    device.doped_flibe_blanket = doped_flibe_blanket
    
    vcrti_VV = vcrti.clone()
    vcrti_VV.name = "VV"
    device.vcrti_VV = vcrti_VV

    vcrti_BI = vcrti.clone()
    vcrti_BI.name = "tank inner"
    device.vcrti_BI = vcrti_BI

    vcrti_BO = vcrti.clone()
    vcrti_BO.name = "tank outer"
    device.vcrti_BO = vcrti_BO

    device.plasma = openmc.Cell(region=plasma, fill=None, name='plasma')
    device.pfc = openmc.Cell(region=pfc, fill=tungsten, name='PFC')
    device.vv = openmc.Cell(region=vv, fill=vcrti_VV, name='VV')
    device.channel = openmc.Cell(region=channel, fill=doped_flibe_channels, name='channel')
    device.multiplier = openmc.Cell(region = multiplier, fill = multiplier_material, name = 'multiplier')
    device.tank_inner = openmc.Cell(region=tank_inner, fill=vcrti_BI, name='tank inner')
    device.blanket = openmc.Cell(region=salt, fill=doped_flibe_blanket, name='blanket')
    device.reflector = openmc.Cell(region = reflector, fill = reflector_material, name = 'reflector')
    device.tank_outer = openmc.Cell(region=tank_outer, fill=vcrti_BO, name='tank outer')
    device.domain.region = device.domain.region & outside

    # ==============================================================================
    # Settings
    # ==============================================================================

    """ Source Definition """
    # source = openmc.Source()
    # source.space = openmc.stats.CylindricalIndependent(openmc.stats.Discrete(400, 1), openmc.stats.Uniform(a=0, b=2*np.pi), openmc.stats.Discrete(0, 1))
    # source.angles = openmc.stats.Isotropic()
    # source.energy = openmc.stats.Discrete([14.1E6], [1.0])

    # device.settings.source = source

    return device
