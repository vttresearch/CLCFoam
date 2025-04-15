from os.path import join, split, isfile
import pandas as pd

from specie_properties.common import AW, MW
from foam.common import read_case_params
from foam.foamcase import FoamCase

rhon = dict( 
    Ar     = 1.633864,
    air    = 1.184318,
    CO2    = 1.807970,
    syngas = 0.61348,
)
""" Densities at normal conditions form coolprop (see mathcad)"""

Yair = dict( # Air mass fractions
    N2  = 0.7551,
    O2  = 0.2314,
    Ar  = 0.0130,
    CO2 = 0.0005,
)
""" Air mass fractions """

# Ysyngas = dict(
#     CO = 0.7523,
#     H2O = 0.1935,
#     H2 = 0.0542,
# )
# """Syngas mass fractions"""

def MW_fluid(fluid: str): # g/mol!
    if fluid == 'syngas':
        return 0.5*MW('CO')+0.5*MW('H2')
    elif fluid == 'syngas_diluted':
        return (5.0/6.0)*MWs['syngas'] + (1.0/6.0)*MW('H2O')
    elif fluid == 'air':
        return 28.946e-3
    else:
        return MW(fluid)

def Vn_to_M(fluid, Vn):
    """
    Vn - normalized volumetric flow [L/min]
    """
    # T_NTP = 273.15 + 25 # [K]  - "normalized" temperature. According to Leion2008 and Moldenhauer2009 it is 25°C.
    # p_NTP = 100000;     # [Pa] - "normalized" pressure. According to Leion2008 it is 1 bar.
    T_NTP = 273.15 + 20 # [K]  - temperature I used, difference is small
    p_NTP = 101325;     # [Pa] - pressure I used
    R     = 8.31446; # [J/K/mol]
    VM    = R*T_NTP/p_NTP # [m^3/mol] Molar volume of an ideal gas at NTP
    massFlow = (Vn/1000.0/60.0)/VM*(MW_fluid(fluid)) # [kg/s]
    return massFlow

class Pilot300W(FoamCase):
    def __init__(self, path, offset = 0, 
            species=dict(
                gas = ("N2", "O2", "CO", "CO2", "H2O", "H2", "CH4"),
                solids = ("Fe2O3", "FeO", "inerts")
            )) -> None:
        super().__init__(path)
        self.case_folder = split(self.path)[-1]
        self.offset = offset
        self.species = species # can be taken from case, but if includes are used it is complicated
        
        self.reaction_species = ['H2', 'CO', 'CH4', 'O2']
        # 1 Fe2O3 + H2  = 2 FeO + H2O
        # 1 Fe2O3 + CO  = 2 FeO + CO2
        # 4 Fe2O3 + CH4 = 8 FeO + CO2 + 2 H2O
        # 4 FeO   + O2  = 2 Fe2O3
        
        # self.patch_names = ['inlet_' +loc for loc in ['AR', 'FR', 'SL', 'DC2']] \
        #                   + ['outlet_'+loc for loc in ['AR', 'FR']]
        walls_name = 'walls' if 'walls' in self.patches.keys() else 'wall'
        self.patch_names = list(self.patches.keys()); self.patch_names.remove(walls_name)

        # self.faceZone_names = ['controlSurface_'+loc for loc in ['DC_AR', 'DC_FR', 'SL_AR', 'SL_FR']]
        self.faceZone_names = list(self.faceZones.keys())

        # self.cellSet_names = ['controlVolume_'+loc for loc in ['AR', 'FR', 'SL', 'DC']]
        self.cellSet_names = list(self.cellSets.keys()); self.cellSet_names.remove('solids_loading'); self.cellSet_names.remove('inlet_DC')

        self._read_case_params() # from caseParamsDict if reacting simuation and from ICBC 

        for patch_name in self.patch_names:
            self.patches[patch_name].register_custom_function_object('mdot_sp', 'surfaceSumSpecieFluxesPatchFn', )
            if self.sim_type == 'leakage':
                self.patches[patch_name].register_custom_function_object('mdot', 'surfaceSumFn', )

        for faceZone_name in self.faceZone_names:
            self.faceZones[faceZone_name].register_custom_function_object('mdot_sp', 'surfaceSumSpecieFluxesFaceZoneFn', )
            self.faceZones[faceZone_name].register_custom_function_object('mdot_solids_sp', 'surfaceSumSolidSpecieFluxesFaceZoneFn', )
            if self.sim_type == 'leakage':
                self.faceZones[faceZone_name].register_custom_function_object('mdot', 'surfaceSumFn', )
        for cellSet_name in self.cellSet_names:
            # TODO: make m_sp(phase='solids')
            self.cellSets[cellSet_name].register_custom_function_object('m_solid_sp', 'volIntegrateSolidSpeciesFn', )
        self.domain.register_custom_function_object('m_solid_sp', 'volIntegrateSolidSpeciesFn', )
        for reaction_specie in self.reaction_species:
            self.domain.register_custom_function_object(f'NRR_{reaction_specie}', 'volIntegrate', {'fields':f'(NRR_Abad_{reaction_specie})'})

        self.patches[walls_name].register_custom_function_object('wallHeatFlux', 'wallHeatFlux',{'patch':None, 'phase':'gas'})
        
        self.domain.register_custom_function_object('pressure', 'pressureProbesFn')

        # hardcoded, because I cannot parse yet if-else in caseParamsDict
        if self.chemistry == 'act':
            self.Ro = 0.033
        elif self.chemistry == 'pre':
            self.Ro = 0.04
        else:
            print(f'Carrier type {self.chemistry} is not "pre" or "act"')
            # exit

        # TODO wallHeatFlux
    
    def _read_case_params(self):
        if isfile(join(self.path,'caseParamsDict')):
            self.__dict__.update(
                read_case_params(self.path, 
                    ['simulationType', 'chemistry', 'solidsMass', 'domainVolume', 
                    'initialBedVolume', 'particleD', 'initialConversion', 'T', 'rhoOxidized'])
            )
            self.sim_type = 'reacting'
        elif isfile(join(self.path,'ICBC')):
            self.ICBC = self._readICBC()
            self.__dict__.update(self.ICBC)
            self.sim_type = 'leakage'
        else:
            print('Warning! Unable to read case parameters!')
        
        if self.sim_type == 'leakage':
            self.inlet_fluid_name = dict(
                inlet_AR  = 'air',
                inlet_FR  = 'Ar',
                inlet_DC1 = 'CO2',
                inlet_DC2 = 'CO2',
                inlet_SL  = 'CO2',
            )
        elif self.sim_type == 'reacting':
            self.inlet_fluid_name = dict(
                inlet_AR  = 'air',
                inlet_FR  = 'syngas',
                inlet_DC1 = 'Ar',
                inlet_DC2 = 'Ar',
                inlet_SL  = 'Ar',
            )
        else:
            raise ValueError('Wrong test_type')

    def _readICBC(self):
        with open(join(self.path,'ICBC'), 'r') as f:
            d = dict()
            for line in f:
                line = line.strip()
                if len(line) == 0 or line.startswith("#"):
                    continue
                param, value = line.split('=')
                d[param] = float(value)
        self.chemistry = 'none'
        return d
    
    def mdot(self, location):
        # TODO: note the sign!
        if location.startswith('outlet') or location.startswith('inlet'):
            if not location in self.patches.keys() or self.patches[location].mdot() is None:
                print('Warning: using inlet value from ICBC, since function object was not found')
                assert location.startswith('inlet')
                return Vn_to_M(self.inlet_fluid_name[location], self.ICBC[location])
            else:
                return self.patches[location].mdot()['sum(alphaRhoPhi.gas)']
        elif location.startswith('controlSurface'):
            return self.faceZones[location].mdot()
    
    def mdot_sp(self, location, sp):
        if location.startswith('outlet') or location.startswith('inlet'):
            return self.patches[location].mdot_sp()[f'sum({sp}.gas)']
        elif location.startswith('controlSurface'):
            return self.faceZones[location].mdot_sp()[f'sum(surfaceInterpolate({sp}.gas))']
    
    def _mdot_air(self, location):
        """Air flow under assumption that N2 is always coming from air.
        Verification:
        (p._mdot_air('outlet_AR') + p._mdot_air('outlet_AR')).plot(); plt.axhline(p.mdot('inlet_AR'))"""
        return self.mdot_sp(location, 'N2') / Yair['N2']

    def _dilution_outlet(self):
        """Calculate dilution, i.e. vol. flow of air going from the inlet of AR to 
        the outlet of FR, normalized by vol. flow at the inlet of FR"""
        fg = self.inlet_fluid_name['inlet_DC1'] # fluidization gas name
        if self.sim_type == 'leakage':
            # mdot_fluidiztion_FR = 
            #     self.mdot_sp('outlet_FR', fg) 
            #   - self._mdot_air('outlet_FR')*Yair[fg]
            Vdot_dilution = self._mdot_air('outlet_FR') / rhon['air']
            return Vdot_dilution/self.mdot('inlet_FR')
    
    def _dilution_PL_outlet(self):
        """Calculate dilution of FR stream by fluidization gases, coming from
        particle locks"""
        fg = self.inlet_fluid_name['inlet_DC1'] # fluidization gas name
        if self.sim_type == 'leakage':
            mdot_fluidiztion_FR = \
                self.mdot_sp('outlet_FR', fg) \
              - self._mdot_air('outlet_FR')*Yair[fg]
            Vdot_dilution_PL = mdot_fluidiztion_FR/rhon[fg]
            return Vdot_dilution_PL/self.mdot('inlet_FR')

    def _total_dilution_outlet(self):
        return self._dilution_outlet + self._dilution_PL_outlet

    
    def _leakage_outlet(self):
        """Calculate leakage, i.e. vol. flow of fuel going from the inlet of FR to
        the outlet of AR, normalized by vol. flow at the inlet of FR"""
        g = self.inlet_fluid_name['inlet_FR']
        if self.sim_type == 'leakage':
            mdot_fuel_AR = self.mdot_sp('outlet_AR', g.upper()) \
                         - self._mdot_air('outlet_FR')*Yair[g]
            Vdot_fuel_AR = mdot_fuel_AR/rhon[g]
            return Vdot_fuel_AR/self.mdot('inlet_FR')

    def Y(self, reactor, sp):
        """Mass fraction of specie at the outlet of a reactor, -"""
        # This kind of averaging is dangerous, since not all species might be
        # in the species['gas']
        sumY = 0.0
        for spi in self.species['gas']:
            # sumY += self._df_sp[f'average({spi}.gas)']
            sumY += self.patches['outlet_'+reactor].mdot_sp()[f'sum({spi}.gas)']
        return self.patches['outlet_'+reactor].mdot_sp()[f'sum({sp}.gas)']/sumY
    

    def X(self, reactor, sp):
        """Mole fractions of a specie at the outlet of a reactor, -"""
        sumN = 0.0
        for spi in self.species['gas']:
            sumN += self.Y(reactor, spi)/MW(spi)
        return self.Y(reactor, sp)/MW(sp)/sumN