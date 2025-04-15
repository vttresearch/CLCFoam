from os.path import join, split, exists, isdir
import pandas as pd
import numpy as np
# import scipy.integrate as integrate
from functools import cache, lru_cache

from specie_properties.common import AW, MW
from foam.common import read_case_params, integrate_time_series
from foam.foamcase import FoamCase


class LabReactor(FoamCase):
    def __init__(self, path, offset = 0, 
                species=dict(
                    gas = ("N2", "O2", "CO", "CO2", "H2O", "H2", "CH4"),
                    solids = ("Fe2O3", "FeO", "inerts")
                )) -> None:
        super().__init__(path)
        # self.path = path
        self.case_folder = split(self.path)[-1]
        self.offset = offset
        self.species = species # can be taken from case
        self.reaction_species = ['H2', 'CO', 'CH4', 'O2']
        # 1 Fe2O3 + H2  = 2 FeO + H2O
        # 1 Fe2O3 + CO  = 2 FeO + CO2
        # 4 Fe2O3 + CH4 = 8 FeO + CO2 + 2 H2O
        # 4 FeO   + O2  = 2 Fe2O3

        self._read_case_params() # from caseParamsDict
        # hardcoded, because I cannot parse yet if-else in caseParamsDict
        if self.chemistry == 'act':
            # self.chemistry_label = 'Activated'
            self.Ro = 0.033
        elif self.chemistry == 'pre':
            # self.chemistry_label = 'Calcined'
            self.Ro = 0.04
        else:
            print(f'Carrier type {self.chemistry} is not "pre" or "act"')
            exit
        
        for patch_name in ['inlet', 'outlet']:
            self.patches[patch_name].register_custom_function_object(
                'mdot', 
                'surfaceSumFn', 
                {
                    'fields':'(alphaRhoPhi.solidsalphaRhoPhi.gas)', 
                    'select':'patch', 
                    'patch': patch_name
                })
            self.patches[patch_name].register_custom_function_object(
                'hdot_gas', 
                'surfaceSumFn', 
                {
                    'fields':'(alphaRhoPhi.gas)', 
                    'weightFields':'(h.gas)', 
                    'select':'patch', 
                    'patch': patch_name
                })
            self.patches[patch_name].register_custom_function_object(
                'hdot_solids', 
                'surfaceSumFn', 
                {
                    'fields':'(alphaRhoPhi.solids)', 
                    'weightFields':'(h.solids)', 
                    'select':'patch', 
                    'patch': patch_name
                })
            self.patches[patch_name].register_custom_function_object(
                'mdot_solids_sp', 
                'surfaceSumSolidSpecieFluxesPatchFn')
            self.patches[patch_name].register_custom_function_object(
                'mdot_sp', 
                'surfaceSumSpecieFluxesPatchFn')
        
        self.patches['outlet'].register_custom_function_object(
            'Y_sp', 
            'surfaceAverageSpecieFluxesPatchFn',
            {
                'patch': 'outlet'
            })
        self.patches['wall'].register_custom_function_object(
        'Qdot', 
        'wallHeatFlux', 
        {
            'phase': 'gas',
            'patch': None,
            'objects': None # quick fix for my setup mistakes
        })

        self.domain.register_custom_function_object(
            'Tavgest', 
            'volAverageFn', 
            {
                'weightField':'alpha.solids',
                'fields':'(T.solids)'
            })
        for phase in ['gas','solids']:
            self.domain.register_custom_function_object(
                f'm_{phase}', 
                'volIntegrateFn', 
                {
                    'weightField': f'alpha.{phase}',
                    'fields': f'(rho.{phase})'
                })
            titled = phase.title() if phase != 'solids' else 'Solid'
            self.domain.register_custom_function_object(
                f'm_{phase}_sp', 
                f'volIntegrate{titled}SpeciesFn', 
                {
                    'weightFields':f'(alpha.{phase}rho.{phase})'
                })
        for reaction_specie in self.reaction_species:
            self.domain.register_custom_function_object(
                f'NRR_{reaction_specie}', 
                'volIntegrateFn', 
                {
                    'fields':f'(NRR_Abad_{reaction_specie})'
                })
        
        self._find_corresponding_experiment()

    
    def _read_case_params(self):
        self.__dict__.update(
            read_case_params(self.path, 
                ['gas', 'chemistry', 'simTime', 'timeStep', 'solidsMass', 
                'particleD', 'initialConversion', 'T', 'nVolRate_mLn_min', 
                'kN', 'defaultSpecie', 'rhoOxidized', ])
        )
    
    def _find_corresponding_experiment(self):
        if self.gas == 'CH4':
            self.timefigname = 'fig5'
            self.eff_fig_name = 'fig8'
            self.plot_species=['CO', 'CO2', 'CH4']
        elif self.gas == 'syngas':
            self.timefigname = 'fig6'
            self.eff_fig_name = 'fig9'
            self.plot_species=['CO', 'CO2']
        elif self.gas == 'CO':
            self.timefigname = None
            self.eff_fig_name = None
            self.plot_species=['CO', 'CO2']
        elif self.gas == 'H2':
            self.timefigname = None
            self.eff_fig_name = None
            self.plot_species=['H2', 'H2O']
        elif self.gas == 'O2':
            print('O2 selected')
            self.timefigname = None # TODO: really?
            self.eff_fig_name = None
            self.plot_species=['O2']
        else:
            print('Error: case folder name should be in the form blabla_OCTYPE_FUEL')
            print('    OCTYPE: ' + self.chemistry)
            print('    FUEL: ' + self.gas)
            exit
    
    # ---------------- access functions --------------------

    def mdot_in(self):
        """Gas mass flow rate at the inlet, kg/s"""
        # return - self._df_mdot_in['sum(alphaRhoPhi.gas)']
        return - self.patches['inlet'].mdot()['sum(alphaRhoPhi.gas)']

    def mdot_out(self):
        """Gas mass flow rate at the outlet, kg/s"""
        return self.patches['outlet'].mdot()['sum(alphaRhoPhi.gas)']
        # return self._df_mdot_out['sum(alphaRhoPhi.gas)']
    
    def hdot_in(self):
        """Gas enthalpy flow rate at the inlet, J/s"""
        return self.patches['inlet'].hdot_gas()['sum(alphaRhoPhi.gas)']
        # return - self._df_hdot_in['sum(alphaRhoPhi.gas)']

    def hdot_out(self):
        """Gas enthalpy flow rate at the outlet, J/s"""
        return self.patches['outlet'].hdot_gas()['sum(alphaRhoPhi.gas)']
        # return self._df_hdot_out['sum(alphaRhoPhi.gas)']
    
    def Qdot_wall(self):
        """Heat flux through the walls, W"""
        # raise NotImplementedError
        return self.patches['wall'].Qdot()['Q[W]']
        # return self._df_Qdot_wall['Q [W]']

    def mdot_in_sp(self, specie):
        """Gas specie mass flow rate at the inlet, kg/s"""
        return -self.patches['inlet'].mdot_sp()[f'sum({specie}.gas)'] # NOTE: added for O2
        # return - self._df_mdot_sp_in[f'sum({specie}.gas)']

    def mdot_out_sp(self, specie):
        """Gas specie mass flow rate at the outlet, kg/s"""
        return self.patches['outlet'].mdot_sp()[f'sum({specie}.gas)']
        # return self._df_mdot_sp_out[f'sum({specie}.gas)']
    
    def m_solid_sp(self, specie):
        """Total mass of solid species in reactor, kg"""
        return self.domain.m_solids_sp()[f'volIntegrate(alpha.solids,rho.solids,{specie}.solids)']
        # return self._df_m_solid_sp[f'volIntegrate(alpha.solids,rho.solids,{specie}.solids)']
    
    def NRR(self, reaction_specie):
        """Net reaction rate of a gaseous specie conversion reaction, mol/s.
        4 Fe2O3 + CH4 = 8 FeO + CO2 + 2 H2O
        1 Fe2O3 + CO = 2 FeO + CO2
        1 Fe2O3 + H2 = 2 FeO + H2O
        4 FeO + O2 = 2 Fe2O3

        Note that OpenFOAM uses kmol
        https://github.com/OpenFOAM/OpenFOAM-dev/blob/4d77261ddf27e756d2a49cc57bd9fa28d0576a24/etc/controlDict#L187C9-L187C43
        """
        # func = getattr(self, f'NRR_{reaction_specie}') # doesn't work for some reason
        # return 1000*func()
        if reaction_specie == 'CH4':
            return 1000*self.domain.NRR_CH4()[f'volIntegrate(NRR_Abad_{reaction_specie})']
        elif reaction_specie == 'H2':
            return 1000*self.domain.NRR_H2()[f'volIntegrate(NRR_Abad_{reaction_specie})']
        elif reaction_specie == 'CO':
            return 1000*self.domain.NRR_CO()[f'volIntegrate(NRR_Abad_{reaction_specie})']
        elif reaction_specie == 'O2':
            return 1000*self.domain.NRR_O2()[f'volIntegrate(NRR_Abad_{reaction_specie})']
        # return 1000*self._dfs_NRR[reaction_specie][f'volIntegrate(NRR_Abad_{reaction_specie})']

    # # TODO: remove?
    # def time(self):
    #     return self._df_sp[self._df_sp.columns[0]]

    def Y(self, sp):
        """Mass fraction of specie at the outlet, -"""
        # This kind of averaging is dangerous, since not all species might be
        # in the species['gas']
        sumY = 0.0
        for spi in self.species['gas']:
            # sumY += self._df_sp[f'average({spi}.gas)']
            sumY += self.patches['outlet'].Y_sp()[f'average({spi}.gas)']
        return self.patches['outlet'].Y_sp()[f'average({sp}.gas)']/sumY
    
    def Yin(self, sp):
        """Mass fraction of a specie at the inlet, -"""
        # sumY = 0.0
        # for spi in self.species['gas']:
        #     sumY += self._df_sp_in[f'average({spi}.gas)']
        # return self._df_sp_in[f'average({sp}.gas)']/sumY
        raise NotImplementedError

    def X(self, sp):
        """Mole fractions of a specie at the outlet"""
        sumN = 0.0
        for spi in self.species['gas']:
            sumN += self.Y(spi)/MW(spi)
        return self.Y(sp)/MW(sp)/sumN
    
    def temperature(self, method='probe'):
        """Reactor temperature. Methods:
        probe - default, a probe at the same location as in Leion et al.
        average - an estimate of solids average temperature"""
        if method == 'probe':
            # return self._df_Tprobe['T.gas']
            return self.probes['probes'].get_data()['0']
        elif method == 'average':
            # return self._df_Tavg_est['volAverage(alpha.solids,T.solids)']
            return self.domain.Tavgest()['volAverage(alpha.solids,T.solids)']
        else:
            raise NotImplemented

    # -------------------- calculated ---------------------

    def _ndot_dO_out(self):
        """Amount of O atoms transfered according to reactions, mol/s"""
        ndot_H2O = self.mdot_out_sp('H2O')/MW('H2O')
        ndot_CO2 = self.mdot_out_sp('CO2')/MW('CO2')

        if self.gas == 'H2': 
            return ndot_H2O
        elif self.gas == 'O2': 
            # mdot_O2 = self._df_mdot_sp_in['sum(O2.gas)'] - self._df_mdot_sp_out['sum(O2.gas)']
            mdot_O2 = self.mdot_in_sp('O2') - self.mdot_out_sp('O2')
            ndot_O2 = mdot_O2/MW('O2') # mol/s
            return ndot_O2*2
        elif self.gas == 'CH4': 
            return 2*ndot_CO2 + ndot_H2O
        elif self.gas == 'CO':
            return ndot_CO2
        elif self.gas == 'syngas':
            print(f'Warning: fuel is {self.gas} -> conversion is based on both H2 and CO')
            return ndot_H2O+ndot_CO2
        else:
            raise NotImplemented(f'Fuel {self.gas} is not known')
    
    def _gamma_outlet(self, noH2O=False):
        """Gas yield (mole fraction of gas converted) based on outlet flow, -"""
        if self.gas == 'O2':
            mdot_consumed = self.mdot_in_sp('O2')-self.mdot_out_sp('O2')
            return mdot_consumed / self.mdot_in_sp('O2')
        elif self.gas == 'H2':
            return self.X('H2O')
        else:
            # Based on equation 3 in Leion et al. 2008
            # "is the fraction of component i in the outgoing gases after the
            # water has been removed (!!!)"
            # However, since mole fraction without H2O is X/(1-X_H2O), and it
            # is present for all species both in nominator and denominator, it
            # cancels out
            # def X_without_H2O(sp):
            #     """
            #     Mole fraction of a specie at the outlet without H2O"""
            #     print('WARNING: X_without_H2O is used for ' + sp)
            #     return self.X(sp) / (1-self.X('H2O'))
            # if noH2O:
            #     return (X_without_H2O('CO2')) \
            #         / (X_without_H2O('CO')+X_without_H2O('CH4')+X_without_H2O('CO2'))
            # else:
            return (self.X('CO2')) \
                / (self.X('CO')+self.X('CH4')+self.X('CO2'))
            # return (self.X('H2O')+self.X('CO2')) \
            #         / (self.X('H2')+self.X('CO')+self.X('CH4')+self.X('H2O')+self.X('CO2'))
    
    def _gamma_NRR(self):
        """Gas yield (mole fraction of gas converted) based on NRR, -
        Very unstable, especially in the beginning, when it gives very large values"""
        RR_all_gas_sp = np.sum([self.NRR(r) for r in self.reaction_species], axis=0) # amount of molecules reacted per unit time
        mdot_in_all_gas_sp = np.sum([self.mdot_in_sp(r)/MW(r) for r in self.reaction_species], axis=0)
        sign = -1 if self.gas == 'O2' else 1
        return RR_all_gas_sp/mdot_in_all_gas_sp
        
    # Conversion is defined differently for Leion et al. 2008 and Berguerand et al. 2011
    # In Leion et al. 2008 water is removed and conversion is calculated for remaining gas
    # In Berguerand et al. 2011 it is calculated for the whole mixture
    @lru_cache
    def gamma(self, method='outlet'):
        """Gas yield (mole fraction of gas converted), -
        Methods:
        outlet - default, stable, consistent with Leion et al.
        NRR - unstable"""
        if self.gas == 'CH4' and method == 'outlet':
            print('WARNING: since gas is CH4, gamma is calculated without H2O, as defined in eq.3 [Leion et al. 2008]')
            return self._gamma_outlet(noH2O=True)
        else:
            f_gamma = getattr(self, f'_gamma_{method}')
            return f_gamma()
    
    def _oxidation_mass(self):
        """Oxidation, defined as fraction of reaction-capable oxygen left in 
        ilmenite"""
        # TODO: connect to scm module (but extra dependency...)
        rho_ratio_ = MW('Fe2O3')/(2*MW('FeO'))
        return self.m_solid_sp('Fe2O3')/(self.m_solid_sp('Fe2O3')+rho_ratio_*self.m_solid_sp('FeO'))

    def _oxidation_NRR(self): # TODO: initial conversion is missing
        """Oxidation based on NRR calculation"""
        mass_O_max = self.Ro*self.solidsMass
        if self.gas == 'syngas':
            NRR_tot = integrate_time_series(self.NRR('CO')) \
                    + integrate_time_series(self.NRR('H2'))
            n_O_per_mol_NRR = 1
        else:
            # NRR_tot = integrate.cumtrapz(self.NRR(self.gas), self.NRR(self.gas).index, initial=0)
            NRR_tot = integrate_time_series(self.NRR(self.gas))
            if self.gas == 'O2':
                n_O_per_mol_NRR = 2
                return n_O_per_mol_NRR*AW('O')*NRR_tot/mass_O_max
            elif self.gas == 'CH4':
                n_O_per_mol_NRR = 4
            else:
                n_O_per_mol_NRR = 1
        return (mass_O_max-n_O_per_mol_NRR*AW('O')*NRR_tot)/mass_O_max
    
    def _oxidation_outlet(self):
        """Oxidation based on outlet mass fractions"""
        mtot_O = integrate_time_series(self._ndot_dO_out()*AW('O'))
        if self.gas == 'O2':
            return mtot_O/(self.Ro*self.solidsMass) # mass-based conversion
        else:
            return (self.Ro*self.solidsMass-mtot_O)/(self.Ro*self.solidsMass) # mass-based conversion

    def oxidation(self, method='mass'):
        """Oxidation of oxygen carrier (0 - reduced, 1 - oxidized), -
        Methods:
        mass - default, based on the mass of oxygen carrier's species
        NRR - based on net reactiion rates
        outlet - basedon outlet compositions, has time lag, but consistent with Leion et al. 2008"""
        f_gamma = getattr(self, f'_oxidation_{method}')
        return f_gamma()
    
    def omega(self, method='mass'):
        """Mass-based oxidation or fraction of maximum possible mass of oxygen carrier
        For methods see oxidation method"""
        return 1 - self.Ro*(1-self.oxidation(method=method))
    
    def done_time(self, fraction=1e-9):
        if self.gas == 'O2':
            done_ind = np.argmax(self.omega()>(np.max(self.omega())-fraction*(1-np.min(self.omega()))))
        else:
            done_ind = np.argmax(self.omega()<(np.min(self.omega())+fraction*(1-np.min(self.omega()))))
        # print(done_ind)
        return self.omega().index[done_ind]
    
    # def analytical_Qdot(self):
        # pass
        # lr.solidsMass*lr.Ro/AW('O')
    
def LabReactorAnalytical():
    pass