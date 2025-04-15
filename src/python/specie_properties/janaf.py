import requests
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class janaf:
    """
    Base class for JANAF polynomials.
    Stores lists of polynomials per temperature intervals (self._coeff_dicts).
    Each dict has 'Tlow', 'Thigh' keys, storing temperature limits (floats).
    
    Child classes need to define three functions (STD - 298.15 K, 1 bar):
    Cp0 - isobaric heat capacity at STD
    H0  - enthalpy (sensible+formation) at STD
    S0  - entropy at STD
    0 stands for standard pressure.
    """
    def __init__(self) -> None:
        self._R = 8.314463      # [J/K/mol] - universal gas constant
        # Note that e.g. NASA polynomials in McBride's 1993 report are
        # calculated for R=8.314510
        self._T_STD = 298.15    # [K]
        self._p_STD = 1e5       # [Pa]
        self._coeff_dicts = list()
    
    def add_coeffs(self, coeffs: list, Tlow: float, Thigh: float):
        """Add coeffs for temperature interval [Tlow, Thigh]"""
        for coeff_dict in self._coeff_dicts:
            if coeff_dict['Tlow'] <  Tlow  < coeff_dict['Thigh'] \
            or coeff_dict['Tlow'] <  Thigh < coeff_dict['Thigh']:
                raise ValueError(f'Error: supplied range {Tlow}-{Thigh} conflicts with'
                      ' existing {coeff_dict["Tlow"]}-{coeff_dict["Thigh"]}')
        self._coeff_dicts.append(dict(coeffs=coeffs, Tlow=Tlow, Thigh=Thigh))
        
    def get_coeffs(self, T: float):
        """Get list of coefficients for temperature T"""
        for coeff_dict in self._coeff_dicts:
            if coeff_dict['Tlow'] <= T < coeff_dict['Thigh']:
                return coeff_dict['coeffs']
        return None
    
    def dH0f(self):
        """Enthalpy of formation at standard pressure and temperature [J/mol]"""
        return self.H0(self._T_STD)

    def H0s(self, T):
        """Sensible enthalpy at standard pressure [J/mol]"""
        return self.H0(T) - self.Hf


class janaf7(janaf):
    """
    NASA7 format polynomials. As defined in:
    McBride, Bonnie J. Coefficients for calculating thermodynamic and transport 
    properties of individual species. Vol. 4513. National Aeronautics and Space 
    Administration, Office of Management, Scientific and Technical Information 
    Program, 1993.

    """
    def Cp0(self, T: float):
        """Isobaric heat capacity at standard pressure [J/K/mol]"""
        t = T
        c = self.get_coeffs(T)
        return (c[0] \
            + c[1]*t \
            + c[2]*t*t \
            + c[3]*t*t*t \
            + c[4]*t*t*t*t)*self._R

    def H0(self, T: float):
        """Enthalpy (sensible+formation) at standard pressure [J/mol]"""
        t = T
        c = self.get_coeffs(T)
        return (c[0]*t \
              + c[1]*t*t/2 \
              + c[2]*t*t*t/3 \
              + c[3]*t*t*t*t/4 \
              + c[4]*t*t*t*t*t/5 \
              + c[5])*self._R
    
    def S0(self, T: float):
        """Entropy at standard pressure [J/mol/K]"""
        t = T
        c = self.get_coeffs(T)
        return (c[0]*np.log(t) \
             +  c[1]*t \
             +  c[2]*t*t/2 \
             +  c[3]*t*t*t/3 \
             +  c[4]*t*t*t*t/4 \
             +  c[6])*self._R
    
    def print_coeffs_openfoam(self):
        """Print coefficients to standard output in format of thermodynamics
        subdictionary in physicalProperties dictionary of OpenFOAM"""
        assert len(self._coeff_dicts) == 2
        if self._coeff_dicts[0]['Tlow'] < self._coeff_dicts[1]['Tlow']:
            lowCoeffs = self._coeff_dicts[0]
            highCoeffs = self._coeff_dicts[1]
        else:
            highCoeffs = self._coeff_dicts[0]
            lowCoeffs = self._coeff_dicts[1]
        assert lowCoeffs['Thigh'] == highCoeffs['Tlow']
        print(f"Tlow            {lowCoeffs['Tlow']};")
        print(f"Thigh           {highCoeffs['Thigh']};")
        print(f"Tcommon         {lowCoeffs['Thigh']};")
        print(f"highCpCoeffs    ( {' '.join([f'{val:15.8e}'.replace('e','E') for val in highCoeffs['coeffs']])} );")
        print(f"lowCpCoeffs     ( {' '.join([f'{val:15.8e}'.replace('e','E') for val in lowCoeffs['coeffs']])} );")

class janaf7Web(janaf7):
    """
    NASA7 coefficients from url with Chemkin format thermo file.
    """
    def __init__(self, formula, 
    thermourl='http://combustion.berkeley.edu/gri-mech/version30/files30/thermo30.dat',
    TmidDefault=1000
    ) -> None:
        super().__init__()
        self.formula = formula
        # self.coeffs = list()
        self.url = thermourl
        self.TmidDefault = TmidDefault
        self.thermo_page = self._retrieve_thermo(self.url)
        self.lines = self._get_specie_lines(self.thermo_page.decode(), formula)
        from cantera import ck2yaml
        p = ck2yaml.Parser()
        species, thermo, composition = p.read_NASA7_entry(self.lines, 
                                                          self.TmidDefault, '')
        if thermo.Tmid is None:
            thermo.Tmid = self.TmidDefault
        self.add_coeffs(thermo.low_coeffs, thermo.Tmin, thermo.Tmid)
        if thermo.high_coeffs is None:
            self.add_coeffs(thermo.low_coeffs, thermo.Tmid, thermo.Tmax)
        else:
            self.add_coeffs(thermo.high_coeffs, thermo.Tmid, thermo.Tmax)

    def _retrieve_thermo(url: str) -> str:
        """Load contents of url"""
        import requests
        return requests.get(url).content

    def _get_specie_lines(thermo_txt, specie):
        """Extract lines related to a specie"""
        lines = thermo_txt.splitlines()
        TmidDefault = float(lines[1].split()[1])
        for i, line in enumerate(lines):
            if line.startswith(specie+' '):
                break
        if i >= len(lines)-1:
            print("Specie not found")
            return None
        else:
            return lines[i:i+4]

class janaf7Cantera(janaf7):
    """
    NASA7 coefficients from local path to Cantera mechanism.
    Default mechanism - gri30, built in Cantera.
    """
    def __init__(self, formula, mech='gri30.yaml') -> None:
        super().__init__()
        self.formula = formula
        self.coeffs = list()
        self.mech_path = mech

        import cantera as ct
        gas = ct.Solution(self.mech_path)
        ind = gas.species_index(formula)
        
        # https://github.com/Cantera/cantera/blob/03f5b82352f47ceec7cc0c044295546b30b2be2d/interfaces/cython/cantera/speciesthermo.pyx#L181
        thermo = gas.species()[ind].thermo
        Tmid = thermo.coeffs[0]
        high_coeffs = thermo.coeffs[1:8]
        low_coeffs  = thermo.coeffs[8:15]

        self.add_coeffs(low_coeffs, thermo.min_temp, Tmid)
        self.add_coeffs(high_coeffs, Tmid, thermo.max_temp)

class janafNist(janaf):
    """
    NASA7 coefficients from local path to Cantera mechanism.
    Default mechanism - gri30, built in Cantera.
    """
    def Cp0(self, T):  # [J/K/mol] - heat capacity
        t = T/1000
        c = self.get_coeffs(T)
        return c[0] \
            + c[1]*t \
            + c[2]*t*t \
            + c[3]*t*t*t \
            + c[4]/(t*t)
    
    def H0(self,T): # [J/mol] - standard enthalpy
        t = T/1000
        c = self.get_coeffs(T)
        return (c[0]*t \
            + c[1]*t*t/2 \
            + c[2]*t*t*t/3 \
            + c[3]*t*t*t*t/4 \
            - c[4]/t \
            + c[5])*1000
    
    def S0(self,T):
        t = T/1000
        c = self.get_coeffs(T)
        return c[0]*np.log(t) \
            + c[1]*t \
            + c[2]*t*t/2 \
            + c[3]*t*t*t/3 \
            - c[4]/(2*t*t) \
            + c[6]

class janafNistWeb(janafNist):
    def search_nist_urls(self, verbose=False):
        import nistchempy as nist
        formula = self.formula
        # nist.print_search_parameters()
        search = nist.Search(MatchIso = True, AllowOther = False, AllowExtra = False, Units = 'SI', NoIon = True)
        search.find_compounds(identifier = formula, search_type = 'formula')
        # print(search.success, search.lost, search.IDs, search.compounds)
        search.load_found_compounds()
        compound = search.compounds[0]
        
        if verbose:
            print(f'Found {len(search.compounds)} compounds. Picking {compound.name}')
        if 'cTC' in compound.data_refs.keys() and 'cTG' in compound.data_refs.keys():
            print('WARNING: Both condensed and gaseous phase data found!')
        if 'cTC' in compound.data_refs.keys():
            self.nist_url_cond = compound.data_refs['cTC'][0]
        else: 
            self.nist_url_cond = None
        if 'cTG' in compound.data_refs.keys():
            self.nist_url_gas  = compound.data_refs['cTG'][0]
        else:
            self.nist_url_gas = None
            
    def exctract_janaf_from_table(self, nist_html_table, verbose=True):
        # returns a list of dictionaries with coeffs and range
        rows = nist_html_table.find_all('tr')
        
        # Header row
        r = rows[0]
        rheader = r.find_all('th')
        assert len(rheader)==1
        if verbose:
            print(rheader[0].contents[0], end='\t')
        cols = rows[0].find_all('td')
        nranges = len(cols)
        coeffs = [dict(coeffs=[0]*8) for i in range(nranges)]
        for i, col in enumerate(cols):
            # l = col.contents[0].split(' - ') # used to be this in 2023
            l = col.contents[0].split(' to ')
            # coeffs[i]['range'] = (float(l[0]), float(l[1]))
            coeffs[i]['Tlow'] = float(l[0])
            coeffs[i]['Thigh'] = float(l[1])
            if verbose:
                print(coeffs[i]['Tlow'], ' ', end='\t')
                print(coeffs[i]['Thigh'], end='\t')
        if verbose:
            print()
        
        # Parsing data
        df = pd.DataFrame(columns=['Name', '1', '2', '3'])
        for i in range(1, 9): # read all 7 JANAF coeffs + standard enthalpy of formation
            row_title = rows[i].find_all('th')[0].contents[0]
            if verbose:
                print(row_title, end='\t')
            cols = rows[i].find_all('td')
            for j in range(nranges):
                coeffs[j]['coeffs'][i-1] = float(cols[j].text.replace('×10', 'e')) # e.g. to fix FeO https://webbook.nist.gov/cgi/cbook.cgi?ID=C1345251&Units=SI&Mask=2#Thermo-Condensed
                if verbose:
                    print(coeffs[j]['coeffs'][i-1], end='\t')
            # df = df.append({'Name': row_title}, ignore_index=True)

            if verbose:
                print()
        return coeffs
    
    def retrieve_nist(self, url, verbose=False):        
        from bs4 import BeautifulSoup

        page = requests.get(url) # TODO
        
        soup = BeautifulSoup(page.content, 'html.parser')
        tables = soup.find_all("table", attrs={'aria-label': lambda L: L and 'Shomate Equation' in L})
        
        # assert len(els) == 1
        if len(tables) > 1:
            print(f'WARNING: multiple tables found ({len(tables)}) for {self.formula}')
        
        for table in tables:
            self._coeff_dicts.extend( self.exctract_janaf_from_table(table, verbose=verbose) )
            # TODO check ranges?

    def __init__(self, formula) -> None:
        super().__init__()
        self.formula = formula
        self.verbose = False
        self.coeffs = list()
        self.search_nist_urls(self.verbose)
        if not self.nist_url_cond is None:
            self.retrieve_nist(self.nist_url_cond, self.verbose)
        if not self.nist_url_gas is None:
            self.retrieve_nist(self.nist_url_gas, self.verbose)
        # print(self.coeffs)
        # self.add_coeffs(self, coeffs, Tlow, Thigh)




# def janafHf(T,c):
#     if len(c) >= 8:
#         return c[7]*1000
#     else:
#         return janafH(298, c)*1000

# def janafdH(T, c): # J/mol
#     return janafH(T, c) - janafHf(T, c)


# def janafG(T, c):
#     t = T/1000
#     return janafH(T, c) - T*janafS(T, c)
