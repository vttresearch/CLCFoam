"""Note that all Tlow are lowered to 280 K in order to be able to calculate
the enthalpy of formation at 298.15 K"""
# janafs_ilmenite = {'Fe2O3': ..., 'FeO': ..., 'TiO2': ...}
# janafs_ilmenite_orig = {'Fe2O3': ..., 'FeO': ..., 'TiO2': ...}


from specie_properties import janaf
R_McBride = 8.314463


def _get_McBride_polynomial_coeffs(specie):
    """Returns a tuple of high temperature (first), low temperature (second) 
    JANAF coefficient lists, and a tuple with low and high temperature limits (third)
    Taken from McBride, Bonnie J. Coefficients for calculating thermodynamic 
    and transport properties of individual species"""
    j7 = janaf.janaf7()

    if specie == 'Fe2O3': # p. 59
        return [  4.04975300E+01, -4.61315960E-02,  3.18264060E-05, 
                 -8.92263310E-09,  8.46554170E-13, -1.13176270E+05, 
                 -2.16350880E+02], \
               [ -7.70378430E+00,  1.36474710E-01, -3.29056550E-04,  
                  3.81504780E-07, -1.63102850E-10, -1.00800760E+05,  
                  2.52920850E+01], \
                (300, 2500)
        # dH0foverR_Fe2O3_table = -9.92620367E+04 # last entry in coeffs
    elif specie == 'FeO': # p. 58
        return [  5.83164890E+00,  1.42751560E-03, -9.32081430E-08, 
                 -6.59977630E-12, -2.25121430E-14, -3.45669020E+04, 
                 -2.64469900E+01], \
               [  5.31954750E+00,  2.20965910E-03,  1.07217750E-06, 
                 -2.79297290E-09,  1.33207330E-12, -3.44071650E+04, 
                 -2.36860340E+01], \
                (300, 1650)
        # dH0foverR_FeO_table = -3.27183475E+04
    elif specie == 'Fe3O4': # p. 59
        return [  2.41337200E+01,  4.15922260E-05, -2.63314920E-08,  
                  6.60350940E-12, -5.69246800E-16, -1.41210520E+05, 
                 -1.20064120E+02], \
               [  3.61981480E+01, -1.74379760E-01,  5.24756730E-04, 
                 -5.42382190E-07,  1.79962020E-10, -1.41387300E+05, 
                 -1.55566830E+02], \
                (300, 5000)
        # dH0foverR_Fe3O4_table = -1.34696136E+05
    elif specie == 'TiO2':
        return [  6.84891510E+00,  4.24634610E-03, -3.00889840E-06,
                  1.06025190E-09, -1.43795970E-13, -1.15992460E+05,
                 -3.45141060E+01], \
               [ -1.61175170E-01,  3.79666600E-02, -6.51547500E-05,
                  5.25521360E-08, -1.62000510E-11, -1.14788970E+05,
                 -1.88740350E+00], \
                (300, 2130)
        # dH0foverR_TiO2_table = -1.13628959+05
    elif specie == 'CaS': # p. 56
        return [  5.65305190E+00,  1.36258740E-03, -7.27811760E-07,
                  2.49897630E-10, -3.09681260E-14, -5.87103410E+04,
                 -2.59063950E+01 ], \
               [  4.64755580E+00,  4.93155160E-03, -5.53089030E-06,
                  3.06639590E-09, -6.07856100E-13, -5.84770470E+04,
                 -2.09227100E+01 ], \
                (300, 3000)
        # dH0foverR_CaS_table = -5.69152785E+04
    elif specie == 'CaSO4': # p. 56
        return [  8.44419050E+00,  1.18762150E-02,  0.00000000E+00,
                  0.00000000E+00,  0.00000000E+00, -1.75532420E+05,
                 -3.88201340E+01 ], \
               [  8.44419050E+00,  1.18762150E-02,  0.00000000E+00,
                  0.00000000E+00,  0.00000000E+00, -1.75532420E+05,
                 -3.88201340E+01 ], \
                (300, 5000)
        # dH0foverR_CaSO4_table = -1.72486926E+05
    # elif specie == 'MnO': # p.
    #     return [   ], \
    #            [   ], \
    #             (, )
    #     # dH0foverR_MnO_table = 
    # elif specie == 'Mn3O4': # p.
    #     return [   ], \
    #            [   ], \
    #             (, )
    #     # dH0foverR_Mn3O4_table = 
    # elif specie == 'NiO': # p. 42 WRONG! GAS
    #     return [  4.10461140E+00,  4.86591600E-04, -1.87867840E-07,
    #               3.55318550E-11, -2.47151660E-15,  3.64456450E+04,
    #               4.07692910E+00 ], \
    #            [  2.99196820E+00,  3.33092080E-03, -1.53524710E-06,
    #              -1.56408330E-09,  1.21285010E-12,  3.67420940E+04,
    #               9.82153990E+00 ], \
    #             (300, 5000)
    #     # dH0foverR_NiO_table = 3.77657540E+04
    # elif specie == 'Ni': # p. 42 WRONG! GAS
    #     return [  3.20614900E+00, -2.09699230E-04, -2.28364480E-08, 
    #               1.50852110E-11, -1.00044450E-15,  5.07081260E+04,
    #               3.53171623E+00 ], \
    #            [  2.77666540E+00, -7.52206380E-04,  4.32561130E-06,
    #              -5.47312870E-09,  2.11075650E-12,  5.09090830E+04,
    #               6.16823253E+00 ], \
    #             (300, 5000)
    #     # dH0foverR_Ni_table = 5.17319098E+04
    else:
        raise NotImplementedError('Specie ' + specie + ' is not implemented')
    
def _get_McBride_janaf(specie):
    highCpCoeffs, lowCpCoeffs, (Tlow, Thigh) = _get_McBride_polynomial_coeffs(specie)
    if Tlow > 280:
        # print(f"Warning: changing Tlow for {specie} from {Tlow} K to 280 K")
        Tlow = 280
    j7 = janaf.janaf7()
    j7.add_coeffs(lowCpCoeffs,  Tlow, 1000)
    j7.add_coeffs(highCpCoeffs, 1000, Thigh)
    return j7

janafs_ilmenite_orig = \
    {sp: _get_McBride_janaf(sp) for sp in ('Fe2O3', 'FeO', 'TiO2', 'Fe3O4')}
"""JANAF objects for ilmenite species from """

def _get_corrected_Fe2O3_for_ilmenite(T_paper = 925 + 273.15):
    #   O2 + 4 FeO -> 2 Fe2O3
    # a O2 + b FeO -> d Fe2O3
    a=1; b=4; d=2

    # Hallberg et al. - A method for determination of reaction enthalpy of 
    # oxygen carriers for chemical looping combustion – Application to ilmenite, 
    # 2011
    dHRm_ox_paper = -468500 # [J/mol O2]
    # T_paper = 298.15 # [K] = I once made a mistake here and used this value

    dHRm_ox_orig =  d*janafs_ilmenite_orig['Fe2O3'].H0(T_paper) \
                  - b*janafs_ilmenite_orig['FeO'].H0(T_paper) \
                  - a*janaf.janaf7Cantera('O2').H0(T_paper) # [kJ/mol O2]
    delta_Hf_Fe2O3 = (dHRm_ox_paper-dHRm_ox_orig)/(d/a)

    highCpCoeffs, lowCpCoeffs, (Tlow, Thigh) = _get_McBride_polynomial_coeffs('Fe2O3')
    lowCpCoeffs[5]  = lowCpCoeffs[5]  + delta_Hf_Fe2O3/R_McBride
    highCpCoeffs[5] = highCpCoeffs[5] + delta_Hf_Fe2O3/R_McBride
    j7 = janaf.janaf7()
    if Tlow > 280:
        # print(f"Warning: changing Tlow for Fe2O3 from {Tlow} K to 280 K")
        Tlow = 280
    j7.add_coeffs(lowCpCoeffs,  Tlow, 1000)
    j7.add_coeffs(highCpCoeffs, 1000, Thigh)
    return j7

janafs_ilmenite = \
    {sp: _get_McBride_janaf(sp) for sp in ('FeO', 'TiO2', 'Fe3O4')}
"""JANAF objects for ilmenite species with correction to Fe2O3 fixing the 
reaction heat of oxidation"""

janafs_ilmenite['Fe2O3'] = _get_corrected_Fe2O3_for_ilmenite()

# Reaction heat
