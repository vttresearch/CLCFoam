import numpy as np
from specie_properties.common import AW, MW
from specie_properties.janaf import janaf
from scipy.optimize import root_scalar

# import cantera as ct
# # https://github.com/Cantera/cantera/blob/03f5b82352f47ceec7cc0c044295546b30b2be2d/src/thermo/Elements.cpp#L67
# MW_FeO   = ct.Element('Fe')._weight   + ct.Element('O')._weight   # 71.844
# MW_Fe2O3 = ct.Element('Fe')._weight*2 + ct.Element('O')._weight*3 # 159.687
# MW_TiO2  = ct.Element('Ti')._weight   + ct.Element('O')._weight*2 # 79.865

const_R = 8.314463
const_MW_Fe2O3 = 159.687
const_MW_FeO = 71.844

def molar_vol_ideal_gas(T, p=101325):
    return const_R*T/p

def concentration_from_molar_fraction(T, x, p=101325):
    return x/molar_vol_ideal_gas(T, p)

# class oxygen_carrier:
#     def __init__(self, reactions: list):
#         self.reactions = reactions

def splitCoeffSpecie(s):
    for i, ch in enumerate(s):
        if not ch.isdigit():
            break
    if i > 0:
        return s[:i], s[i:]
    else:
        return 1, s

def parseReaction(rstr): # "H2 + Fe2O3 -> H2O + 2FeO"
    lhs, rhs = rstr.split('->')
    (a,A), (b,B) = [splitCoeffSpecie(el.strip()) for el in lhs.split('+')]
    rhsL = [splitCoeffSpecie(el.strip()) for el in rhs.split('+')]
    d, D = rhsL[-1]
    Cs = rhsL[:-1]
    return (int(a), A), (int(b), B), [(int(c), C)for (c, C) in Cs], (int(d), D)


class Abad_SCR:
    # a A(g) + b B(s) -> c_i C_i(g) + d D(s)
    def __init__(self, rstr, rho_m, r_g, bavg, ks0, Echr, n) -> None:
        self.rstr = rstr
        self.rho_m = rho_m
        self.r_g = r_g
        self.bavg = bavg
        self.ks0 = ks0
        self.Echr = Echr
        self.n = n

        (a, A), (b, B), Cs, (d, D) = parseReaction(rstr)
        self.f_OF = (b*MW(B) - d*MW(D)) / (b*AW('O'))
        if B == 'Fe2O3':
            self.reduction = True
        elif B == 'FeO':
            self.f_OF = - self.f_OF
            self.reduction = False
        else: 
            raise NotImplementedError

        self.a = a
        self.b = b
        self.d = d
        self.A = A
        self.B = B
        self.D = D
        self.Cs = Cs
            
    def k(self, T):
        return self.ks0*np.exp(-self.Echr/const_R/T)

    def _tau_reac(self, T, C):
        return self.rho_m*self.r_g/(self.bavg*self.k(T)*(C**self.n))

    def tau(self, T, C):
        return self._tau_reac(T, C)
    
    def _t_over_tau(self, X):
        return 1 - (1-X)**(1./3.)
    
    def _conversion(self, t_over_tau):
        """Inverse of _t_over_tau(X)"""
        return 1-(1-t_over_tau)**3
    

    def conversion(self, T, C, t):
        return self._conversion(t/self._tau_reac(T, C))
        # return 1-(1-t/self.tau(T, C))**3
    
    def heat_release_rate(self, h):
        """
        Heat release rate enthalpy function h(sp), where sp is a string with 
        a specie name.
        Example usage:
        def getj7(specie):
            if specie in janafs_ilmenite.keys():
                return janafs_ilmenite[specie]
            else:
                return janaf.janaf7Cantera(specie)

            def h(specie, T=None):
                if T is None:
                    return getj7(specie).dH0f()
                else:
                    return getj7(specie).H0(T)
        """
        return - self.a*h(self.A) \
               - self.b*h(self.B) \
               + np.sum([c[0]*h(c[1]) for c in self.Cs]) \
               + self.d*h(self.D)

    
class Abad_SCRd(Abad_SCR):
    def __init__(self, rstr, rho_m, r_g, bavg, ks0, Echr, n, De0, Edif, Xcrit, r_p) -> None:
        super().__init__(rstr, rho_m, r_g, bavg, ks0, Echr, n)
        self.De0 = De0
        self.Edif = Edif
        self.Xcrit = Xcrit # lower -> reaction, higher -> diffusion
        self.r_p = r_p # diffusion is controlled by particle radius, not grain
    
    def De(self, T):
        return self.De0*np.exp(-self.Edif/const_R/T)

    def _tau_diff(self, T, C):
        return self.rho_m*self.r_p*self.r_p/(6*self.bavg*self.De(T)*(C**self.n))
    
    def _t_over_tau_diff(self, X):
        """For pure diffusion-limited reaction"""
        # assert 0 <= X <= 1
        return ( 1 - 3*(1-X)**(2./3.) + 2*(1-X) )
    
    def _conversion_diff(self, t_over_tau):
        """For pure diffusion-limited reaction
        Inverse of _t_over_tau_diff(X)"""
        # Vectorized version using numpy
        t_over_tau = np.asarray(t_over_tau)
        roots = np.zeros_like(t_over_tau)

        def func(X, t_over_tau):
            return self._t_over_tau_diff(X) - t_over_tau

        for i, t in np.ndenumerate(t_over_tau):
            sol = root_scalar(func, args=(t,), bracket=[0, 1], method='brentq')
            roots[i] = sol.root

        return roots

    def tau(self, T, C):
        tau_r = super().tau(T, C)
        tau_d = self._tau_diff(T, C)
        tau_critr = super()._t_over_tau(self.Xcrit)*tau_r # time to reach Xcrit
        tau_critd = self._t_over_tau_diff(self.Xcrit)*tau_d

        return tau_critr + (tau_d - tau_critd)

    def conversion(self, T, C, t):
        t = np.asarray(t)
        X = super().conversion(T, C, t)
        Xcrit_mask = X < self.Xcrit

        tau_r = super().tau(T, C)
        tau_d = self._tau_diff(T, C)
        tau_critr = super()._t_over_tau(self.Xcrit) * tau_r
        tau_critd = self._t_over_tau_diff(self.Xcrit) * tau_d
        print(f'tau_critr = {tau_critr}, tau_critd = {tau_critd}')

        t_in_diff = t - tau_critr
        t_in_diff[t_in_diff < 0] = 0

        X_diff = self._conversion_diff((tau_critd + t_in_diff) / tau_d)

        return np.where(Xcrit_mask, X, X_diff)

class OCphase():
    """Density of inerts is assumed rhoi = rhoox"""
    def __init__(self, rhoox, rhored, Ro) -> None:
        self.rhoox = rhoox
        self.rhored = rhored
        self.Ro = Ro

        self.rhoi = rhoox
        self.gamma = rhoox/rhored
        self.Ro_max = 1 - 1/self.gamma
    
    def volFractions(self, X):
        """Returns vol fractions for ox, red, and inert"""
        vi = 1.0 - self.Ro / self.Ro_max
        vox = X*(1-vi)
        vred = 1-vi-vox
        return vox, vred, vi
    
    def density(self, X):
        vox, vred, vi = self.volFractions(X)
        return vox*self.rhoox + vred*self.rhored + vi*self.rhoi
    
    def massFractions(self, X):
        """Returns mass fractions for ox, red, and inert"""
        vox, vred, vi = self.volFractions(X)
        rhos = self.density(X)
        return vox*self.rhoox/rhos, vred*self.rhored/rhos, vi*self.rhoi/rhos

class IlmenitePhase(OCphase):
    def __init__(self, rhoox, Ro) -> None:
        rhored = 2*const_MW_FeO/const_MW_Fe2O3 *rhoox
        super().__init__(rhoox, rhored, Ro)


rstr_H2  = "H2 + Fe2O3 -> H2O + 2FeO"
rstr_CO  = "CO + Fe2O3 -> CO2 + 2FeO"
rstr_CH4 = "CH4 + 4Fe2O3 -> CO2 + 2H2O + 8FeO"
rstr_O2  = "O2 + 4FeO -> 2Fe2O3"

# not corrected values ! TODO: implement correction
ilmenite_orig = dict(
    pre = dict(
        H2  = Abad_SCR (rstr=rstr_H2,  rho_m=13590, r_g=0.5e-6,  bavg=1.19, ks0=0.51, Echr=109200, n=1),
        CO  = Abad_SCR (rstr=rstr_CO,  rho_m=13590, r_g=0.5e-6,  bavg=1.19, ks0=0.21, Echr=113300, n=1),
        CH4 = Abad_SCR (rstr=rstr_CH4, rho_m=13590, r_g=0.5e-6,  bavg=4.74, ks0=8.8 , Echr=165200, n=1),
        O2  = Abad_SCRd(rstr=rstr_O2,  rho_m=31100, r_g=0.48e-6, bavg=4,    ks0=8e-5, Echr= 11800, n=1, De0=1.37e-5, Edif=77400, Xcrit=0.25, r_p=160e-6),
    ),
    act = dict(
        H2  = Abad_SCR(rstr=rstr_H2,  rho_m=13590, r_g=1.25e-6, bavg=1.45, ks0=6.2e-2, Echr= 65000, n=1),
        CO  = Abad_SCR(rstr=rstr_CO,  rho_m=13590, r_g=1.25e-6, bavg=1.45, ks0=0.1,    Echr= 80700, n=0.8),
        CH4 = Abad_SCR(rstr=rstr_CH4, rho_m=13590, r_g=1.25e-6, bavg=5.78, ks0=9.8,    Echr=135200, n=1),
        O2  = Abad_SCR(rstr=rstr_O2,  rho_m=31100, r_g=1.20e-6, bavg=4,    ks0=1.9e-3, Echr= 25500, n=1),
    )
)

ilmenitePhase_Abad = dict(
    pre = IlmenitePhase(4100, 0.040),
    act = IlmenitePhase(4250, 0.033)
)

def correct_reaction(r, ip):
    rho_m_new = ip.Ro * ip.rhoox / r.f_OF / AW('O')
    ks0_new = r.ks0 * r.bavg / r.b * rho_m_new / r.rho_m

    if type(r) == Abad_SCR:
        return Abad_SCR ( rstr    = r.rstr, 
                        rho_m   = rho_m_new, # modified
                        r_g     = r.r_g,  
                        bavg    = r.b,      # modified
                        ks0     = ks0_new,  # modified
                        Echr    = r.Echr, 
                        n       = r.n)
    elif type(r) == Abad_SCRd:
        De0_new = r.De0 * r.bavg / r.b * rho_m_new / r.rho_m
        return Abad_SCRd( rstr    = r.rstr, 
                        rho_m   = rho_m_new, # modified
                        r_g     = r.r_g,  
                        bavg    = r.b,      # modified
                        ks0     = ks0_new,  # modified
                        Echr    = r.Echr, 
                        n       = r.n,
                        De0     = De0_new, # modified
                        Edif    = r.Edif, 
                        Xcrit   = r.Xcrit,
                        r_p     = r.r_p)
    else:
        raise NotImplementedError

ilmenite = {
    chem: {gas: correct_reaction(ilmenite_orig[chem][gas], ilmenitePhase_Abad[chem]) 
           for gas in ('H2', 'CO', 'CH4', 'O2')} 
           for chem in ('pre', 'act')
    }

# # validation of correction: tau must be the same
# T = 950 + 273.15
# x = 1.0
# p = 101325
# for chem in ['pre', 'act']:
#     for gas in ('H2', 'CO', 'CH4', 'O2'):
#         if chem =='pre' and gas == 'O2':
#             continue
#         print(ilmenite[chem][gas].tau(T, concentration_from_molar_fraction(T, x=x, p=p)))
#         print(ilmenite_orig[chem][gas].tau(T, concentration_from_molar_fraction(T, x=x, p=p)))