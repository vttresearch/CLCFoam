"""Analytical model for 1d batch reactor.
Read the paper for details.
In short, it assumes quasi-steady state with oxidation of solids and temperature 
being constant and calculates gas distributions accross the height of the 
reactor.

Main methods to use are:
- oned_Ca_Vdot_vs_hbar(...) - gives concentrations and volumetric flows versus 
non-dimensionalized reactor height.
- gas_yield_Leion(...) - gives gas yield, as defined in Eq. (3) of 
Leion et al. (2008).

See also analytical_model.ipynb for examples."""

import numpy as np
from scm.scm import IlmenitePhase
from scipy.integrate import odeint

# See analytical_model

def core_total_area(phase, reaction, ms, oxidation):
    """Total area of unreacted cores of all grain of all oxygen carrier in the
    reactor. [m^2]"""
    conversion = 1.-oxidation if reaction.reduction else oxidation
    return 3*np.power(1.-conversion, 2./3.) * ms / (reaction.r_g * phase.density(oxidation))

def derivatives(w, hbar, params):
    """
    Returns derivatives of concentration of reactant specie cA and volumetric
    gas flow Vdot with respect to non-dimensionized reactor height hbar.

    Arguments:
        w :  vector of the state variables:
                  w = [cA, Vdot]
        hbar :  normalized height (0..1)
        p :  vector of the parameters, defining oxygen carrier state, reaction,
        and reactor conditions:
                  p = [phase, reaction, ms, oxidation, T, p]
            where 
                phase       - an instance of OCphase
                reaction    - instance of Abad_SCR
                ms          - mass of oxygen carrier in reactor
                oxidation   - oxidation of oxygen carrier (0..1)
                T           - temperatures
    Note: for O2 reaction, don't forger to set cA lower than 1/Vm (in other 
    words, dilute it)
    """
    cA, Vdot = w
    phase, reaction, ms, oxidation, T, p = params
    Vm = 8.314*T/p
    rhos = phase.density(oxidation)
    Acore = core_total_area(phase, reaction, ms, oxidation)
    csstoich = np.sum([c[0] for c in reaction.Cs])
    dVdotdhbar = Vm*(csstoich-reaction.a)*Acore*reaction.k(T)*np.power(cA,reaction.n)

    return [
        -reaction.a*Acore*reaction.k(T)*np.power(cA,reaction.n)/Vdot - cA/Vdot*dVdotdhbar, # = dcA/dhbar
        dVdotdhbar
    ]

def oned_cA_Vdot_vs_hbar_firstord(phase, reaction, ms, oxidation, T, p, xA0, Vdot0, hbar):
    """
    Calculate concentration of reactant specie cA and volumetric
    gas flow Vdot with respect to non-dimensionized reactor height hbar.
    Formulation for a first-order reaction with no molar flow change.
    """
    assert hasattr(hbar, '__iter__')
    assert np.isclose(reaction.n, 1.)
    assert np.isclose(np.sum([c[0] for c in reaction.Cs]), reaction.a)
    cA0 = xA0/(8.314*T/p)
    alpha = reaction.k(T) * core_total_area(phase, reaction, ms, oxidation) / Vdot0
    return cA0*np.exp(-alpha*np.array(hbar)), np.array([Vdot0]*len(hbar))

def oned_cA_Vdot_vs_hbar_nthord(phase, reaction, ms, oxidation, T, p, xA0, Vdot0, hbar):
    """
    Calculate concentration of reactant specie cA and volumetric
    gas flow Vdot with respect to non-dimensionized reactor height hbar.
    Formulation for a n-order reaction with no molar flow change.
    """
    assert hasattr(hbar, '__iter__')
    assert np.isclose(np.sum([c[0] for c in reaction.Cs]), reaction.a)
    cA0 = xA0/(8.314*T/p)
    alpha = reaction.k(T) * core_total_area(phase, reaction, ms, oxidation) / Vdot0
    z = 1 - reaction.n
    cA = np.power(z, 1./z)*np.power(np.clip(np.power(cA0, z)/z-alpha*np.array(hbar), a_min=0, a_max=None), 1/z)
    cA[cA < 0] = 0
    # cA = np.nan_to_num(cA)
    return cA, np.array([Vdot0]*len(hbar))


def oned_cA_Vdot_vs_hbar_odes(phase, reaction, ms, oxidation, T, p, xA0, Vdot0, hbar):
    """
    Calculate concentration of reactant specie cA and volumetric
    gas flow Vdot with respect to non-dimensionized reactor height hbar.
    System of two ODEs is used, the most general solution.
    """
    assert hasattr(hbar, '__iter__')
    cA0 = xA0/(8.314*T/p)
    w0 = [cA0, Vdot0]
    params = [phase, reaction, ms, oxidation, T, p]
    wsol = odeint(derivatives, w0, hbar, args=(params,), atol=1.0e-8, rtol=1.0e-6)
    return wsol.T[0], wsol.T[1]

def oned_Ca_Vdot_vs_hbar(phase, reaction, ms, oxidation, T, p, xA0, Vdot0, hbar):
    if np.isclose(np.sum([c[0] for c in reaction.Cs]), reaction.a):
        if np.isclose(reaction.n, 1.):
            return oned_cA_Vdot_vs_hbar_firstord(phase, reaction, ms, oxidation, T, p, xA0, Vdot0, hbar)
        else:
            return oned_cA_Vdot_vs_hbar_nthord(phase, reaction, ms, oxidation, T, p, xA0, Vdot0, hbar)
    else:
        return oned_cA_Vdot_vs_hbar_odes(phase, reaction, ms, oxidation, T, p, xA0, Vdot0, hbar)

def gas_yield_Leion(phase, reaction, ms, oxidation, T, p, xA0, Vdot0):
    hbar = np.linspace(0, 1, 1001)
    cA0 = xA0/(8.314*T/p)
    cA, Vdot = oned_Ca_Vdot_vs_hbar(phase, reaction, ms, oxidation, T, p, xA0, Vdot0, hbar)
    if reaction.A in ['CH4', 'H2', 'CO']:
        # if np.isclose(xA0, 1): # 100% of inlet gas is CH4
        # in case of CH4, the molar amount of CO2 at the outlet is equal to the one of consumed CH4
        # in case of H2, amount of H2O is equal to the amount of H2 consumed
        # in case of CO, amount of CO2 is equal to the amount of CO consumed
        return 1 - cA[-1]*Vdot[-1]/(cA0*Vdot[0])
        # else:
        #     raise NotImplementedError('Inlet gas contains not only '+reaction.A)
    raise NotImplementedError('Reaction fuel ' + reaction.A + ' is not supported')
