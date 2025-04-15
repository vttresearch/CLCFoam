# 1. Leion et al. The use of ilmenite 2008

from os.path import join
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


repo_root_path = join('..','..','..')
python_modules_path = join(repo_root_path, 'src', 'python')
if python_modules_path not in sys.path:
    sys.path.append(python_modules_path)

from specie_properties.common import AW
from specie_properties.common import MW

species = ("N2", "O2", "CO", "CO2", "H2O", "H2", "CH4")

case_fuel = os.path.split(os.getcwd())[-1].split('_')[2]
case_carrier = os.path.split(os.getcwd())[-1].split('_')[1]
print(f"Fuel is {case_fuel}, carrier is {case_carrier}")

# =====================
# Conversions plot
# =====================

case_path = '.'
# path_species = os.path.join(case_path, 'postProcessing/outletSpecies/1/surfaceFieldValue.dat')
# path_mass_flow = os.path.join(case_path, 'postProcessing/outlet_massFlowRates/1/surfaceFieldValue.dat')
path_species = os.path.join(case_path, 'postProcessing/surfaceAverageSpecieFluxesPatchFn(patch=outlet)/0/surfaceFieldValue.dat')
path_mass_flow = os.path.join(case_path, 'postProcessing/surfaceSumFn(select=patch,patch=outlet,fields=(alphaRhoPhi.solidsalphaRhoPhi.gas))/0/surfaceFieldValue.dat')

df_sp = pd.read_csv(path_species, skiprows=3, delimiter='\t')[1:]
df_mw = pd.read_csv(path_mass_flow, skiprows=3, delimiter='\t')[1:]

def Y(sp):
    sumY = 0.0
    for spi in species:
        sumY += df_sp[f'average({spi}.gas)']
    return df_sp[f'average({sp}.gas)']/sumY

def X(sp):
    sumN = 0.0
    for spi in species:
        sumN += Y(spi)/MW(spi)
    return Y(sp)/MW(sp)/sumN


# Conversion calculation similar to one used in the paper
# how much oxygen we took from 15g of ilmenite
# 1 mol CO2 -> 4 mol O taken
m_ilmenite = 0.015 # kg
# mdot_CO2 = df_sp['average(CO2.gas)']*df_mw['sum(alphaRhoPhi.gas)'] # kg/s #!!!
mdot_CO2 = Y('CO2')*df_mw['sum(alphaRhoPhi.gas)'] # kg/s #!!!

ndot_CO = mdot_CO2/MW('CO2') # mol/s
ndot_O = 4*ndot_CO # amount of O atoms transfered according to reaction   4 Fe2O3 + CH4 = 8 FeO + CO2 + 2 H2O
mdot_O = ndot_O*AW('O') # kg/s O atom mass transfer rate at time step
dm_O = df_sp[df_sp.columns[0]].diff() * mdot_O # kg O atom transported through outlet at time step
mtot_O = dm_O.cumsum() # kg O atom total transported since t=0
omega = 1 - mtot_O/m_ilmenite # mass-based conversion


# gamma calculation, equation 3 [1]
gamma = X('CO2')/(X('CO2')+X('CO')+X('CH4'))

fig, ax = plt.subplots()
indent = 600
ax.plot(omega[indent:], gamma[indent:])
ax.invert_xaxis()
ax.set_xlim([1,None])
ax.set_ylim([0,None])
plt.savefig('gamma_omega.png', dpi=300, transparent=True)




# =====================
# Outlet concentrations
# =====================

def plot_validation_species(ax, path, style={}):
    dat = np.loadtxt(path, unpack=True)
    ax.plot(dat[0], dat[1], **style)

def plot_outlet_species(ax, df_sp, specie, style={}, shift=0):
    correction = (1-X('H2O'))

    col_name = f'average(X_{specie}.gas)'
    ind = df_sp[df_sp.columns[0]]
    vals = X(specie)
    ax.plot(ind+shift, vals/correction, label=specie, **style)

    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Mole fraction [-]')


fig, ax = plt.subplots()

if case_fuel == 'CH4':
    figname = 'fig5'
    plot_species=['CO', 'CO2', 'CH4']
elif case_fuel == 'syngas':
    figname = 'fig6'
    plot_species=['CO', 'CO2']
elif case_fuel == 'O2':
    figname = None
    plot_species = ['O2']
elif case_fuel == 'CO':
    figname = None
    plot_species = ['CO', 'CO2']
else:
    print('Error: case folder name should be in the form blabla_OCTYPE_FUEL')
    exit


colors = ['k', 'r', 'b']

# Old version without error zones
for i, specie in enumerate(plot_species):
    plot_outlet_species(ax, df_sp, specie=specie, style=dict(color=colors[i]), shift=30)
    if figname is not None:
        plot_validation_species(ax, f'../validation/Leion2008_{figname}_{specie}.dat', style=dict(color=colors[i], linestyle='--'))


# New version with error zones
# validation_dats = []
# for i, specie in enumerate(plot_species):
#     plot_outlet_species(ax, df_sp, specie=specie, style=dict(color=colors[i]), shift=25)
    
#     validation_path = os.path.join(case_path, 'validation', f'Leion2008_{figname}_{specie}.dat')
#     dat = np.loadtxt(validation_path, unpack=True)
#     validation_dats.append(dat)
#     # ax.plot(dat[0], dat[1], color=colors[i], linestyle='--')
#     # plot_validation_species(ax, , style=dict())

# ts = np.linspace(0,np.max([dat[0][-1] for dat in validation_dats]),101)
# Xtot = np.sum([np.interp(ts, dat[0],dat[1]) for dat in validation_dats], 0)
# # ax.plot(ts, Xtot)
# for i, specie in enumerate(plot_species):
#     validation_path = os.path.join(case_path, 'validation', f'Leion2008_{figname}_{specie}.dat')
#     dat = np.loadtxt(validation_path, unpack=True)
#     Xtot = np.sum([np.interp(dat[0], datt[0],datt[1]) for datt in validation_dats], 0)
#     Xtot[dat[0]<25]=1
#     Xtot[dat[0]>120]=1
#     pltdat = ax.plot(dat[0], dat[1]/Xtot, color=colors[i], linestyle='--')
#     ax.fill_between(dat[0], dat[1]/Xtot-(1-Xtot)/2, dat[1]/Xtot+(1-Xtot)/2, alpha=0.1, color=pltdat[0]._color, linestyle='--')
#     ax.set_ylim([0,1])
#     # ax.set_xlim([0,120])
    

ax.legend(frameon=False)
ax.set_title('solid - present work\ndashed - Leion et al. 2008')
plt.savefig('outlet_species.png')
