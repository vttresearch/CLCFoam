from os.path import join, split
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from labreactor.labreactor import LabReactor
from foam.common import get_wpd_columns



# def plot_validation_species(ax, path, style={}):
#     dat = np.loadtxt(path, unpack=True)
#     ax.plot(dat[0], dat[1], **style)

# def plot_outlet_species(ax, labReactor, specie, style={}, shift=0):
#     print('Warning: plot_outlet_species funtion removes H2O from outlet composition!')
#     correction = (1-labReactor.X('H2O'))
#     vals = labReactor.X(specie)
#     ind = vals.index
#     ax.plot(ind+shift, vals/correction, label=specie, **style)


# def plot_labReactors_temporal(case_paths: list, format='png', common_style = dict(linewidth=1.5)):
#     cases = []
#     for case_path in case_paths:
#         cases.append(LabReactor(case_path))

#     timefigname = cases[0].timefigname
#     plot_species=cases[0].plot_species
#     styles_specie = [
#         dict(color='r'),
#         dict(color='g'),
#         dict(color='b'),
#     ]

#     specie_annotation_xy = {
#         'CH4': [120, 0.6],
#         'CO2': [62, 0.39],
#         'CO':  [80, 0.07], 
#     }

#     specie_labels = {
#         'CH4':'$\mathrm{CH}_4$',
#         'CO2':'$\mathrm{CO}_2$',
#         'CO':'$\mathrm{CO}$', 
#     }
#     sim_line_styles = ['--', ':', (0, (1, 10))]

#     fig, ax = plt.subplots()
#     for i, specie in enumerate(plot_species):
#         style = common_style | styles_specie[i]
#         validation_path = join(cases[0].path,'validation',f'Leion2008_{timefigname}_{specie}.dat')
#         for j, labReactor in enumerate(cases):
#             plot_outlet_species(ax, labReactor, specie=specie, style=dict(linestyle=sim_line_styles[j]) | style, shift=25)
#         plot_validation_species(ax, validation_path,   style=dict(linestyle='-')  | style)

#         ax.annotate(specie_labels[specie], specie_annotation_xy[specie], size=16, color=styles_specie[i]['color'])

#     ax.set(
#         xlim = [0,160],
#         ylim = [0,1.0],
#         xlabel = 'Time, s',
#         ylabel = 'Mole fraction, -',
#     )

#     exp_legend = mlines.Line2D([], [], color='k', linestyle='-', label='Leion et al. 2008')
#     sim_legends = []
#     for j, labReactor in enumerate(cases):
#         sim_legends.append(mlines.Line2D([], [], color='k', linestyle=sim_line_styles[j], label=labReactor.carrier_label))
#     ax.legend(handles=[exp_legend] + sim_legends, loc='upper right', framealpha=0.9)

#     fig.savefig('temporal_'+'_'.join(case_paths)+'.'+format, dpi=300)

def plot_labReactors_exp_eff(ax, eff_fig_name, **kwargs):
    if eff_fig_name is not None:
        exp_path = join('validation',f'Leion2008_{eff_fig_name}.csv')
        df = pd.read_csv(exp_path, skiprows=2, names=get_wpd_columns(exp_path))
        for i in range(1,9):
            ax.plot(df[f'cycle_{i}_X'], df[f'cycle_{i}_Y'], color='k', **kwargs)

def plot_labReactor_eff(ax, case, first_ind_treshold = 0.98, **kwargs):
    omega = case.omega()
    gamma = case.gamma()
    first_ind = gamma.loc[gamma > first_ind_treshold*gamma.max()].index[0]
    min_ind = min(omega.shape[0], gamma.shape[0])
    # ax.plot(omega.iloc[:min_ind], gamma.iloc[:min_ind], **kwargs)
    ax.plot(omega.iloc[:min_ind].loc[first_ind:], gamma.iloc[:min_ind].loc[first_ind:], **kwargs)


def plot_labReactors_eff(ax, cases: list, label_func=None, indent = 0, **kwargs):
    for i, case in enumerate(cases):
        label = label_func(case) if label_func is not None else None
        plot_labReactor_eff(ax, case, label = label, **kwargs)