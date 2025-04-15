import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate

def read_IC(field):
    with open(os.path.join("0",field)) as f:
         for line in f:
            if "internalField" in line:
                elements = line.replace(";"," ").split()
                return float(elements[elements.index("uniform")+1])

def read_var(s, varname):
    elements = s.replace(";"," ").split()
    return float(elements[elements.index(varname)+1])

with open(os.path.join("caseDictExpanded")) as f:
    for line in f:
        # if line.startswith("rho_ox"):
        #     rho_ox = read_var(line, "rho_ox")
        if line.startswith("solidsMass"):
            m_ilm = float(read_var(line, "solidsMass"))
        elif line.startswith("initialConversion"):
            initialConversion = float(read_var(line, "initialConversion"))
        elif line.startswith("rhoOxidized"):
            rho_ox = float(read_var(line, "rhoOxidized"))
        elif line.startswith("Ro"):
            Ro = float(read_var(line, "Ro"))
        elif line.startswith("domainVolume"):
            v_domain = float(read_var(line, "domainVolume"))
        elif line.startswith("alpha_solids"):
            alpha_s = float(read_var(line, "alpha_solids"))

# m_ilm = 0.000050 # [kg] - mass of ilmenite
W_O = 15.9994e-3
MW_FeO = 71.844
MW_Fe2O3 = 159.687
# v_domain = 0.001
# alpha_s = 1e-6
g_species = ("N2", "O2", "CO", "CO2", "H2O", "H2", "CH4")
s_species = ("Fe2O3", "FeO", "inerts")

case_specie = os.path.split(os.getcwd())[-1].split('_')[-1]
case_carrier = os.path.split(os.getcwd())[-1].split('_')[-2]

path_s = os.path.join("postProcessing","volIntegrateFn(fields=(rho.solids),weightField=alpha.solids)","0","volFieldValue.dat")
path_g = os.path.join("postProcessing","volIntegrateFn(fields=(rho.gas),weightField=alpha.gas)","0","volFieldValue.dat")
path_rhos = os.path.join("postProcessing","volIntegrateFn(fields=(rho.solids),operation=sum)","0","volFieldValue.dat")
path_g_sp  = os.path.join("postProcessing","volIntegrateFn(fields=({}.gas),weightFields=(alpha.gasrho.gas))","0","volFieldValue.dat")
path_s_sp = os.path.join("postProcessing",'volIntegrateFn(fields=({}.solids),weightFields=(alpha.solidsrho.solids))',"0","volFieldValue.dat")
path_probes = os.path.join("postProcessing","probesFn","0")
exp_dat = np.loadtxt(f'../validation/Abad2011_fig5_{case_carrier}_{case_specie}.dat', unpack=True)
dat_g_sp = dict()
for sp in g_species:
    dat_g_sp[sp] = np.loadtxt(path_g_sp.format(sp), unpack=True) # [kg]
dat_s_sp = dict()
for sp in s_species:
    dat_s_sp[sp] = np.loadtxt(path_s_sp.format(sp), unpack=True) # [kg]

dat_conversion = np.loadtxt(os.path.join(path_probes, "conversion"), unpack=True)
dat_s = np.loadtxt(path_s, unpack=True)
dat_g = np.loadtxt(path_g, unpack=True)
dat_rhos = np.loadtxt(path_rhos, unpack=True)
threshold = 1e-3
if case_specie == 'O2':
    ind = (dat_conversion[1][1:]>1-threshold).argmax()+1
else:
    ind = (dat_conversion[1][1:]<threshold).argmax()+1
time = dat_conversion[0][ind]
# if ind>1:
#     print(f"Full conversion time = {time}s, experimental is {exp_dat[0][-1]}s, diff={100*(time-exp_dat[0][-1])/exp_dat[0][-1]:5.2f}%")
# else:
#     print(f"Full conversion time > {dat_conversion[0][-1]}s, experimental is {exp_dat[0][-1]}s")

# This does not work anymore, since a C++ code has beed added
# def read_density(specie):
#     with open(os.path.join("constant","physicalProperties.solids")) as f:
#         lineFound = False
#         for line in f:
#             if line.startswith(specie):
#                 lineFound = True
#             if lineFound:
#                 if "rho" in line:
#                     elements = line.replace(";"," ").split()
#                     return float(elements[elements.index("rho")+1])
# def read_var(s, varname):
#     elements = s.replace(";"," ").split()
#     return float(elements[elements.index(varname)+1])
# with open(os.path.join("constant","physicalProperties.solids")) as f:
#     for line in f:
#         # if line.startswith("rho_ox"):
#         #     rho_ox = read_var(line, "rho_ox")
#         if line.startswith("MW_Fe2O3"):
#             MW_Fe2O3 = read_var(line, "MW_Fe2O3")
#         elif line.startswith("MW_FeO"):
#             MW_FeO = read_var(line, "MW_FeO")

rho_red = 2*MW_FeO/MW_Fe2O3*rho_ox
rho_i = rho_ox

fig, axs = plt.subplots(figsize=[15,11], nrows=3, ncols=3)


ax = axs[0][0]
if case_specie == 'O2':
    ax.plot(dat_conversion[0][1:], dat_conversion[1][1:], label='sim')
else:
    ax.plot(dat_conversion[0][1:], 1-dat_conversion[1][1:], label='sim')
ax.plot(exp_dat[0], exp_dat[1], 'kx', label='exp')
ax.legend()
ax.set_ylim([0,1])
ax.set_title('conversion')


ax = axs[1][0]
datafiles=['Fe2O3.solids', 'FeO.solids', 'inerts.solids']
sum = 0
for filename in datafiles:
    path = os.path.join(path_probes, filename)
    dat = np.loadtxt(path, unpack=True)
    ax.plot(dat[0], dat[1], label=filename)
    sum += dat[1]
if not np.allclose(sum,1):
    print('Warning: mass fractions does not always sum up to 1!')
    print('t=', dat[0])
    print('sumY=', sum)
ax.set_ylim([0,1])
ax.set_xlabel('time, s')
ax.set_ylabel('Solids')
ax.set_title('Mass fractions')
ax.legend()


ax = axs[2][0]
datafiles=['CO2.gas', 'CO.gas', 'H2.gas', 'H2O.gas', 'O2.gas', 'CH4.gas', 'N2.gas']
for filename in datafiles:
    path = os.path.join(path_probes, filename)
    dat = np.loadtxt(path, unpack=True)
    ax.plot(dat[0], dat[1], label=filename)
# ax.axhline(y=read_IC(case_specie+'.gas'), color='black', linestyle='--', label=f'exp {case_carrier} {case_specie}') # removed until (if?) read_IC is fixed
# ax.set_ylim([0,1])
ax.set_xlabel('time, s')
ax.set_ylabel('Gases')
ax.legend()


ax = axs[1][1]
datafiles=['X_Fe2O3.solids', 'X_FeO.solids', 'X_inerts.solids']
sum = 0
for filename in datafiles:
    path = os.path.join(path_probes, filename)
    dat = np.loadtxt(path, unpack=True)
    ax.plot(dat[0], dat[1], label=filename)
    sum += dat[1]
if not np.allclose(sum,1):
    print('Warning: mole fractions does not always sum up to 1!')
    print('t=', dat[0])
    print('sumX=', sum)
ax.set_ylim([0,1])
ax.set_xlabel('time, s')
ax.set_ylabel('Solids')
ax.set_title('Molar fractions')
ax.legend()


ax = axs[2][1]
datafiles=['X_CO2.gas', 'X_CO.gas', 'X_H2.gas', 'X_H2O.gas', 'X_O2.gas', 'X_CH4.gas', 'X_N2.gas' ]
for filename in datafiles:
    path = os.path.join(path_probes, filename)
    dat = np.loadtxt(path, unpack=True)
    ax.plot(dat[0], dat[1], label=filename)
ax.set_xlabel('time, s')
ax.legend()

ax = axs[1][2]
ax.plot(dat_rhos[0], np.loadtxt(os.path.join(path_probes, "Fe2O3.solids"),  unpack=True)[1]*dat_rhos[1]/rho_ox,  label='Fe2O3')
ax.plot(dat_rhos[0], np.loadtxt(os.path.join(path_probes, "FeO.solids"),    unpack=True)[1]*dat_rhos[1]/rho_red, label='FeO')
ax.plot(dat_rhos[0], np.loadtxt(os.path.join(path_probes, "inerts.solids"), unpack=True)[1]*dat_rhos[1]/rho_i,   label='inerts')
Y_Fe2O3 = np.loadtxt(os.path.join(path_probes, "Fe2O3.solids"),  unpack=True)[1]
Y_FeO   = np.loadtxt(os.path.join(path_probes, "FeO.solids"),    unpack=True)[1]
Y_i     = np.loadtxt(os.path.join(path_probes, "inerts.solids"), unpack=True)[1]
v_Fe2O3 = Y_Fe2O3*dat_rhos[1]/rho_ox
v_FeO   = Y_FeO  *dat_rhos[1]/rho_red
v_i     = Y_i    *dat_rhos[1]/rho_i
sum = v_Fe2O3 + v_FeO + v_i
# print(f"v_Fe2O3={v_Fe2O3[0]}={Y_Fe2O3[0]}*{dat_rhos[1][0]}/{rho_ox}")
# print(f"v_FeO  ={v_FeO[0]  }={Y_FeO[0]  }*{dat_rhos[1][0]}/{rho_red}")
# print(f"v_i    ={v_i[0]    }={Y_i[0]    }*{dat_rhos[1][0]}/{rho_i}")

if not np.allclose(sum,1):
    print('Warning: volume fractions does not always sum up to 1!')
    print('t=', dat[0])
    print('sumv=', sum)
# print(np.loadtxt(os.path.join(path_probes, "Fe2O3.solids"),  unpack=True)[1]*dat_rhos[1]/rho_ox)
# print(np.loadtxt(os.path.join(path_probes, "FeO.solids"),    unpack=True)[1]*dat_rhos[1]/rho_red)
# print(np.loadtxt(os.path.join(path_probes, "inerts.solids"), unpack=True)[1]*dat_rhos[1]/rho_i)
ax.set_title('Volume fractions')
ax.set_ylim([0,1])
ax.set_xlabel('time, s')
ax.legend()


ax = axs[2][2]
# datafiles=['CO2.gas', 'CO.gas', 'H2.gas', 'H2O.gas', 'O2.gas', 'CH4.gas', 'N2.gas']
for sp in g_species:
    # path = os.path.join(path_probes, filename)
    path = path_g_sp.format(sp)
    dat = np.loadtxt(path, unpack=True)
    ax.plot(dat[0], 1e6*(dat[1]-dat[1][0]), label=sp)
# ax.axhline(y=read_IC(case_specie+'.gas'), color='black', linestyle='--', label=f'exp {case_carrier} {case_specie}')
# ax.set_ylim([0,1])
ax.set_xlabel('time, s')
ax.set_ylabel('mass change, mkg')
ax.legend()



ax = axs[0][1]
# if case_specie == 'O2':
w = 1-(rho_ox-dat_rhos[1])/rho_ox
# else:
#     w = 1+(dat_s[1]-rho_ox)/rho_ox
ax.plot(dat_s[0], w)
ax.set_title('mass-based conversion')


ax = axs[0][2]
ax.plot(dat_s[0], 1e6*(dat_s[1]-dat_s[1][0]), label='solids')
ax.plot(dat_g[0], -1e6*(dat_g[1]-dat_g[1][0]), '--', label='gas')
ax.set_title('mass balance')
ax.set_ylabel('mass change, mg')
ax.legend()
ax2 = ax.twinx()
ax2.plot(dat_s[0], 1e6*dat_s[1], label='ilmenite mass, mg')
ax2.set_ylabel('total mass, mg')

fig.savefig("pp.png")

dir_path = os.path.dirname(os.path.realpath(__file__))

if case_specie == 'O2':
    f = interpolate.interp1d(dat_conversion[0], (dat_conversion[1]), kind='cubic', fill_value="extrapolate")
else:
    f = interpolate.interp1d(dat_conversion[0], (1-dat_conversion[1]), kind='cubic', fill_value="extrapolate")



def calc_reaction_heat_per_O_atom(carrier, reaction, case_path = None):
    if case_path is None:
            case_path = '_'.join(['0D', carrier, reaction])
    Ro = 0.033 if carrier == 'act' else 0.04
    
    rhomO = Ro/W_O # [mol/kg] - molar density of reacting O atoms
    
    path_Hs = os.path.join("postProcessing","volIntegrateFn(fields=(h.solids),weightFields=(alpha.solidsrho.solids))","0","volFieldValue.dat")
    path_Hg = os.path.join("postProcessing","volIntegrateFn(fields=(h.gas),weightFields=(alpha.gasrho.gas))","0","volFieldValue.dat")
    path_p  = os.path.join("postProcessing","probesFn","0","p")

    dat_Hs = np.loadtxt(path_Hs, unpack=True)
    dat_Hg = np.loadtxt(path_Hg, unpack=True)
    dat_p   = np.loadtxt(path_p, unpack=True)

    
    change_g = dat_Hg[1][-1]-dat_Hg[1][0] - v_domain*(dat_p[1][-1]-dat_p[1][0])
    change_s = dat_Hs[1][-1]-dat_Hs[1][0] - alpha_s*v_domain*(dat_p[1][-1]-dat_p[1][0])

    change_total = change_g+change_s

    m = m_ilm if reaction != 'O2' else m_ilm/(1-Ro)
    # print(f"{case_path:10s}: {change_total/m:10.2f} J/kg oxidized ilmenite")
    # print(f"{case_path:10s}: {change_total/m/rhomO/1000:10.2f} kJ/mol O atom")
    return change_total/m/rhomO

dQ = calc_reaction_heat_per_O_atom(case_carrier, case_specie) # compare it to the data in the presentation or Ilmenite_Abad_2010.mcdx
L0 = np.sum(exp_dat[1]-f(exp_dat[0]))/len(exp_dat[0])
L1 = np.sum(np.abs(exp_dat[1]-f(exp_dat[0])))/len(exp_dat[0])
L2 = np.sqrt(np.sum((exp_dat[1]-f(exp_dat[0]))**2))/len(exp_dat[0])
# print(f"L0 = {L0:.6f} - mean deviation from exp")
# print(f"L1 = {L1:.6f} - mean absolute deviation from exp")
# print(f"L2 = {L2:.6f} - root mean square deviation from exp")

cwd = os.path.basename(os.getcwd())
print(f"{cwd:20}{exp_dat[0][-1]:-10.1f}{time:-10.1f}{L0:-10.6f}{L1:-10.6f}{L2:-10.6f}{dQ/1000:-10.2f}")