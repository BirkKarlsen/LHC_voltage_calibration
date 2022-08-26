'''
File to analyse error in inferred voltage as a function of total voltage.

Author: Birk Emil Karlsen-Bæck
'''

# Import
import numpy as np
import matplotlib.pyplot as plt

import utility_files.data_utilities as dut
import utility_files.mathematical_relations as mre

plt.rcParams.update({
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{fourier}',
        'font.family': 'serif',
        'font.size': 16
    })

# Options
PLT_FWHM = True
PLT_COM = True
PLT_TOT = False

# Parameters
scanmode = 2
emittance = 'nominal'

dEs = np.array([10, 50, 100])    # [MeV]
intensity = 90                  # e8 [#]
turns = 30000

exclude = 3999

data_dir = f'../data_files/{emittance}_emittance_scanmode_{scanmode}/'

if scanmode == 1:
    Vs = np.array([0.5, 1.0, 1.5])
else:
    Vs = np.array([4, 5, 6, 7, 8, 9, 10, 11, 12])

mV_fwhm = np.zeros((dEs.shape[0], Vs.shape[0]))
mV_com = np.zeros((dEs.shape[0], Vs.shape[0]))
bli = np.zeros((dEs.shape[0], Vs.shape[0]))
blf = np.zeros((dEs.shape[0], Vs.shape[0]))

for j in range(dEs.shape[0]):
    for i in range(Vs.shape[0]):
        fit_dit = dut.get_sim_fit(emittance, intensity, Vs[i] * 1e3, dEs[j], turns, data_dir, mode='fwhm')
        mV_fwhm[j, i] = mre.RF_voltage_from_synchrotron_frequency(fit_dit['freq'])

        fit_dit = dut.get_sim_fit(emittance, intensity, Vs[i] * 1e3, dEs[j], turns, data_dir, mode='com')
        mV_com[j, i] = mre.RF_voltage_from_synchrotron_frequency(fit_dit['freq'])

        bli[j, i], blf[j, i] = dut.get_sim_init_and_final_bunch_lengths(emittance, intensity,
                                                                        Vs[i] * 1e3, dEs[j], turns, data_dir)

Vfrac_fwhm = mV_fwhm / Vs
Vfrac_com = mV_com / Vs

if PLT_TOT:
    fig, ax = plt.subplots()
    fig.suptitle('Voltage error as a function of RF Voltage')
    mark = ['--', '-.', '-']
    for i in range(dEs.shape[0]):
        ax.plot(Vs, Vfrac_fwhm[i, :] * 1e-6, linestyle=mark[i], label=f'FWHM {dEs[i]}', color='black')
        ax.plot(Vs, Vfrac_com[i, :] * 1e-6, linestyle=mark[i], label=f'COM {dEs[i]}', color='g')

    ax.set_xlabel('Set point [MV]')
    ax.set_ylabel(r'$V_\textrm{ant}/V_\textrm{set}$ [-]')
    ax.legend()


plt.show()