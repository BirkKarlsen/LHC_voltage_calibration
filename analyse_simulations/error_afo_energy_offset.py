'''
File to study the difference in inferred RF voltage as a function of injection energy error.

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

# Parameters
scanmode = 2
emittance = 'nominal'

V = 12000                       # [kV]
intensity = 90                  # e8 [#]
turns = 30000

exclude = 3999

data_dir = f'../data_files/{emittance}_emittance_scanmode_{scanmode}/'

if scanmode == 1:
    dEs = np.linspace(0, 50, 11)
else:
    dEs = np.linspace(0, 100, 21)

mV_fwhm = np.zeros(dEs.shape)
mV_com = np.zeros(dEs.shape)

for i in range(dEs.shape[0]):
    sim_str_i = dut.get_sim_name(emittance, intensity, V, dEs[i], turns)

    bp_str = 'bunch_position_' + sim_str_i + '.npy'
    time_str = 'time_since_injection_' + sim_str_i + '.npy'
    bp = np.load(data_dir + bp_str)[:-exclude]
    time = np.load(data_dir + time_str)[:-exclude]
    fit_dict = dut.fit_sin(time, bp)
    mV_fwhm[i] = mre.RF_voltage_from_synchrotron_frequency(fit_dict['freq']) * 1e-6

    bp_str = 'bunch_position_com_' + sim_str_i + '.npy'
    time_str = 'time_since_injection_' + sim_str_i + '.npy'
    bp = np.load(data_dir + bp_str)[:-exclude] * 1e9
    time = np.load(data_dir + time_str)[:-exclude]
    fit_dict = dut.fit_sin(time, bp)
    mV_com[i] = mre.RF_voltage_from_synchrotron_frequency(fit_dict['freq']) * 1e-6



fig, ax = plt.subplots()
fig.suptitle('Voltage as a function of injection energy error')
ax.plot(dEs, mV_fwhm, 'x', label='FWHM', color='black')
ax.plot(dEs, mV_com, 'x', label='COM', color='g')

ax.set_xlabel('Injection Error [MeV]')
ax.set_ylabel('Measured Voltage [MV]')
ax.legend()



plt.show()