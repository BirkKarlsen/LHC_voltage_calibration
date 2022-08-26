'''
File to analyse a single simulation.

Author Birk Emil Karlsen-Bæck
'''

# Imports -------------------------------------------------------------------------------------------------------------
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


# Options -------------------------------------------------------------------------------------------------------------
PLT_BP = True
PLT_BL = True

# Parameters ----------------------------------------------------------------------------------------------------------
scanmode = 2
emittance = 'nominal'

injection_error = 70        # [MeV]
V = 4000                     # [kV]
intensity = 90              # e8 [#]
turns = 30000

exclude = 3999

data_dir = f'../data_files/{emittance}_emittance_scanmode_{scanmode}/'

# Analysis ------------------------------------------------------------------------------------------------------------
if PLT_BP:
    sim_str = dut.get_sim_name(emittance, intensity, V, injection_error, turns)

    bp_str = 'bunch_position_' + sim_str + '.npy'
    time_str = 'time_since_injection_' + sim_str + '.npy'
    bp = np.load(data_dir + bp_str)[:-exclude]
    time = np.load(data_dir + time_str)[:-exclude]
    fit_dict = dut.fit_sin(time, bp)
    Vs = mre.RF_voltage_from_synchrotron_frequency(fit_dict['freq']) * 1e-6

    print(f'FWHM Freq {fit_dict["freq"]:.3f} Hz which is {Vs:.3f} MV')
    fig, ax = plt.subplots(2, 1, sharex=True)
    fig.suptitle('Bunch Position')
    ax[0].plot(time, bp)
    ax[0].plot(time, fit_dict['fitfunc'](time), linestyle='--')

    plt.setp(ax[0].get_xticklabels(), visible=False)

    sim_str = dut.get_sim_name(emittance, intensity, V, injection_error, turns)

    bp_str = 'bunch_position_com_' + sim_str + '.npy'
    time_str = 'time_since_injection_' + sim_str + '.npy'
    bp = np.load(data_dir + bp_str)[:-exclude] * 1e9
    time = np.load(data_dir + time_str)[:-exclude]
    fit_dict = dut.fit_sin(time, bp)
    Vs = mre.RF_voltage_from_synchrotron_frequency(fit_dict['freq']) * 1e-6

    print(f'COM Freq {fit_dict["freq"]:.3f} Hz which is {Vs:.3f} MV')
    ax[1].plot(time, bp)
    ax[1].plot(time, fit_dict['fitfunc'](time), linestyle='--')
    ax[1].set_xlabel('Time [s]')


if PLT_BL:
    sim_str = dut.get_sim_name(emittance, intensity, V, injection_error, turns)

    bl_str = 'bunch_length_' + sim_str + '.npy'
    time_str = 'time_since_injection_' + sim_str + '.npy'

    bl = np.load(data_dir + bl_str)[:-exclude]
    time = np.load(data_dir + time_str)[:-exclude]

    fig, ax = plt.subplots()
    fig.suptitle('Bunch Position')
    ax.plot(time, bl)
    ax.set_xlabel('Time [s]')

plt.show()