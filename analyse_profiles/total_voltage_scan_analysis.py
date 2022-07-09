'''
File to analyse the synchrotron frequency for the total voltage scan that was performed as a part of the
LHC voltage calibration MD. The voltage scan was performed after the variation in voltage was compensated for.

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
PLT_DATA = True

# Parameters
fdir = f'../data_files/voltage_calibration_sorted/'
T_rev = 8.892465516509656e-05
turn_constant = 2

init_osc_length = 2000
final_osc_start = 0

voltages = np.linspace(4, 12, 9)

init_freqsB1 = np.zeros(len(voltages))
init_freqsB2 = np.zeros(len(voltages))

final_freqsB1 = np.zeros(len(voltages))
final_freqsB2 = np.zeros(len(voltages))

# Retrieve synchrotron frequencies
for i in range(len(voltages)):
    freqs_init, freqs_final = dut.analyse_synchrotron_frequency_with_all_cavities(V=voltages[i], QL=60,
                                                                                  emittance='nominal', fdir=fdir,
                                                                                  T_rev=T_rev,
                                                                                  turn_constant=turn_constant,
                                                                                  init_osc_length=init_osc_length,
                                                                                  final_osc_start=final_osc_start,
                                                                                  add='_corr')

    init_freqsB1[i] = freqs_init[0]
    init_freqsB2[i] = freqs_init[1]

    final_freqsB1[i] = freqs_final[0]
    final_freqsB2[i] = freqs_final[1]

init_VB1 = mre.RF_voltage_from_synchrotron_frequency(init_freqsB1)
init_VB2 = mre.RF_voltage_from_synchrotron_frequency(init_freqsB2, eta=mre.eta2)
final_VB1 = mre.RF_voltage_from_synchrotron_frequency(final_freqsB1)
final_VB2 = mre.RF_voltage_from_synchrotron_frequency(final_freqsB2, eta=mre.eta2)


if PLT_DATA:
    plt.figure()
    plt.title('Total Voltage After Correction')
    plt.plot(voltages, init_freqsB1, color='r')
    plt.plot(voltages, final_freqsB1, color='r', linestyle='--')
    plt.plot(voltages, init_freqsB2, color='b')
    plt.plot(voltages, final_freqsB2, color='b', linestyle='--')

    plt.xlabel('RF Voltage [MV]')
    plt.ylabel('Synchrotron Frequency [Hz]')

    V_s = 1e-6
    fig, ax1 = plt.subplots()
    ax1.set_title('Total Voltage After Correction')
    ax1.plot(voltages, init_VB1 * V_s/voltages, color='r')
    ax1.plot(voltages, final_VB1 * V_s/voltages, color='r', linestyle='--')
    ax1.plot(voltages, init_VB2 * V_s/voltages, color='b')
    ax1.plot(voltages, final_VB2 * V_s/voltages, color='b', linestyle='--')
    ax1.set_xlabel('RF Voltage [MV]')
    ax1.set_ylabel(r'$V_\textrm{ant}/V_\textrm{set}$ [-]')


plt.show()