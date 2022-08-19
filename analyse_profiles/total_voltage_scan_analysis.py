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

init_freqs = np.zeros((len(voltages), 2))
final_freqs = np.zeros((len(voltages), 2))

init_bl = np.zeros((len(voltages), 2))
final_bl = np.zeros((len(voltages), 2))

plt_volt = [13]

# Retrieve synchrotron frequencies
for i in range(len(voltages)):
    print(f'Analysing V = {voltages[i]} MV...')
    freqs_init, freqs_final, ibl, fbl = dut.analyse_profiles_all_cavities(V=voltages[i], QL=60, emittance='nominal',
                                                                          fdir=fdir, T_rev=T_rev,
                                                                          turn_constant=turn_constant,
                                                                          init_osc_length=init_osc_length,
                                                                          final_osc_start=final_osc_start,
                                                                          add='_corr',
                                                                          plt1=voltages[i] in plt_volt,
                                                                          plt2=voltages[i] in plt_volt)

    init_freqs[i, :] = freqs_init
    final_freqs[i, :] = freqs_final
    init_bl[i, :] = ibl
    final_bl[i, :] = fbl

init_VB1 = mre.RF_voltage_from_synchrotron_frequency(init_freqs[:, 0])
init_VB2 = mre.RF_voltage_from_synchrotron_frequency(init_freqs[:, 1], eta=mre.eta2)
final_VB1 = mre.RF_voltage_from_synchrotron_frequency(final_freqs[:, 0])
final_VB2 = mre.RF_voltage_from_synchrotron_frequency(final_freqs[:, 1], eta=mre.eta2)


if PLT_DATA:
    plt.figure()
    plt.title('Total Voltage After Correction')
    plt.plot(voltages, init_freqs[:, 0], color='b')
    plt.plot(voltages, final_freqs[:, 0], color='b', linestyle='--')
    plt.plot(voltages, init_freqs[:, 1], color='r')
    plt.plot(voltages, final_freqs[:, 1], color='r', linestyle='--')

    plt.xlabel('RF Voltage [MV]')
    plt.ylabel('Synchrotron Frequency [Hz]')

    V_s = 1e-6
    fig, ax1 = plt.subplots()
    ax1.set_title('Total Voltage After Correction')
    ax1.plot(voltages, init_VB1 * V_s/voltages, color='b')
    ax1.plot(voltages, final_VB1 * V_s/voltages, color='b', linestyle='--')
    ax1.plot(voltages, init_VB2 * V_s/voltages, color='r')
    ax1.plot(voltages, final_VB2 * V_s/voltages, color='r', linestyle='--')
    ax1.set_xlabel('RF Voltage [MV]')
    ax1.set_ylabel(r'$V_\textrm{ant}/V_\textrm{set}$ [-]')
    ax1.grid()
    ax1.set_xlim((voltages[0], voltages[-1]))
    dummy_lines = []
    linestyles = ['-', '--']
    for i in range(2):
        dummy_lines.append(ax1.plot([], [], c="black", ls=linestyles[i])[0])
    lines = ax1.get_lines()
    legend2 = plt.legend([dummy_lines[i] for i in [0, 1]], ["Initial", "Final"])


    fig, ax1 = plt.subplots()
    ax1.set_title('Total Voltage Scan, Bunch Length')
    ax1.plot(voltages, init_bl[:, 0], color='b')
    ax1.plot(voltages, final_bl[:, 0], color='b', linestyle='--')
    ax1.plot(voltages, init_bl[:, 1], color='r')
    ax1.plot(voltages, final_bl[:, 1], color='r', linestyle='--')
    ax1.set_xlabel('RF Voltage [MV]')
    ax1.set_ylabel(r'$\tau_b$ [ns]')
    ax1.grid()
    ax1.set_xlim((voltages[0], voltages[-1]))

    dummy_lines = []
    linestyles = ['-', '--']
    for i in range(2):
        dummy_lines.append(ax1.plot([], [], c="black", ls=linestyles[i])[0])
    lines = ax1.get_lines()
    legend2 = plt.legend([dummy_lines[i] for i in [0, 1]], ["Initial", "Final"])

plt.show()