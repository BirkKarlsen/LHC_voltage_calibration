'''
File to analyze the small emittance bunches from the LHC Voltage Calibration MD.

Author: Birk Emil Karlsen-Bæck
'''

# Import
import numpy as np
import matplotlib.pyplot as plt

import utility_files.data_utilities as dut
import utility_files.mathematical_relations as mat

plt.rcParams.update({
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{fourier}',
        'font.family': 'serif',
        'font.size': 16
    })

# Options
PLT_DATA = True
VOLTAGE_05MV = False
VOLTAGE_10MV = False
VOLTAGE_15MV = True
VOLTAGE_15MV_CORR = True
NOMINAL_EMITTANCE = False

# Parameters
fdir = f'../data_files/voltage_calibration_sorted/'
cavities = np.linspace(1, 8, 8, dtype=int)
T_rev = 8.892465516509656e-05
turn_constant = 2

# Function to analyze the different cases

# Analysis of 0.5 MV
if VOLTAGE_05MV:
    init_osc_length = 2000
    final_osc_start = 10000

    freqs_init, freqs_final = dut.analyse_synchrotron_frequency_cavity_by_cavity(V=0.5, QL=20, cavitiesB1=cavities,
                                                                                 cavitiesB2=cavities,
                                                                                 emittance='short', fdir=fdir,
                                                                                 T_rev=T_rev,
                                                                                 turn_constant=turn_constant,
                                                                                 init_osc_length=init_osc_length,
                                                                                 final_osc_start=final_osc_start)

    V_init1 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 0])
    V_init2 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 1], eta=mat.eta2)
    V_final1 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 0])
    V_final2 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 1], eta=mat.eta2)

    f_s1 = mat.synchrotron_frequency(0.5e6, eta=mat.eta1)
    f_s2 = mat.synchrotron_frequency(0.5e6, eta=mat.eta2)

    if PLT_DATA:
        fig, ax1 = plt.subplots()

        ax1.set_title(f'$V$ = 0.5 MV, $Q_L$ = 20k')
        ax1.plot(cavities, freqs_init[:, 0], 'x', color='r')
        ax1.plot(cavities, freqs_final[:, 0], 'D', color='r')

        ax1.plot(cavities, freqs_init[:, 1], 'x', color='b')
        ax1.plot(cavities, freqs_final[:, 1], 'D', color='b')

        ax1.set_xlabel('Cavity Number [-]')
        ax1.set_ylabel('Synchrotron Frequency [Hz]')
        ax1.grid()
        ax1.set_xticks(cavities)

        ax2 = ax1.twinx()
        mn, mx = ax1.get_ylim()
        mn = mat.RF_voltage_from_synchrotron_frequency(mn)
        mx = mat.RF_voltage_from_synchrotron_frequency(mx)
        V_s = 1e-6
        ax2.set_ylim(mn * V_s, mx * V_s)
        ax2.set_ylabel('RF Voltage [MV]')


        V_s = 1e-6
        plt.figure()
        plt.title(f'$V$ = 0.5 MV, $Q_L$ = 20k')
        plt.plot(cavities, V_init1 * V_s, 'x', color='r')
        plt.plot(cavities, V_final1 * V_s, 'D', color='r')

        plt.plot(cavities, V_init2 * V_s, 'x', color='b')
        plt.plot(cavities, V_final2 * V_s, 'D', color='b')

        plt.xlabel('Cavity Number [-]')
        plt.ylabel('RF Voltage [MV]')
        plt.xticks(cavities)
        plt.grid()

    print(f'With 0.5 MV we expect f_s = {f_s1} Hz in B1 and f_s = {f_s2} Hz in B2.')


# Analysis of 1.0 MV
if VOLTAGE_10MV:
    init_osc_length = 2000
    final_osc_start = 10000

    freqs_init, freqs_final = dut.analyse_synchrotron_frequency_cavity_by_cavity(V=1.0, QL=60, cavitiesB1=cavities,
                                                                                 cavitiesB2=cavities,
                                                                                 emittance='short', fdir=fdir,
                                                                                 T_rev=T_rev,
                                                                                 turn_constant=turn_constant,
                                                                                 init_osc_length=init_osc_length,
                                                                                 final_osc_start=final_osc_start)

    V_init1 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 0])
    V_init2 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 1], eta=mat.eta2)
    V_final1 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 0])
    V_final2 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 1], eta=mat.eta2)

    f_s1 = mat.synchrotron_frequency(1.0e6, eta=mat.eta1)
    f_s2 = mat.synchrotron_frequency(1.0e6, eta=mat.eta2)

    if PLT_DATA:
        fig, ax1 = plt.subplots()

        ax1.set_title(f'$V$ = 1.0 MV, $Q_L$ = 60k')
        ax1.plot(cavities, freqs_init[:, 0], 'x', color='r')
        ax1.plot(cavities, freqs_final[:, 0], 'D', color='r')

        ax1.plot(cavities, freqs_init[:, 1], 'x', color='b')
        ax1.plot(cavities, freqs_final[:, 1], 'D', color='b')

        ax1.set_xlabel('Cavity Number [-]')
        ax1.set_ylabel('Synchrotron Frequency [Hz]')
        ax1.grid()
        ax1.set_xticks(cavities)

        ax2 = ax1.twinx()
        mn, mx = ax1.get_ylim()
        mn = mat.RF_voltage_from_synchrotron_frequency(mn)
        mx = mat.RF_voltage_from_synchrotron_frequency(mx)
        V_s = 1e-6
        ax2.set_ylim(mn * V_s, mx * V_s)
        ax2.set_ylabel('RF Voltage [MV]')

        V_s = 1e-6
        plt.figure()
        plt.title(f'$V$ = 1.0 MV, $Q_L$ = 60k')
        plt.plot(cavities, V_init1 * V_s, 'x', color='r')
        plt.plot(cavities, V_final1 * V_s, 'D', color='r')

        plt.plot(cavities, V_init2 * V_s, 'x', color='b')
        plt.plot(cavities, V_final2 * V_s, 'D', color='b')

        plt.xlabel('Cavity Number [-]')
        plt.ylabel('RF Voltage [MV]')
        plt.xticks(cavities)
        plt.grid()

    print(f'With 1.0 MV we expect f_s = {f_s1} Hz in B1 and f_s = {f_s2} Hz in B2.')

# Analysis of 1.5 MV
if VOLTAGE_15MV:
    init_osc_length = 2000
    final_osc_start = 10000

    freqs_init, freqs_final = dut.analyse_synchrotron_frequency_cavity_by_cavity(V=1.5, QL=60, cavitiesB1=cavities,
                                                                                 cavitiesB2=cavities,
                                                                                 emittance='short', fdir=fdir,
                                                                                 T_rev=T_rev,
                                                                                 turn_constant=turn_constant,
                                                                                 init_osc_length=init_osc_length,
                                                                                 final_osc_start=final_osc_start)

    V_init1 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 0])
    V_init2 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 1], eta=mat.eta2)
    V_final1 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 0])
    V_final2 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 1], eta=mat.eta2)

    f_s1 = mat.synchrotron_frequency(1.5e6, eta=mat.eta1)
    f_s2 = mat.synchrotron_frequency(1.5e6, eta=mat.eta2)

    if PLT_DATA:
        fig, ax1 = plt.subplots()

        ax1.set_title(f'$V$ = 1.5 MV, $Q_L$ = 60k')
        ax1.plot(cavities, freqs_init[:, 0], 'x', color='r')
        ax1.plot(cavities, freqs_final[:, 0], 'D', color='r')

        ax1.plot(cavities, freqs_init[:, 1], 'x', color='b')
        ax1.plot(cavities, freqs_final[:, 1], 'D', color='b')

        ax1.set_xlabel('Cavity Number [-]')
        ax1.set_ylabel('Synchrotron Frequency [Hz]')
        ax1.grid()
        ax1.set_xticks(cavities)

        ax2 = ax1.twinx()
        mn, mx = ax1.get_ylim()
        mn = mat.RF_voltage_from_synchrotron_frequency(mn)
        mx = mat.RF_voltage_from_synchrotron_frequency(mx)
        V_s = 1e-6
        ax2.set_ylim(mn * V_s, mx * V_s)
        ax2.set_ylabel('RF Voltage [MV]')

        V_s = 1e-6
        plt.figure()
        plt.title(f'$V$ = 1.5 MV, $Q_L$ = 60k')
        plt.plot(cavities, V_init1 * V_s, 'x', color='r')
        plt.plot(cavities, V_final1 * V_s, 'D', color='r')

        plt.plot(cavities, V_init2 * V_s, 'x', color='b')
        plt.plot(cavities, V_final2 * V_s, 'D', color='b')

        plt.xlabel('Cavity Number [-]')
        plt.ylabel('RF Voltage [MV]')
        plt.xticks(cavities)
        plt.grid()

    print(f'With 1.5 MV we expect f_s = {f_s1} Hz in B1 and f_s = {f_s2} Hz in B2.')

# Analysis of 1.5 MV with corrections
if VOLTAGE_15MV_CORR:
    init_osc_length = 2000
    final_osc_start = 10000

    freqs_init, freqs_final = dut.analyse_synchrotron_frequency_cavity_by_cavity(V=1.5, QL=60, cavitiesB1=cavities,
                                                                                 cavitiesB2=cavities,
                                                                                 emittance='short', fdir=fdir,
                                                                                 T_rev=T_rev,
                                                                                 turn_constant=turn_constant,
                                                                                 init_osc_length=init_osc_length,
                                                                                 final_osc_start=final_osc_start,
                                                                                 add='_corr')

    V_init1 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 0])
    V_init2 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 1], eta=mat.eta2)
    V_final1 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 0])
    V_final2 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 1], eta=mat.eta2)

    f_s1 = mat.synchrotron_frequency(1.5e6, eta=mat.eta1)
    f_s2 = mat.synchrotron_frequency(1.5e6, eta=mat.eta2)

    if PLT_DATA:
        fig, ax1 = plt.subplots()

        ax1.set_title(f'$V$ = 1.5 MV, $Q_L$ = 60k, Corrected')
        ax1.plot(cavities, freqs_init[:, 0], 'x', color='r')
        ax1.plot(cavities, freqs_final[:, 0], 'D', color='r')

        ax1.plot(cavities, freqs_init[:, 1], 'x', color='b')
        ax1.plot(cavities, freqs_final[:, 1], 'D', color='b')

        ax1.set_xlabel('Cavity Number [-]')
        ax1.set_ylabel('Synchrotron Frequency [Hz]')
        ax1.grid()
        ax1.set_xticks(cavities)

        ax2 = ax1.twinx()
        mn, mx = ax1.get_ylim()
        mn = mat.RF_voltage_from_synchrotron_frequency(mn)
        mx = mat.RF_voltage_from_synchrotron_frequency(mx)
        V_s = 1e-6
        ax2.set_ylim(mn * V_s, mx * V_s)
        ax2.set_ylabel('RF Voltage [MV]')

        V_s = 1e-6
        plt.figure()
        plt.title(f'$V$ = 1.5 MV, $Q_L$ = 60k, Corrected')
        plt.plot(cavities, V_init1 * V_s, 'x', color='r')
        plt.plot(cavities, V_final1 * V_s, 'D', color='r')

        plt.plot(cavities, V_init2 * V_s, 'x', color='b')
        plt.plot(cavities, V_final2 * V_s, 'D', color='b')

        plt.xlabel('Cavity Number [-]')
        plt.ylabel('RF Voltage [MV]')
        plt.xticks(cavities)
        plt.grid()

    print(f'With 1.5 MV we expect f_s = {f_s1} Hz in B1 and f_s = {f_s2} Hz in B2.')

# Analysis of 1.5 MV with nominal emittance
if NOMINAL_EMITTANCE:
    init_osc_length = 2000
    final_osc_start = 10000

    freqs_init, freqs_final = dut.analyse_synchrotron_frequency_cavity_by_cavity(V=1.5, QL=60, cavitiesB1=cavities,
                                                                                 cavitiesB2=cavities,
                                                                                 emittance='nominal', fdir=fdir,
                                                                                 T_rev=T_rev,
                                                                                 turn_constant=turn_constant,
                                                                                 init_osc_length=init_osc_length,
                                                                                 final_osc_start=final_osc_start)

    V_init1 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 0])
    V_init2 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 1], eta=mat.eta2)
    V_final1 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 0])
    V_final2 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 1], eta=mat.eta2)

    f_s1 = mat.synchrotron_frequency(1.5e6, eta=mat.eta1)
    f_s2 = mat.synchrotron_frequency(1.5e6, eta=mat.eta2)

    if PLT_DATA:
        fig, ax1 = plt.subplots()

        ax1.set_title(f'$V$ = 1.5 MV, $Q_L$ = 60k, Nominal Emittance')
        ax1.plot(cavities, freqs_init[:, 0], 'x', color='r')
        ax1.plot(cavities, freqs_final[:, 0], 'D', color='r')

        ax1.plot(cavities, freqs_init[:, 1], 'x', color='b')
        ax1.plot(cavities, freqs_final[:, 1], 'D', color='b')

        ax1.set_xlabel('Cavity Number [-]')
        ax1.set_ylabel('Synchrotron Frequency [Hz]')
        ax1.grid()
        ax1.set_xticks(cavities)

        ax2 = ax1.twinx()
        mn, mx = ax1.get_ylim()
        mn = mat.RF_voltage_from_synchrotron_frequency(mn)
        mx = mat.RF_voltage_from_synchrotron_frequency(mx)
        V_s = 1e-6
        ax2.set_ylim(mn * V_s, mx * V_s)
        ax2.set_ylabel('RF Voltage [MV]')

        V_s = 1e-6
        plt.figure()
        plt.title(f'$V$ = 1.5 MV, $Q_L$ = 60k, Nominal Emittance')
        plt.plot(cavities, V_init1 * V_s, 'x', color='r')
        plt.plot(cavities, V_final1 * V_s, 'D', color='r')

        plt.plot(cavities, V_init2 * V_s, 'x', color='b')
        plt.plot(cavities, V_final2 * V_s, 'D', color='b')

        plt.xlabel('Cavity Number [-]')
        plt.ylabel('RF Voltage [MV]')
        plt.xticks(cavities)
        plt.grid()

    print(f'With 1.5 MV we expect f_s = {f_s1} Hz in B1 and f_s = {f_s2} Hz in B2.')

plt.show()