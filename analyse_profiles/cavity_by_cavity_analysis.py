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
VOLTAGE_15MV = False
VOLTAGE_15MV_CORR = False
NOMINAL_EMITTANCE = True

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

    plt_cav1 = None
    plt_cav2 = None

    freqs_init, freqs_final, ibl, fbl = dut.analyze_profiles_cavity_by_cavity(V=0.5, QL=20, cavitiesB1=cavities,
                                                                              cavitiesB2=cavities,
                                                                              emittance='short', fdir=fdir,
                                                                              T_rev=T_rev,
                                                                              turn_constant=turn_constant,
                                                                              init_osc_length=init_osc_length,
                                                                              final_osc_start=final_osc_start,
                                                                              plt_cav1=plt_cav1, plt_cav2=plt_cav2,
                                                                              fbl_mean=1000)

    V_init1 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 0])
    V_init2 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 1], eta=mat.eta2)
    V_final1 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 0])
    V_final2 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 1], eta=mat.eta2)

    f_s1 = mat.synchrotron_frequency(0.5e6, eta=mat.eta1)
    f_s2 = mat.synchrotron_frequency(0.5e6, eta=mat.eta2)

    if PLT_DATA:
        dut.plot_cavity_by_cavity_voltage(cavities, freqs_init, freqs_final, 0.5, 20)

        dut.plot_cavity_by_cavity(cavities, f'Bunch Length, $V$ = 0.5 MV, $Q_L$ = 20k',
                                  r'$\tau_b$ [ns]', ibl, fbl)

    print(f'With 0.5 MV we expect f_s = {f_s1} Hz in B1 and f_s = {f_s2} Hz in B2.')


# Analysis of 1.0 MV
if VOLTAGE_10MV:
    init_osc_length = 2000
    final_osc_start = 10000

    plt_cav1 = None
    plt_cav2 = None

    freqs_init, freqs_final, ibl, fbl = dut.analyze_profiles_cavity_by_cavity(V=1.0, QL=60, cavitiesB1=cavities,
                                                                              cavitiesB2=cavities,
                                                                              emittance='short', fdir=fdir,
                                                                              T_rev=T_rev,
                                                                              turn_constant=turn_constant,
                                                                              init_osc_length=init_osc_length,
                                                                              final_osc_start=final_osc_start,
                                                                              plt_cav1=plt_cav1, plt_cav2=plt_cav2)

    V_init1 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 0])
    V_init2 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 1], eta=mat.eta2)
    V_final1 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 0])
    V_final2 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 1], eta=mat.eta2)

    f_s1 = mat.synchrotron_frequency(1.0e6, eta=mat.eta1)
    f_s2 = mat.synchrotron_frequency(1.0e6, eta=mat.eta2)

    if PLT_DATA:
        dut.plot_cavity_by_cavity_voltage(cavities, freqs_init, freqs_final, 1.0, 60)
        dut.plot_cavity_by_cavity(cavities, f'Bunch Length, $V$ = 1.0 MV, $Q_L$ = 60k',
                                  r'$\tau_b$ [ns]', ibl, fbl)

    print(f'With 1.0 MV we expect f_s = {f_s1} Hz in B1 and f_s = {f_s2} Hz in B2.')

# Analysis of 1.5 MV
if VOLTAGE_15MV:
    init_osc_length = 2000
    final_osc_start = 10000

    plt_cav1 = None
    plt_cav2 = None

    freqs_init, freqs_final, ibl, fbl = dut.analyze_profiles_cavity_by_cavity(V=1.5, QL=60, cavitiesB1=cavities,
                                                                              cavitiesB2=cavities,
                                                                              emittance='short', fdir=fdir,
                                                                              T_rev=T_rev,
                                                                              turn_constant=turn_constant,
                                                                              init_osc_length=init_osc_length,
                                                                              final_osc_start=final_osc_start,
                                                                              plt_cav1=plt_cav1, plt_cav2=plt_cav2)

    V_init1 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 0])
    V_init2 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 1], eta=mat.eta2)
    V_final1 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 0])
    V_final2 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 1], eta=mat.eta2)

    f_s1 = mat.synchrotron_frequency(1.5e6, eta=mat.eta1)
    f_s2 = mat.synchrotron_frequency(1.5e6, eta=mat.eta2)

    if PLT_DATA:
        dut.plot_cavity_by_cavity_voltage(cavities, freqs_init, freqs_final, 1.5, 60)
        dut.plot_cavity_by_cavity(cavities, f'Bunch Length, $V$ = 1.5 MV, $Q_L$ = 60k',
                                  r'$\tau_b$ [ns]', ibl, fbl)

    print(f'With 1.5 MV we expect f_s = {f_s1} Hz in B1 and f_s = {f_s2} Hz in B2.')

# Analysis of 1.5 MV with corrections
if VOLTAGE_15MV_CORR:
    init_osc_length = 2000
    final_osc_start = 10000

    plt_cav1 = None
    plt_cav2 = None

    freqs_init, freqs_final, ibl, fbl = dut.analyze_profiles_cavity_by_cavity(V=1.5, QL=60, cavitiesB1=cavities,
                                                                              cavitiesB2=cavities,
                                                                              emittance='short', fdir=fdir,
                                                                              T_rev=T_rev,
                                                                              turn_constant=turn_constant,
                                                                              init_osc_length=init_osc_length,
                                                                              final_osc_start=final_osc_start,
                                                                              plt_cav1=plt_cav1, plt_cav2=plt_cav2,
                                                                              add='_corr')

    V_init1 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 0])
    V_init2 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 1], eta=mat.eta2)
    V_final1 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 0])
    V_final2 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 1], eta=mat.eta2)

    f_s1 = mat.synchrotron_frequency(1.5e6, eta=mat.eta1)
    f_s2 = mat.synchrotron_frequency(1.5e6, eta=mat.eta2)

    if PLT_DATA:
        dut.plot_cavity_by_cavity_voltage(cavities, freqs_init, freqs_final, 1.5, 60, add_str=', Corrected')
        dut.plot_cavity_by_cavity(cavities, f'Bunch Length, $V$ = 1.5 MV, $Q_L$ = 60k, Corrected',
                                  r'$\tau_b$ [ns]', ibl, fbl)

    print(f'With 1.5 MV we expect f_s = {f_s1} Hz in B1 and f_s = {f_s2} Hz in B2.')

# Analysis of 1.5 MV with nominal emittance
if NOMINAL_EMITTANCE:
    init_osc_length = 2000
    final_osc_start = 10000

    plt_cav1 = None
    plt_cav2 = None

    freqs_init, freqs_final, ibl, fbl = dut.analyze_profiles_cavity_by_cavity(V=1.5, QL=60, cavitiesB1=cavities,
                                                                              cavitiesB2=cavities,
                                                                              emittance='nominal', fdir=fdir,
                                                                              T_rev=T_rev,
                                                                              turn_constant=turn_constant,
                                                                              init_osc_length=init_osc_length,
                                                                              final_osc_start=final_osc_start,
                                                                              plt_cav1=plt_cav1, plt_cav2=plt_cav2)

    V_init1 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 0])
    V_init2 = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 1], eta=mat.eta2)
    V_final1 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 0])
    V_final2 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 1], eta=mat.eta2)

    f_s1 = mat.synchrotron_frequency(1.5e6, eta=mat.eta1)
    f_s2 = mat.synchrotron_frequency(1.5e6, eta=mat.eta2)

    if PLT_DATA:
        dut.plot_cavity_by_cavity_voltage(cavities, freqs_init, freqs_final, 1.5, 60, add_str=', Nom. Em.')
        dut.plot_cavity_by_cavity(cavities, f'Bunch Length, $V$ = 1.5 MV, $Q_L$ = 60k, Nom. Em.',
                                  r'$\tau_b$ [ns]', ibl, fbl)

    print(f'With 1.5 MV we expect f_s = {f_s1} Hz in B1 and f_s = {f_s2} Hz in B2.')

plt.show()