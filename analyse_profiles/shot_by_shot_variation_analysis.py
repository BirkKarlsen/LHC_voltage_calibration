'''
File to compute the shot-by-shot variation for the cavities that had the best and worst match with computed
synchrotron frequency.

Author: Birk Emil Karlsen-Bæck
'''

# Import
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot

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
cavitiesB1 = np.array([4, 6])
cavitiesB2 = np.array([7, 8])
T_rev = 8.892465516509656e-05
turn_constant = 2

init_osc_length = 2000
final_osc_start = 0

acquisitons = np.linspace(1, 10, 10, dtype=int)

init_freqsB1 = np.zeros((len(acquisitons), len(cavitiesB1)))
init_freqsB2 = np.zeros((len(acquisitons), len(cavitiesB2)))

final_freqsB1 = np.zeros((len(acquisitons), len(cavitiesB1)))
final_freqsB2 = np.zeros((len(acquisitons), len(cavitiesB2)))

all_data_init = np.concatenate((init_freqsB2, init_freqsB1), axis=1)
means_init = np.mean(all_data_init, axis=0)
print(all_data_init.shape)
print(means_init.shape)

# Find synchrotron frequencies for each acquisition

for i in range(len(acquisitons)):
    add_str = '_Acq' + str(round(0.1 * acquisitons[i], 1)).replace('.', '')
    print('Acquisition:', acquisitons[i])

    freqs_init, freqs_final = dut.analyse_synchrotron_frequency_cavity_by_cavity(V=1.5, QL=60, cavitiesB1=cavitiesB1,
                                                                                 cavitiesB2=cavitiesB2,
                                                                                 emittance='nominal', fdir=fdir,
                                                                                 T_rev=T_rev,
                                                                                 turn_constant=turn_constant,
                                                                                 init_osc_length=init_osc_length,
                                                                                 final_osc_start=final_osc_start,
                                                                                 add=add_str)

    init_freqsB1[i, :] = freqs_init[:, 0]
    init_freqsB2[i, :] = freqs_init[:, 1]

    final_freqsB1[i, :] = freqs_final[:, 0]
    final_freqsB2[i, :] = freqs_final[:, 1]


if PLT_DATA:
    plt.figure()
    plt.title('Voltage Shot-by-shot Variation')

    plt.plot(acquisitons, init_freqsB1[:, 0], color='r')
    plt.plot(acquisitons, init_freqsB1[:, 1], color='g')

    plt.plot(acquisitons, final_freqsB1[:, 0], color='r', linestyle='--')
    plt.plot(acquisitons, final_freqsB1[:, 1], color='g', linestyle='--')


    plt.plot(acquisitons, init_freqsB2[:, 0], color='b')
    plt.plot(acquisitons, init_freqsB2[:, 1], color='black')

    plt.plot(acquisitons, final_freqsB2[:, 0], color='b', linestyle='--')
    plt.plot(acquisitons, final_freqsB2[:, 1], color='black', linestyle='--')

    plt.xlabel('Acquisition [-]')
    plt.ylabel('Synchrotron Frequency [Hz]')

    init_VB1 = mre.RF_voltage_from_synchrotron_frequency(init_freqsB1)
    final_VB1 = mre.RF_voltage_from_synchrotron_frequency(final_freqsB1)
    init_VB2 = mre.RF_voltage_from_synchrotron_frequency(init_freqsB2, eta=mre.eta2)
    final_VB2 = mre.RF_voltage_from_synchrotron_frequency(final_freqsB2, eta=mre.eta2)

    V_s = 1e-6
    plt.figure()
    plt.title('Voltage Shot-by-shot Variation')
    cavities = ['4B1', '6B1', '7B2', '8B2']
    all_data_init = np.concatenate((init_VB1, init_VB2), axis=1)
    means_init = np.mean(all_data_init, axis=0)
    std_init = np.std(all_data_init, axis=0)

    all_data_final = np.concatenate((final_VB1, final_VB2), axis=1)
    means_final = np.mean(all_data_final, axis=0)
    std_final = np.std(all_data_final, axis=0)

    #plt.errorbar(cavities, means_init * V_s, yerr=std_init * V_s, fmt='o', color='r',
    #             label=r'Initial', capsize=6)
    plt.errorbar(cavities, means_final * V_s, yerr=std_final * V_s, fmt='o', color='b',
                 label=r'Final', capsize=6)

    plt.grid()

    plt.xlabel('Cavity [-]')
    plt.ylabel(r'$V_\textrm{ant}$ [MV]')



plt.show()