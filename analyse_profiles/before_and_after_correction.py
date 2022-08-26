'''
File to plot the comparison of offset between cavities before and after correction.

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

# Parameters
fdir = f'../data_files/voltage_calibration_sorted/'
cavities = np.linspace(1, 8, 8, dtype=int)
T_rev = 8.892465516509656e-05
turn_constant = 2


# Without correction

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

V_init1_wo = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 0])
V_init2_wo = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 1], eta=mat.eta2)
V_final1 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 0])
V_final2 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 1], eta=mat.eta2)

f_s1 = mat.synchrotron_frequency(1.5e6, eta=mat.eta1)
f_s2 = mat.synchrotron_frequency(1.5e6, eta=mat.eta2)


# With correction
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

V_init1_w = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 0])
V_init2_w = mat.RF_voltage_from_synchrotron_frequency(freqs_init[:, 1], eta=mat.eta2)
V_final1 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 0])
V_final2 = mat.RF_voltage_from_synchrotron_frequency(freqs_final[:, 1], eta=mat.eta2)

f_s1 = mat.synchrotron_frequency(1.5e6, eta=mat.eta1)
f_s2 = mat.synchrotron_frequency(1.5e6, eta=mat.eta2)



fig, ax1 = plt.subplots()

ax1.set_title(f'$V$ = {1.5} MV, $Q_L$ = {60}k')
ax1.plot(cavities, freqs_init[:, 0], 'x', color='b')
ax1.plot(cavities, freqs_final[:, 0], 'D', color='b')

ax1.plot(cavities, freqs_init[:, 1], 'x', color='r')
ax1.plot(cavities, freqs_final[:, 1], 'D', color='r')

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
fig, ax1 = plt.subplots()
ax1.set_title(f'Measured Voltage, $V$ = {1.5} MV, $Q_L$ = {60}k')
ax1.plot(cavities, V_init1_wo * V_s, 'x', color='b')
ax1.plot(cavities, V_init1_w * V_s, 'D', color='b')

ax1.plot(cavities, V_init2_wo * V_s, 'x', color='r')
ax1.plot(cavities, V_init2_w * V_s, 'D', color='r')

ax1.set_xlabel('Cavity Number [-]')
ax1.set_ylabel('RF Voltage [MV]')
ax1.set_xticks(cavities)
ax1.grid()

dummy_lines = []
linestyles = ['x', 'D']
for i in range(2):
    dummy_lines.append(ax1.plot([], [], linestyles[i], c="black")[0])
lines = ax1.get_lines()
legend2 = ax1.legend([dummy_lines[i] for i in [0, 1]], ["Before", "After"])

plt.show()