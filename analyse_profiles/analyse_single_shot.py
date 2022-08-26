'''
File to analyse single shots from the MD.

Author: Birk Emil Karlsen-Bæck
'''

# Import --------------------------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

import utility_files.data_utilities as dut
import utility_files.mathematical_relations as mre

plt.rcParams.update({
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{fourier}',
        'font.family': 'serif',
        'font.size': 16
    })

# Parameters
data_dir = '../data_files/voltage_calibration_sorted/'

# Shot configuration --------------------------------------------------------------------------------------------------
V = 4.0                 # [MV]
QL = 60                 # k [-]
cavity = 'all'              # Cavity number
beam = 2                # Beamline
emittance = 'nominal'     # Bunch Emittance mode
add = '_corr'                # Additional configurations

if emittance == 'short':
    emit_str = 'small'
else:
    emit_str = 'nominal'


PLT_SHOT = True
PLT_BUNCH_POS = True

# Plot the entire shot ------------------------------------------------------------------------------------------------
file_name = dut.get_sorted_files(V, QL, cavity, beam, emittance, add)
profile, t = dut.get_profile_data(file_name, data_dir)
profile = dut.set_profile_reference(profile, new_reference=0, sample=25)

if PLT_SHOT:
    n_start = 10000
    n_end = 16000
    n_jump = 500

    print('Profile dimensions:', profile.shape)

    fig, ax = plt.subplots()
    fig.suptitle(f'Shot in {cavity}B{beam}, V = {V} MV, {emit_str} emit.')

    profiles_cut = profile[:, n_start//2:n_end//2:n_jump//2]

    #plt.plot(t, profile[:, n_start//2:n_end//2:n_jump//2])

    n_p = profiles_cut.shape[1]
    c = np.linspace(n_start, n_end - n_jump, n_p, dtype=int)
    cc = np.arange(0, n_p)
    cmap = mpl.cm.get_cmap('jet', n_p)
    dummie_cax = ax.scatter(c, c, c=cc, cmap=cmap)

    ax.cla()
    for i in range(n_p):
        ax.plot((t - t[-1]/2) * 1e9, profiles_cut[:, i], color=cmap(i))

    ax.set_xlabel(r'$\Delta t$ [ns]')
    ax.set_ylabel(r'Arb. Units [-]')
    plt.setp(ax.get_yticklabels(), visible=False)

    cbar = fig.colorbar(dummie_cax, ticks=c)
    cbar.set_label(r'Turn', rotation=270, labelpad=60)

    cbar.ax.get_yaxis().set_ticks([])
    for j, lab in enumerate(c):
        cbar.ax.text(30, 1.85 * (2 * j + 1) / 4.0, lab, ha='center', va='center')



if PLT_BUNCH_POS:
    init_osc_length = 3000
    final_osc_start = 6000

    dict_init, dict_final, bpos, blen, ts = dut.analyze_profile(profile, t, T_rev=mre.T_rev, turn_constant=2,
                                                                init_osc_length=init_osc_length,
                                                                final_osc_start=final_osc_start,
                                                                mode='fwhm', wind_len=2.5)
    ts_turns = np.linspace(0, 2 * len(ts), len(ts))

    fig, ax = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    fig.suptitle(f'Synchro. Osc. in {cavity}B{beam}, V = {V} MV, {emit_str} emit.')

    ax[0].plot(ts_turns, bpos, color='black', label='Data')
    ax[0].plot(ts_turns, dict_init['fitfunc'](ts), color='g', label='Fit 1', linestyle='-.')
    ax[0].plot(ts_turns, dict_final['fitfunc'](ts), color='magenta', label='Fit 2', linestyle='--')
    ax[0].set_xlim((0, 10000))
    ax[0].set_xlabel('Turns [-]')
    ax[0].set_ylabel(r'$\Delta t_s$ [ns]')

    ax[1].plot(ts_turns, bpos, color='black', label='Data')
    ax[1].plot(ts_turns, dict_init['fitfunc'](ts), color='g', label='Fit 1', linestyle='-.')
    ax[1].plot(ts_turns, dict_final['fitfunc'](ts), color='magenta', label='Fit 2', linestyle='--')
    plt.setp(ax[1].get_yticklabels(), visible=False)
    ax[1].set_xlim((50000, ts_turns[-1]))
    ax[1].set_xlabel('Turns [-]')

    handles, labels = ax[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='center', bbox_to_anchor=(0.5, 0.05), ncol=3)


    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(7, 7))
    fig.suptitle(f'FWHM Fit to {cavity}B{beam}, V = {V} MV, {emit_str} emit.')

    ax[0].plot(ts_turns, bpos, color='black')
    ax[0].set_ylabel(r'$\Delta t_s$ [ns]')
    ax[0].set_xlim((ts_turns[0], ts_turns[-1]))
    plt.setp(ax[0].get_xticklabels(), visible=False)

    ax[1].plot(ts_turns, blen, color='black')
    ax[1].set_ylabel(r'$\tau_b$ [ns]')
    ax[1].set_xlabel('Turns [-]')



plt.show()