'''
File to study the nominal bunches in greater detail.

Author: Birk Emil Karlsen-Bæck
'''

# Imports
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
PLT_HIST = True

# Parameters
fdir = f'../data_files/voltage_calibration_sorted/'
emittance = 'short'
n_samples = 250

# Analysis

# Small Bunches -------------------------------------------------------------------------------------------------------
profiles, ts, prof_id = dut.get_first_profiles(fdir, emittance, n_samples)

N_bunches, Bunch_positions, Bunch_peaks, Bunch_lengths, Bunch_intensities, Bunch_positionsFit, \
    Bunch_peaksFit, Bunch_Exponent, Goodness_of_fit = dut.getBeamPattern(ts[:, 0], profiles, heightFactor=30,
                                                                         wind_len=5, fit_option='fwhm')

profiles = dut.set_profile_reference(profiles, new_reference=0, sample=25)
profiles = dut.center_profiles(profiles, ts, mre.T_rev/mre.h/2)
profiles = dut.renormalize_profiles(profiles, ts)

small_bl = Bunch_lengths

# Nominal Bunches -----------------------------------------------------------------------------------------------------
profiles, ts, prof_id = dut.get_first_profiles(fdir, 'nominal', n_samples)

ids = dut.find_weird_bunches(profiles, ts, PLOT=True)
ts = np.delete(ts, ids, axis=1)
profiles = np.delete(profiles, ids, axis=1)

N_bunches, Bunch_positions, Bunch_peaks, Bunch_lengths, Bunch_intensities, Bunch_positionsFit, \
    Bunch_peaksFit, Bunch_Exponent, Goodness_of_fit = dut.getBeamPattern(ts[:, 0], profiles, heightFactor=30,
                                                                         wind_len=5, fit_option='fwhm')

profiles = dut.set_profile_reference(profiles, new_reference=0, sample=25)
profiles = dut.center_profiles(profiles, ts, mre.T_rev/mre.h/2)
profiles = dut.renormalize_profiles(profiles, ts)

nominal_bl = Bunch_lengths



# Weird Bunches -------------------------------------------------------------------------------------------------------
times_b1 = np.array([2025, 2028, 2031, 2037, 2040])
times_b2 = np.array([2026, 2028, 2032, 2037, 2041])
weird_original_files = []

for i in range(len(times_b1)):
    weird_original_files.append(f'PROFILE_B1_b1_20220625{times_b1[i]}')
    weird_original_files.append(f'PROFILE_B2_b321_20220625{times_b2[i]}')

orig_dir = f'../data_files/2022-06-25_voltageCalibration/'
wprofiles, wts, wids = dut.retrieve_profile_measurements_based_on_file_names(weird_original_files,
                                                                             fdir=orig_dir, )

N_bunches, Bunch_positions, Bunch_peaks, Bunch_lengths, Bunch_intensities, Bunch_positionsFit, \
    Bunch_peaksFit, Bunch_Exponent, Goodness_of_fit = dut.getBeamPattern(wts[:, 0], wprofiles, heightFactor=30,
                                                                         wind_len=5, fit_option='fwhm')

wprofiles = dut.set_profile_reference(wprofiles, new_reference=0, sample=25)
wprofiles = dut.center_profiles(wprofiles, wts, mre.T_rev/mre.h/2)
wprofiles = dut.renormalize_profiles(wprofiles, wts)

wbl = Bunch_lengths


fig, ax = plt.subplots()
fig.suptitle('Distribution of Bunch Lengths')
ax.hist(small_bl, facecolor='black', alpha=0.7, label='Small')
ax.hist(nominal_bl, facecolor='g', alpha=0.7, label='Nominal')
ax.hist(wbl, facecolor='pink', alpha=0.7, label='Anomalous')
ax.legend()
ax.set_xlabel(r'$\tau_b$ [ns]')
ax.set_ylabel(r'Num. Bunches [-]')

plt.show()

