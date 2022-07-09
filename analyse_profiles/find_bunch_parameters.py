'''
File to find the correct parameters of the bunches from the first turn of the LHC Voltage Calibration MD.

Author: Birk Emil Karlsen-Bæck
'''

# Import
import numpy as np
import matplotlib.pyplot as plt

import utility_files.data_utilities as dut
import utility_files.mathematical_relations as mre

# Parameters
fdir = f'../data_files/voltage_calibration_sorted/'
emittance = 'short'
n_samples = 250

profiles, ts = dut.get_first_profiles(fdir, emittance, n_samples)

plt.figure()
plt.plot(ts[:, 0], profiles[:, 0])


N_bunches, Bunch_positions, Bunch_peaks, Bunch_lengths, Bunch_intensities, Bunch_positionsFit, \
           Bunch_peaksFit, Bunch_Exponent, Goodness_of_fit = dut.getBeamPattern(ts[:, 0], profiles, heightFactor=30,
                                                                                wind_len=5, fit_option='fwhm')

bl_avg = np.mean(Bunch_lengths)
bl_std = np.std(Bunch_lengths)
mu_avg = np.mean(Bunch_Exponent)
mu_std = np.std(Bunch_Exponent)
print(bl_avg, bl_std)
print(mu_avg, mu_std)







plt.show()