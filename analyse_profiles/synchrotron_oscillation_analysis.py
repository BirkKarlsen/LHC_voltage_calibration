'''
File to find and analyse the synchrotron oscillation from measurements of profiles.

Author: Birk Emil Karlsen-Bæck
'''

# Import
import numpy as np
import matplotlib.pyplot as plt
import h5py

import utility_files.data_utilities as dut


# Options
PLT_PROF = True

# Parameters
fdir = f'../data_files/voltage_calibration_sorted/'
f = f'profile_V05MV_QL20k_C1B1_short_emittance.h5'  # File names on form: profile_V05MV_QL20k_C1B1_short_emittance.h5
N_shots = 10
shot_size = 2*6
sample_rate = 1

# Retrieve profiles
profiles = np.zeros((N_shots, shot_size))



for i in range(N_shots):
    profiles, _ = dut.get_profile_data(f, fdir)

# Analyse each profile
f_s = np.zeros(N_shots)
dipole_osc = np.zeros((N_shots, shot_size))

for i in range(N_shots):
    pass


plt.show()