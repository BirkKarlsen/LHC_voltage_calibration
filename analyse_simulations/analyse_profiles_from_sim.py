'''
File to analyse the profiles from simulations.

Author: Birk Emil Karlsen-Bæck
'''

# Import
import numpy as np
import matplotlib.pyplot as plt
import os

import utility_files.data_utilities as dut
import utility_files.mathematical_relations as mre


ddir = f'../data_files/Aug-26-2022/small_emittance_int90e8_v500kV_dE0MeV_30000turns/'
fn = f'bunch_profiles.npy'

profiles = np.load(ddir + fn)
print(profiles.shape)
plt.plot(profiles)
plt.show()