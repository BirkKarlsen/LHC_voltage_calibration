'''
File to relabel the MD data to include the voltage and emittance of the bunch and in general sort the data.

Author: Birk Emil Karlsen-Bæck
'''

# Imports
import numpy as np
import matplotlib.pyplot as plt
import os
import h5py

import utility_files.data_utilities
import import_data as id

# Options
PLT_DATA = False


# Parameters
data_directory = f'../data_files/'
original_folder = data_directory + f'2022-06-25_voltageCalibration/'
new_folder = data_directory + f'voltage_calibration_sorted/'

date_string = f'20220625'

beam_line = 2
time_string = f'2041'

voltage = '120MV'
cavity = f'allB{beam_line}_corr'
emittance = 'nominal'
QL = '60k'

if beam_line == 1:
    bucket = 1
else:
    bucket = 321


# File names
old_file = f'PROFILE_B{beam_line}_b{bucket}_{date_string}{time_string}'
old_file = id.find_file_in_folder(old_file, original_folder)

new_file = f'profile_V{voltage}_QL{QL}_C{cavity}_{emittance}_emittance.h5'

# Sort file
os.system(f'cp {original_folder + old_file} {new_folder + new_file}')

print(f'Retrieved {old_file}')
print(f'Saved successfully as {new_file}')