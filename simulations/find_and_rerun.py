'''
File to find and rerun simulations that did not run properly during some parameter scan on some day.

Author: Birk Emil Karlsen-Bæck
'''

# Import --------------------------------------------------------------------------------------------------------------
import numpy as np
import os
import argparse


# Parse arguments -----------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="File to retrieve files from parameter scans.")

# Arguments related to retrieving and saving
parser.add_argument("--date_str", "-ds", type=str, default='Aug-10-2022/',
                    help="Option to specify the date of the simulations, default is 'Mar-17-2022/'.")
parser.add_argument("--rerun", '-rr', type=int, default=0,
                    help="Option to rerun simulations that were not found, default is not (0)")

# Arguments related to finding the specific set of simulations
parser.add_argument("--emittance", '-em', type=int, default=1,
                    help="Option to either have small emittance (True) or nominal (False).")
parser.add_argument("--intensity", '-in', type=float, default=9.0,
                    help="Intensity of the injected bunch, default is 9.0e9.")
parser.add_argument("--scan_mode", '-sm', type=int, default=1,
                    help="Different parameter scan modes.")

args = parser.parse_args()

# Define parameters for script ----------------------------------------------------------------------------------------
N_t = 30000
EMIT = args.emittance
N_p = args.intensity
IMP = 0

if args.scan_mode == 1:
    Vs = np.array([0.5, 1.0, 1.5])
    Es = np.linspace(0, 50, 11)
elif args.scan_mode == 2:
    Vs = np.array([1.5, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    Es = np.linspace(0, 100, 21)

if EMIT:
    emit_str = 'small'
else:
    emit_str = 'nominal'

mst_dir = os.getcwd()[:-len('simulations')]
sub_dir = os.getcwd()[:-len('Simulation_Files/LHC_voltage_calibration/simulations')] \
          + 'Submittion_Files/LHC_voltage_calibration/'

TBT_PARAM = args.tbt_param

data_dir = mst_dir + 'data_files/'


# Implement the changes made by the parser
data_dir += args.date_str

# Find data -----------------------------------------------------------------------------------------------------------

# Search for the data within the given directory
for i in range(len(Vs)):
    for j in range(len(Es)):
        sim_folder_i = f'{emit_str}_emittance_int{N_p * 10:.0f}e8' \
                       f'_v{Vs[i] * 1e3:.0f}kV_dE{Es[j]:.0f}MeV_{N_t}turns/'
        sim_dir_i = data_dir + sim_folder_i

        file_i = None

        if not os.path.isdir(sim_dir_i[:-1]):
            print(f'{sim_folder_i} was not found!\n')

            if args.rerun:
                print('Rerunning this simulation')
                sub_file_names = sim_folder_i[:-1] + '.sub'
                os.system(f'condor_submit {sub_dir}{sub_file_names}')