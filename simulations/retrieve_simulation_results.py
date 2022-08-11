'''
File to retrieve all simulations results from a parameter scan and put them into a folder.

Author: Birk Emil Karlsen-Bæck
'''

# Import --------------------------------------------------------------------------------------------------------------
import numpy as np
import os
import argparse


# Parse arguments -----------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="File to retrieve files from parameter scans.")

# Arguments related to retrieving and saving
parser.add_argument("--signal_type", "-st", type=str, default='profile_data_',
                    help="The signal that is being retrived from parameter scan; default is profiles")
parser.add_argument("--save_dir", "-sd", type=str,
                    help="Name of directory and file to save the signals into.")
parser.add_argument("--date_str", "-ds", type=str, default='Aug-10-2022/',
                    help="Option to specify the date of the simulations„ default is 'Mar-17-2022/'.")
parser.add_argument("--tbt_param", "-tb", type=int, default=1,
                    help="Option to determine whether or not it is a turn-by-turn signal, default is true (1)")

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

signal_name = args.signal_type
mst_dir = os.getcwd()[:-len('utility_files')]
TBT_PARAM = args.tbt_param

save_dir = mst_dir + 'data_files/'
data_dir = mst_dir + 'data_files/'

save_name = f'{emit_str}_emittance_scanmode_{args.scan_mode}/'
if args.save_dir is not None:
    save_name = args.save_dir

# Implement the changes made by the parser
data_dir += args.date_str

# Retrieve data -------------------------------------------------------------------------------------------------------

# Make the directory that the files will be saved to
if not os.path.exists(save_dir + save_name):
    os.mkdir(save_dir + save_name)

# Search for the data within the given directory
for i in range(len(Vs)):
    for j in range(len(Es)):
        sim_folder_i = f'{emit_str}_emittance_int{N_p * 10:.0f}e8' \
                       f'_v{Vs[i] * 1e3:.0f}kV_dE{Es[j]:.0f}MeV_{N_t}turns/'
        sim_dir_i = data_dir + sim_folder_i

        if TBT_PARAM:
            # Find all files for the given signal type prefix
            file_list = []
            turns = []
            for file in os.listdir(sim_dir_i[:-1]):
                if file.startswith(signal_name):
                    file_list.append(file)
                    turns.append(file[len(signal_name):-4])
            turns = np.array(turns, dtype=int)

            # Find latest turn that the signals was recorded and saved
            final_index = np.where(turns == np.amax(turns))[0][0]
            file_i = file_list[final_index]
        else:
            for file in os.listdir(sim_dir_i[:-1]):
                if file.startswith(signal_name):
                    print(f'Found {file} in {sim_folder_i}!\n')
                    file_i = file

        os.system(f"cp {sim_dir_i}{file_i} {save_dir + save_name}{file_i[:-4]}"
                  f"_{emit_str}_emittance_int{N_p * 10:.0f}e8"
                  f"_v{Vs[i] * 1e3:.0f}kV_dE{Es[j]:.0f}MeV_{N_t}turns.npy")
