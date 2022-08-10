'''
File to launch a simulation using the lhc_md_simulation.py script in lxplus.

Author: Birk Emil Karlsen-Bæck
'''

# Import --------------------------------------------------------------------------------------------------------------
import argparse
import os

# Arguments -----------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="This file launches simulations using the lhc_md_simulation.py script "
                                             "in lxplus.")

# File managment
parser.add_argument("--sim_name", "--sn", type=str,
                    help="Option to give custom name to the simulation. If none is given, then a default name is given "
                         "based on the simulation configurations.")
parser.add_argument("--save_dir", "-sd", type=str,
                    help="Name of directory to save the results to.")

# Arguments for LHC simulation file
parser.add_argument("--n_turns", '-nt', type=int, default=30000,
                    help="The number of turns to simulates, default is 30000")
parser.add_argument("--e_err", '-ee', type=float, default=0.0,
                    help="Injection energy error in MeV into the LHC.")
parser.add_argument("--emittance", '-em', type=int, default=1,
                    help="Option to either have small emittance (1) or nominal (0).")
parser.add_argument("--rf_voltage", '-rv', type=float, default=0.5,
                    help="RF voltage in MV, default is 0.5 MV")
parser.add_argument("--intensity", '-in', type=float, default=9.0,
                    help="Intensity of the injected bunch, default is 9.0e9.")
parser.add_argument("--impedance", '-im', type=int, default=0,
                    help="Option to include a LHC flat-bottom impedance model.")

args = parser.parse_args()

# Parameter Values ----------------------------------------------------------------------------------------------------
N_p = args.intensity
V = args.rf_voltage
E_err = args.e_err
N_t = args.n_turns
IMP = args.impedance
EMIT = args.emittance

imp_str = ''
bash_dir = '/afs/cern.ch/work/b/bkarlsen/Simulation_Files/LHC_voltage_calibration/bash_files/'
sub_dir = '/afs/cern.ch/work/b/bkarlsen/Submittion_Files/LHC_voltage_calibration/'

if args.emittance:
    emitt_str = f'small'
else:
    emitt_str = f'nominal'

if IMP:
    imp_str = '_with_impedance'
else:
    imp_str = ''

# Make necessary preparations for Sims --------------------------------------------------------------------------------
print('\nMaking shell scripts...')

bash_file_names = f'{emitt_str}_emittance_int{N_p * 10:.0f}e8' \
                  f'_v{V * 1e3:.0f}kV_dE{E_err:.0f}MeV{imp_str}_{N_t}turns.sh'
sub_file_names = f'{emitt_str}_emittance_int{N_p * 10:.0f}e8' \
                 f'_v{V * 1e3:.0f}kV_dE{E_err:.0f}MeV{imp_str}_{N_t}turns.sub'
file_names = f'{emitt_str}_emittance_int{N_p * 10:.0f}e8' \
             f'_v{V * 1e3:.0f}kV_dE{E_err:.0f}MeV{imp_str}_{N_t}turns'

save_dir = f'{emitt_str}_emittance_int{N_p * 10:.0f}e8' \
           f'_v{V * 1e3:.0f}kV_dE{E_err:.0f}MeV{imp_str}_{N_t}turns/'

if args.sim_name is not None:
    bash_file_names = args.sim_name + '.sh'
    sub_file_names = args.sim_name + '.sub'
    file_names = args.sim_name

if args.save_dir is not None:
    save_dir = args.save_dir


# Make bash file
print(f'Launching {file_names}...')
os.system(f'touch {bash_dir}{bash_file_names}')

bash_content = f'#!/bin/bash\n' \
               f'source /afs/cern.ch/user/b/bkarlsen/.bashrc\n' \
               f'python /afs/cern.ch/work/b/bkarlsen/Simulation_Files/LHC_voltage_calibration/simulations/' \
               f'lhc_md_simulation.py ' \
               f'-nt {N_t} -ee {E_err} -em {EMIT} -rv {V} -in {N_p} -im {IMP} -st {save_dir}'

os.system(f'echo "{bash_content}" > {bash_dir}{bash_file_names}')
os.system(f'chmod a+x {bash_dir}{bash_file_names}')


print('\nMaking and submitting simulations...')

# Make submission file
os.system(f'touch {sub_dir}{sub_file_names}')

sub_content = f'executable = {bash_dir}{bash_file_names}\n' \
              f'arguments = \$(ClusterId)\$(ProcId)\n' \
              f'output = {bash_dir}{file_names}.\$(ClusterId)\$(ProcId).out\n' \
              f'error = {bash_dir}{file_names}.\$(ClusterId)\$(ProcId).err\n' \
              f'log = {bash_dir}{file_names}.\$(ClusterId)\$(ProcId).log\n' \
              f'+JobFlavour = \\"testmatch\\"\n' \
              f'queue'

os.system(f'echo "{sub_content}" > {sub_dir}{sub_file_names}')
os.system(f'chmod a+x {sub_dir}{sub_file_names}')

os.system(f'condor_submit {sub_dir}{sub_file_names}')



