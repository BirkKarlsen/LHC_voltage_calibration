'''
File to launch the voltage and injection error scans.

Author: Birk Emil Karlsen-Bæck
'''

# Imports -------------------------------------------------------------------------------------------------------------
import numpy as np
import argparse
import os

# Arguments -----------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="This file launches parameter scans in RF voltage and injection energy "
                                             "error.")

parser.add_argument("--emittance", '-em', type=int, default=1,
                    help="Option to either have small emittance (True) or nominal (False).")
parser.add_argument("--intensity", '-in', type=float, default=9.0,
                    help="Intensity of the injected bunch, default is 9.0e9.")
parser.add_argument("--scan_mode", '-sm', type=int, default=1,
                    help="Different parameter scan modes.")

args = parser.parse_args()

# Parameters ----------------------------------------------------------------------------------------------------------
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

print(f'Launching {len(Vs) * len(Es)} simulations...')
print(f'Injection energy error: Max = {Es[-1]:.1f} MeV, Min = {Es[0]:.1f} MeV, Step = {Es[2] - Es[1]:.1f} MeV')
print(f'RF voltage: Max = {Vs[-1]:.1f} MV, Min = {Vs[0]:.1f} MV, Step = {Vs[2] - Vs[1]:.1f} MV')
print()

current_working_directory = os.getcwd() + '/'
sim_file = 'launch_simulation.py'

for i in range(len(Vs)):
    for j in range(len(Es)):
        argument_str_i = f'-nt {N_t} -em {EMIT} -im {IMP} -in {N_p} -ee {Es[j]} -rv {Vs[i]}'
        os.system(f'python {current_working_directory}{sim_file} {argument_str_i}')

print('\nDone!')

