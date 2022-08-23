'''
File to simulate the LHC conditions during the Voltage Calibration MD.

Author: Birk Emil Karlsen-Bæck
'''

# Arguments -----------------------------------------------------------------------------------------------------------
import argparse

parser = argparse.ArgumentParser(description="This file simulates the LHC during the voltage calibration MD of 2022.")

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
parser.add_argument("--save_to", '-st', type=str,
                    help="Option for custom name to the save-to-folder.")

args = parser.parse_args()

# Options -------------------------------------------------------------------------------------------------------------
LXPLUS = False
IMP = bool(args.impedance)
if IMP:
    imp_str = '_with_impedance'
else:
    imp_str = ''

# Imports -------------------------------------------------------------------------------------------------------------
import numpy as np
import os
from datetime import date
from scipy.constants import c
import matplotlib.pyplot as plt

import utility_files.data_utilities as dut
import utility_files.mathematical_relations as mre

from blond.input_parameters.rf_parameters import RFStation
from blond.input_parameters.ring import Ring
from blond.beam.beam import Beam, Proton
from blond.beam.profile import Profile, CutOptions
from blond.trackers.tracker import RingAndRFTracker
from blond.impedances.impedance_sources import InputTable
from blond.impedances.impedance import InducedVoltageFreq, TotalInducedVoltage
from blond.trackers.utilities import separatrix


# Parameters ----------------------------------------------------------------------------------------------------------
# Accelerator parameters
C = 26658.883                       # Machine circumference [m]
p_s = 450e9                         # Synchronous momentum [eV/c]
h = 35640                           # Harmonic number [-]
gamma_t = 53.606713                 # Transition gamma [-]
alpha = 1./gamma_t/gamma_t          # First order mom. comp. factor [-]
V = args.rf_voltage * 1e6           # RF voltage [V]
dphi = 0                            # Phase modulation/offset [rad]

# Beam parameters
N_p = args.intensity * 1e9          # Bunch intensity [p/b]
E_err = args.e_err * 1e6            # Injection energy error [eV]

# Simulation parameters
N_t = args.n_turns                  # Number of turns
N_m = int(1e6)                      # Number of macroparticles
freqRes = 200e3                     # Frequency resolution [Hz]

if LXPLUS:
    lxdir = "/afs/cern.ch/work/b/bkarlsen/Simulation_Files/LHC_voltage_calibration/"
else:
    lxdir = '../'

data_files_dir = f'data_files/'

if args.emittance:
    emitt_str = f'small'
else:
    emitt_str = f'nominal'

if args.save_to is not None:
    save_to = args.save_to
else:
    save_to = f'{emitt_str}_emittance_int{N_p * 1e-8:.0f}e8' \
              f'_v{V * 1e-3:.0f}kV_dE{E_err * 1e-6:.0f}MeV{imp_str}_{N_t}turns/'

# Objects for simulation ----------------------------------------------------------------------------------------------
print('Initializing Objects...\n')
# LHC ring
ring = Ring(C, alpha, p_s, Proton(), n_turns=N_t)

# 400MHz RF station
rfstation = RFStation(ring, [h], [V], [dphi])

print(rfstation.t_rev)

# Beam
beam = Beam(ring, N_m, N_p)

# Beam Profile
profile = Profile(beam, CutOptions(cut_left=-1.5 * rfstation.t_rf[0, 0],
                                   cut_right=2.5 * rfstation.t_rf[0, 0],
                                   n_slices=4 * 2**7))

BQM = Profile(beam, CutOptions(cut_left=-0.75 * rfstation.t_rf[0, 0],
                               cut_right=1.75 * rfstation.t_rf[0, 0],
                               n_slices=250))

# Fetching the beam
for file in os.listdir(lxdir + data_files_dir + 'generated_beams/'):
    if emitt_str in file:
        beam_str = file

beam_data = np.load(lxdir + data_files_dir + 'generated_beams/' + beam_str)
# Delta t to compensate for difference in bucket length between the LHC and SPS
Dt = (((2 * np.pi * mre.R_SPS)/(mre.h_SPS * c * mre.beta)) - rfstation.t_rf[0, 0])/2

beam.dE = beam_data[0, :] + E_err
beam.dt = beam_data[1, :] - Dt
profile.track()

# Impedance
if IMP:
    n_necessary = 57418             # Necessary indices to keep when we want to resolve up to 50 GHz
    imp_data = np.loadtxt(lxdir + 'impedance_models/Zlong_Allthemachine_450GeV_B1_LHC_inj_450GeV_B1.dat', skiprows=1)
    imp_table = InputTable(imp_data[:n_necessary, 0], imp_data[:n_necessary, 1], imp_data[:n_necessary, 2])

    ind_volt_freq = InducedVoltageFreq(beam, profile, [imp_table])
    total_Vind = TotalInducedVoltage(beam, profile, [ind_volt_freq])
else:
    total_Vind = None

rftracker = RingAndRFTracker(rfstation, beam, Profile=profile,
                             interpolation=True, TotalInducedVoltage=total_Vind)

# Set up directories for saving results -------------------------------------------------------------------------------
today = date.today()
sim_dir = f'data_files/{today.strftime("%b-%d-%Y")}/{save_to}'
if not os.path.exists(lxdir + sim_dir):
    os.makedirs(lxdir + sim_dir)


# Simulate ------------------------------------------------------------------------------------------------------------
dt_track = 2
dt_plot = 1000

bunch_pos = np.zeros(N_t//2)
bunch_pos_COM = np.zeros(N_t//2)
bunch_length = np.zeros(N_t//2)
time_since_injection = np.zeros(N_t//2)
j = 0

for i in range(N_t):
    if i == 0:
        print('Entered For-Loop!\n')

    rftracker.track()
    profile.track()
    BQM.track()
    if IMP:
        total_Vind.induced_voltage_sum()

    if i % dt_track == 0:
        bunch_pos[j], bunch_length[j] = dut.extract_bunch_position(BQM.bin_centers, BQM.n_macroparticles,
                                                                   wind_len=2.5)
        bunch_pos_COM[j] = dut.bunch_position_from_COM(BQM.bin_centers, BQM.n_macroparticles)
        time_since_injection[j] = np.sum(rfstation.t_rev[:i])
        j += 1

    if i % dt_plot == 0:
        print(f'Turn {i}')
        dut.save_profile(BQM, i, lxdir + sim_dir)
        dut.plot_profile(BQM, i, lxdir + sim_dir)


        dts = np.linspace(-1.25e-9, (2.5 + 1.25) * 1e-9, 1000)
        des = separatrix(ring, rfstation, dts)
        dut.plot_phase_space(beam, des, dts)
        plt.savefig(lxdir + sim_dir + f'phase_space_{i}.png')


        dut.plot_bunch_position(bunch_pos, time_since_injection, j - 1, lxdir + sim_dir)
        dut.plot_bunch_position(bunch_pos_COM, time_since_injection, j - 1, lxdir + sim_dir, COM=True)
        dut.plot_bunch_length(bunch_length, time_since_injection, j - 1, lxdir + sim_dir)

        dut.save_array(bunch_pos, 'bunch_position', lxdir + sim_dir)
        dut.save_array(bunch_length, 'bunch_length', lxdir + sim_dir)
        dut.save_array(bunch_pos_COM, 'bunch_position_com', lxdir + sim_dir)
        dut.save_array(time_since_injection, 'time_since_injection', lxdir + sim_dir)





