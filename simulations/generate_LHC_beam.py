'''
File to generate the beam in the LHC using the measured profiles.

Author: Birk Emil Karlsen-Bæck
'''


# Imports -------------------------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.fft as spf

import utility_files.data_utilities as dut
import utility_files.mathematical_relations as mre

from blond.input_parameters.rf_parameters import RFStation
from blond.input_parameters.ring import Ring
from blond.beam.beam import Beam, Proton
from blond.beam.profile import Profile, CutOptions
from blond.beam.distributions import matched_from_line_density
from blond.trackers.tracker import RingAndRFTracker, FullRingAndRF
from blond.impedances.impedance import TotalInducedVoltage, InducedVoltageFreq
from blond.impedances.impedance_sources import InputTable
from blond.trackers.utilities import separatrix

from SPS.impedance_scenario import scenario, impedance2blond

# Options
SAVE = True

# Parameters ----------------------------------------------------------------------------------------------------------
# Accelerator parameters
C = 2 * np.pi * 1100.009            # Machine circumference [m]
p_s = 450e9                         # Synchronous momentum [eV/c]
h = 4620                            # Harmonic number [-]
gamma_t = 18.0                      # Transition gamma [-]
alpha = 1./gamma_t/gamma_t          # First order mom. comp. factor [-]
V = 4.9e6                           # 200 MHz RF voltage [V]
V_800 = 0.10 * V                    # 800 MHz RF voltage [V]
dphi = 0                            # 200 MHz Phase modulation/offset [rad]
dphi_800 = np.pi                    # 800 MHz Phase modulation/offset [rad]


# Beam parameters
N_p = 9e9                           # Bunch intensity [p/b]
mu = 1.5

# Simulation parameters
N_t = 1                             # Number of turns
N_m = int(1e6)                      # Number of macroparticles

# Parameters for the SPS Impedance Model
freqRes = 43.3e3                                # Frequency resolution [Hz]
modelStr = "futurePostLS2_SPS_f1.txt"           # Name of Impedance Model

# Options -------------------------------------------------------------------------------------------------------------
fdir = f'../data_files/voltage_calibration_sorted/'
emittance = 'short'
n_samples = 250

bl, bl_std = dut.find_bunch_length(fdir=fdir, emittance=emittance, n_samples=n_samples)
bl *= 1e-9
bl_std *= 1e-9

profiles, ts, prof_id = dut.get_first_profiles(fdir, emittance, n_samples)

if emittance == 'nominal':
    em = 0.7
    em_str = 'nominal'
else:
    em = 0.026
    em_str = 'small'


# Objects -------------------------------------------------------------------------------------------------------------

# SPS Ring
ring = Ring(C, alpha, p_s, Proton(), n_turns=1)


# RF Station
rfstation = RFStation(ring, [h, 4 * h], [V, V_800], [dphi, dphi_800], n_rf=2)

# Beam
beam = Beam(ring, N_m, N_p)

# Profile
profile = Profile(beam, CutOptions=CutOptions(cut_left=rfstation.t_rf[0, 0] * (-1.5),
            cut_right=rfstation.t_rf[0, 0] * (1.5),
            n_slices=int(round(2 ** 7 * (3)))))

# SPS Impedance Model
impScenario = scenario(modelStr)
impModel = impedance2blond(impScenario.table_impedance)

impFreq = InducedVoltageFreq(beam, profile, impModel.impedanceList, freqRes)
SPSimpedance_table = InputTable(impFreq.freq,impFreq.total_impedance.real*profile.bin_size,
                                    impFreq.total_impedance.imag*profile.bin_size)
impedance_freq = InducedVoltageFreq(beam, profile, [SPSimpedance_table],
                                       frequency_resolution=freqRes)
total_imp = TotalInducedVoltage(beam, profile, [impedance_freq])



# Tracker Object without SPS OTFB
SPS_rf_tracker = RingAndRFTracker(rfstation, beam, TotalInducedVoltage=None,
                                  CavityFeedback=None, Profile=profile, interpolation=True)
SPS_tracker = FullRingAndRF([SPS_rf_tracker])

ids = dut.find_weird_bunches(profiles, ts, PLOT=True)
profiles = dut.set_profile_reference(profiles, new_reference=0, sample=25)
profiles = dut.center_profiles(profiles, ts, rfstation.t_rf[0, 0]/2)
profiles = dut.renormalize_profiles(profiles, ts)
print(ts.shape)
n_pr = 74
n_pr1 = 0
plt.figure()
plt.title(f'Profiles')
plt.plot(ts[:, n_pr1:n_pr], profiles[:, n_pr1:n_pr])

if emittance == 'nominal':
    plt.figure()
    plt.title('Only nominal emittance bunches')
    plt.plot(np.delete(ts, ids, axis=1), np.delete(profiles, ids, axis=1))

    ts = np.delete(ts, ids, axis=1)
    profiles = np.delete(profiles, ids, axis=1)

mean_profile = np.mean(profiles, axis=1)
std_profile = np.std(profiles, axis=1)

plt.figure()
plt.title('Average Bunch')
plt.plot(ts, profiles, alpha=0.3, color='black')
plt.plot(ts[:, 0], mean_profile, color='r', linestyle='--')
plt.fill_between(ts[:, 0], mean_profile, mean_profile + std_profile, color='r', alpha=0.5)
plt.fill_between(ts[:, 0], mean_profile, mean_profile - std_profile, color='r', alpha=0.5)

line_density_input = {'time_line_den': ts[:, 0],
                      'line_density': mean_profile}
matched_from_line_density(beam, SPS_tracker, line_density_type='user_input',
                          line_density_input=line_density_input, TotalInducedVoltage=total_imp)



# The Weird Profiles --------------------------------------------------------------------------------------------------
times_b1 = np.array([2025, 2028, 2031, 2037, 2040])
times_b2 = np.array([2026, 2028, 2032, 2037, 2041])
weird_original_files = []

for i in range(len(times_b1)):
    weird_original_files.append(f'PROFILE_B1_b1_20220625{times_b1[i]}')
    weird_original_files.append(f'PROFILE_B2_b321_20220625{times_b2[i]}')

orig_dir = f'../data_files/2022-06-25_voltageCalibration/'
wprofiles, wts, wids = dut.retrieve_profile_measurements_based_on_file_names(weird_original_files,
                                                                             fdir=orig_dir, )
wprofiles = dut.set_profile_reference(wprofiles, new_reference=0, sample=25)
wprofiles = dut.center_profiles(wprofiles, wts, rfstation.t_rf[0, 0]/2)
wprofiles = dut.renormalize_profiles(wprofiles, wts)

plt.figure()
plt.plot(ts, profiles, color='r', alpha=0.5)
plt.plot(wts, wprofiles, color='b', alpha=0.5)
# 40 Gsamples per second


profile.track()

plt.figure()
plt.plot(profile.bin_centers, profile.n_macroparticles)
plt.plot(ts[:, 0], mean_profile * np.max(profile.n_macroparticles)/np.max(mean_profile))


file_name = f'../data_files/generated_beams/' \
            f'generated_beam_{em_str}_emittance_bl_{bl * 1e12:.0f}ps_mu{mu * 10:.0f}d.npy'
generated_beam = np.zeros((2, N_m))
generated_beam[0, :] = beam.dE
generated_beam[1, :] = beam.dt
if SAVE:
    np.save(file_name, generated_beam)

dts = np.linspace(-1.25e-9, (2.5 + 1.25) * 1e-9, 1000)
des = separatrix(ring, rfstation, dts)
dut.plot_phase_space(beam, des, dts)

plt.show()