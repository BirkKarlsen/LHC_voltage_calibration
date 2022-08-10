'''
File to generate the MD bunches in the SPS.

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
from blond.beam.distributions import matched_from_distribution_function, distribution_function
from blond.beam.distributions_multibunch import matched_from_distribution_density_multibunch
from blond.trackers.tracker import RingAndRFTracker, FullRingAndRF
from blond.impedances.impedance import TotalInducedVoltage, InducedVoltageFreq
from blond.impedances.impedance_sources import InputTable
from blond.trackers.utilities import separatrix

from SPS.impedance_scenario import scenario, impedance2blond

# Parameters ----------------------------------------------------------------------------------------------------------
# Accelerator parameters
C = 2 * np.pi * 1100.009            # Machine circumference [m]
p_s = 450e9                         # Synchronous momentum [eV/c]
h = 4620                            # Harmonic number [-]
gamma_t = 18.0                      # Transition gamma [-]
alpha = 1./gamma_t/gamma_t          # First order mom. comp. factor [-]
V = 6.7e6                           # 200 MHz RF voltage [V]
V_800 = 0.19 * V                    # 800 MHz RF voltage [V]
dphi = 0                            # 200 MHz Phase modulation/offset [rad]
dphi_800 = np.pi                    # 800 MHz Phase modulation/offset [rad]

# Beam parameters
N_p = 6e9                           # Bunch intensity [p/b]         # TODO: Find value
mu = 1.5                              # Binomial exponent

# Simulation parameters
N_t = 1                             # Number of turns
N_m = int(1e6)                      # Number of macroparticles

# Parameters for the SPS Impedance Model
freqRes = 43.3e3                                # Frequency resolution [Hz]
modelStr = "futurePostLS2_SPS_noMain200TWC.txt" # Name of Impedance Model

# Options -------------------------------------------------------------------------------------------------------------
fdir = f'../data_files/voltage_calibration_sorted/'
emittance = 'nominal'
n_samples = 250

bl, bl_std = dut.find_bunch_length(fdir=fdir, emittance=emittance, n_samples=n_samples)
bl *= 1e-9
bl_std *= 1e-9

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
SPS_rf_tracker = RingAndRFTracker(rfstation, beam, TotalInducedVoltage=total_imp,
                                  CavityFeedback=None, Profile=profile, interpolation=True)
SPS_tracker = FullRingAndRF([SPS_rf_tracker])


distribution_options_list = {'bunch_length': [bl],
                             'type': 'binomial',
                             'density_variable': 'Hamiltonian',
                             'bunch_length_fit': 'fwhm',
                             'exponent': [mu],
                             'emittance': [em]}

matched_from_distribution_density_multibunch(beam, ring, SPS_tracker, distribution_options_list,
                                             1, np.zeros(1),
                                             intensity_list=[N_p],
                                             n_iterations=4, TotalInducedVoltage=total_imp)

profile.track()

plt.figure()
plt.plot(profile.bin_centers, profile.n_macroparticles)

file_name = f'../data_files/generated_beams/' \
            f'generated_beam_{em_str}_emittance_bl_{bl * 1e12:.0f}ps_mu{mu * 10:.0f}d.npy'
generated_beam = np.zeros((2, N_m))
generated_beam[0, :] = beam.dE
generated_beam[1, :] = beam.dt

np.save(file_name, generated_beam)

dts = np.linspace(-1.25e-9, (2.5 + 1.25) * 1e-9, 1000)
des = separatrix(ring, rfstation, dts)
dut.plot_phase_space(beam, des, dts)

plt.show()