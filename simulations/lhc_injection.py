'''
File to simulate LHC flat bottom with a single pilot bunch with a certain energy offset.

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
from blond.trackers.tracker import RingAndRFTracker, FullRingAndRF
from blond.trackers.utilities import separatrix

# Parameters ----------------------------------------------------------------------------------------------------------
# Accelerator parameters
C = 26658.883                       # Machine circumference [m]
p_s = 450e9                         # Synchronous momentum [eV/c]
h = 35640                           # Harmonic number [-]
gamma_t = 53.606713                 # Transition gamma [-]
alpha = 1./gamma_t/gamma_t          # First order mom. comp. factor [-]
V = 0.5e6                            # RF voltage [V]
dphi = 0                            # Phase modulation/offset [rad]

# Beam parameters
bl = 0.75e-9                        # Bunch length [s]              # TODO: Find value
N_p = 6e9                           # Bunch intensity [p/b]         # TODO: Find value
E_err = 20e6                        # Injection energy error [eV]
mu = 2                              # Binomial exponent

# Simulation parameters
N_t = 5000                          # Number of turns
N_m = int(1e6)                      # Number of macroparticles

# Objects for simulation ----------------------------------------------------------------------------------------------
print('Initializing Objects...\n')
# LHC ring
ring = Ring(C, alpha, p_s, Proton(), n_turns=N_t)

# 400MHz RF station
rfstation = RFStation(ring, [h], [V], [dphi])

# Beam
beam = Beam(ring, N_m, N_p)

# Beam Profile
profile = Profile(beam, CutOptions(cut_left=-0.5 * rfstation.t_rf[0, 0],
                                   cut_right=1.5 * rfstation.t_rf[0, 0],
                                   n_slices=2 * 2**7))


rftracker = RingAndRFTracker(rfstation, beam, Profile=profile)
full_rftracker = FullRingAndRF([rftracker])

distribution_options_list = {'bunch_length': bl,
                             'type': 'binomial',
                             'density_variable': 'Hamiltonian',
                             'bunch_length_fit': 'fwhm',
                             'exponent': mu}

matched_from_distribution_function(beam, full_rftracker, distribution_type='binomial', distribution_exponent=mu,
                                   bunch_length_fit=bl, bunch_length=bl)

beam.dE += E_err

# Simulate ------------------------------------------------------------------------------------------------------------
dts = np.linspace(-1.25e-9, (2.5 + 1.25) * 1e-9, 1000)
des = separatrix(ring, rfstation, dts)
plt_sns = False
dt_plot = 100
dt_track = 100

bunch_pos = np.zeros(N_t)
bunch_pos_COM = np.zeros(N_t)
bunch_length = np.zeros(N_t)

for i in range(N_t):
    rftracker.track()
    profile.track()

    bunch_pos[i], bunch_length[i] = dut.extract_bunch_position(profile.bin_centers, profile.n_macroparticles)
    bunch_pos_COM[i] = dut.bunch_position_from_COM(profile.bin_centers, profile.n_macroparticles)

    if i % dt_track == 0:
        print(f'Turn {i}')

    if i % dt_plot == 0 and plt_sns:
        plt.figure()
        plt.plot(profile.bin_centers, profile.n_macroparticles)

        data = {r'$\Delta E$': beam.dE * 1e-6, r'$\Delta t$': beam.dt * 1e9}
        cp = sns.color_palette('coolwarm', as_cmap=True)
        sns.displot(data, x=r'$\Delta t$', y=r'$\Delta E$', cbar=True, cmap=cp, vmin=0, vmax=150)
        plt.xlabel(r'$\Delta t$ [ns]')
        plt.ylabel(r'$\Delta E$ [MeV]')
        plt.xlim((-1.25, 2.5 + 1.25))
        plt.ylim((-600, 600))

        plt.plot(dts * 1e9, des * 1e-6, color='black')
        plt.plot(dts * 1e9, -des * 1e-6, color='black')

        plt.show()

plt.figure()
t = np.linspace(0, N_t, N_t)
bunch_pos -= np.mean(bunch_pos)
bunch_pos_COM *= 1e9
bunch_pos_COM -= np.mean(bunch_pos_COM)

plt.plot(t * rfstation.t_rev[0], bunch_pos, label='FWHM')
plt.plot(t * rfstation.t_rev[0], bunch_pos_COM, label='COM')
popt = dut.fit_sine_curve_fit(bunch_pos, t)

#plt.plot(t, popt[1] * np.sin(popt[0] * t + popt[2]))
#plt.plot(t, 0.15 * np.sin(0.009 * t))
fit_param = dut.fit_sin(t, bunch_pos)
fit_param_COM = dut.fit_sin(t, bunch_pos_COM)

plt.plot(t * rfstation.t_rev[0], fit_param['fitfunc'](t), label='fit_sin')

plt.legend()
freq_theory = mre.synchrotron_frequency(V)
print(f'Expected from theory {freq_theory:.4f} Hz')
print(f"FWHM position {fit_param['freq'] / rfstation.t_rev[0]:.4f} Hz")
print(f"COM position {fit_param_COM['freq'] / rfstation.t_rev[0]:.4f} Hz")

plt.show()