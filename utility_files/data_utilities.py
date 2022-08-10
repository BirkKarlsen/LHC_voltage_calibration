'''
Functions for data analysis

Author: Birk Emil Karlsen-Bæck
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from blond_common.fitting.profile import binomial_amplitudeN_fit, FitOptions
from blond_common.interfaces.beam.analytic_distribution import binomialAmplitudeN
import os
import h5py
import seaborn as sns

from scipy.stats import linregress
from scipy.optimize import curve_fit


def fwhm(x, y, level=0.5):
    offset_level = np.mean(y[0:5])
    amp = np.max(y) - offset_level
    t1, t2 = interp_f(x, y, level)
    mu = (t1 + t2) / 2.0
    sigma = (t2 - t1) / 2.35482
    popt = (mu, sigma, amp)

    return popt


def interp_f(time, bunch, level):
    bunch_th = level * bunch.max()
    time_bet_points = time[1] - time[0]
    taux = np.where(bunch >= bunch_th)
    taux1, taux2 = taux[0][0], taux[0][-1]
    t1 = time[taux1] - (bunch[taux1] - bunch_th) / (bunch[taux1] - bunch[taux1 - 1]) * time_bet_points
    t2 = time[taux2] + (bunch[taux2] - bunch_th) / (bunch[taux2] - bunch[taux2 + 1]) * time_bet_points

    return t1, t2


def getBeamPattern(timeScale, frames, heightFactor=0.015, distance=500, N_bunch_max=3564,
                     fit_option='fwhm', plot_fit=False, baseline_length=1, BASE=False,
                     wind_len=10):
    dt = timeScale[1] - timeScale[0]
    fit_window = int(round(wind_len * 1e-9 / dt / 2))
    N_frames = frames.shape[1]
    N_bunches = np.zeros((N_frames,), dtype=int)
    Bunch_positions = np.zeros((N_frames, N_bunch_max))
    Bunch_lengths = np.zeros((N_frames, N_bunch_max))
    Bunch_peaks = np.zeros((N_frames, N_bunch_max))
    Bunch_intensities = np.zeros((N_frames, N_bunch_max))
    Bunch_positionsFit = np.zeros((N_frames, N_bunch_max))
    Bunch_peaksFit = np.zeros((N_frames, N_bunch_max))
    Bunch_Exponent = np.zeros((N_frames, N_bunch_max))
    Goodness_of_fit = np.zeros((N_frames, N_bunch_max))

    for i in np.arange(N_frames):
        frame = frames[:, i]

        # pos, _ = find_peaks(frame,height=np.max(frames[:,i])*heightFactor,distance=distance)
        pos, _ = find_peaks(frame, height=heightFactor, distance=distance)
        N_bunches[i] = len(pos)
        Bunch_positions[i, 0:N_bunches[i]] = timeScale[pos]
        Bunch_peaks[i, 0:N_bunches[i]] = frame[pos]

        for j, v in enumerate(pos):
            x = 1e9 * timeScale[v - fit_window:v + fit_window]
            y = frame[v - fit_window:v + fit_window]
            if BASE:
                baseline = np.mean(y[:baseline_length])
                y = y - baseline

            try:
                if fit_option == 'fwhm':
                    (mu, sigma, amp) = fwhm(x, y, level=0.5)

                #                    (mu2, sigma2, amp2) = fwhm(x,y,level=0.95)
                else:
                    (amp, mu, sigma, exponent) = binomial_amplitudeN_fit(x, y)
                    y_fit = binomialAmplitudeN(x, *[amp, mu, sigma, exponent])


                    if plot_fit: #or exponent > 5:
                        print(amp, mu, sigma, exponent)

                        plt.plot(x, y, label='measurement')
                        plt.plot(x, y_fit, label='fit')
                        plt.vlines(x[baseline_length], np.min(y), np.max(y), linestyle='--')
                        plt.legend()
                        plt.show()

                    sigma /= 4
            except:
                print(len(x), len(y))
                print(x, y)
                plt.figure()
                plt.plot(x, y)
                plt.show()
                x_71 = x
                y_71 = y


            Bunch_lengths[i, j] = 4 * sigma
            Bunch_intensities[i, j] = np.sum(y)
            Bunch_positionsFit[i, j] = mu
            Bunch_peaksFit[i, j] = amp
            if fit_option != 'fwhm':
                Bunch_Exponent[i, j] = exponent
                Goodness_of_fit[i, j] = np.mean(np.abs(y - y_fit)/np.max(y)) * 100

    N_bunches_max = np.max(N_bunches)
    Bunch_positions = Bunch_positions[:, 0:N_bunches_max]
    Bunch_peaks = Bunch_peaks[:, 0:N_bunches_max]
    Bunch_lengths = Bunch_lengths[:, 0:N_bunches_max]
    Bunch_intensities = Bunch_intensities[:, 0:N_bunches_max]
    Bunch_positionsFit = Bunch_positionsFit[:, 0:N_bunches_max]
    Bunch_peaksFit = Bunch_peaksFit[:, 0:N_bunches_max]
    Bunch_Exponent = Bunch_Exponent[:, 0:N_bunches_max]
    Goodness_of_fit = Goodness_of_fit[:, 0:N_bunches_max]

    return N_bunches, Bunch_positions, Bunch_peaks, Bunch_lengths, Bunch_intensities, Bunch_positionsFit, \
           Bunch_peaksFit, Bunch_Exponent, Goodness_of_fit


def extract_bunch_position(time, profile, heighFactor=0.015, wind_len=10):
    N_bunches, Bunch_positions, Bunch_peaks, Bunch_lengths, Bunch_intensities, Bunch_positionsFit, \
    Bunch_peaksFit, Bunch_Exponent, Goodness_of_fit = getBeamPattern(time, np.array([profile]).T,
                                                                     heightFactor=heighFactor, wind_len=wind_len)
    return Bunch_positionsFit[0, 0], Bunch_lengths[0, 0]

def bunch_position_from_COM(time, profile):
    M = np.trapz(profile, time)
    return np.trapz(profile * time, time) / M


def naive_fit_sine(signal, t):
    signal = np.arcsin(signal)

    slope, intercept, r, p, se = linregress(t, signal)
    return slope

def sine_wave(t, omega, A, phi):
    return A * np.sin(omega * t + phi)

def fit_sine_curve_fit(signal, t):
    popt, pcov = curve_fit(sine_wave, t, signal, p0=[0.008, 0.15, 0])

    return popt


import numpy, scipy.optimize

def fit_sin(tt, yy):
    '''Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'''
    tt = numpy.array(tt)
    yy = numpy.array(yy)
    ff = numpy.fft.fftfreq(len(tt), (tt[1]-tt[0]))   # assume uniform spacing
    Fyy = abs(numpy.fft.fft(yy))
    guess_freq = abs(ff[numpy.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
    guess_amp = numpy.std(yy) * 2.**0.5
    guess_offset = numpy.mean(yy)
    guess = numpy.array([guess_amp, 2.*numpy.pi*guess_freq, 0., guess_offset, 0])

    def sinfunc(t, A, w, p, c, alpha):  return A * numpy.sin(w*t + p) * np.exp(-alpha * t) + c
    try:
        popt, pcov = scipy.optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    except:
        plt.figure()
        plt.plot(tt, yy)
        plt.show()

    A, w, p, c, alpha = popt
    f = w/(2.*numpy.pi)
    fitfunc = lambda t: A * numpy.sin(w*t + p) * np.exp(-alpha * t) + c
    return {"amp": A, "omega": w, "phase": p, "offset": c, "alpha": alpha, "freq": f,
            "period": 1./f, "fitfunc": fitfunc, "maxcov": numpy.max(pcov), "rawres": (guess,popt,pcov)}


def reshape_data(data, t, T):
    N_turns = t[-1] // T
    n_samples = T // (t[1] - t[0])

    rdata = np.zeros((N_turns, n_samples))

    for i in range(N_turns):
        data_i = data[n_samples * i: n_samples * (i + 1)]

        rdata[i, :] = data_i

    return rdata

def analyse_profile(profile, sample_rate):
    t = np.linspace(0, len(profile) / sample_rate, len(profile))


def get_profile_data(f, fdir):
    data = h5py.File(fdir + f, 'r')
    profile = data['Profile']['profile'][:]
    t = np.linspace(0, data['Profile']['profile'][:].shape[0] / data['Profile']['samplerate'][0],
                    data['Profile']['profile'][:].shape[0])

    return profile, t

def find_file_in_folder(f, fdir):
    file_name = None
    for file in os.listdir(fdir):
        if file.startswith(f):
            file_name = file

    return file_name

def get_sorted_files(V, QL, cavity, beam, emittance, add=''):
    V = str(V).replace('.', '')
    return f'profile_V{V}MV_QL{QL}k_C{cavity}B{beam}{add}_{emittance}_emittance.h5'


def find_synchrotron_frequency_from_profile(profile, t, T_rev, turn_constant, init_osc_length, final_osc_start):
    N_bunches, Bunch_positions, Bunch_peaks, Bunch_lengths, Bunch_intensities, Bunch_positionsFit, \
    Bunch_peaksFit, Bunch_Exponent, Goodness_of_fit = getBeamPattern(t, profile, heightFactor=30,
                                                                         wind_len=4)

    bpos = Bunch_positionsFit[:, 0]
    t = np.linspace(0, len(bpos), len(bpos)) * T_rev * turn_constant

    fit_dict_init = fit_sin(t[:init_osc_length], bpos[:init_osc_length])
    fit_dict_final = fit_sin(t[final_osc_start:], bpos[final_osc_start:])

    return fit_dict_init, fit_dict_final

def analyse_synchrotron_frequency_cavity_by_cavity(V, QL, cavitiesB1, cavitiesB2, emittance, fdir, T_rev, turn_constant,
                                                   init_osc_length, final_osc_start, add=''):
    B1_files = []
    B2_files = []
    for i in range(len(cavitiesB1)):
        B1_files.append(get_sorted_files(V=V, QL=QL, cavity=cavitiesB1[i],
                                             beam=1, emittance=emittance, add=add))
        B2_files.append(get_sorted_files(V=V, QL=QL, cavity=cavitiesB2[i],
                                             beam=2, emittance=emittance, add=add))

    freqs_init = np.zeros((len(cavitiesB1), 2))
    freqs_final = np.zeros((len(cavitiesB1), 2))

    for i in range(len(cavitiesB1)):
        # Beam 1
        B1_profiles, t = get_profile_data(B1_files[i], fdir)

        fit_dict_init, fit_dict_final = find_synchrotron_frequency_from_profile(B1_profiles, t, T_rev,
                                                                                turn_constant,
                                                                                init_osc_length, final_osc_start)
        freqs_init[i, 0] = fit_dict_init['freq']
        freqs_final[i, 0] = fit_dict_final['freq']

        # Beam 2
        B2_profiles, t = get_profile_data(B2_files[i], fdir)

        fit_dict_init, fit_dict_final = find_synchrotron_frequency_from_profile(B2_profiles, t, T_rev,
                                                                                turn_constant,
                                                                                init_osc_length, final_osc_start)
        freqs_init[i, 1] = fit_dict_init['freq']
        freqs_final[i, 1] = fit_dict_final['freq']

    return freqs_init, freqs_final


def analyse_synchrotron_frequency_with_all_cavities(V, QL, emittance, fdir, T_rev, turn_constant,
                                                    init_osc_length, final_osc_start, add=''):

    B1_file = get_sorted_files(V=V, QL=QL, cavity='all', beam=1, emittance=emittance, add=add)
    B2_file = get_sorted_files(V=V, QL=QL, cavity='all', beam=2, emittance=emittance, add=add)

    freqs_init = np.zeros(2)
    freqs_final = np.zeros(2)


    # Beam 1
    B1_profiles, t = get_profile_data(B1_file, fdir)

    fit_dict_init, fit_dict_final = find_synchrotron_frequency_from_profile(B1_profiles, t, T_rev,
                                                                            turn_constant,
                                                                            init_osc_length, final_osc_start)
    freqs_init[0] = fit_dict_init['freq']
    freqs_final[0] = fit_dict_final['freq']

    # Beam 2
    B2_profiles, t = get_profile_data(B2_file, fdir)

    fit_dict_init, fit_dict_final = find_synchrotron_frequency_from_profile(B2_profiles, t, T_rev,
                                                                            turn_constant,
                                                                            init_osc_length, final_osc_start)
    freqs_init[1] = fit_dict_init['freq']
    freqs_final[1] = fit_dict_final['freq']

    return freqs_init, freqs_final

def get_first_profiles(fdir, rev_str, profile_length):
    '''
    File to get first profiles from a folder fdir and containing rev_string in the filename.

    :param fdir:
    :param rev_str:
    :param profile_length:
    :return:
    '''
    file_names = []
    for file in os.listdir(fdir):
        if rev_str in file:
            file_names.append(file)

    data = np.zeros((profile_length, len(file_names)))
    ts = np.zeros((profile_length, len(file_names)))

    for i in range(len(file_names)):
        data_i, ti = get_profile_data(file_names[i], fdir)

        data[:, i] = data_i[:, 0]
        ts[:, i] = ti

    return data, ts


def find_bunch_length(fdir, emittance, n_samples=250):
    r'''
    File to get the average and standard deviation of the first turn bunch in the LHC with a small or
    nominal emittance.

    :param fdir:
    :param emittance:
    :param n_samples:
    :return:
    '''
    profiles, ts = get_first_profiles(fdir, emittance, n_samples)

    N_bunches, Bunch_positions, Bunch_peaks, Bunch_lengths, Bunch_intensities, Bunch_positionsFit, \
    Bunch_peaksFit, Bunch_Exponent, Goodness_of_fit = getBeamPattern(ts[:, 0], profiles, heightFactor=30,
                                                                         wind_len=5, fit_option='fwhm')

    return np.mean(Bunch_lengths), np.std(Bunch_lengths)



def plot_profile(Profile, turn, save_to):
    fig, ax = plt.subplots()

    ax.set_title(f'Profile at turn {turn}')
    ax.plot(Profile.bin_centers * 1e9, Profile.n_macroparticles)
    ax.set_xlabel(r'$\Delta t$ [ns]')
    ax.set_ylabel(r'$N_m$ [-]')

    fig.savefig(save_to + f'profile_{turn}.png')


def save_profile(Profile, turn, save_to):
    profile_data = np.zeros((2, len(Profile.n_macroparticles)))
    profile_data[0, :] = Profile.n_macroparticles
    profile_data[1, :] = Profile.bin_centers

    np.save(save_to + f'profile_data_{turn}', profile_data)


def plot_bunch_position(bp, time, j, save_to, COM=False):
    fig, ax = plt.subplots()

    if COM:
        ax.set_title('Bunch Position COM')
    else:
        ax.set_title('Bunch Position')

    ax.plot(time[:j] * 1e-6, bp[:j])
    ax.set_xlabel(r'Time since injection [$\mu$s]')
    ax.set_ylabel(r'Bunch position')

    if COM:
        fig.savefig(save_to + 'bunch_position_com.png')
    else:
        fig.savefig(save_to + 'bunch_position.png')


def plot_bunch_length(bl, time, j, save_to):
    fig, ax = plt.subplots()

    ax.set_title('Bunch Length')

    ax.plot(time[:j] * 1e-6, bl[:j])
    ax.set_xlabel(r'Time since injection [$\mu$s]')
    ax.set_ylabel(r'Bunch length')

    fig.savefig(save_to + 'bunch_length.png')


def save_array(arr, filename, save_to):
    np.save(save_to + filename, arr)


def plot_phase_space(Beam, des, dts):
    data = {r'$\Delta E$': Beam.dE * 1e-6, r'$\Delta t$': Beam.dt * 1e9}
    cp = sns.color_palette('coolwarm', as_cmap=True)
    sns.displot(data, x=r'$\Delta t$', y=r'$\Delta E$', cbar=True, cmap=cp, vmin=0, vmax=150)
    plt.xlabel(r'$\Delta t$ [ns]')
    plt.ylabel(r'$\Delta E$ [MeV]')
    plt.xlim((-1.25, 2.5 + 1.25))
    plt.ylim((-600, 600))

    plt.plot(dts * 1e9, des * 1e-6, color='black')
    plt.plot(dts * 1e9, -des * 1e-6, color='black')