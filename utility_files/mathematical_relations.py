'''
File with mathematical relations from beam dynamics.

Author: Birk Emil Karlsen-Bæck
'''

# Import
import numpy as np
import scipy.constants as spc

# Physical Parameters for the MD
q = 1                           # [e]
phi_s = 0                       # [rad]
p_s = 450e9                     # [eV]
h = 35640                       # [-]
gamma_t1 = 53.606713            # [-]
gamma_t2 = 53.674152            # [-]
T_rev = 8.892465516509656e-05   # [s]
m_p = spc.physical_constants['proton mass energy equivalent in MeV'][0] * 1e6
E_s = np.sqrt(p_s**2 + m_p**2)
gamma_s = E_s / m_p
eta1 = (1/gamma_t1**2) - (1/gamma_s**2)
eta2 = (1/gamma_t2**2) - (1/gamma_s**2)
beta = np.sqrt(1 - 1/gamma_s**2)

# Functions

def synchrotron_frequency(V, h=h, eta=eta1, q=q, phi_s=phi_s, beta=beta, T_rev=T_rev, E_s=E_s):
    r'''Relation for synchrotron frequency when assuming small amplitude oscillations in
    longitudinal phase-space.
    :param h:
    :param eta:
    :param q:
    :param V:
    :param phi_s:
    :param beta:
    :param T_rev:
    :param E_s:
    :return:
    '''
    return np.sqrt((2 * np.pi * h * eta * q * V * np.cos(phi_s))/(beta**2 * T_rev**2 * E_s)) / (2 * np.pi)


def RF_voltage_from_synchrotron_frequency(f_s, E_s=E_s, beta=beta, T_rev=T_rev, h=h, eta=eta1, q=q, phi_s=phi_s):
    r'''
    Relation between the RF voltage and synchrotron frequency when assuming small amplitude oscillations in
    longitudinal phase-space.
    :param E_s:
    :param omega_s:
    :param beta:
    :param T_rev:
    :param h:
    :param eta:
    :param q:
    :param phi_s:
    :return:
    '''
    omega_s = 2 * np.pi * f_s
    return (E_s * omega_s**2 * beta**2 * T_rev**2)/(2 * np.pi * h * eta * q * np.cos(phi_s))