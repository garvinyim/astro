        #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 23:14:50 2023

@author: garvinyim
"""

# =============================================================================
# Preamble
# =============================================================================

import numpy
import math
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font',**{'family':'serif','serif':['Times']})
import matplotlib.ticker as ticker
import csv

import JKS

# Global variables
solar_mass = 1.9885e30 # kg
kiloparsec = 3.086e19 # m
G = 6.67428e-11 # kg^-1 m^3 s^-2
c = 2.99792458e8 # m s^-1
pi = numpy.pi

# =============================================================================
# Functions
# =============================================================================


def plot_sqrt_AT_over_sqrt_ATplusB(T_obs_list, detectors, alpha1, delta1, psi1, iota1, Omega_r, phi_r):
    
    hanford_sqrt_AT_over_sqrt_ATplusB = []
    livingston_sqrt_AT_over_sqrt_ATplusB = []
    virgo_sqrt_AT_over_sqrt_ATplusB = []
    
    # Calculation of A_k and B_k
    for detector in detectors:
        lambda_det = detector[0]
        gamma_det = detector[1]
        # zeta_det = detector[2]
        
        detector_sqrt_AT_over_sqrt_ATplusB = []
        
        for duration in T_obs_list:
            A = JKS.A(2, delta1, psi1, iota1, lambda_det, gamma_det)
            B = JKS.B(2, alpha1, delta1, duration, psi1, iota1, lambda_det, gamma_det, Omega_r, phi_r)
            numerator = numpy.sqrt(A * duration)
            denominator = numpy.sqrt(A * duration + B)
            # print(numerator/denominator - 1.0)
            detector_sqrt_AT_over_sqrt_ATplusB.append(numerator/denominator - 1.0)
            
        if detector == hanford:
            hanford_sqrt_AT_over_sqrt_ATplusB  = detector_sqrt_AT_over_sqrt_ATplusB
        elif detector == livingston:
            livingston_sqrt_AT_over_sqrt_ATplusB  = detector_sqrt_AT_over_sqrt_ATplusB
        elif detector == virgo:
            virgo_sqrt_AT_over_sqrt_ATplusB  = detector_sqrt_AT_over_sqrt_ATplusB
     
    # Plotting figure
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_xlabel(r'$T_\mathrm{GW}$ [d]', fontsize = 36, labelpad = 10)
    ax.set_ylabel(r'$\frac{\rho_\mathrm{approx} - \rho_\mathrm{exact}}{\rho_\mathrm{exact}}$', fontsize = 36, labelpad = 6)
    ax.set_xscale('linear')
    ax.set_yscale('linear')
    ax.set_xlim([-1.0,30.99])
    ax.set_ylim([-0.32,0.219])
    # ax.set_ylim([-0.24,0.235])
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(axis='both', which='major', labelsize=26, width = 3, length = 12, direction = 'in', top = True, right = True, pad = 12)
    ax.tick_params(axis='both', which='minor', labelsize=26, width = 1.5, length = 8, direction = 'in', top = True, right = True)
    
    ax.plot(T_obs_list, hanford_sqrt_AT_over_sqrt_ATplusB, label="Hanford", linestyle=":", color = "tab:blue")
    ax.plot(T_obs_list, livingston_sqrt_AT_over_sqrt_ATplusB, label="Livingston", linestyle=":", color = "tab:red")
    ax.plot(T_obs_list, virgo_sqrt_AT_over_sqrt_ATplusB, label="Virgo", linestyle=":", color = "dimgrey")
    ax.text(29.0, 0.16, r'$\alpha = {:.1f}$\textdegree, $\delta = {:.1f}$\textdegree'.format(alpha1*(180.0/pi), delta1*(180.0/pi)), size=22, horizontalalignment='right')
    ax.text(26.0, 0.12, r'$\iota = {:.1f}$\textdegree'.format(iota1*(180.0/pi)), size=22, horizontalalignment='right')
    # ax.text(29.0, 0.18, r'$\alpha = {:.1f}$\textdegree, $\delta = {:.1f}$\textdegree'.format(alpha1*(180.0/pi), delta1*(180.0/pi)), size=22, horizontalalignment='right')
    # ax.text(26.0, 0.14, r'$\iota = {:.1f}$\textdegree'.format(iota1*(180.0/pi)), size=22, horizontalalignment='right')
    
    ax.legend(fontsize=20, loc='lower right', bbox_to_anchor=(0.94, 0.08))
    fig.tight_layout()
    # plt.savefig('../figures/new/alpha_{:.0f}_delta_{:.0f}.eps'.format(alpha1*(180.0/pi), delta1*(180.0/pi)))
    # plt.savefig('../figures/new/alpha_{:.0f}_delta_{:.0f}.pdf'.format(alpha1*(180.0/pi), delta1*(180.0/pi)))
    plt.show() 
    return None


# =============================================================================
# Main code
# =============================================================================

Omega_r = 2.0*pi # day^-1
phi_r = 0.0

hanford = [46.45*(pi/180.0), 171.8*(pi/180.0)] 
livingston = [30.56*(pi/180.0), 243.0*(pi/180.0)] 
virgo = [43.63*(pi/180.0), 116.5*(pi/180.0)]

detectors = [hanford, livingston, virgo]
T_GW = numpy.linspace(0.0, 30, 3001)
# T_GW = [i * 60*60*24 for i in T_GW_days]

alpha, delta = 40.0*(pi/180.0), 30.0*(pi/180.0)
# alpha, delta = 83.6*(pi/180.0), 22.0*(pi/180.0) # Crab
# alpha, delta = 128.8*(pi/180.0), -45.2*(pi/180.0) # Vela
# alpha, delta = 84.4*(pi/180.0), -69.2*(pi/180.0) # Crab


plot_sqrt_AT_over_sqrt_ATplusB(T_GW, detectors, alpha, delta, 0.0, 0.0, Omega_r, phi_r)





