#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 23:23:18 2023

@author: garvinyim
"""

import numpy
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font',**{'family':'serif','serif':['Times']})

# Global variables
solar_mass = 1.9885e30 # kg
kiloparsec = 3.086e19 # m
G = 6.67428e-11 # kg^-1 m^3 s^-2
c = 2.99792458e8 # m s^-1
pi = numpy.pi

# Definitions (mostly from Appendix B of Jaranowski, Krolak & Schutz, 1998)

def j_1(lambda_det, gamma_det):
    return (1.0/256.0) * (4.0 - 20.0 * pow(numpy.cos(lambda_det), 2) + 35.0 * pow(numpy.sin(2.0*gamma_det), 2) * pow(numpy.cos(lambda_det), 4))

def j_2(lambda_det, gamma_det):
    return (1.0/1024.0) * (68.0 - 20.0 * pow(numpy.cos(lambda_det), 2) - 13.0 * pow(numpy.sin(2.0*gamma_det), 2) * pow(numpy.cos(lambda_det), 4))

def j_3(lambda_det, gamma_det):
    return (1.0/128.0) * (28.0 - 44.0 * pow(numpy.cos(lambda_det), 2) + 5.0 * pow(numpy.sin(2.0*gamma_det), 2) * pow(numpy.cos(lambda_det), 4))

def j_4(lambda_det, gamma_det):
    return (1.0/32.0) * (2.0 - 7.0 * pow(numpy.sin(2.0*gamma_det), 2) * pow(numpy.cos(lambda_det), 2)) * numpy.sin(2.0*lambda_det)

def j_5(lambda_det, gamma_det):
    return (1.0/32.0) * (3.0 - 7.0 * numpy.cos(4.0*gamma_det) - 7.0 * pow(numpy.sin(2.0*gamma_det), 2) * pow(numpy.cos(lambda_det), 2)) * pow(numpy.cos(lambda_det), 2)
                         
def j_6(lambda_det, gamma_det):
    return (1.0/96.0) * (2.0 * numpy.cos(4.0*gamma_det) + pow(numpy.sin(2.0*gamma_det), 2) * pow(numpy.cos(lambda_det), 2)) * numpy.sin(2.0*lambda_det)
                         
def j_7(lambda_det, gamma_det):
    return (1.0/1024.0) * (4.0 * numpy.cos(4.0*gamma_det) * pow(numpy.sin(lambda_det), 2) - pow(numpy.sin(2.0*gamma_det), 2) * pow(numpy.cos(lambda_det), 4))

def j_8(lambda_det, gamma_det):
    return (1.0/32.0) * numpy.sin(4.0*gamma_det) * pow(numpy.cos(lambda_det), 3)

def j_9(lambda_det, gamma_det):
    return (1.0/32.0) * numpy.sin(4.0*gamma_det) * pow(numpy.cos(lambda_det), 2) * numpy.sin(lambda_det)

def j_10(lambda_det, gamma_det):
    return (1.0/196.0) * numpy.sin(4.0*gamma_det) * (5.0 - 3.0 * numpy.cos(2.0*lambda_det)) * numpy.cos(lambda_det)

def j_11(lambda_det, gamma_det):
    return (1.0/1024.0) * numpy.sin(4.0*gamma_det) * (3.0 - numpy.cos(2.0*lambda_det)) * numpy.sin(lambda_det)

def j_12(lambda_det, gamma_det):
    return (1.0/32.0) * (14.0 - pow(numpy.sin(2.0*gamma_det), 2) * pow(numpy.cos(lambda_det), 2)) * numpy.sin(2.0*lambda_det)

def j_13(lambda_det, gamma_det):
    return (1.0/32.0) * (9.0 - 5.0 * numpy.cos(4.0*gamma_det) - 5.0 * pow(numpy.sin(2.0*gamma_det), 2) * pow(numpy.cos(lambda_det), 2)) * pow(numpy.cos(lambda_det), 2)



def e_1(delta, lambda_det, gamma_det):
    return 4.0 * j_1(lambda_det, gamma_det) * pow(numpy.cos(delta), 4)

def e_2(delta, lambda_det, gamma_det):
    return 4.0 * j_2(lambda_det, gamma_det) - j_3(lambda_det, gamma_det) * numpy.cos(2.0*delta) + j_1(lambda_det, gamma_det) * pow(numpy.cos(2.0*delta), 2)

def f_11(delta, lambda_det, gamma_det):
    return - 4.0 * j_4(lambda_det, gamma_det) * pow(numpy.cos(delta), 3) * numpy.sin(delta)

def f_12(delta, lambda_det, gamma_det):
    return j_5(lambda_det, gamma_det) * pow(numpy.cos(delta), 2) * (3.0 - numpy.cos(2.0*delta))

def f_13(delta, lambda_det, gamma_det):
    return - j_6(lambda_det, gamma_det) * (7.0 - numpy.cos(2.0*delta)) * numpy.sin(2.0*delta)

def f_14(delta, lambda_det, gamma_det):
    return - j_7(lambda_det, gamma_det) * (35.0 - 28.0 * numpy.cos(2.0*delta) + numpy.cos(4.0*delta))

def f_21(delta, lambda_det, gamma_det):
    return - 28.0 * j_8(lambda_det, gamma_det) * pow(numpy.cos(delta), 3) * numpy.sin(delta)

def f_22(delta, lambda_det, gamma_det):
    return - 7.0 * j_9(lambda_det, gamma_det) * (3.0 - numpy.cos(2.0*delta)) * pow(numpy.cos(delta), 2) 

def f_23(delta, lambda_det, gamma_det):
    return - j_10(lambda_det, gamma_det) * (7.0 - numpy.cos(2.0*delta)) * numpy.sin(2.0*delta)

def f_24(delta, lambda_det, gamma_det):
    return - j_11(lambda_det, gamma_det) * (35.0 - 28.0 * numpy.cos(2.0*delta) + numpy.cos(4.0*delta))

def g_11(delta, lambda_det, gamma_det):
    return 28.0 * j_8(lambda_det, gamma_det) * pow(numpy.cos(delta), 3)

def g_12(delta, lambda_det, gamma_det):
    return 28.0 * j_9(lambda_det, gamma_det) * pow(numpy.cos(delta), 2) * numpy.sin(delta)

def g_13(delta, lambda_det, gamma_det):
    return 2.0 * j_10(lambda_det, gamma_det) * (5.0 - 3.0 * numpy.cos(2.0*delta)) * numpy.cos(delta) 

def g_14(delta, lambda_det, gamma_det):
    return 16.0 * j_11(lambda_det, gamma_det) * (3.0 - numpy.cos(2.0*delta)) * numpy.sin(delta) 

def g_21(delta, lambda_det, gamma_det):
    return - 4.0 * j_4(lambda_det, gamma_det) * pow(numpy.cos(delta), 3)

def g_22(delta, lambda_det, gamma_det):
    return 4.0 * j_5(lambda_det, gamma_det) * pow(numpy.cos(delta), 2) * numpy.sin(delta)

def g_23(delta, lambda_det, gamma_det):
    return - 2.0 * j_6(lambda_det, gamma_det) * (5.0 - 3.0 * numpy.cos(2.0*delta)) * numpy.cos(delta) 

def g_24(delta, lambda_det, gamma_det):
    return - 16.0 * j_7(lambda_det, gamma_det) * (3.0 - numpy.cos(2.0*delta)) * numpy.sin(delta) 

def h_11(delta, lambda_det, gamma_det):
    return (j_12(lambda_det, gamma_det) - j_4(lambda_det, gamma_det) * numpy.cos(2.0*delta)) * numpy.sin(2.0*delta) 

def h_12(delta, lambda_det, gamma_det):
    return (j_13(lambda_det, gamma_det) - j_5(lambda_det, gamma_det) * numpy.cos(2.0*delta)) * pow(numpy.cos(delta), 2)

def h_13(delta, lambda_det, gamma_det):
    return 4.0 * j_6(lambda_det, gamma_det) * pow(numpy.cos(delta), 3) * numpy.sin(delta)

def h_14(delta, lambda_det, gamma_det):
    return - 8.0 * j_7(lambda_det, gamma_det) * pow(numpy.cos(delta), 4) 

def h_21(delta, lambda_det, gamma_det):
    return j_8(lambda_det, gamma_det) * (1.0 - 7.0 * numpy.cos(2.0*delta)) * numpy.sin(2.0*delta)

def h_22(delta, lambda_det, gamma_det):
    return - j_9(lambda_det, gamma_det) * (5.0 - 7.0 * numpy.cos(2.0*delta)) * pow(numpy.cos(delta), 2) 

def h_23(delta, lambda_det, gamma_det):
    return 4.0 * j_10(lambda_det, gamma_det) * pow(numpy.cos(delta), 3) * numpy.sin(delta)

def h_24(delta, lambda_det, gamma_det):
    return - 8.0 * j_11(lambda_det, gamma_det) * pow(numpy.cos(delta), 4) 

for k in [1, 2]:
    for n in [1, 2, 3, 4]:
        def f(k, n, delta, lambda_det, gamma_det):
            subscript = "{}{}".format(k, n)
            if subscript == "11":
                return f_11(delta, lambda_det, gamma_det)
            elif subscript == "12":
                return f_12(delta, lambda_det, gamma_det)
            elif subscript == "13":
                return f_13(delta, lambda_det, gamma_det)
            elif subscript == "14":
                return f_14(delta, lambda_det, gamma_det)
            elif subscript == "21":
                return f_21(delta, lambda_det, gamma_det)
            elif subscript == "22":
                return f_22(delta, lambda_det, gamma_det)
            elif subscript == "23":
                return f_23(delta, lambda_det, gamma_det)
            elif subscript == "24":
                return f_24(delta, lambda_det, gamma_det)
            
for k in [1, 2]:
    for n in [1, 2, 3, 4]:
        def g(k, n, delta, lambda_det, gamma_det):
            subscript = "{}{}".format(k, n)
            if subscript == "11":
                return g_11(delta, lambda_det, gamma_det)
            elif subscript == "12":
                return g_12(delta, lambda_det, gamma_det)
            elif subscript == "13":
                return g_13(delta, lambda_det, gamma_det)
            elif subscript == "14":
                return g_14(delta, lambda_det, gamma_det)
            elif subscript == "21":
                return g_21(delta, lambda_det, gamma_det)
            elif subscript == "22":
                return g_22(delta, lambda_det, gamma_det)
            elif subscript == "23":
                return g_23(delta, lambda_det, gamma_det)
            elif subscript == "24":
                return g_24(delta, lambda_det, gamma_det)                

for k in [1, 2]:
    for n in [1, 2, 3, 4]:
        def h(k, n, delta, lambda_det, gamma_det):
            subscript = "{}{}".format(k, n)
            if subscript == "11":
                return h_11(delta, lambda_det, gamma_det)
            elif subscript == "12":
                return h_12(delta, lambda_det, gamma_det)
            elif subscript == "13":
                return h_13(delta, lambda_det, gamma_det)
            elif subscript == "14":
                return h_14(delta, lambda_det, gamma_det)
            elif subscript == "21":
                return h_21(delta, lambda_det, gamma_det)
            elif subscript == "22":
                return h_22(delta, lambda_det, gamma_det)
            elif subscript == "23":
                return h_23(delta, lambda_det, gamma_det)
            elif subscript == "24":
                return h_24(delta, lambda_det, gamma_det)  



def F_1(iota):
    return - (1.0/16.0) * pow(numpy.sin(iota), 4)

def F_2(iota):
    return (1.0/4.0) * pow(numpy.sin(iota), 4)

def G_1(iota):
    return (1.0/16.0) * pow(numpy.sin(iota), 2) (1.0 + pow(numpy.cos(iota), 2))

def G_2(iota):
    return (1.0/4.0) * (1.0 + 6.0 * pow(numpy.cos(iota), 2) + pow(numpy.cos(iota), 4))

for k in [1, 2]:
    def F(k, iota):
        if k == 1:
            return F_1(iota)
        elif k == 2:
            return F_2(iota)
        
for k in [1, 2]:
    def G(k, iota):
        if k == 1:
            return G_1(iota)
        elif k == 2:
            return G_2(iota)

def C(k, n, delta, psi, iota, lambda_det, gamma_det):
    return F(k, iota) * (f(1, n, delta, lambda_det, gamma_det) * numpy.cos(4.0*psi) + g(1, n, delta, lambda_det, gamma_det) * numpy.sin(4.0*psi)) + G(k, iota) * h(1, n, delta, lambda_det, gamma_det)
    
def D(k, n, delta, psi, iota, lambda_det, gamma_det):
    return F(k, iota) * (f(2, n, delta, lambda_det, gamma_det) * numpy.cos(4.0*psi) + g(2, n, delta, lambda_det, gamma_det) * numpy.sin(4.0*psi)) + G(k, iota) * h(2, n, delta, lambda_det, gamma_det)

def A(k, delta, psi, iota, lambda_det, gamma_det):
    return F(k, iota) * e_1(delta, lambda_det, gamma_det) * numpy.cos(4.0*psi) + G(k, iota) * e_2(delta, lambda_det, gamma_det)

def B(k, alpha, delta, T_obs, psi, iota, lambda_det, gamma_det, Omega_r, phi_r):
    summation = 0
    for n in [1, 2, 3, 4]:
        contribution = numpy.sin(n * Omega_r * T_obs / 2.0) * (C(k, n, delta, psi, iota, lambda_det, gamma_det) * numpy.cos(n*(alpha - phi_r)) + D(k, n, delta, psi, iota, lambda_det, gamma_det) * numpy.sin(n*(alpha - phi_r)))
        summation = summation + contribution
    return summation / Omega_r