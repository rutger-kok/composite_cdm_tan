# -*- coding: utf-8 -*-
'''
Script to determine the in-situ properties of a composite material.
Based on:

[1] G. Catalanotti, P. P. Camanho, and A. T. Marques
Three-dimensional failure criteria for fiber-reinforced laminates
Compos. Struct., vol. 95, pp. 63-79, 2013
http://doi.org/10.1016/j.compstruct.2012.07.016

[2] P. P. Camanho, C. G. Davila, S. T. Pinho, L. Iannucci, and P. Robinson
Prediction of in situ strengths and matrix cracking in composites under
transverse tension and in-plane shear
Compos. Part A Appl. Sci. Manuf., vol. 37, no. 2, pp. 165-176, 2006.
http://doi.org/10.1016/j.compositesa.2005.04.023

IM7-8552 properties from:
[3] P. P. Camanho, P. Maimi, and C. G. Dávila
Prediction of size effects in notched laminates using continuum damage
mechanics
Compos. Sci. Technol., vol. 67, no. 13, pp. 1715-1727, 2007.
http://doi.org/10.1016/j.compscitech.2007.02.005

(c) Rutger Kok 2021
'''
import numpy as np
from catalanotti import catalanotti
from scipy.optimize import fsolve

# # Material Properties SHD Composites VTC401-UD300-T700-24K
# XT = 2.354  # tensile strength fiber direction
# XC = 1.102  # compressive strength fiber direction
# beta = 2.98e-8  # shear response factor
# g12 = 6470  # shear modulus (in MPa)
# t = 0.207  # cured ply thickness (in mm)
# GIc = 0.30  # mode I fracture toughness (in kJ/m2)
# GIIc = 0.80  # mode II fracture toughness (in kJ/m2)
# E11 = 116600  # longitudinal modulus (MPa)
# E22 = 7700  # transverse modulus (MPa)
# nu_21 = 0.30 * (E22 / E11)  # Poisson's ratio 21-direction
# alpha0 = np.radians(53.0)

# Material Properties IM7-8552
XT = 2.3235  # tensile strength fiber direction
XC = 1.2001  # compressive strength fiber direction
beta = 2.98e-8  # shear response factor
g12 = 5290  # shear modulus (in MPa)
t = 0.131  # cured ply thickness (in mm)
GIc = 0.2774  # mode I fracture toughness (in kJ/m2)
GIIc = 0.7879  # mode II fracture toughness (in kJ/m2)
E11 = 171420  # longitudinal modulus (MPa)
E22 = 9080  # transverse modulus (MPa)
nu_21 = 0.32 * (E22 / E11)  # Poisson's ratio 21-direction
alpha0 = np.radians(53.0)

# Calculate in-situ properties
# in-situ longitudinal shear strength
phi = (48.0 * GIIc) / (np.pi * t)  # assume thin embedded ply (Eq. 49 [2])
SL = np.sqrt((np.sqrt(1.0 + beta * phi * g12**2.0) - 1.0) /
             (3.0 * beta * g12)) / 1000.0  # Eq. 48 [2] (in GPa)
# in-situ transverse tensile strength
A022 = 2.0 * ((1.0 / E22) - (nu_21**2.0 / E11))  # Eq. 27 [2]
YT = np.sqrt((8.0 * GIc) / (np.pi * t * A022)) / 1000.0  # Eq. 43 [2] (in GPa)


# calculate etaL iteratively to satisfy the condition
# fi_mc == 1.0 @ trial_stress[1] = -YC
def func(etaLTest):
    # in-situ transverse shear strength (Eq. 12 [1])
    ST = (0.5 * (((2.0 * np.sin(alpha0)**2.0) - 1.0) * SL) /
          (((1.0 - np.sin(alpha0)**2.0)**0.5) * np.sin(alpha0) * etaLTest))
    # in-situ transverse compressive strength (Eq. 11 [1])
    YC = (-(SL * (2.0 * np.cos(alpha0)**2.0 - 1.0)) /
          (etaLTest * np.cos(alpha0)**2.0))
    etaT = (etaLTest * ST) / SL  # Eq.10 [1]
    kappa = (ST**2.0 - YT**2.0) / (ST * YT)  # Eq.43 [1]
    lmbda = ((2.0 * etaLTest * ST) / SL) - kappa  # Eq.45 [1]

    trial_stress = np.zeros((6, 1))
    trial_stress[1] = -YC  # set sigma_22 to equal

    # calculate failure indices
    fi_mt, fi_mc = catalanotti(trial_stress, ST, SL, etaLTest, etaT, lmbda,
                               kappa)
    # subtract 1.0 from fi_mc for root finding algorithm
    result = fi_mc - 1.0
    return result


# use scipy root finding algorithm to determine etaL (initial guess = 0.5)
etaL = fsolve(func, [0.5])
# calculate final value of in-situ transverse shear strength
ST = (0.5 * (((2.0 * np.sin(alpha0)**2.0) - 1.0) * SL) /
      (((1.0 - np.sin(alpha0)**2.0)**0.5) * np.sin(alpha0) * etaL))
# calculate final value of in-situ transverse compressive strength
YC = (-(SL * (2.0 * np.cos(alpha0)**2.0 - 1.0)) /
      (etaL * np.cos(alpha0)**2.0))

# print results
print 'etaL = ', etaL
print 'ST = ', ST
print 'SL = ', SL
print 'YT = ', YT
print 'YC = ', YC
