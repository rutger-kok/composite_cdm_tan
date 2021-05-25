# -*- coding: utf-8 -*-
'''
This script determine the in-situ properties of a composite material.
Copyright (C) 2021 Rutger Wouter Kok

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301 USA

In-situ material property calculations are based on [1] and [2]. The etaL
parameter is normally experimentally determined but can be estimated using this
script by iteratively determining the value of etaL that satisfies the
condition: FI_M = 1.0 when sigma_22 = -YC.

[1] G. Catalanotti, P. P. Camanho, and A. T. Marques
Three-dimensional failure criteria for fiber-reinforced laminates
Compos. Struct., vol. 95, pp. 63-79, 2013
http://doi.org/10.1016/j.compstruct.2012.07.016

[2] P. P. Camanho, C. G. Davila, S. T. Pinho, L. Iannucci, and P. Robinson
Prediction of in situ strengths and matrix cracking in composites under
transverse tension and in-plane shear
Compos. Part A Appl. Sci. Manuf., vol. 37, no. 2, pp. 165-176, 2006.
http://doi.org/10.1016/j.compositesa.2005.04.023

IM7-8552 properties are taken from:
[3] P. P. Camanho, P. Maimi, and C. G. Dávila
Prediction of size effects in notched laminates using continuum damage
mechanics
Compos. Sci. Technol., vol. 67, no. 13, pp. 1715-1727, 2007.
http://doi.org/10.1016/j.compscitech.2007.02.005

'''
import numpy as np
from catalanotti import fail_initiation
from scipy.optimize import fsolve

# Material Properties IM7-8552
# XT = 2.3235  # tensile strength fiber direction
# XC = 1.2001  # compressive strength fiber direction
# beta = 2.98e-8  # shear response factor
# g12 = 5290  # shear modulus (in MPa)
# t = 0.131  # cured ply thickness (in mm)
# GIc = 0.2774  # mode I fracture toughness (in kJ/m2)
# GIIc = 0.7879  # mode II fracture toughness (in kJ/m2)
# E11 = 171420  # longitudinal modulus (MPa)
# E22 = 9080  # transverse modulus (MPa)
# nu_21 = 0.32 * (E22 / E11)  # Poisson's ratio 21-direction
# alpha0 = np.radians(53.0)

# Material Properties SHD Composites VTC401
XT = 2.180  # tensile strength fiber direction
XC = 0.812  # compressive strength fiber direction
beta = 4.72e-8  # shear response factor
g12 = 3267  # shear modulus (in MPa)
t = 0.205  # cured ply thickness (in mm)
GIc = 0.380  # mode I fracture toughness (in kJ/m2)
GIIc = 1.62  # mode II fracture toughness (in kJ/m2)
E11 = 116600  # longitudinal modulus (MPa)
E22 = 7231  # transverse modulus (MPa)
nu_21 = 0.339 * (E22 / E11)  # Poisson's ratio 21-direction
alpha0 = np.radians(53.0)

# Calculate in-situ properties
# in-situ longitudinal shear strength
phi = (48.0 * GIIc) / (np.pi * t)  # assume thin embedded ply Eq. 49 [2]
SL = np.sqrt((np.sqrt(1.0 + beta * phi * g12**2.0) - 1.0) /
             (3.0 * beta * g12)) / 1000.0  # Eq. 48 [2] (in GPa)
# in-situ transverse tensile strength
A022 = 2.0 * ((1.0 / E22) - (nu_21**2.0 / E11))  # Eq. 27 [2]
YT = np.sqrt((8.0 * GIc) / (np.pi * t * A022)) / 1000.0  # Eq. 43 [2] (in GPa)


# calculate etaL iteratively to satisfy the condition
# fi_m == 1.0 @ trial_stress[1] = -YC
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

    trial_stress = np.zeros((6, 1))  # initialize trial_stress array
    trial_stress[1] = -YC  # set sigma_22 = transverse compressive strength

    # calculate failure index
    fi_m, a_fail = fail_initiation(trial_stress, ST, SL, etaLTest, etaT, lmbda,
                                   kappa)
    # subtract 1.0 from fi_m for root finding algorithm
    result = fi_m - 1.0
    return result


# use scipy root finding algorithm to determine etaL (initial guess = 0.5)
etaL = fsolve(func, [0.5])
# calculate final value of in-situ transverse shear strength
ST = (0.5 * (((2.0 * np.sin(alpha0)**2.0) - 1.0) * SL) /
      (((1.0 - np.sin(alpha0)**2.0)**0.5) * np.sin(alpha0) * etaL))
# calculate final value of in-situ transverse compressive strength
YC = (-(SL * (2.0 * np.cos(alpha0)**2.0 - 1.0)) /
      (etaL * np.cos(alpha0)**2.0))
etaT = (etaL * ST) / SL

# print results
print 'etaL = ', etaL
print 'etaT = ', etaT
print 'ST = ', ST
print 'SL = ', SL
print 'YT = ', YT
print 'YC = ', YC
