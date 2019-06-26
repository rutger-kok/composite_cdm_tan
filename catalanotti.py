import matplotlib.pyplot as plt
import numpy as np 
from math import sin, cos, atan, pi, radians, tan

# Python implementation of the Catalanotti failure criteria
# Used to plot failure envelopes

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Material properties

#IM78552 taken from Catalanotti (2013)
# 'is' denotes in-situ properties
Xt = 2323.5
Xc = 1200.1
Yt_is = 160.2
Yc_is = 198.0
Sl_is = 130.2
etaL = 0.5
G12 = 5290.0
alpha0 = radians(53.0)
St_is = (0.5*(((2*sin(alpha0)**2.0)-1.0)*Sl_is)/
        (((1-sin(alpha0)**2.0)**0.5)*sin(alpha0)*etaL))
phic = (atan((1.0-(1.0-4.0*(Sl_is/Xc)*((Sl_is/Xc)+etaL))**0.5)/
        (2.0*((Sl_is/Xc)+etaL))))
etaT = (etaL*St_is)/Sl_is  # Eq 10 catalanotti

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

# Function to determine  f ailure in transverse tension and compression
# and (with rotated coord system) in longitudinal compression 
def fail(tStress,St_is,Yt_is,Sl_is,etaL,etaT):
    
    angleList = [stress1*(pi/200) for stress1 in range(0,101)] # array of angles 0 to pi/2

    phiMC = 0.0 # initialize failure criteria
    phiMT = 0.0

    for angle in angleList:
        # tractions: eq 3. see also eqs 59-61
        tN = (tStress[1]*cos(angle)**2.0 + 2.0*tStress[4]
            *sin(angle)*cos(angle) + tStress[2]*sin(angle)**2.0)
        tT = (-1.0*tStress[1]*cos(angle)*sin(angle)+ 
            tStress[2]*sin(angle)*cos(angle) - tStress[4]*
            (cos(angle)**2.0 - sin(angle)**2.0)) 
        tL = tStress[3]*cos(angle) + tStress[5]*sin(angle)  
  
        kappa = (St_is**2.0-Yt_is**2.0)/(St_is*Yt_is)  # eq 43
        lmbda = ((2.0*etaL*St_is)/Sl_is)-kappa  # eq 45
  
        trialPhiMC = (tL/(Sl_is-etaL*tN))**2.0 + (tT/(St_is-etaT*tN))**2.0 # eq. 5
        trialPhiMT = ((tN/St_is)**2.0 + (tL/Sl_is)**2.0 + (tT/St_is)**2.0 +
            lmbda*(tN/St_is)*(tL/Sl_is)**2.0 + kappa*(tN/St_is)) # eq. 42

        if trialPhiMC > phiMC: # update if criteria at angle is max 
            phiMC = trialPhiMC
            aFailC = angle # record failure plane at max fail criteria
        if trialPhiMT > phiMT:
            phiMT = trialPhiMT
            aFailT = angle
    return phiMC, phiMT

# ----------------------------------------------------------------------------
# Function to determine the misalignment angle phi using the Newton-
# Raphson method. 
# NOTE: not converging well, using a different method for now see line 160
# def newtonRaphson(initial, tStress, X):
#     tol = 1.0*10**-5.0
#     x0 = initial
#     delta = abs(0.0-f(x0, tStress, X))
#     n = 0
#     while delta > tol:
#         x0 = x0 - f(x0, tStress, X)/df(x0, tStress, X)
#         delta = abs(0.0-f(x0, tStress, X))
#         n += 1
#     return x0

# # Definition of function used in the Newton Raphson iteration method
# def f(gamma, tStress, X):
#     val = (X*gamma + 0.5*(tStress[0] - tStress[1]*
#         sin(2.0*gamma)) - abs(tStress[3])*cos(2.0*gamma)) # eq. 88
#     return val     
      
# # Definiton of the derivative function used in NR method
# def df(gamma, tStress, X):
#     val = (X + (tStress[0] - tStress[1])*cos(2.0*gamma)+
#         2.0*abs(tStress[3])*sin(2.0*gamma)) # eq. 89
#     return val

# ----------------------------------------------------------------------------
# Function to rotate stresses into misalignment frame for compressive failure
def rotateStresses(tStress, tht, p):
    tStressTheta = np.zeros(6)
    tStressPhi = np.zeros(6)

    # Rotate stresses by angle tht
    tStressTheta[0] = tStress[0]
    tStressTheta[3] = tStress[3]*cos(tht)+tStress[5]*sin(tht)
    tStressTheta[5] = tStress[5]*cos(tht)-tStress[3]*sin(tht)
    tStressTheta[1] = (tStress[1]*cos(tht)**2.0 + 
                        2.0*tStress[4]*cos(tht)*sin(tht) + 
                        tStress[2]*sin(tht)**2)
    tStressTheta[4] = (tStress[4]*(cos(tht)**2 - sin(tht)**2) -
                        tStress[1]*sin(tht)*cos(tht) +
                        tStress[2]*sin(tht)*cos(tht))
    tStressTheta[2] = (tStress[2]*cos(tht)**2.0 - 
                        2.0*tStress[4]*cos(tht)*sin(tht) +
                        tStress[1]*sin(tht)**2.0)

    # Rotate stresses by angle p
    tStressPhi[0] = (tStressTheta[0]*cos(p)**2 +
                        2.0*tStressTheta[3]*cos(p)*sin(p) +
                        tStressTheta[1]*sin(p)**2)
    tStressPhi[3] = (tStressTheta[3]*(cos(p)**2 - 
                        sin(p)**2) + tStressTheta[1]*sin(p)*cos(p) -
                        tStressTheta[0]*sin(p)*cos(p))
    tStressPhi[5] = tStressTheta[5]*cos(p) + tStressTheta[4]*sin(p)
    tStressPhi[1] = (tStressTheta[1]*cos(p)**2 -
                        2.0*tStressTheta[3]*sin(p)*cos(p) +
                        tStressTheta[0]*sin(p)**2)
    tStressPhi[4] = tStressTheta[4]*cos(p) - tStressTheta[5]*sin(p)
    tStressPhi[2] = tStressTheta[2]
    return tStressPhi 

# ----------------------------------------------------------------------------
# Create failure envelopes

# Xt = 2323.5, Xc = 1200.1, Yt_is = 160.2, Yc_is = 198.0, Sl_is = 130.2

# trialStress specified by Abaqus in the format: (11, 22, 33, 12, 23, 31)

# to plot envelope choose two stresses. stress1 is plotted on x-axis, stress2
# is plotted on the y-axis. e.g. currently sigma_11 plotted against tau_12
envelope = []
n = 250  # number of points along axes
for stress1 in np.linspace(-200, 180, n):
    for stress2 in np.linspace(0.0, 180, n):

        # define loading
        trialStress = np.array([0.0, stress1, 0.0, stress2, 0.0, 0.0])

        # --------------------------------------------------------------------
        # determine fracture plane angle theta (for fiber kinking)   
        if trialStress[3] == 0 and trialStress[5] == 0:
            # eq 55
            theta = 0.5*atan(
                (2.0*trialStress[4])/(trialStress[1]-trialStress[2]))
        else: 
            # eq 56
            theta = atan(trialStress[5]/trialStress[3])
        if np.isnan(theta):
            # if theta is NaN (because of division by zero, set theta to 0)
            theta = 0.0

        # --------------------------------------------------------------------
        # determine the misalignment angle phi
        # X = (sin(2.0*phic)*Xc)/(2.0*phic) # eq. 86 - unused, only needed for 
        # newton-raphson method
        gammaMC = (phic*Xc)/G12  # eq 74
        phi0 = phic - gammaMC  # eq 75
        gammaM = ((phi0*G12 + abs(trialStress[3]))/
                (G12+trialStress[0]-trialStress[1]) - phi0)  # eq 81
        # gammaM = newtonRaphson(3.07, trialStress, X) # NR method - unused
        if trialStress[3] >= 0:  # eq 77
            phi = phi0 + gammaM
        else:
            phi = -1.0*(phi0+gammaM)

        # --------------------------------------------------------------------
        # determine failure indices
        phiLT = 0
        phiLC = 0
        phiMC2 = 0
        phiMT2 = 0
        phiMC3 = 0
        phiMT3 = 0
        
        # longitudinal failure
        if trialStress[0] > 0.0:
            # longitudinal tensile failure critetion
            phiLT = trialStress[0]/Xt # eq. 54
        elif trialStress[0] < 0.0:
            trialStressRot = rotateStresses(trialStress, theta, phi)
            phiLC_MC, phiLC_MT = fail(trialStressRot,St_is,Yt_is,Sl_is,etaL,etaT)
            phiLC = max(phiLC_MC, phiLC_MT)

        # transverse failure (in the 2-direction)
        if trialStress[1] >= 0:
            _, phiMT2 = fail(trialStress,St_is,Yt_is,Sl_is,etaL,etaT)
        elif trialStress[1] < 0:
            phiMC2, _ = fail(trialStress,St_is,Yt_is,Sl_is,etaL,etaT)

        # transverse failure (in the 3 direction)
        if trialStress[2] > 0:
            _, phiMT3 = fail(trialStress,St_is,Yt_is,Sl_is,etaL,etaT)
        elif trialStress[2] < 0:
            phiMC3, _ = fail(trialStress,St_is,Yt_is,Sl_is,etaL,etaT)

        # max failure index
        failI = max(phiLT, phiLC, phiMC2, phiMT2, phiMC3, phiMT3)
        # if failure index approx. equal to 1 -> add to envelope line coords
        if failI < 0.99: continue
        elif fail >= 0.99 and failI < 1.001: 
            envelope.append((stress1,stress2))
            break
        else: break

# plot the failure envelope
x,y = zip(*envelope)
plt.plot(x,y)
plt.grid(True)
plt.xlabel(r'$\sigma_{22}$')
plt.ylabel(r'$\tau_{12}$')
plt.show()

