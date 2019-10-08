import matplotlib.pyplot as plt
import numpy as np 
from math import sin, cos, atan, pi, radians, tan

# Python implementation of the Catalanotti failure criteria
# Used to plot failure envelopes

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Material properties

#IM78552
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

# Initial values
phic = (atan((1.0-(1.0-4.0*(Sl_is/Xc)*((Sl_is/Xc)+etaL))**0.5)/
        (2.0*((Sl_is/Xc)+etaL))))
etaT = (etaL*St_is)/Sl_is  # Eq 10 catalanotti

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

# Function to determine  f ailure in transverse tension and compression
# and in longitudinal com p ression (with rotated coord system)
def fail(trialStress,St_is,Yt_is,Sl_is,etaL,etaT):
    
    angleList = [s*(pi/200) for s in range(0,101)] # array of angles 0 to pi/2

    phiMC = 0.0 # initialize max failure criteria
    phiMT = 0.0

    for angle in angleList:
        tN = (trialStress[1]*cos(angle)**2.0 + 2.0*trialStress[4]
            *sin(angle)*cos(angle) + trialStress[2]*sin(angle)**2.0)
        tT = (-1.0*trialStress[1]*cos(angle)*sin(angle)+ 
            trialStress[2]*sin(angle)*cos(angle) - trialStress[4]*
            (cos(angle)**2.0 - sin(angle)**2.0)) 
        tL = trialStress[3]*cos(angle) + trialStress[5]*sin(angle)  
  
        kappa = (St_is**2.0-Yt_is**2.0)/(St_is*Yt_is) 
        lmbda = ((2.0*etaL*St_is)/Sl_is)-kappa   # CHANGE!
  
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

# Function to determine the misalignment angle phi using the Newton-
# Raphson method (uses functions FUNC and DFUNC as inputs)
def newtonRaphson(initial, tStress, X):
    tol = 1.0*10**-5.0
    x0 = initial
    delta = abs(0.0-f(x0, tStress, X))
    n = 0
    while delta > tol:
        print x0, delta, n
        x0 = x0 - f(x0, tStress, X)/df(x0, tStress, X)
        delta = abs(0.0-f(x0, tStress, X))
        n += 1
    return x0

# Definition of function used in the Newton Raphson iteration method
def f(gamma, tStress, X):
    # val = (phic - gamma - abs((Xc*sin(2.0*gamma))/(2.0*G12) + 
    #     (B/8.0)*((Xc**3.0)*(sin(2.0*gamma)**3.0))))
    val = (X*gamma + 0.5*(tStress[0] - tStress[1]*
        sin(2.0*gamma)) - abs(tStress[3])*cos(2.0*gamma)) # eq. 88
    return val     
      
# Definiton of the derivative function of FUNC
def df(gamma, tStress, X):
    # val = (-1.0 - abs((Xc*cos(2*gamma))/(G12) + 
    #     (3.0*B/4.0)*((Xc**3.0)*(sin(2*gamma)**2.0)*cos(2.0*gamma))))
    val = (X + (tStress[0] - tStress[1])*cos(2.0*gamma)+
        2.0*abs(tStress[3])*sin(2.0*gamma)) # eq. 89
    return val

# Function to rotate stresses into misalignment frame for com. failure
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
# Loading
# Xt = 2323.5, Xc = 1200.1, Yt_is = 160.2, Yc_is = 198.0, Sl_is = 130.2
# (11, 22, 33, 12, 23, 31)

# NOTE envelope s22-tau12: PASSED
# NOTE envelope s22-tau13: PASSED
# NOTE envelope tau12-tau13: PASSED
# NOTE envelope tau12-tau23: PASSED
# NOTE envelope s11-tau12: PASSED (note: noise at tensile failure)

env = []
for s in np.linspace(-198.0, 160.2, 250):
    for t in np.linspace(0.0,200.0,250):
        trialStress = np.array([0.0, s, 0.0, t, 0.0, 0.0])

        # Evaluation of the damage activation functions

        # determine fracture plane angle theta (for fiber kinking)   
        if trialStress[3] == 0 and trialStress[5] == 0:
            # eq 55
            theta = 0.5*atan((2.0*trialStress[4])/(trialStress[1]-trialStress[2]))
        else: 
            # eq 56
            theta = atan(trialStress[5]/trialStress[3])
        if np.isnan(theta):
            theta = 0.0

        # determine the misalignment angle phi
        X = (sin(2.0*phic)*Xc)/(2.0*phic) # eq. 86
        gammaMC = (phic*Xc)/G12
        phi0 = phic - gammaMC
        gammaM = (phi0*G12 + abs(trialStress[3]))/(G12+trialStress[0]-trialStress[1]) - phi0
        # gammaM = newtonRaphson(3.07, trialStress, X) # initial value of 3.0
        if trialStress[3] >= 0:
            phi = phi0 + gammaM
        else:
            phi = -1.0*(phi0+gammaM)

        phiLT = 0
        phiK = 0
        phiMC = 0
        phiMT = 0
        # call subroutine to determine longitudinal compressive failure criteria
        if trialStress[0] >= 0.0:
            # longitudinal tensile failure critetion
            phiLT = trialStress[0]/Xt # eq. 54
        elif trialStress[0] < 0:
            trialStressRot = rotateStresses(trialStress, theta, phi)
            phiKMC, phiKMT = fail(trialStressRot,St_is,Yt_is,Sl_is,etaL,etaT)
            phiK = max(phiKMC, phiKMT)

        if trialStress[1] or trialStress[2] >= 0:
            _, phiMT = fail(trialStress,St_is,Yt_is,Sl_is,etaL,etaT)
        elif trialStress[1] or trialStress[2] < 0:
            phiMC, _ = fail(trialStress,St_is,Yt_is,Sl_is,etaL,etaT)
        # print trialStress
        failI = max(phiLT, phiK,phiMC,phiMT)
        if failI < 0.99: continue
        elif fail >= 0.99 and failI < 1.001: 
            # print '{} | {} | {} | {}'.format(phiLT, phiK, phiMC, phiMT)
            env.append((s,t))
            # break
        else: break

x,y = zip(*env)
plt.plot(x,y)
plt.grid(True)
plt.show()

# print env