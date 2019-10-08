import matplotlib.pyplot as plt
import numpy as np 
from math import sin, cos, atan, pi, radians, tan
import pandas as pd

# Python implementation of the Catalanotti failure criteria
# Used to plot failure envelopes
# NOTE: Python 2.7.XX, requires numpy and matplotlib libraries

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
        (((1-sin(alpha0)**2.0)**0.5)*sin(alpha0)*etaL))  # eq 12
phic = (atan((1.0-(1.0-4.0*(Sl_is/Xc)*((Sl_is/Xc)+etaL))**0.5)/
        (2.0*((Sl_is/Xc)+etaL))))  # eq 72
etaT = (etaL*St_is)/Sl_is  # eq 10
kappa = (St_is**2.0-Yt_is**2.0)/(St_is*Yt_is)  # eq 43
lmbda = ((2.0*etaL*St_is)/Sl_is)-kappa  # eq 45 

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

fille = open(r'C:\Workspace\FI.txt','w') 

# Function to determine  f ailure in transverse tension and compression
# and (with rotated coord system) in longitudinal compression 
def fail(tStress,St_is,Yt_is,Sl_is,etaL,etaT):
    # array of angles 0 to pi/2
    angleList = [k*(pi/60) for k in range(0,31)] 
    phiM = 0
    for angle in angleList:
        # tractions: eq 3. see also eqs 59-61
        tN = (tStress[1]*cos(angle)**2.0 + 2.0*tStress[4]
            *sin(angle)*cos(angle) + tStress[2]*sin(angle)**2.0)
        tT = ((-1.0*cos(angle)*sin(angle)*(tStress[1]-tStress[2])) + 
            (tStress[4]*(cos(angle)**2.0 - sin(angle)**2.0)))
        tL = tStress[3]*cos(angle) + tStress[5]*sin(angle) 
        
        if tN >= 0:  # tension
            trialPhiM = ((tN/St_is)**2.0 + (tL/Sl_is)**2.0 + (tT/St_is)**2.0 +
            lmbda*(tN/St_is)*(tL/Sl_is)**2.0 + kappa*(tN/St_is)) # eq. 42 
        else:  # compression
            trialPhiM = ((tL/(Sl_is-etaL*tN))**2.0) + ((tT/(St_is-etaT*tN))**2.0)
        if trialPhiM > phiM: # update if criteria at angle is max 
            phiM = trialPhiM
            aFailC = angle # record failure plane at max fail criteria
    return phiM

# ----------------------------------------------------------------------------
# Function to determine the misalignment angle phi using the Newton-
# Raphson method. 
# NOTE: not converging well, using a different method for now see line 161
def newtonRaphson(initial, tStress, X):
    tol = 1.0*10**-4.0
    x0 = initial
    delta = abs(0.0-f(x0, tStress, X))
    n = 0
    while delta > tol:
        x0 = x0 - f(x0, tStress, X)/df(x0, tStress, X)
        delta = abs(0.0-f(x0, tStress, X))
        n += 1
    return x0

# Definition of function used in the Newton Raphson iteration method
def f(gamma, tStress, X):
    val = (X*gamma + 0.5*(tStress[0] - tStress[1])*
        sin(2.0*gamma) - abs(tStress[3])*cos(2.0*gamma)) # eq. 88
    return val     
      
# Definiton of the derivative function used in NR method
def df(gamma, tStress, X):
    val = (X + (tStress[0] - tStress[1])*cos(2.0*gamma)+
        2.0*abs(tStress[3])*sin(2.0*gamma)) # eq. 89
    return val

# ----------------------------------------------------------------------------
# Function to rotate stresses into misalignment frame for compressive failure
def rotateStresses(tStress, tht, p):
    tStressTheta = np.zeros(6)
    tStressPhi = np.zeros(6)

    # Rotate stresses by angle tht
    tStressTheta[0] = tStress[0] #11
    tStressTheta[1] = (tStress[1]*cos(tht)**2.0 + 
                        2.0*tStress[4]*cos(tht)*sin(tht) + 
                        tStress[2]*sin(tht)**2) #22
    tStressTheta[2] = (tStress[2]*cos(tht)**2.0 - 
                        2.0*tStress[4]*cos(tht)*sin(tht) +
                        tStress[1]*sin(tht)**2.0) #33
    tStressTheta[3] = tStress[3]*cos(tht)+tStress[5]*sin(tht) #12
    tStressTheta[4] = (tStress[4]*(cos(tht)**2 - sin(tht)**2) -
                        tStress[1]*sin(tht)*cos(tht) +
                        tStress[2]*sin(tht)*cos(tht)) #23
    tStressTheta[5] = tStress[5]*cos(tht)-tStress[3]*sin(tht) #13

    # Rotate stresses by angle p
    tStressPhi[0] = (tStressTheta[0]*cos(p)**2 +
                        2.0*tStressTheta[3]*cos(p)*sin(p) +
                        tStressTheta[1]*sin(p)**2)
    tStressPhi[1] = (tStressTheta[1]*cos(p)**2 -
                        2.0*tStressTheta[3]*sin(p)*cos(p) +
                        tStressTheta[0]*sin(p)**2)
    tStressPhi[2] = tStressTheta[2]
    tStressPhi[3] = (tStressTheta[3]*(cos(p)**2 - 
                        sin(p)**2) + tStressTheta[1]*sin(p)*cos(p) -
                        tStressTheta[0]*sin(p)*cos(p))
    tStressPhi[4] = tStressTheta[4]*cos(p) - tStressTheta[5]*sin(p)
    tStressPhi[5] = tStressTheta[5]*cos(p) + tStressTheta[4]*sin(p)

    return tStressPhi, tStressTheta

# ----------------------------------------------------------------------------
# Create failure envelopes

# Xt = 2323.5, Xc = 1200.1, Yt_is = 160.2, Yc_is = 198.0, Sl_is = 130.2

# trialStress specified by Abaqus in the format: (11, 22, 33, 12, 23, 31)

# to plot envelope choose two stresses. stress1 is plotted on x-axis, stress2
# is plotted on the y-axis. e.g. currently sigma_11 plotted against tau_12

envelope = {0:[],1:[],2:[],3:[],4:[],5:[]}
n = 250  # number of points along axes
for stress1 in np.linspace(-1200.1, 2323.5, n):
    for stress2 in np.linspace(0,160.2, n):
        # define loading
        trialStress = np.array([stress1, stress2, 0.0, 0.0, 0.0, 0.0])

        # ----------------------------------------------------------------
        # determine fracture plane angle theta (for fiber kinking)   
        if trialStress[3] == 0 and trialStress[5] == 0:
            if (trialStress[1]-trialStress[2]) == 0:
                theta = pi/4.0
            else:
                # eq 55
                theta = 0.5*atan(
                    (2.0*trialStress[4])/(trialStress[1]-trialStress[2]))
        else: 
            if trialStress[3] == 0:
                theta = pi/2.0
            else:
                # eq 56
                theta = atan(trialStress[5]/trialStress[3])

        # ----------------------------------------------------------------
        # misalignment angle phi using the newton-raphson method
        X = (sin(2.0*phic)*Xc)/(2.0*phic) # eq. 86 
        gammaM = newtonRaphson(0.1, trialStress, X)  # initial guess = 0.1
        if trialStress[3] >= 0:  # eq 77
            phi = gammaM
        else:
            phi = -1.0*(gammaM)

        # ----------------------------------------------------------------
        # determine failure indices
        phiLT = 0
        phiLC = 0
        phiM = 0

        # longitudinal failure
        if trialStress[0] > 0.0:
            # longitudinal tensile failure criterion
            phiLT = trialStress[0]/Xt # eq. 54
        elif trialStress[0] < 0.0:
            trialStressRot, _ = rotateStresses(trialStress, theta, phi)
            phiLC = fail(trialStressRot,St_is,Yt_is,Sl_is,etaL,etaT)

        # transverse failure
        phiM = fail(trialStress,St_is,Yt_is,Sl_is,etaL,etaT)

        # max failure index
        failCriteriaList = (phiLT, phiLC, phiM)
        maxFailCriteria = max(failCriteriaList)
        maxFailIndex = failCriteriaList.index(maxFailCriteria)

        

        # if failure index approx. = 1 -> add to envelope line coords
        if maxFailCriteria < 0.98: continue
        elif maxFailCriteria >= 0.98 and maxFailCriteria < 1.02: 
            envelope[maxFailIndex].append((stress1,stress2))
            fille.write('{},{},{},{}\n'.format(trialStress[0],trialStress[1],theta,phi))
            break
        else: break

# plot the failure envelope
# (red = phiLT, blue = phiLC, green = phiM)
colors = ['r','b','g']
for idd in envelope.keys():
    try:
        x,y = zip(*envelope[idd])
        plt.plot(x,y, color=colors[idd], label='me')
    except ValueError:
        continue

filePath = r"C:\Users\s1342398\OneDrive - University of Edinburgh\Rutger Kok PhD\s11_s22_envelope(pos).csv"
data = pd.read_csv(filePath, skiprows=[0])
s11cln = np.array(data['X'])  # units: kN
s22cln = np.array(data['Y'])  # units: kN
plt.plot(s11cln,s22cln, label='catalanotti')
plt.legend()
plt.grid(True)
plt.xlabel(r'$\sigma_{11}$')
plt.ylabel(r'$\sigma_{22}$')
plt.show()
