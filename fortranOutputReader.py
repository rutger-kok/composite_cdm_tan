import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import os

filePath = r"C:\Workspace\catalanotti_failure_criteria\FI_data _c.txt"
data = pd.read_csv(filePath)
s11 = np.array(data['trialStress1'])  # units: kN
s22 = np.array(data['trialStress2'])  # units: kN
s33 = np.array(data['trialStress3'])  # units: kN
s12 = np.array(data['trialStress4'])  # units: kN
s23 = np.array(data['trialStress5'])  # units: kN
s13 = np.array(data['trialStress6'])  # units: kN
e11 = np.array(data['trialStrain1'])  # units: kN
e22 = np.array(data['trialStrain2'])  # units: kN
e33 = np.array(data['trialStrain3'])  # units: kN
e12 = np.array(data['trialStrain4'])  # units: kN
e23 = np.array(data['trialStrain5'])  # units: kN
e13 = np.array(data['trialStrain6'])  # units: kN
d1 = np.array(data['xOmega1'])  # units: kN
d2 = np.array(data['xOmega2'])  # units: kN
unused1 = np.array(data['xOmega3'])  # units: kN
d6 = np.array(data['xOmega4'])  # units: kN
unused2 = np.array(data['xOmega5'])  # units: kN
unused3 = np.array(data['xOmega6'])  # units: kN
# fi_LT = np.array(data['FI_LT'])  # units: -
# fi_LC = np.array(data['FI_LC'])  # units: -
# fi_MT = np.array(data['FI_MT'])  # units: -
# fi_MC = np.array(data['FI_MC'])  # units: -

fig = plt.figure()
ax = fig.add_subplot(111)
# ax.plot(fi_LC, label='FI_LC')
# ax.plot(fi_LT, label='FI_LT')
# ax.plot(fi_MC, label='FI_MC')
# ax.plot(fi_MT, label='FI_MT')
# ax.legend()
ax.plot(e22,s22)
ax.plot(e22,d2)
# ax.plot(s22)
# ax.plot(s33)
# ax.plot(s12)
# ax.plot(s23)
# ax.plot(s13)
plt.grid(True)
plt.xlabel('le11')
plt.ylabel('s11')
plt.show()
print len(s11)