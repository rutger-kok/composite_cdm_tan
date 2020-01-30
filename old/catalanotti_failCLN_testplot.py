import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# fortran data
filePath = r"C:\Workspace\fail.txt"
data = pd.read_csv(filePath)
s11 = np.array(data['trialStress(1)'])  # units: kN
s22 = np.array(data['trialStress(2)'])  # units: kN
s33 = np.array(data['trialStress(3)'])  # units: kN
s12 = np.array(data['trialStress(4)'])  # units: kN
s23 = np.array(data['trialStress(5)'])  # units: kN
s13 = np.array(data['trialStress(6)'])  # units: kN
fi_LT = np.array(data['FI_LT'])  # units: -
fi_LC = np.array(data['FI_LC'])  # units: -
fi_MT = np.array(data['FI_MT'])  # units: -
fi_MC = np.array(data['FI_MC'])  # units: -
theta = np.array(data['theta'])  # units: -
phi = np.array(data['phi'])  # units: -

# # python data
# filePath2 = r"C:\Workspace\FI.txt"
# data2 = pd.read_csv(filePath2)
# s11_2 = np.array(data2['trialStress1'])/1000.0  # units: kN
# s22_2 = np.array(data2['trialStress2'])/1000.0  # units: kN
# theta_2 = np.array(data2['theta'])  # units: -
# phi_2 = np.array(data2['phi'])  # units: -

# plot
fig = plt.figure()
ax = fig.add_subplot(111)
# ax.plot(fi_LC, label='FI_LC')
ax.scatter(s11,s12)
# ax.scatter(s11,theta, label='fortran theta')
# ax.plot(s11_2,phi_2, label='python phi')
# ax.plot(s11,phi, label='fortran phi')
# ax.plot(s11_2,s22_2, label='python')
# ax.plot(s11,s22, label='fortran')
plt.legend()
plt.grid(True)
plt.xlabel('s22')
plt.ylabel('s11')
plt.show()
