import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

filePath = r"C:\Workspace\catalanotti_failure_criteria\fi_lc_test.txt"
data = pd.read_csv(filePath)
s11 = np.array(data['trialStress1'])  # units: kN
theta = np.array(data['theta'])  # units: kN
phi = np.array(data['phi'])  # units: kN
fi_LC = np.array(data['FI_LC'])  # units: -
fi_LT = np.array(data['FI_LT'])

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# scatterplot = ax.scatter(theta,phi,s11)

# plt.show(fig)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
plot = ax2.plot(s11,fi_LC)
plot = ax2.plot(s11,fi_LT)
plt.show(fig2)