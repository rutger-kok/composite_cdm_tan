import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

filePath = r"C:\Workspace\fail2.txt"
data = pd.read_csv(filePath)
a = np.array(data['angle'])  # units: kN alpha
fi = np.array(data['func'])  # units: kN FI

fig = plt.figure()
ax = fig.add_subplot(111)
plot = ax.plot(a,fi)

plt.show(fig)
