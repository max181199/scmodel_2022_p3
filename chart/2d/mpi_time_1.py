import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

fig = plt.figure(figsize=(32, 32))
ax = fig.add_subplot(1,1,1)
a = [1, 4, 8, 16, 32]
ax.set_xticks(a)
ax.xaxis.set_tick_params(labelsize=45)
ax.yaxis.set_tick_params(labelsize=45)
y  = [1.97, 0.72, 0.41, 0.34, 0.34]
y2 = [11.19, 3.11, 2.30, 1.63, 1.02]
y3 = [91.26, 23.27, 12.7, 6.77, 5.13]
scat_plot  = ax.plot(a, y, linewidth=5, label='N = 128')
scat_plot  = ax.plot(a, y2, linewidth=5, label='N = 256')
scat_plot  = ax.plot(a, y3, linewidth=5, label='N = 512')
ax.legend(fontsize=45)
plt.savefig('mpi_time_1.png')


