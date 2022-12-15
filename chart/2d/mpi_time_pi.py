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
y  = [1.43, 0.57, 0.34, 0.26, 0.28]
y2 = [11.57, 3.26, 1.91, 1.35, 0.97]
y3 = [93.08, 24.96, 12.87, 8.09, 4.56]
scat_plot  = ax.plot(a, y, linewidth=5, label='N = 128')
scat_plot  = ax.plot(a, y2, linewidth=5, label='N = 256')
scat_plot  = ax.plot(a, y3, linewidth=5, label='N = 512')
ax.legend(fontsize=45)
plt.savefig('mpi_time_pi.png')


