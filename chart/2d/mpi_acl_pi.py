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
y  = [1, 2.509, 4.205, 5.500, 5.100]
y2 = [1, 3.549, 6.057, 8.570, 11.927]
y3 = [1, 3.729, 7.232, 11.505, 20.412]
scat_plot  = ax.plot(a, y, linewidth=5, label='N = 128')
scat_plot  = ax.plot(a, y2, linewidth=5, label='N = 256')
scat_plot  = ax.plot(a, y3, linewidth=5, label='N = 512')
ax.legend(fontsize=45)
plt.savefig('mpi_acl_pi.png')


