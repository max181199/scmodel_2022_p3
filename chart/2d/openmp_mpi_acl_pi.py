import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

fig = plt.figure(figsize=(32, 32))
ax = fig.add_subplot(1,1,1)
a = [1, 2, 4, 8]
ax.set_xticks(a)
ax.xaxis.set_tick_params(labelsize=45)
ax.yaxis.set_tick_params(labelsize=45)
y  = [1, 1.925, 3.225, 4.961]
y2 = [1, 1.969, 3.683, 6.748]
y3 = [1, 2.020, 4.040, 7.493]
scat_plot  = ax.plot(a, y, linewidth=5, label='N = 128')
scat_plot  = ax.plot(a, y2, linewidth=5, label='N = 256')
scat_plot  = ax.plot(a, y3, linewidth=5, label='N = 512')
ax.legend(fontsize=45)
plt.savefig('openmp_mpi_acl_pi.png')


