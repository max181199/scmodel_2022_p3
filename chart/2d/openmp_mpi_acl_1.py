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
y  = [1, 1.797, 2.755, 3.351]
y2 = [1, 2.090, 4.042, 7.187]
y3 = [1, 2.007, 3.989, 7.312]
scat_plot  = ax.plot(a, y, linewidth=5, label='N = 128')
scat_plot  = ax.plot(a, y2, linewidth=5, label='N = 256')
scat_plot  = ax.plot(a, y3, linewidth=5, label='N = 512')
ax.legend(fontsize=45)
plt.savefig('openmp_mpi_acl_1.png')


