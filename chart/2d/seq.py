import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

fig = plt.figure(figsize=(32, 32))
ax = fig.add_subplot(1,1,1)
a = [128, 256, 512]
ax.set_xticks(a)
ax.xaxis.set_tick_params(labelsize=45)
ax.yaxis.set_tick_params(labelsize=45)
# transform = interp1d([-1, 1],[-5, 100])
y = [2.13, 17.454, 140.62]
y2 = [2.098, 18.686, 148.195]
scat_plot  = ax.plot(a, y, linewidth=5, label='Lx=Ly=Lz = 1.0')
scat_plot  = ax.plot(a, y2, linewidth=5, label='Lx=Ly=Lz = Pi')
ax.legend(fontsize=45)
plt.savefig('seq_time.png')


