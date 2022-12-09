import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

size = 511
step = 16
k = 20
g = 1



for t in range(0, k + 1):
    if (t == 0 or t == 20 or t == k): 
        print('data' + str(t) + '.json' + '\n')
        fig = plt.figure(figsize=(128, 128))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title('Lx=Ly=Lz=1, N=' + str(size) + ', T=1, K=' + str(886) + ', t=' + str(t) + ', s=' + str(step), fontsize = 120)
        a = list(map(lambda x: int(x*step), range(0, 1 + int(size / step))))
        ra = a
        ra.reverse()
        ax.set_xticks(ra)
        ax.set_yticks(a)
        ax.set_zticks(a)
        ax.xaxis.set_tick_params(labelsize=60)
        ax.yaxis.set_tick_params(labelsize=60)
        ax.zaxis.set_tick_params(labelsize=60)
        f = open('data' + str(t) + '.json')
        data = json.load(f)
        transform = interp1d([0, 1],[0, 100])
        x = list(map(lambda x: x['i'], data['data']))
        y = list(map(lambda x: x['j'], data['data']))
        z = list(map(lambda x: x['k'], data['data']))
        c = list(map(lambda x: float(transform(x['err'] * 10000)), data['data']))
        s = list(map(lambda x: 600, data['data']))
        scat_plot  = ax.scatter(x, y, z, s=s, c=c, cmap='hot')
        # cb = plt.colorbar(scat_plot, pad=0.04)
        # cb.set_ticks([-data['max_val'],data['max_val']])
        # cb.set_ticklabels([-data['max_val'], data['max_val']], fontsize = 120)
        f.close()
        plt.savefig('err_' + str(size) + '_' +  str(step) + '_' + str(t) + '_' + str(g) + '.png')
        fig.clear()

