# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

L=751
grid=np.load("grid.npy")
piece=200
location=grid[piece]
for i in range(10):
	location+=grid[piece+i+1]
mean=np.sum(location)/(L*L)
location=(location-mean)/mean
print location

plt.matshow(location,cmap=plt.cm.Reds,vmin=-1,vmax=0.05*mean)
plt.savefig("/home/yhwu/pic/slice.png",dpi=L)

plt.matshow(location,cmap=plt.cm.Reds,vmin=-1,vmax=0.1*mean)
plt.savefig("/home/yhwu/pic/slice1.png",dpi=L)

plt.matshow(location,cmap=plt.cm.Reds,vmin=-1,vmax=0.2*mean)
plt.savefig("/home/yhwu/pic/slice2.png",dpi=L)


