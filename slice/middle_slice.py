# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
data=open("/data/s5/yhwu/ELUCID/tree/outputs/consistent_99.list")
L=500

location=np.zeros((L,L))
for i in range(80):
	line=data.readline()

while(line):
	info=line.split()
	m=float(info[2])
	x=int(float(info[8]))
	y=int(float(info[9]))
	z=int(float(info[10]))
	if m>math.pow(10,10.7) and m<math.pow(10,10.72):
		location[x][y]+=1
	line=data.readline()

mean=1.0*np.sum(location)/(L*L)

location=(location-mean)/mean


plt.matshow(location,cmap=plt.cm.Reds,vmin=-1,vmax=0.05*mean)
plt.savefig("/home/yhwu/pic/slice.png",dpi=L)

plt.matshow(location,cmap=plt.cm.Reds,vmin=-1,vmax=0.1*mean)
plt.savefig("/home/yhwu/pic/slice1.png",dpi=L)

plt.matshow(location,cmap=plt.cm.Reds,vmin=-1,vmax=0.2*mean)
plt.savefig("/home/yhwu/pic/slice2.png",dpi=L)


