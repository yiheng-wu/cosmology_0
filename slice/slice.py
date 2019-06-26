# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

data=open("/data/s5/yhwu/data/SLICE")
L=3000

location=np.zeros((L,L))
a=0
for i in range(0,L):
	for j in range(0,L):
		location[i][j]=int(data.readline())
		a+=location[i][j]
mean=a*1.0/L/L
print (mean)
for i in range(0,L):
	for j in range(0,L):
		location[i][j]=(location[i][j]-mean)/mean



plt.matshow(location,cmap=plt.cm.Reds,vmin=-1,vmax=0.05*mean)
plt.savefig("/home/yhwu/pic/slice.png",dpi=L)

plt.matshow(location,cmap=plt.cm.Reds,vmin=-1,vmax=0.1*mean)
plt.savefig("/home/yhwu/pic/slice1.png",dpi=L)

plt.matshow(location,cmap=plt.cm.Reds,vmin=-1,vmax=0.2*mean)
plt.savefig("/home/yhwu/pic/slice2.png",dpi=L)


