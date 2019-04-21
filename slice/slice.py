# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

data=open("/data/s5/yhwu/work/Jing/SLICE")
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

plt.matshow(location,cmap=plt.cm.Reds,vmin=-1,vmax=10)
plt.savefig("/home/yhwu/pic/slice.png",dpi=L)




