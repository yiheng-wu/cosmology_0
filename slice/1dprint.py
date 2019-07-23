# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

data=open("/data/s5/yhwu/work/Jing/2")
#data=open("/data/s5/yhwu/data/SLICE")
L=600
line=data.readline()
x=np.zeros(L)
y=np.zeros(L)

for i in range(0,L):
	x[i]=i

while(line):
	info=line.split()
	y[int(float(info[0]))]+=1
	line=data.readline()

plt.plot(x,y,'.')
plt.semilogy()
plt.savefig("/home/yhwu/pic/1D1.png")


