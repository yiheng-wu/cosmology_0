# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time

print time.time()

L=751
grid=np.zeros((L,L,L))

for i in range(64):
	print time.asctime(time.localtime(time.time()))
	dir="/data/s5/yhwu/data/TNG/snap/position0/snap_099."+str(i)+".hdf5_position"
	data=open(dir)
	line=data.readline()
	info=line.split()

	while line:
		info=line.split()
		x=float(info[0])*10
		y=float(info[1])*10
		z=float(info[2])*10

		grid[int(x)][int(y)][int(z)]+=1	

		line=data.readline()
	data.close()

np.save("grid",grid)
mean=np.sum(grid)*1.0/L/L
grid=(grid-mean)/mean

plt.matshow(grid,cmap=plt.cm.Reds,vmin=-1,vmax=0.05*mean)
plt.savefig("/home/yhwu/pic/slice.png",dpi=L)

plt.matshow(grid,cmap=plt.cm.Reds,vmin=-1,vmax=0.1*mean)
plt.savefig("/home/yhwu/pic/slice1.png",dpi=L)

plt.matshow(grid,cmap=plt.cm.Reds,vmin=-1,vmax=0.2*mean)
plt.savefig("/home/yhwu/pic/slice2.png",dpi=L)


quit()










L=3000

grid=np.zeros((L,L))
a=0
for i in range(0,L):
	for j in range(0,L):
		grid[i][j]=int(data.readline())
		a+=grid[i][j]
mean=a*1.0/L/L
print (mean)
for i in range(0,L):
	for j in range(0,L):
		grid[i][j]=(grid[i][j]-mean)/mean



plt.matshow(grid,cmap=plt.cm.Reds,vmin=-1,vmax=0.05*mean)
plt.savefig("/home/yhwu/pic/slice.png",dpi=L)

plt.matshow(grid,cmap=plt.cm.Reds,vmin=-1,vmax=0.1*mean)
plt.savefig("/home/yhwu/pic/slice1.png",dpi=L)

plt.matshow(grid,cmap=plt.cm.Reds,vmin=-1,vmax=0.2*mean)
plt.savefig("/home/yhwu/pic/slice2.png",dpi=L)


