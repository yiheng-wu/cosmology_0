import numpy as np
import time

N=1024
position="position.txt"
l=300000/N



grid=np.zeros((N,N,N),dtype=np.float32)
data=open(position,'rb')
line=data.readline()
while line:
	pos=line.split()
	x=int(float(pos[1])/l)
	y=int(float(pos[2])/l)
	z=int(float(pos[3])/l)
	grid[x][y][z]+=1
	line=data.readline()
grid=grid/(np.sum(grid)/N**3)-1
np.save(str(N),grid)

