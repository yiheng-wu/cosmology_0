# -*- coding: utf-8 -*-
import numpy as np
import time
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
#import MyLib as ml
#import ReadSnap as rs

def bins(xmin,xmax):
	x=[]
	base=math.pow(xmax/xmin,1.0/Nbin)

	for i in range(Nbin):
		x.append(xmin*math.pow(base,i))
	return x


def projection(PK):
	xmin=1.0
	xmax=N/2.0
	x=bins(xmin,xmax)
	block=len(x)
	y=[0]*block
	y_n=[0]*block
	base=math.pow(xmax/xmin,1.0/Nbin)

	for i in range(N//2):
		for j in range(N//2):
			for k in range(N//2):
				dis=math.sqrt((i-N/2-1)**2+(j-N/2-1)**2+(k-N/2-1)**2)
				if dis <N/2 and dis>0:
					y[int(math.log(dis/xmin,base))]+=PK[i][j][k]
					y_n[int(math.log(dis/xmin,base))]+=1
	for i in range(block):
		if y_n[i]!=0:
			y[i]=y[i]/y_n[i]
	return x,y

def projection_k3(PK):
	xmin=1.0
	xmax=N/2.0
	x=bins(xmin,xmax)
	block=len(x)
	y=[0]*block
	y_n=[0]*block
	base=math.pow(xmax/xmin,1.0/Nbin)

	for i in range(N//2):
		for j in range(N//2):
			for k in range(N//2):
				dis=math.sqrt((i-N/2-1)**2+(j-N/2-1)**2+(k-N/2-1)**2)
				if dis <N/2 and dis>0:
					y[int(math.log(dis/xmin,base))]+=PK[i][j][k]*(dis*2*math.pi/L)**3
					y_n[int(math.log(dis/xmin,base))]+=1
	for i in range(block):
		if y_n[i]!=0:
			y[i]=y[i]/y_n[i]
	return x,y

#load the file and calculate power spectrum
def power_spectrum(file1,file2=None):
	if (file2==None):
		d=np.load(file1)
		d=(d*1.0/(np.sum(d)/N**3)-1)
		d=np.fft.fftn(d,norm='ortho')
		d=np.multiply(d,d.conjugate())
		d=np.fft.fftshift(d)
		x,y=projection(d)
		del d
		x=list(np.array(x)*2*math.pi/L)
		y=list(np.array(y)*(L)**3/pow(N,3))
		return x,y
	else:
		d1=np.load(file1)
		d1=(d1*1.0/(np.sum(d1)/N**3)-1)
		d1=np.fft.fftn(d1,norm='ortho')
		d2=np.load(file2)
		d2=(d2*1.0/(np.sum(d2)/N**3)-1)
		d2=np.fft.fftn(d2,norm='ortho')
		p=np.multiply(d1,d2.conjugate())
		del d1,d2
		p=np.fft.fftshift(p)
		x,y=projection(p)
		del p
		x=list(np.array(x)*2*math.pi/L)
		y=list(np.array(y)*L**3/pow(N,3))
		return x,y

def dimensionless_power_spectrum(file1,file2=None):
	if (file2==None):
		d=np.load(file1)
		d=(d*1.0/(np.sum(d)/N**3)-1)
		d=np.fft.fftn(d,norm='ortho')
		d=np.multiply(d,d.conjugate())
		d=np.fft.fftshift(d)
		x,y=projection_k3(d)
		del d
		x=list(np.array(x)*2*math.pi/L)
		y=list(np.array(y)*(L)**3/pow(N,3))
		return x,y
	else:
		d1=np.load(file1)
		d1=(d1*1.0/(np.sum(d1)/N**3)-1)
		d1=np.fft.fftn(d1,norm='ortho')
		d2=np.load(file2)
		d2=(d2*1.0/(np.sum(d2)/N**3)-1)
		d2=np.fft.fftn(d2,norm='ortho')
		p=np.multiply(d1,d2.conjugate())
		del d1,d2
		p=np.fft.fftshift(p)
		x,y=projection_k3(p)
		del p
		x=x*2*math.pi/L
		y=y*L**3/pow(N,3)
		return x,y

def correlation_function(file1,file2=None):
	if (file2==None):
		d=np.load(file1)
		d=(d*1.0/(np.sum(d)/N**3)-1)
		d=np.fft.fftn(d,norm='ortho')
		d=np.multiply(d,d.conjugate())
		d=np.fft.fftn(d,norm='ortho')
		d=np.fft.fftshift(d)
		x,y=projection(d)
		del d
		#x=x*L/N
		#y=y*L**3/pow(N,3)
		return x,y
	else:
		d1=np.load(file1)
		d1=(d1*1.0/(np.sum(d1)/N**3)-1)
		d1=np.fft.fftn(d1,norm='ortho')
		d2=np.load(file2)
		d2=(d2*1.0/(np.sum(d2)/N**3)-1)
		d2=np.fft.fftn(d2,norm='ortho')
		p=np.multiply(d1,d2.conjugate())
		del d1,d2
		p=np.fft.fftshift(p)
		x,y=projection(p)
		del p
		x=list(np.array(x)*2*math.pi/L)
		y=list(np.array(y)*L**3/pow(N,3))
		return x,y

def bias_auto(file1,file2):
	d1=np.load(file1)
	d1=(d1*1.0/(np.sum(d1)/N**3)-1)
	d1=np.fft.fftn(d1,norm='ortho')
	d1=np.multiply(d1,d1.conjugate())
	d1=np.fft.fftshift(d1)
	x1,y1=projection(d1)
	del d1
	d2=np.load(file2)
	d2=(d2*1.0/(np.sum(d2)/N**3)-1)
	d2=np.fft.fftn(d2,norm='ortho')
	d2=np.multiply(d2,d2.conjugate())
	d2=np.fft.fftshift(d2)
	x2,y2=projection(d2)
	del d2
	x=[];y=[]
	for i in range(len(y1)):
		if y1[i]!=0 and y2[i]!=0:
			x.append(x1[i])
			y.append(y2[i]/y1[i])
	x=np.array(x)*2*math.pi/L
	y=np.sqrt(np.abs(y))
	return x,y

def bias_cross(file1,file2):
	d1=np.load(file1)
	d2=np.load(file2)
	d1=(d1*1.0/(np.sum(d1)/N**3)-1)
	d2=(d2*1.0/(np.sum(d2)/N**3)-1)
	d1=np.fft.fftn(d1,norm='ortho')
	d2=np.fft.fftn(d2,norm='ortho')
	d2=np.multiply(d1,d2.conjugate())
	d2=np.fft.fftshift(d2)
	x2,y2=projection(d2)
	del d2
	d1=np.multiply(d1,d1.conjugate())
	d1=np.fft.fftshift(d1)
	x1,y1=projection(d1)
	del d1
	x=[];y=[]
	for i in range(len(y1)):
		if y1[i]!=0 and y2[i]!=0:
			x.append(x1[i])
			y.append(y2[i]/y1[i])
	x=np.array(x)*2*math.pi/L
	y=np.array(y)
	return x,y


###########################################################################
#Power Spectrum
def ps_mm():
	print ("Begin: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	px_mm,py_mm=power_spectrum(dm_dir)
	print ("mmdone: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

	np.save("/home/yhwu/pic/x_mm",np.array(px_mm))
	np.save("/home/yhwu/pic/y_mm",np.array(py_mm))

def ps_mgh():
	

	print ("begin: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	px_mm,py_mm=power_spectrum(dm_dir)
	print ("mm: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	px_mgas,py_mgas=power_spectrum(dm_dir,gas_dir)
	print ("mgas: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	px_mh,py_mh=power_spectrum(dm_dir,halo_dir)
	print ("mh: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

	plt.plot(px_mm,py_mm,'.-')
	plt.plot(px_mgas,py_mgas,'.-')
	plt.plot(px_mh,py_mh,'.-')
	plt.loglog()

	plt.savefig("/home/yhwu/pic/ps1024.png",dpi=1000)

	x=np.vstack((px_mm,px_mgas,px_mh))
	y=np.vstack((py_mm,py_mgas,py_mh))
	np.save("/home/yhwu/pic/x",x)
	np.save("/home/yhwu/pic/y",y)

###########################################################################
#dimensionless Power Spectrum
def ps_mm_dimensionless():
	print ("Begin: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	px_mm,py_mm=dimensionless_power_spectrum(dm_dir)
	print ("mmdone: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

	np.save("/home/yhwu/pic/x_mm",np.array(px_mm))
	np.save("/home/yhwu/pic/y_mm",np.array(py_mm))

def ps_mgh_dimensionless():
	

	print ("begin: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	px_mm,py_mm=dimensionless_power_spectrum(dm_dir)
	print ("mm: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	px_mgas,py_mgas=dimensionless_power_spectrum(dm_dir,gas_dir)
	print ("mgas: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	px_mh,py_mh=dimensionless_power_spectrum(dm_dir,halo_dir)
	print ("mh: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

	plt.plot(px_mm,py_mm,'.-')
	plt.plot(px_mgas,py_mgas,'.-')
	plt.plot(px_mh,py_mh,'.-')
	plt.loglog()

	plt.savefig("/home/yhwu/pic/dimensionless_ps"+str(N)+".png",dpi=1000)

	x=np.vstack((px_mm,px_mgas,px_mh))
	y=np.vstack((py_mm,py_mgas,py_mh))
#	np.save("/home/yhwu/pic/x",x)
#	np.save("/home/yhwu/pic/y",y)

###########################################################################
#Correlation function
def cf_mm():
	print ("Begin: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	cfx_mm,cfy_mm=correlation_function(dm_dir)
	print ("mmdone: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

	np.save("/home/yhwu/pic/cf_x_mm",np.array(cfx_mm))
	np.save("/home/yhwu/pic/cf_y_mm",np.array(cfy_mm))

def cf_mgh():
	

	print ("begin: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	px_mm,py_mm=power_spectrum(dm_dir)
	print ("mm: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	px_mgas,py_mgas=power_spectrum(dm_dir,gas_dir)
	print ("mgas: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	px_mh,py_mh=power_spectrum(dm_dir,halo_dir)
	print ("mh: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

	plt.plot(px_mm,py_mm,'.-')
	plt.plot(px_mgas,py_mgas,'.-')
	plt.plot(px_mh,py_mh,'.-')
	plt.loglog()

	plt.savefig("/home/yhwu/pic/ps1024.png",dpi=1000)

	x=np.vstack((px_mm,px_mgas,px_mh))
	y=np.vstack((py_mm,py_mgas,py_mh))
	np.save("/home/yhwu/pic/x",x)
	np.save("/home/yhwu/pic/y",y)
###########################################################################
#bias
def b_mg():
	print ("Begin: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	b_a_x,b_a_y=bias_auto(dm_dir,halo_dir)
	print ("autodone: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	b_c_x,b_c_y=bias_cross(dm_dir,halo_dir)
	print ("crossdone: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	plt.plot(b_a_x,b_a_y,".-")
	plt.plot(b_c_x,b_c_y,".-")
	plt.legend(['auto','cross'],loc='lower left')
	plt.loglog()
	
	plt.savefig("/home/yhwu/pic/bias"+str(N)+".png",dpi=1000)
#	np.save("/home/yhwu/pic/x_mm",np.array(px_mm))
#	np.save("/home/yhwu/pic/y_mm",np.array(py_mm))

def b_mgh():
	

	print ("begin: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	px_mm,py_mm=dimensionless_power_spectrum(dm_dir)
	print ("mm: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	px_mgas,py_mgas=dimensionless_power_spectrum(dm_dir,gas_dir)
	print ("mgas: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	px_mh,py_mh=dimensionless_power_spectrum(dm_dir,halo_dir)
	print ("mh: "+time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

	plt.plot(px_mm,py_mm,'.-')
	plt.plot(px_mgas,py_mgas,'.-')
	plt.plot(px_mh,py_mh,'.-')
	plt.loglog()

	plt.savefig("/home/yhwu/pic/dimensionless_ps"+str(N)+".png",dpi=1000)

	x=np.vstack((px_mm,px_mgas,px_mh))
	y=np.vstack((py_mm,py_mgas,py_mh))
#	np.save("/home/yhwu/pic/x",x)
###########################################################################
#data
N=1024
L=75 	#Mpc
Nbin=50

dm_dir="/data/s5/yhwu/data/TNG/dm"+str(N)+"_99.npy"
gas_dir="/data/s5/yhwu/data/TNG/gas"+str(N)+"_99.npy"
halo_dir="/data/s5/yhwu/data/TNG/halo"+str(N)+"_99.npy"

###########################################################################
#main

#mgh()
#mm()
#mgh_dimensionless()
b_mg()
#ps_mm()
'''
d1="128.npy"
d2="256.npy"
d3="512.npy"
d4="1024.npy"

N=128
x1,y1=power_spectrum(d1)



#y1=y1*x1*x1*x1
N=256
x2,y2=power_spectrum(d2)
#y2=y2*x2*x2*x2
N=512
x3,y3=power_spectrum(d3)
#y3=y3*x3*x3*x3
N=1024
x4,y4=power_spectrum(d4)

x=np.vstack((x1,x2,x3,x4))
y=np.vstack((y1,y2,y3,y4))
np.save("/home/yhwu/pic/x",x)
np.save("/home/yhwu/pic/y",y)


print (x1,y1)
print (x2,y2)
print (x3,y3)
print (x4,y4)

plt.plot(x1,y1,'.')
plt.plot(x2,y2,'.')
plt.plot(x3,y3,'.')
plt.plot(x4,y4,'.')
plt.loglog()
plt.savefig("/home/yhwu/pic/N3.png")
'''


