# -*- coding: utf-8 -*-
import numpy as np
import time
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
import MyLib as ml
import ReadSnap as rs


halo_file="/data/s5/yhwu/ELUCID/rockstar/halos_76"
#halo_file="/home/yhwu/work/rockstar/OUTBASE0/halos_0"
halo_file_amount=16
particle_pos_file="p_pos/p_pos"
particle_pos_file_amount=8

h_length=500.0#Mpc
p_length=300000.0#kpc
number=512
h_side_length=h_length/number
p_side_length=p_length/number



class Common_Use:
	def __init__(self):
		self.grid=1

common_use=Common_Use()
grid=common_use.grid

def particle_pos_list(p_pos_file=particle_pos_file,files=particle_pos_file_amount):
	p_pos_list=[]
	for i in range(files):
		p_pos_list.append(p_pos_file+str(i))
	return p_pos_list

def halo_file_list(halofile,files):
	halo_file_list=[]
	for i in range(files):
		halo_file_list.append(halofile+"."+str(i)+".ascii")
	return halo_file_list

# read rockstar file
# get halo position,index means how many pieces do 1E10 needs to be divided by

def Make_Grids(length,number):
	x_grid=[0]
	x_pos=[]
	for i in range(1,number+1):
		x_grid.append(i*length/number)
		x_pos.append(x_grid[i-1]+0.5*length/number)

	grid=[]
	grid.append(x_grid)
	grid.append(x_grid)
	grid.append(x_grid)

	pos=[]
	pos.append(x_pos)
	pos.append(x_pos)
	pos.append(x_pos)

	volume=(length/number)**3

	return grid,pos,volume

def Particle_Number_Density_Sub(particle_pos_file,p_side_length=p_side_length,number=number):
	p_grid=np.zeros(shape=(number,number,number))
	data=open(particle_pos_file,'rb')
	line=data.readline()
	amount=0
	while line:
		amount+=1
		p_info=line.split()
		line=data.readline()

		x=int(float(p_info[1])/p_side_length)
		y=int(float(p_info[2])/p_side_length)
		z=int(float(p_info[3])/p_side_length)
		p_grid[x][y][z]+=1
	return p_grid,amount
			
def Particle_Delta_N(particle_pos_list,number=number):
	p_grid=np.zeros(shape=(number,number,number))

	pool=ProcessPoolExecutor(max_workers=8)
	p_grid_list=list(pool.map(Particle_Number_Density_Sub,particle_pos_list))
	p_total=0
	for i in range(len(p_grid_list)):
		p_grid+=p_grid_list[i][0]
		p_total+=p_grid_list[i][1]
	p_mean_density=p_total*1.0/number**3#mean density per grid
	p_delta_n=(p_grid-p_mean_density)/p_mean_density
	return p_delta_n

def Halo_Number_Density_Sub(halo_file):
	h_grid=np.zeros(shape=(number,number,number))
	data=open(halo_file,'rb')
	for i in range(21):
		line=data.readline()
	amount=0
	while line:
		amount+=1
		h_info=line.split()
		line=data.readline()

		x=int(float(h_info[8])/h_side_length)
		y=int(float(h_info[9])/h_side_length)
		z=int(float(h_info[10])/h_side_length)
		if x>=number:
			x=0
		if y>=number:
			y=0
		if z>=number:
			z=0
		h_grid[x][y][z]+=1
	return h_grid,amount
	
	file_read=open(halo_file,'rb')
	# throw the header
	for i in range(16):
		file_read.readline()
	line=file_read.readline()

	while(line):
		halo_data=line.split()
		line=file_read.readline()
		hm=float(halo_data[2])#halo mass

def Halo_Delta_N(halo_file_list,*halo_bin):
	h_grid=np.zeros(shape=(number,number,number))
	pool=ProcessPoolExecutor(max_workers=16)

	if len(halo_bin)==0:#All halo in one bins
		h_grid_list=list(pool.map(Halo_Number_Density_Sub,halo_file_list))
		h_total=0
		for i in range(len(h_grid_list)):
			h_grid+=h_grid_list[i][0]
			h_total+=h_grid_list[i][1]
		h_mean_density=h_total*1.0/number**3
		h_delta_n=(h_grid-h_mean_density)/h_mean_density
		return h_delta_n	

	devides=int(round(math.log(halo_max/halo_min,10)))*index
	hm_list=[None]*(devides+1)
	for i in range(int(math.log(halo_max/halo_min,10))):
		for j in range(index):
			hm_list[i*index+j]=(j+1)*math.pow(10,i)*halo_min
	hm_list[-1]=halo_max

	h_p=[[]]
	h_m=[[]]
	for i in range(devides-1):
		h_p.append([])
		h_m.append([])




	#n_mean=len(h_p_list)/(length**3)
	#for i in range(0,number):
	#	for j in range(0,number):
	#		for k in range(0,number):
	#			n_grid[i][j][k]=(n_grid[i][j][k]/vol-n_mean)/n_mean
	return n_grid,len(pos_list)
def Number_Density(pos_list,length=h_length,number=number):
	side_length=length/number
	n_grid=np.zeros(shape=(number,number,number))

	for i in range(len(pos_list)):
		x=int(pos_list[i][0]/side_length)
		y=int(pos_list[i][1]/side_length)
		z=int(pos_list[i][2]/side_length)
		if x>=number:
			x=number-1
		if y>=number:
			y=number-1
		if z>=number:
			z=number-1
		n_grid[x][y][z]=n_grid[x][y][z]+1
	
	#n_mean=len(h_p_list)/(length**3)
	#for i in range(0,number):
	#	for j in range(0,number):
	#		for k in range(0,number):
	#			n_grid[i][j][k]=(n_grid[i][j][k]/vol-n_mean)/n_mean
	return n_grid,len(pos_list)

def P_k(*args):
	if len(args)==1:
		delta_k=np.fft.fftn(args[0],norm='ortho')
		p_k=np.multiply(delta_k,delta_k.conjugate())
#		for i in range(number):
#			for j in range(number):
#				for k in range(number):
#					p_k[i][j][k]=p_k[i][j][k]*2/math.pi/math.pi*math.pow((i-number/2)**2+(j-number/2)**2+(k-number/2)**2,1.5)
		return p_k
	else:
		delta_k1=np.fft.fftn(args[0],norm='ortho')
		delta_k2=np.fft.fftn(args[1],norm='ortho')
		p_k=delta_k1*delta_k2.conjugate()
		return p_k

def BDelta2_k(p_k):
	for i in range(number):
		for j in range(number):
			for k in range(number):
				p_k[i][j][k]=p_k[i][j][k]/2/math.pi/math.pi*math.pow((i-number/2)**2+(j-number/2)**2+(k-number/2)**2,1.5)
	return p_k#big delta square

def R2k_r(k_r):
	for i in range(number):
		for j in range(number):
			for k in range(number):
				k_r[i][j][k]=k_r[i][j][k]*((i-number/2)**2+(j-number/2)**2+(k-number/2)**2)
	return k_r#r2*kesi_r
	

def TD_to_1D(pk_or_kr):
	if pk_or_kr[0]==0:#Normal plot
		partition=range(number/2)
		block=len(partition)
		y=[0]*block
		y_n=[0]*block
		for i in range(number):
			for j in range(number):
				for k in range(number):
					distance=math.sqrt(i**2+j**2+k**2)
					distance=math.sqrt((i-number/2)**2+(j-number/2)**2+(k-number/2)**2)
					if distance<number/2:
						y[int(distance)]+=pk_or_kr[1][i][j][k]
						y_n[int(distance)]+=1
		for i in range(block):
			if y_n[i]!=0:
				y[i]=y[i]/y_n[i]
		return partition,y

	elif pk_or_kr[0]==1:#Log plot
		partition=ml.Create_Partition(number/2,0)
		block=len(partition)-1
		y=[0]*block
		y_n=[0]*block
		for i in range(number):
			for j in range(number):
				for k in range(number):
					distance=math.sqrt((i-number/2)**2+(j-number/2)**2+(k-number/2)**2)
					for l in range(block):
						if distance>=partition[l] and distance<partition[l+1]:
							y_n[l]+=1
							y[l]+=pk_or_kr[1][i][j][k]
		del partition[-1]
		for i in range(block):
			if y_n[i]!=0:
				y[i]=y[i]/y_n[i]
		return partition,y


def TD_to_1D_All_Point(pk_or_kr):
	x=[]
	y=[]
	for i in range(number):
		for j in range(number):
			for k in range(number):
				x.append(math.sqrt((i-number/2)**2+(j-number/2)**2+(k-number/2)**2))
				y.append(pk_or_kr[i][j][k])
	xx=list(set(x))
	xx.sort()
	stop=ml.Fib_Search(number/2,xx)
	del xx[stop+1:]
	del yy[stop+1:]

	block=len(xx)
	yy_n=[0]*block
	yy=[0]*block
	for i in range(len(x)):
		j=ml.Fib_Search(x[i],xx)
		if j!=-1:
			yy_n[j]+=1
			yy[j]+=y[i]
	for j in range(block):
		yy[j]=yy[j]/yy_n[j]
	return xx,yy

def main():
	print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
	p_list=particle_pos_list(particle_pos_file,particle_pos_file_amount)
	h_list=halo_file_list(halo_file,halo_file_amount)

	h_delta_n=Halo_Delta_N(h_list)
	h_Pk=P_k(h_delta_n)
	h_bdelta2_k=BDelta2_k(h_Pk)
	h_Kr=np.fft.fftn(h_Pk,norm='ortho')
	h_r2Kr=R2k_r(h_Kr)
	h_Pk=np.fft.fftshift(h_Pk)	
	h_Kr=np.fft.fftshift(h_Kr)
	h_bdelta2_k=np.fft.fftshift(h_bdelta2_k)
	h_r2Kr=np.fft.fftshift(h_r2Kr)

	in_pk=[1]
	in_pk.append(h_Pk)
	in_kr=[1]
	in_kr.append(h_Kr)
	in_bdk=[1]
	in_bdk.append(h_bdelta2_k)
	in_r2kr=[0]
	in_r2kr.append(h_r2Kr)

	in_list=[]
	in_list.append(in_pk)
	in_list.append(in_kr)
	in_list.append(in_bdk)
	in_list.append(in_r2kr)

	pool=ProcessPoolExecutor(max_workers=4)
	out_list=list(pool.map(TD_to_1D,in_list))
	h_pkxy=out_list[0]
	h_krxy=out_list[1]
	h_bdkxy=out_list[2]
	h_r2krxy=out_list[3]

	h_pkxy_x=np.array(h_pkxy[0])*h_side_length*2*math.pi/h_length
#	plt.plot(h_pkxy_x,h_pkxy[1],'*')
#	plt.loglog()
#	plt.savefig('pic_h_pk.png')
#	plt.clf()
#
	h_krxy_x=np.array(h_krxy[0])*h_side_length
#	plt.plot(h_krxy_x,h_krxy[1],'*')
#	plt.semilogx()
#	plt.savefig('pic_h_kr_semi.png')
#	plt.loglog()
#	plt.savefig('pic_h_kr.png')
#	plt.clf()
#
	h_bdkxy_x=np.array(h_bdkxy[0])*h_side_length*2*math.pi/h_length
#	plt.plot(h_bdkxy_x,h_bdkxy[1],'*')
#	plt.loglog()
#	plt.savefig('pic_h_bdk.png')
#	plt.clf()
#
	h_r2krxy_x=np.array(h_r2krxy[0])*h_side_length
#	plt.plot(h_r2krxy_x,h_r2krxy[1],'*')
#	plt.savefig('pic_h_r2kr.png')
#	plt.semilogy()
#	plt.savefig('pic_h_r2kr_semi.png')
#	plt.clf()

	sys.exit()
#
#	p_delta_n=Particle_Delta_N(p_list)
#	p_Pk=P_k(p_delta_n)
#	p_Kr=(np.fft.ifftn(p_Pk,norm='ortho'))
#
#	p_pk_x,p_pk_y=TD_to_1D(p_Pk)
#	plt.plot(p_pk_x,p_pk_y,'*')
#	plt.savefig('pic_pk.png')
#	plt.semilogx()
#	plt.savefig('pic_pk_x.png')
#	plt.loglog()
#	plt.savefig('pic_pk_xy.png')
#	plt.clf()
#	
##	print p_Kr[0][0][0]
##	for i in range(1,number):
##		print p_Kr[i][i][i].real," ",p_Kr[number-i][number-i][number-i].real
##	for i in range(number):
##		print p_Kr[i][i][i]
##	sys.exit()
#
#
#	p_kr_x,p_kr_y=TD_to_1D_r2k_r(p_Kr.real)
#	plt.plot(p_kr_x,p_kr_y,'*')
#	plt.savefig('pic_kr.png')
#	plt.semilogx()
#	plt.savefig('pic_kr_x.png')
#	plt.loglog()
#	plt.savefig('pic_kr_xy.png')
#	plt.clf()



	# every piece of h_pos_list consist of
	#halo_positon list of halo mass bins and relevant halo mass
	h_pos_list=list(pool.map(Halo_Pos,h_list))


	

if __name__=='__main__':
	main()


