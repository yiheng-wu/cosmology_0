# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from concurrent.futures import ProcessPoolExecutor

#halo数密度分布
#记得改变box的g大小
def sim_hmf(mass,box):
	#volume
	sim_volum = box**3
	#sorted mass
	masses_sim = np.sort(mass)
	#0,1,2...len(mass)
	n_cumulative_sim = np.arange(len(mass),0,-1)
	#去重后的质量，各元素在原masses_sim中的位置
	masses_sim, unique_indices = np.unique(masses_sim,return_index=True)
	#
	n_cumulative_sim = n_cumulative_sim[unique_indices]/sim_volum
	hmf = np.vstack((masses_sim,n_cumulative_sim))
	print (len(mass))
	print (len(unique_indices))

	return hmf

def dm_hmf(mass,box):
	devide=100
	volum=box**3
	mass.sort()
	dm=(mass[-1]-mass[0])/devide
	m_dm=np.zeros(devide)
	n_dm=np.zeros(devide)

	for i in range(0,devide):
		m_dm[i]=mass[0]+dm*i

	for i in range(0,devide-1):
		for j in range(0,len(mass)):
			if mass[j]>m_dm[i] and mass[j]<m_dm[i+1]:
				n_dm[i]+=1
	for i in range(0,devide):
		n_dm[i]=n_dm[i]/dm/volum
	
	hmf=np.vstack((m_dm,n_dm))
	return hmf	
	
def ascii_mass(name):
	read=open(name)
	for i in range(0,20):
		useless=read.readline()
	halo=[]
	line=read.readline()
	while line:
		data=line.split()
		halo.append(float(data[2]))
		line=read.readline()
	read.close()
	return halo


def out_mass(name):
	read=open(name)
	for i in range(0,16):
		useless=read.readline()
	halo=[]
	line=read.readline()
	while line:
		data=line.split()
		halo.append(float(data[2]))
		line=read.readline()
	read.close()
	return halo


def get_names(name,n):
	names=[]
	for i in range(0,n):
		strn=str(i)
		names.append(name+strn+".ascii")
	return names




name="/data/s5/yhwu/LJing/halos_5000."
#name="/data/s5/yhwu/P512/rockstar/halos_0."
n_block=80

names=get_names(name,n_block)

pool=ProcessPoolExecutor(max_workers=n_block)
masses=[]
for mass in pool.map(ascii_mass,names):
	masses+=mass

hmf=sim_hmf(masses,1200)
plt.plot(hmf[0],hmf[1])

name1="/data/s5/yhwu/ELUCID/rockstar_restart/halos_100."
names1=get_names(name1,n_block)
masses1=[]
for mass in pool.map(ascii_mass,names1):
	masses1+=mass

hmf1=sim_hmf(masses1,500)
plt.plot(hmf1[0],hmf1[1])




plt.loglog()
plt.savefig("/home/yhwu/pic/hmf.png")
		
	