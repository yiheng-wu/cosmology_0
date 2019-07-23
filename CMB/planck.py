nside=2048
nside_dg=512
LMAX=1200



import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import math


def example():
	mask=hp.read_map("/data/s5/yhwu/work/CMB/example/mask.fits")
	map=hp.read_map("/data/s5/yhwu/work/CMB/example/wmap.fits")

	
	plt.figure()
	hp.mollview(map,norm="hist",min=-1,max=1)
	plt.savefig("/home/yhwu/pic/example_unmasked.png",dpi=1000)
	plt.clf()
	
	map_masked=hp.ma(map)
	map_masked.mask=np.logical_not(mask)

	
	plt.figure()
	hp.mollview(map_masked.filled())
	plt.savefig("/home/yhwu/pic/example_masked.png",dpi=1000)
	plt.clf()


	Cl=hp.anafast(map_masked.filled(),lmax=LMAX)
	ell=np.arange(len(Cl))


	plt.plot(ell,Cl*l*(ell+el1))
	plt.savefig("/home/yhwu/pic/Cl.png",dpi=1000)
	plt.clf()

def map_masked(mask_path,data_path):
	mask=hp.read_map(mask_path)
	map=hp.read_map(data_path)

	map_masked=hp.ma(map)
	map_masked.mask=np.logical_not(mask)

	plt.figure()
	hp.mollview(map_masked.filled())
	plt.savefig("/home/yhwu/pic/map_masked.png",dpi=1000)
	plt.clf()

def alm2map_masked(mask_path,data_path):
	mask=hp.read_map(mask_path)
	alm=hp.read_alm(data_path)
	map=hp.sphtfunc.alm2map(alm,nside)

	map_masked=hp.ma(map)
	map_masked.mask=np.logical_not(mask)

	plt.figure()
	hp.mollview(map_masked)
	plt.savefig("/home/yhwu/pic/alm2map_masked.png",dpi=1000)
	plt.clf()

def Cl_from_map_mask(mask_path,data_path):
	mask=hp.read_map(mask_path)
	map=hp.read_map(data_path)
	
	map_masked=hp.ma(map)
	map_masked.mask=np.logical_not(mask)

	Cl=hp.anafast(map_masked.filled(),lmax=LMAX)
	ell=np.arange(len(Cl))

	return ell,Cl

def Cl_from_map(data_path):
	map=hp.read_map(data_path)
	Cl=hp.anafast(map,lmax=LMAX)
	ell=np.arange(len(Cl))
	return ell,Cl

def Cl_from_alm(data_path):
	alm=hp.read_alm(data_path)
	Cl=hp.sphtfunc.alm2cl(alm,lmax=LMAX)
	ell=np.arange(len(Cl))
	return ell,Cl

def Cl_plot(Cl,name):
	ell=np.arange(len(Cl))
	plt.plot(ell,Cl*ell*(ell+1))
	plt.savefig("/home/yhwu/pic/"+name+".png",dpi=1000)
	plt.clf()


def denoise(mask_path,data_path,tmp_path):
	mask=hp.read_map(mask_path)
	klm=hp.read_alm(data_path)
	tmp=np.loadtxt(tmp_path)

	ell=np.array(tmp[:,0])
	N=np.array(tmp[:,1])
	C=np.array(tmp[:,2]-tmp[:,1])
	nlm=hp.sphtfunc.almxfl(klm,N/(C+N))
	klm=hp.sphtfunc.almxfl(klm,C/(C+N))

	ell,m=hp.sphtfunc.Alm.getlm(lmax=LMAX)
	ell+=1
	ls=klm*2.0/np.sqrt(ell*(ell+1))
	lss=klm*2.0/(ell*(ell+1))
	
	kmap=hp.sphtfunc.alm2map(klm,nside)
	nmap=hp.sphtfunc.alm2map(nlm,nside)
	lsmap=hp.sphtfunc.alm2map(ls,nside)
	lssmap=hp.sphtfunc.alm2map(lss,nside)
	
	kappa=hp.sphtfunc.alm2cl(klm,lmax=LMAX)
	kappa_name="kappa_cl"
	#Cl_plot(kappa,kappa_name)
	lensing=hp.sphtfunc.alm2cl(ls,lmax=LMAX)
	lensing_name="lensing_filtered_cl"
	#Cl_plot(lensing,lensing_name)

	
	plt.figure()
	'''
	hp.mollview(kmap)
	plt.savefig("/home/yhwu/pic/kappa.png",dpi=1000)
	plt.clf()
	hp.mollview(nmap)
	plt.savefig("/home/yhwu/pic/noise.png",dpi=1000)
	plt.clf()
	'''
	hp.mollview(lsmap)
	plt.savefig("/home/yhwu/pic/lensing.png",dpi=1000)
	plt.clf()
	hp.mollview(lssmap)
	plt.savefig("/home/yhwu/pic/lss.png",dpi=1000)
	plt.clf()




'''
LMAX=300
isw_path="/data/s5/yhwu/data/planck/COM_CompMap_ISW_0064_R2.00.fits"
isw="isw_Cl"
ell,Cl=Cl_from_map(isw_path)
Cl_plot(ell,Cl,isw)
quit()
'''
LMAX=4096
mask_path_example="/data/s5/yhwu/work/CMB/example/mask.fits"
map_path_example="/data/s5/yhwu/work/CMB/example/wmap.fits"

mask_path='/data/s5/yhwu/data/planck/COM_Lensing_4096_R3.00/mask.fits.gz'
alm_path='/data/s5/yhwu/data/planck/COM_Lensing_4096_R3.00/MV/dat_klm.fits'
tmp_path='/data/s5/yhwu/data/planck/COM_Lensing_4096_R3.00/MV/nlkk.dat'

denoise(mask_path,alm_path,tmp_path)
quit()


LMAX=300
Cl_pic="Cl"
ell,Cl=Cl_from_alm(alm_path)
Cl_plot(ell,Cl,Cl_pic)
quit()


nside=2048
alm2map_masked(mask_path,alm_path)
quit()




masked_map(mask_path,map_path)
ell,Cl=Cl_from_map(mask_path,map_path)

plt.plot(l,Cl*ell*(ell+1))
plt.savefig("/home/yhwu/pic/Cl.png")	
quit()



###########################################################################
#read the alm and transfer to map
klm=hp.read_alm('/data/s5/yhwu/data/planck/COM_Lensing_4096_R3.00/TT/dat_klm.fits')
map=hp.sphtfunc.alm2map(klm,nside)

###########################################################################
#get the Cl from alm
l=np.arange(len(Cl))
#plt.plot(l,Cl)
plt.plot(l,Cl*l*(l+1))
#plt.plot(l,Cl*np.sqrt(l)*(l+1))
plt.savefig("/home/yhwu/pic/Cl.png",dpi=1000)
quit()

###########################################################################
#mask the map with the mask
map_masked=hp.ma(map)
map_masked.mask=np.logical_not(mask)

###########################################################################
#get the Cl from map
Cl=hp.anafast(map_masked.filled(),lmax=LMAX)
l=np.arange(len(Cl))
#plt.plot(l,Cl)
plt.plot(l,Cl*l*(l+1))
plt.savefig("/home/yhwu/pic/Cl.png",dpi=1000)
quit()


###########################################################################
#plot the map
plt.figure()
hp.mollview(map_masked)
#hp.mollview(map_masked.filled())
#hp.mollview(map_masked.filled(),coord=["G","E"])#Galactic coordinates to Ecliptic coordinates
plt.savefig("/home/yhwu/pic/cmb.png",dpi=1000)







tmp=np.loadtxt('/data/s5/yhwu/data/planck/COM_Lensing_4096_R3.00/MV/nlkk.dat')
ell=np.array(tmp[:,0])
N=np.array(tmp[:,1])
C=np.array(tmp[:,2]-tmp[:,1])
N_lm=hp.almxfl(K_lm,N/(C+N))
K_lm=hp.almxfl(K_lm,C/(C+N))
K=hp.alm2map(K_lm,nside=nside)











