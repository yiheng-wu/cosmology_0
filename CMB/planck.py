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
	plt.savefig("/home/yhwu/pic/example_unmasked.png",dpi=3000)
	plt.clf()
	
	map_masked=hp.ma(map)
	map_masked.mask=np.logical_not(mask)

	
	plt.figure()
	hp.mollview(map_masked.filled())
	plt.savefig("/home/yhwu/pic/example_masked.png",dpi=3000)
	plt.clf()


	Cl=hp.anafast(map_masked.filled(),lmax=LMAX)
	ell=np.arange(len(Cl))


	plt.plot(ell,Cl*l*(ell+el1))
	plt.savefig("/home/yhwu/pic/Cl.png",dpi=1000)
	plt.clf()

def masked_map(mask_path,data_path):
	mask=hp.read_map(mask_path)
	map=hp.read_map(data_path)

	map_masked=hp.ma(map)
	map_masked.mask=np.logical_not(mask)

	plt.figure()
	hp.mollview(map_masked.filled())
	plt.savefig("/home/yhwu/pic/masked_map.png",dpi=3000)
	plt.clf()

def Cl_from_map(mask_path,data_path):
	mask=hp.read_map(mask_path)
	map=hp.read_map(data_path)
	
	map_masked=hp.ma(map)
	map_masked.mask=np.logical_not(mask)

	Cl=hp.anafast(map_masked.filled(),lmax=LMAX)
	ell=np.arange(len(Cl))

	return ell,Cl

mask_path="/data/s5/yhwu/work/CMB/example/mask.fits"
map_path="/data/s5/yhwu/work/CMB/example/wmap.fits"

masked_map(mask_path,map_path)
ell,Cl=Cl_from_map(mask_path,map_path)

plt.plot(l,Cl*ell*(ell+1))
:wq
plt.savefig("/home/yhwu/pic/Cl.png")	
quit()


###########################################################################
#read the mask
mask=hp.read_map('/data/s5/yhwu/data/planck/COM_Lensing_4096_R3.00/mask.fits.gz')
#K_mask=hp.ud_grade(mask,nside_dg)
#K_mask[K_mask<1]=0

###########################################################################
#read the alm and transfer to map
klm=hp.read_alm('/data/s5/yhwu/data/planck/COM_Lensing_4096_R3.00/TT/dat_klm.fits')
map=hp.sphtfunc.alm2map(klm,nside)

###########################################################################
#get the Cl from alm
Cl=hp.sphtfunc.alm2cl(klm,lmax=LMAX)
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
plt.savefig("/home/yhwu/pic/cmb.png",dpi=3000)







tmp=np.loadtxt('/data/s5/yhwu/data/planck/COM_Lensing_4096_R3.00/MV/nlkk.dat')
ell=np.array(tmp[:,0])
N=np.array(tmp[:,1])
C=np.array(tmp[:,2]-tmp[:,1])
N_lm=hp.almxfl(K_lm,N/(C+N))
K_lm=hp.almxfl(K_lm,C/(C+N))
K=hp.alm2map(K_lm,nside=nside)











