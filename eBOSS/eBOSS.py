import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import healpy as hp
from pandas import value_counts
import treecorr as tc

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
	potential=klm*2.0/(ell*(ell+1))
	
	kmap=hp.sphtfunc.alm2map(klm,nside)
	nmap=hp.sphtfunc.alm2map(nlm,nside)
	lsmap=hp.sphtfunc.alm2map(ls,nside)
	potentialmap=hp.sphtfunc.alm2map(potential,nside)
	
	kappa=hp.sphtfunc.alm2cl(klm,lmax=LMAX)
	kappa_name="kappa_cl"
	Cl_plot(kappa,kappa_name)
	lensing=hp.sphtfunc.alm2cl(ls,lmax=LMAX)
	lensing_name="lensing_cl"
	Cl_plot(lensing,lensing_name)
	potential=hp.sphtfunc.alm2cl(potential,lmax=LMAX)
	potential_name="potential_cl"
	Cl_plot(potential,potential_name)

	dec,ra=np.rad2deg(hp.pix2ang(nside,range(hp.nside2npix(nside))))
	dec=90.-dec
	ra=ra[mask>0]
	dec=dec[mask>0]
	kmap=kmap[mask>0]

	return ra,dec,kmap
	
	plt.figure()
	hp.mollview(kmap)
	plt.savefig("/home/yhwu/pic/kappa.png",dpi=1000)
	plt.clf()
	hp.mollview(nmap)
	plt.savefig("/home/yhwu/pic/noise.png",dpi=1000)
	plt.clf()
	hp.mollview(lsmap)
	plt.savefig("/home/yhwu/pic/lensing.png",dpi=1000)
	plt.clf()
	hp.mollview(potentialmap)
	plt.savefig("/home/yhwu/pic/potential.png",dpi=1000)
	plt.clf()

def load_galaxy(path):
	data=fits.open(path)
	galaxy=[]
	for i in data[1].data:
		if (i['IMATCH']==1 or i['IMATCH']==2) and i['SPECTYPE']=="GALAXY":
			galaxy.append(np.array(i[:3]))
	data.close()
	return np.array(galaxy)

def ecliptic2galactic(ra,dec):
	galactic_coordinate=SkyCoord(ra=ra*u.degree,dec=dec*u.degree,frame='fk5')
	gal_l=galactic_coordinate.galactic.l.deg
	gal_b=galactic_coordinate.galactic.b.deg
	return gal_l,gal_b

def galactic_pixel(l,b):
	pixel=hp.ang2pix(nside,np.pi/2.-np.deg2rad(b),np.deg2rad(l))
	pixel_vc=value_counts(pixel)
	mean_vc=np.mean(pixel_vc.values)
	pixel_vc-=mean_vc

	[decc,raa]=np.rad2deg(hp.pix2ang(nside,pixel_vc.index))
	decc=90.-decc

	return raa,decc,pixel_vc.values
	
	gal_plot=np.ones(hp.nside2npix(nside))*(-mean_vc)
	gal_plot[pixel_vc.index]=pixel_vc.values
	plt.figure()
	hp.mollview(gal_plot)
	plt.savefig("/home/yhwu/pic/galaxy_density.png",dpi=1000)
	plt.clf()

def load_rand(path):
	data=fits.open(path)
	ra=data[1].data['RA']
	dec=data[1].data['DEC']
	z=data[1].data['Z']
	data.close()
	rand_gal=np.vstack((ra,dec,z))
	rand_gal=rand_gal.T#transpose()
	return np.array(rand_gal)

def NNCorr(galaxy,rand_galaxy):
	#galaxy
	ra=galaxy[:,0]
	dec=galaxy[:,1]
	z=galaxy[:,2]

	l,b=ecliptic2galactic(ra,dec)
	raa,decc,k=galactic_pixel(l,b)

	galaxy_catalogue=tc.Catalog(ra=raa,dec=decc,k=k,ra_units="deg",dec_units="deg")
	gg=tc.NNCorrelation(nbins=NBINS,min_sep=MIN_SEP,max_sep=MAX_SEP,bin_slop=0.01,verbose=0,sep_units='degrees')
	gg.process(galaxy_catalogue)

	#rand_galaxy
	rra=rand_gal[:,0]
	rdec=rand_gal[:,1]
	rz=rand_gal[:,2]

	rl,rb=ecliptic2galactic(rra,rdec)

	rraa,rdecc,rk=galactic_pixel(rl,rb)

	rand_gal_catalogue=tc.Catalog(ra=rraa,dec=rdecc,k=rk,ra_units="deg",dec_units="deg")
	rr=tc.NNCorrelation(nbins=NBINS,min_sep=MIN_SEP,max_sep=MAX_SEP,bin_slop=0.01,verbose=0,sep_units='degrees')
	rr.process(rand_gal_catalogue)

	#xi-rr
	xi,varxi=gg.calculateXi(rr)
	r=np.exp(gg.meanlogr)
	sig=np.sqrt(varxi)

	plt.plot(r,xi,color='blue')
	plt.plot(r,-xi,color='blue',ls=':')
	plt.loglog()
	plt.savefig("/home/yhwu/pic/xi_rr.png",png=1000)
	plt.clf()

	#xi-dr
	dr=tc.NNCorrelation(nbins=NBINS,min_sep=MIN_SEP,max_sep=MAX_SEP,bin_slop=0.01,verbose=0,sep_units='degrees')
	dr.process(galaxy_catalogue,rand_gal_catalogue)

	xi,varxi=gg.calculateXi(rr,dr)
	r=np.exp(gg.meanlogr)
	sig=np.sqrt(varxi)

	plt.plot(r,xi,color='blue')
	plt.plot(r,-xi,color='blue',ls=':')
	plt.loglog()
	plt.savefig("/home/yhwu/pic/xi_dr.png",png=1000)

def NKCorr(galaxy):

	ra=galaxy[:,0]
	dec=galaxy[:,1]
	z=galaxy[:,2]

	l,b=ecliptic2galactic(ra,dec)
	raa,decc,k=galactic_pixel(l,b)

	galaxy_catalogue=tc.Catalog(ra=raa,dec=decc,k=k,ra_units="deg",dec_units="deg")
	k_ra,k_dec,k_k=denoise(mask_path,alm_path,tmp_path)
	CMB_catalogue=tc.Catalog(ra=k_ra,dec=k_dec,k=k_k,ra_units="deg",dec_units="deg")

	nk=tc.NKCorrelation(nbins=NBINS,min_sep=MIN_SEP,max_sep=MAX_SEP,bin_slop=0.01,verbose=0,sep_units='degrees')
	nk.process(galaxy_catalogue,CMB_catalogue)

	xi=nk.xi
	r=np.exp(nk.meanlogr)
	
	plt.plot(r,xi,color='blue')
	plt.plot(r,-xi,color='blue',ls=':')
	plt.loglog()
	plt.savefig("/home/yhwu/pic/xi_cross.png",png=1000)






#CMB parameter
nside=2048
nside_degree=512
LMAX=4096

#Corr parameter
NBINS=20
MIN_SEP=0.05
MAX_SEP=50.0

elg_path="/data/s5/yhwu/data/eBOSS/eBOSS_ELG_full_ALL_v5.dat.fits"
elg_rand_path="/data/s5/yhwu/data/eBOSS/eBOSS_ELG_full_ALL_v5.ran.fits"
lrg_path="/data/s5/yhwu/data/eBOSS/eBOSS_LRG_full_ALL_v5.dat.fits"


mask_path='/data/s5/yhwu/data/planck/COM_Lensing_4096_R3.00/mask.fits.gz'
alm_path='/data/s5/yhwu/data/planck/COM_Lensing_4096_R3.00/MV/dat_klm.fits'
tmp_path='/data/s5/yhwu/data/planck/COM_Lensing_4096_R3.00/MV/nlkk.dat'

#galaxy=load_galaxy(elg_path)
#np.save("tem/galaxy",galaxy)
galaxy=np.load("tem/galaxy.npy")


NKCorr(galaxy)
quit()











rand_gal=load_rand(elg_rand_path)

NNCorr(galaxy,rand_gal)




quit()


elg=fits.open(elg_path)
lrg=fits.open(lrg_path)
redshiftbin()





































