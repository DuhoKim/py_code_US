#####################################################
# Python script that psfmatches SDSS images to Spitzer 
# images 
# written by Duho Kim (1/31/18)        
######################################################
from pyraf import iraf
from astropy.io import ascii
import os

work_dir = '/Users/dhk/work/data/NGC_IC/'
cat_dir = '/Users/dhk/work/cat/NGC_IC/'

sha_cat=ascii.read(cat_dir+'sha_quarry_batch_257_without_prefix.txt')
iraf.daophot()

for i in range(0,len(sha_cat)):
	name	= sha_cat['id'][i]	# NGC/IC name
	g	= work_dir+'SDSS/g/'+name+'-g.fits'
	gi	= work_dir+'SDSS/g/'+name+'-gi.fits'
	gir	= work_dir+'SDSS/g/'+name+'-gir.fits'
	r	= work_dir+'SDSS/r/'+name+'-r.fits'
	ri	= work_dir+'SDSS/r/'+name+'-ri.fits'
	rir	= work_dir+'SDSS/r/'+name+'-rir.fits'
	i1	= work_dir+'SHA/SHA_NGC_IC_LONG/'+name+'.fits'
	g_psf 	= work_dir+'SDSS/g_psf6/'+name+'-g.psf3.fits'
	r_psf 	= work_dir+'SDSS/r_psf6/'+name+'-r.psf3.fits'
	i_psf 	= work_dir+'SHA/psf_master/'+name+'.psf.fits'
	
#	iraf.psfmatch(g,reference=i_psf,psfdata=g_psf,output=gi,convolution='psf',kernel='')
	iraf.wregister(gi,reference=i1,output=gir)
#        iraf.psfmatch(r,reference=i_psf,psfdata=r_psf,output=ri,convolution='psf',kernel='')
	iraf.wregister(ri,reference=i1,output=rir)

