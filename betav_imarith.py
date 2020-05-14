#####################################################
# Python script that makes betaV images which are
# FITS images that has pixel values of the flux ratios
# between M_V = M_g-0.59(M_g-M_r)-0.01 (Jester+05) and 
# Spitzer IRAC 1 band
# written by Duho Kim (1/31/18)        
######################################################
from pyraf import iraf
from astropy.io import fits
from astropy.io import ascii
import numpy as np

work_dir = '/Users/dhk/work/data/NGC_IC/'
cat_dir = '/Users/dhk/work/cat/NGC_IC/'

sha_cat=ascii.read(cat_dir+'sha_quarry_batch_257_without_prefix.txt')

for x in range(0,len(sha_cat)):
#for x in range(0,1):
	name	= sha_cat['id'][x]	# NGC/IC name
	gn 	= work_dir+'SDSS/g/'+name+'-gir.fits'
	rn 	= work_dir+'SDSS/r/'+name+'-rir.fits'
	i1n 	= work_dir+'SHA/SHA_NGC_IC_LONG_mean/'+name+'.fits'
        ##### READ FITS #####
        hdu 	= fits.open(gn)
        g 	= hdu[0].data
        hdu     = fits.open(rn)
        r       = hdu[0].data
        hdu     = fits.open(i1n)
        i       = hdu[0].data

	# Jester+05 All stars R-I < 1.15	
	V = 10**(np.log10(g)+0.59*(np.log10(r/g))+0.004)

	# zeropoint difference factor
	fits.writeto(work_dir+'betav_mean/'+name+'.betaV.fits',V/(i*2.330451129),overwrite=True)
	fits.writeto(work_dir+'betav_mean/'+name+'.betag.fits',g/(i*2.330451129),overwrite=True)
	fits.writeto(work_dir+'betav_mean/'+name+'.betar.fits',r/(i*2.330451129),overwrite=True)

