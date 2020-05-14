from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import random as rd

hgpsf=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/gr.psf.160.out.5.fits')
hipsf=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/i2.psf.160.out.6.fits')

gpsf=hgpsf[0].data
ipsf=hipsf[0].data

x,y = np.ogrid[0:161,0:161]
for r in range(30,162):
	index1= (r-1)**2 < (x-80)**2 + (y-80)**2 
	index2= (x-80)**2 + (y-80)**2	<= r**2
	index = index1 * index2
	anul=gpsf[index]
	gpsf[index]=rd.sample(anul,len(anul))
	anul=ipsf[index]
        ipsf[index]=rd.sample(anul,len(anul))


hgpsf[0].data=gpsf
hipsf[0].data=ipsf

hgpsf.writeto('/Users/dhk/work/data/NGC_IC/pilot4/gr.psf.160.out.5.flat.fits',overwrite=True)
hipsf.writeto('/Users/dhk/work/data/NGC_IC/pilot4/i2.psf.160.out.6.flat.fits',overwrite=True)

hgpsf.close()
hipsf.close()

