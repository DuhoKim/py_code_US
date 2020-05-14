from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import random as rd
from scipy import interpolate

hgpsf=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/gr.psf.160.out.5.60D.fits')
hipsf=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/i2.psf.160.out.6.60D.fits')

gpsf=hgpsf[0].data
ipsf=hipsf[0].data

x,y = np.ogrid[0:61,0:61]
out = 15**2 < (x-30)**2 + (y-30)**2

gadd = gpsf[out].mean()
iadd = ipsf[out].mean()

gsum = gpsf + gadd
isum = ipsf + iadd

gsum[out] = 0
isum[out] = 0

hgpsf[0].data=gsum
hipsf[0].data=isum

hgpsf.writeto('/Users/dhk/work/data/NGC_IC/pilot4/gr.psf.160.out.5.60D.add.flat.fits',overwrite=True)
hipsf.writeto('/Users/dhk/work/data/NGC_IC/pilot4/i2.psf.160.out.6.60D.add.flat.fits',overwrite=True)

hgpsf.close()
hipsf.close()
