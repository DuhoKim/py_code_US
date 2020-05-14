from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import random as rd
from scipy import interpolate

hgpsf=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/gr.psf.160.5.grad2.out.fits')
hipsf=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/i2.psf.160.6.flat44.out.fits')

gpsf=hgpsf[0].data
ipsf=hipsf[0].data

fig = plt.figure(figsize=(10,10))

pix = range(0,161)
ipix = range(0,161,160)
for r in range(2,159):
	h=fig.add_axes([(r % 13)*(1.0/13),(1-1.0/13)-(r / 13)*(1.0/13),1.0/13,1.0/13])	
	vert1=gpsf[:,r]
	(i1) = vert1.nonzero()
	vert2=ipsf[:,r]
	(i2) = vert2.nonzero()
	h.set_ylim([-0.1,0.1])
	h.set_axis_off()
	h.plot(pix[i2],vert2[i2]+0.05)
	f = interpolate.interp1d(pix[i1],vert2[i2])
	h.plot(ipix,f(ipix)+0.05)
	h.text(0,0.025,'%d' % (f(ipix)[1]-f(ipix)[0]))
	h.plot(pix[i1],vert1[i1]*3)
	f = interpolate.interp1d(pix[i1],vert1[i1]*3)
	h.plot(ipix,f(ipix))
	h.text(0,-0.025,'%d' % (f(ipix)[1]-f(ipix)[0]))
	
fig.savefig("/Users/dhk/work/data/NGC_IC/pilot4/psf_bg_xgrad.pdf")




hgpsf.close()
hipsf.close()

plt.show()
