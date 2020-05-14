from astropy.io import fits
import numpy as np
for x in range(10,11):
	hdulist=fits.open('/Users/dhk/work/data/M31/Mayall_V/MESSIER_031F%d-I-V-moh2006.fits' % (x))
	sci=hdulist[0].data
	sci[sci==50000]=np.nan
	hdulist[0].data=sci
	hdulist.writeto('/Users/dhk/work/data/M31/Mayall_V/%d.fits' % (x))
	hdulist.close()
			
