from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

data_dir = '/Users/dhk/work/data/NGC_IC/SHA/'

psf_hdu = fits.open(data_dir+'psf4_rotate/med.fits')
psf=psf_hdu[1].data

fig = plt.figure(figsize=(8,10))

h1=fig.add_axes([0.1,0.6,0.8,0.4])
h2=fig.add_axes([0.1,0.1,0.8,0.4])

#h1.invert_yaxis()
#h2.invert_yaxis()

#h1.set_xlim([0,80])
#h2.set_xlim([0,80])

h1.plot(range(0,61),psf[:,30])
h2.plot(range(0,61),psf[30,:])

#h1.text(40,-3,'SDSS g PSF + 0.012 [pixel value]')
#h2.text(40,-3,'IRAC CH 1 PSF + 0.0126 [pixel value]')

#h1.set_xlabel('radius [pixel]')
#h2.set_xlabel('radius [pixel]')
#h1.set_ylabel('arbitrary magnitude')
#h2.set_ylabel('arbitrary magnitude')

#h1.errorbar(tab1['col1'],tab1['col8'],yerr=[tab1['col10'],tab1['col9']],fmt='-o')
#h2.errorbar(tab2['col1'],tab2['col8'],yerr=[tab2['col10'],tab2['col9']],fmt='-o')

fig.savefig(data_dir+"psf4_rotate/med_prof.pdf")

