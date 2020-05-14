import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from astropy.io import ascii
from matplotlib.ticker import NullFormatter

i1=ascii.read("/Users/dhk/work/data/NGC_IC/pilot/i1.cat")
#r=ascii.read("/Users/dhk/work/data/NGC_IC/pilot/r.cat")
r_r=ascii.read("/Users/dhk/work/data/NGC_IC/pilot/r_r.cat")
r_r_p=ascii.read("/Users/dhk/work/data/NGC_IC/pilot/r_r_p8.cat")

i12=ascii.read("/Users/dhk/work/data/NGC_IC/pilot2/i1.cat")
g_p_r=ascii.read("/Users/dhk/work/data/NGC_IC/pilot2/g.p.r.cat")
g_r=ascii.read("/Users/dhk/work/data/NGC_IC/pilot2/g.r.cat")
g_r_p=ascii.read("/Users/dhk/work/data/NGC_IC/pilot2/g.r.p.cat")

bins=[0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0]

plt.figure(1)

f, axarr = plt.subplots(4,2,sharex=True,sharey=True)


axarr[0,0].hist(i1['FWHM_IMAGE'],bins=bins)
axarr[1,0].hist(r_r['FWHM_IMAGE'],bins=bins)
axarr[2,0].hist(r_r_p['FWHM_IMAGE'],bins=bins)
axarr[3,0].set_xlabel('FWHM [pix(0.6\")]')
axarr[1,0].set_ylabel('N')
axarr[0,0].set_title('IC1613')
axarr[0,0].text(5,2000,'IRAC Ch1')
axarr[1,0].text(5,2000,'SDSS r')
axarr[2,0].text(5,2000,'SDSS r PSF')

axarr[0,1].hist(i12['FWHM_IMAGE'],bins=bins)
axarr[1,1].hist(g_r['FWHM_IMAGE'],bins=bins)
axarr[2,1].hist(g_r_p['FWHM_IMAGE'],bins=bins)
axarr[3,1].hist(g_p_r['FWHM_IMAGE'],bins=bins)
axarr[3,1].set_xlabel('FWHM [pix(0.6\")]')
axarr[1,1].set_ylabel('N')
axarr[0,1].set_title('NGC2403')
axarr[0,1].text(5,2000,'IRAC Ch1')
axarr[1,1].text(5,2000,'SDSS g reg')
axarr[2,1].text(4,2000,'SDSS g reg -> PSF')
axarr[3,1].text(4,2000,'SDSS g PSF -> reg')

plt.savefig('/Users/dhk/work/pyplot/fwhm1.png')

