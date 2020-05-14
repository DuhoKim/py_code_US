#####################################################
# Python script of making a plot of comparison of output file generated during
# building PSF for given FITS images 
# written by Duho Kim (12/01/17)        
######################################################
from astropy.io import ascii
from astropy.table import Table
import os
import numpy as np
import matplotlib.pyplot as plt


data_dir = '/Users/dhk/work/data/NGC_IC/SDSS/'

out1=ascii.read(data_dir+'g_bench1/numstars_fit_func_norm_mag.txt',format='rst')
#out2=ascii.read(data_dir+'bench4/numstars_fit_func_norm_mag.txt',format='ascii.rst')

fig = plt.figure(figsize=(8,8))

h=fig.add_axes([0.1,0.6,0.4,0.4])
h2=fig.add_axes([0.6,0.6,0.4,0.4])
h3=fig.add_axes([0.1,0.1,0.4,0.4])
h4=fig.add_axes([0.6,0.1,0.4,0.4])

h.hist(out1['NUMSTARS'],bins=np.arange(5,20,1),histtype='bar')
h.set_xlabel('# of stars for generating PSFs')
h.set_ylabel('Frequency')

h2.hist(out1['FUNC3'][np.where(out1['FUNC3']!='')],histtype='bar')
h2.set_xlabel('Analytic Functions for fitting PSFs')
h2.set_ylabel('Frequency')

h3.hist(out1['NORM3'][~np.isnan(out1['NORM3'])],bins=np.arange(0,0.2,0.01),histtype='bar')
h3.set_xlabel('Norm scatter for fitting PSFs')
h3.set_ylabel('Frequency')

h4.hist(out1['MAG3'][~np.isnan(out1['MAG3'])],bins=np.arange(10,17,0.5),histtype='bar')
h4.set_xlabel('Mag of brightest star for fitting PSFs')
h4.set_ylabel('Frequency')


fig.savefig("/Users/dhk/work/data/NGC_IC/SDSS/psf_stars_param.pdf")

