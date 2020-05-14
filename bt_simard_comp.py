from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import os
from shutil import copyfile
from pandas import DataFrame, read_csv
import pandas as pd
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from operator import itemgetter
from numpy import exp
from scipy.optimize import curve_fit
import scipy.integrate as integrate

work_dir= '/Users/dhk/work/data/NGC_IC/'


sha_cat=ascii.read("/Users/dhk/work/cat/NGC_IC/sha_quarry_batch_257_without_prefix.txt")

btot=ascii.read('/Users/dhk/work/cat/NGC_IC/btot.txt')

sim_table1=ascii.read('/Users/dhk/work/cat/NGC_IC/Simard2011/table1.dat')
#sim_table2=ascii.read('/Users/dhk/work/cat/NGC_IC/Simard2011/table2.dat')
#sim_table3=ascii.read('/Users/dhk/work/cat/NGC_IC/Simard2011/table3.dat')

hdu1=fits.open('/Users/dhk/work/cat/NGC_IC/SDSS/PhotoObjDR7_1.fits')
hdu2=fits.open('/Users/dhk/work/cat/NGC_IC/SDSS/PhotoObjDR7_2.fits')

sdss1=hdu1[1].data
sdss2=hdu2[1].data

names1=sdss1['name']
names2=sdss2['name']

type1=sdss1['type']
type2=sdss2['type']

dr71=sdss1['dr7objid']
dr72=sdss2['dr7objid']

dr81=sdss1['dr8objid']
dr82=sdss2['dr8objid']

simard=open('/Users/dhk/work/cat/NGC_IC/Simard2011/match.txt','w')

for x in range(0,len(sha_cat)):
#for x in range(0,1):
	name = sha_cat['id'][x]		# NGC/IC name
	if name in names1:
		ind = np.where(names1 == name)
		if type1[ind]!='GALAXY':
			break
		ind2 = np.where(sim_table1['col1']==dr71[ind])
		if len(ind2[0]):
			simard.writelines(name+' '+str(sim_table1['col15'][ind2][0])+' \n')
	if name in names2:
                ind = np.where(names2 == name)
                if type2[ind]!='GALAXY':
                        break
                ind2 = np.where(sim_table1['col1']==dr72[ind])
		if len(ind2[0]):
			simard.writelines(name+' '+str(sim_table1['col15'][ind2][0])+' \n')
simard.close()
