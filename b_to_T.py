from astropy.io import ascii
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

ned1=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_1.txt")
ned2=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_2.txt")
ned3=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_3.txt")
ned4=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_4.txt")
ned5=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_5.txt")

ned_tot=vstack([ned1,ned2,ned3,ned4,ned5])

s257=ascii.read("/Users/dhk/work/cat/NGC_IC/sha_quarry_batch_257_without_prefix.txt")

rc3=ascii.read("/Users/dhk/work/cat/NGC_IC/RC3/myrc3.dat")
file=r'/Users/dhk/work/cat/NGC_IC/ttype.xls'
ttype = pd.read_excel(file)
file=r'/Users/dhk/work/cat/NGC_IC/SDSS/sdss_query_result.xls'
sdss = pd.read_excel(file)
file=r'/Users/dhk/work/cat/NGC_IC/Steinicke/NI2018.xls'
df = pd.read_excel(file)


z=np.zeros(len(s257))
V=np.zeros(len(s257))
T=np.zeros(len(s257))

reff_ell=np.zeros(len(s257))
reff_sdss=np.zeros(len(s257))

hubble=np.zeros(9)	# E S0 Sa Sb Sbc Sc Sd Irr Pec
hubble_label=('E', 'S0', 'Sa' ,'Sb', 'Sbc', 'Sc', 'Sd', 'Irr', 'Pec')

def bulge (r, I0, R0):
	return I0 * exp(-(r/R0)**0.25)

def disk (r, I0, R0):
	return I0 * exp(-r/R0)

def func (r,I0_bulge,R0_bulge,I0_disk,R0_disk):
	return bulge(r,I0_bulge,R0_bulge) + disk(r,I0_disk,R0_disk) 

btot=open('/Users/dhk/work/cat/NGC_IC/btot.txt','w')
hlr=open('/Users/dhk/work/cat/NGC_IC/hlr.txt','w')

i=0

for x in range(0,len(ned_tot)):
	if ned_tot['col2'][x] in s257['id']:
		name = ned_tot['col2'][x]
		if name[0]=='n':
			galnum = name[3:].strip()
			ngc_match = df.loc[(df['N']=='N') & (df['NI']==int(galnum))]
		elif name[0]=='i':
			galnum = name[2:].strip()
			ngc_match = df.loc[(df['N']=='I') & (df['NI']==int(galnum))]

		z[i]=ned_tot['col8'][x]
		T[i]=ttype.iloc[i,5]
		i=i+1

		########## read Reff from ELLIPSE ##########
		ell = ascii.read(work_dir+'ellipse/pa_ellip/'+name+'_I1.txt')
                min_idx=0
                hlr_idx=0

		if ell['col2'][0] != 0:
			min_idx = min(enumerate(ell['col2']),key=itemgetter(1))[0]                                      # find the index where INTENS is minimum
			hlr_idx = min(enumerate(np.abs(ell['col11'][min_idx]/2.0-ell['col11'])),key=itemgetter(1))[0]   # find half-light radius
		else:
			half_light = max(ell['col11'])/2.0
			hlr_idx = min(enumerate(np.abs(ell['col11']-half_light)),key=itemgetter(1))[0]   # find half-light radius

                hlr_pix = ell['col1'][hlr_idx] * np.sqrt(1-ell['col5'][hlr_idx])					# sqrt(ab) = a * sqrt(1-E); E = 1-b/a
		reff_ell[i-1]=hlr_pix*0.6	# pix -> arcsec
		
		hlr.writelines(name+' '+str(reff_ell[i-1])+' \n')

		if ell['col2'][0] == 0:		# skip zero flux at the center for now
			reff_ell[i-1]=0

		try:
			popt, pcov = curve_fit(func, ell['col1'], ell['col2'],sigma=ell['col4'],p0=[ell['col2'][0],hlr_pix/10.0,ell['col2'][0],hlr_pix])
		except:
			print "curve_fit error at "+name
			continue

		chi = 0
		for y in range(0,len(ell['col1'])):
			chi = chi + (ell['col2'][y]-func(ell['col1'][y],*popt))**2/ell['col4'][y]**2

		bulge_int = integrate.quad(bulge,ell['col1'][0],ell['col1'][-1],args=(popt[0],popt[1]))		
		disk_int = integrate.quad(disk,ell['col1'][0],ell['col1'][-1],args=(popt[2],popt[3]))		

		plt.figure(figsize=(5,5))
		plt.plot(ell['col1'],ell['col2'],'k-',label='data')
		try:
			plt.plot(ell['col1'],func(ell['col1'], *popt),'g-',label='bulge + disk, B/T=%4.1f' % tuple([bulge_int[0]/(bulge_int[0]+disk_int[0])]))
		except:
			print 'plot error at '+name
			continue
		plt.plot(ell['col1'],bulge(ell['col1'], popt[0], popt[1]),'r--',label='bulge: %4.1f*exp(-(r/%4.1f)^0.25)' % tuple(popt[0:2]))
		plt.plot(ell['col1'],disk(ell['col1'], popt[2], popt[3]),'b--',label='disk: %4.1f*exp(-(r/%4.1f))' % tuple(popt[2:4]))
		plt.fill_between(ell['col1'],ell['col2']-ell['col3'],ell['col2']+ell['col4'],color='lightgray')
		plt.ylabel('Intensity')
		plt.yscale('log')
		plt.ylim([max(ell['col2'])/1e4,max(ell['col2'])*1e1])
		plt.xlabel('radius')
		plt.tick_params(direction='in')
		plt.legend()
		plt.text(0.8,0.7,r'$\chi^{2}$='+'%4.1f' % (chi/len(ell['col1'])),horizontalalignment='center',verticalalignment='center',transform=plt.gca().transAxes)
		plt.savefig('/Users/dhk/work/data/NGC_IC/ellipse/pa_ellip/'+name+'.png')

		btot.writelines(name+' '+str(bulge_int[0]/(bulge_int[0]+disk_int[0]))+' '+str(chi/len(ell['col1']))+'\n')

btot.close()
hlr.close()
