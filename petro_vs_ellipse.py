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
		ell = ascii.read(work_dir+'ellipse/pa_ellip/'+name+'.txt')
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
		if ell['col2'][0] == 0:		# skip zero flux at the center for now
			reff_ell[i-1]=0

		########## read Reff from SDSS ##############
		sdss_match = sdss.loc[sdss['name']==name]
		if len(sdss_match):
			g = sdss_match.iloc[0,10]
			r = sdss_match.iloc[0,11]
			sdss_V = g-0.59*(g-r)-0.01
			if np.abs(ngc_match.iloc[0,16]-sdss_V) < 1:
				reff_sdss[i-1] = (sdss_match.iloc[0,14]+sdss_match.iloc[0,16])/2.0
		
idx_E = np.where(np.logical_and.reduce((reff_ell > 0.01, reff_sdss > 0.01, T < 0)))
idx_S = np.where(np.logical_and.reduce((reff_ell > 0.01, reff_sdss > 0.01, T >= 0, T < 5)))
idx_L = np.where(np.logical_and.reduce((reff_ell > 0.01, reff_sdss > 0.01, T >= 5)))

plt.figure(figsize=(5,5))
plt.scatter(reff_ell[idx_E],reff_sdss[idx_E],s=3.0,color='red',alpha=0.5,label='E, S0')
plt.scatter(reff_ell[idx_S],reff_sdss[idx_S],s=3.0,color='green',alpha=0.5,label='Sa, Sb, Sbc')
plt.scatter(reff_ell[idx_L],reff_sdss[idx_L],s=3.0,color='blue',alpha=0.5,label='Sc, Sd, Irr, Pec')
plt.plot([0,45],[0,45],':',alpha=0.5,color='black')
plt.xlim(0,45)
plt.ylim(0,45)
plt.ylabel(r'$(r_{50,g}^{P} + r_{50,r}^{P})/2$ from SDSS ["]')
plt.xlabel(r'$r_{50,V}^{HL}$ from ELLIPSE fitting ["]')
#plt.gca().set_aspect('equal')
#plt.axis('scaled')
plt.legend()
plt.tick_params(direction='in')
plt.savefig('/Users/dhk/Documents/publish/ngcic/figs/reff_pa_ellip.pdf')
