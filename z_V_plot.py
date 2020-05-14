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

ned1=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_1.txt")
ned2=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_2.txt")
ned3=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_3.txt")
ned4=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_4.txt")
ned5=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_5.txt")

ned_tot=vstack([ned1,ned2,ned3,ned4,ned5])

s257=ascii.read("/Users/dhk/work/cat/NGC_IC/sha_quarry_batch_257_without_prefix.txt")

rc3=ascii.read("/Users/dhk/work/cat/NGC_IC/RC3/myrc3.dat")
file=r'/Users/dhk/work/cat/NGC_IC/Steinicke/NI2018.xls'
df = pd.read_excel(file)
file=r'/Users/dhk/work/cat/NGC_IC/ttype.xls'
ttype = pd.read_excel(file)

z=np.zeros(len(s257))
V=np.zeros(len(s257))
T=np.zeros(len(s257))
agn=np.zeros(len(s257))
hubble=np.zeros(9)	# E S0 Sa Sb Sbc Sc Sd Irr Pec
agn_sum=np.zeros(6)	# AGN type
hubble_label=('E', 'S0', 'Sa' ,'Sb', 'Sbc', 'Sc', 'Sd', 'Irr', 'Pec')
agn_label=['Not active','Starburst','LIRG','LINER','Sy 2','Sy 1']

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
		V[i]=ngc_match.iloc[0,16]
		T[i]=ttype.iloc[i,5]
		agn[i]=ttype.iloc[i,6]
		i=i+1
hubble[0] = sum(1 for i in T if i < -3)			# E
hubble[1] = sum(1 for i in T if i >= -3 and i < 0)	# S0
hubble[2] = sum(1 for i in T if i >= 0 and i < 2)	# Sa
hubble[3] = sum(1 for i in T if i >= 2 and i < 4)	# Sb
hubble[4] = sum(1 for i in T if i >= 4 and i < 5)	# Sbc
hubble[5] = sum(1 for i in T if i >= 5 and i < 7)	# Sc
hubble[6] = sum(1 for i in T if i >= 7 and i < 9)	# Sd
hubble[7] = sum(1 for i in T if i >= 9 and i < 91)	# Irr
hubble[8] = sum(1 for i in T if i >= 91 and i < 100)	# Pec
for x in range(0,6):
	agn_sum[x] = sum(agn==x)

fig, axs = plt.subplots(2,2,tight_layout=True,figsize=(6,6))
axs[0,0].hist(z,bins=20)
axs[0,0].set_xlabel('z',size=20)
axs[0,0].set_ylabel('N',size=20)
axs[0,0].tick_params(direction='in',top=True,right=True,labelsize='large')
axs[0,0].text(0.75, 0.8, '(a)', size=20, transform=axs[0,0].transAxes)
#axs[0,0].set_ylim(0,70)
axs[0,1].hist(V,bins=20)
#axs[0,1].set_ylim(0,40)
axs[0,1].set_xlabel('V (mag)',size=20)
axs[0,1].tick_params(direction='in',top=True,right=True,labelsize='large')
axs[0,1].text(0.75, 0.8, '(b)', size=20, transform=axs[0,1].transAxes)
axs[1,0].bar(range(len(hubble)),hubble,align='center')
plt.setp(axs[1,0],xticks=range(len(hubble)),xticklabels=hubble_label)
#axs[1,0].set_ylim(0,50)
axs[1,0].set_xlabel('Hubble type',size=20)
axs[1,0].set_ylabel('N',size=20)
axs[1,0].tick_params(direction='in',top=True,right=True)
axs[1,0].text(0.75, 0.8, '(c)', size=20, transform=axs[1,0].transAxes)
axs[1,1].bar(range(len(agn_sum)),agn_sum,align='center')
#axs[1,1].set_ylim(0,160)
plt.setp(axs[1,1],xticks=range(len(agn_label)),xticklabels='')
axs[1,1].set_xlabel('SF and AGN type',size=20)
#axs[1,1].set_ylabel('N',size=20)
axs[1,1].tick_params(direction='in',top=False,right=True,bottom=False)
axs[1,1].text(0.75, 0.8, '(d)', size=20, transform=axs[1,1].transAxes)
axs[1,1].text(-0.2, 90, 'Not active', size=15, rotation=90, color='white')
axs[1,1].text(0.7, 80, 'Starburst', size=15, rotation=90)
axs[1,1].text(1.7, 40, 'LIRG', size=15, rotation=90)
axs[1,1].text(2.7, 80, 'LINER', size=15, rotation=90)
axs[1,1].text(3.7, 55, 'Sy 2', size=15, rotation=90)
axs[1,1].text(4.7, 40, 'Sy 1', size=15, rotation=90)
fig.savefig('/Users/dhk/Documents/publish/ngcic/hist.pdf')

