from astropy.io import ascii
from astropy.table import Table, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import os
from shutil import copyfile
from pandas import DataFrame, read_csv
import pandas as pd
from astropy.wcs import WCS

ned1=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_1.txt")
ned2=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_2.txt")
ned3=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_3.txt")
ned4=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_4.txt")
ned5=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_5.txt")

ned_tot=vstack([ned1,ned2,ned3,ned4,ned5])

sha_cat=ascii.read('/Users/dhk/work/cat/NGC_IC/sha_ned_match_off_lt_1_SDSS_100psfs_b_over_a_gt_05_from_Tom_2_erase_Tamura_match_segmentation.csv') # Sample NGC/IC numbers

rc3=ascii.read("/Users/dhk/work/cat/NGC_IC/RC3/myrc3.dat")
#nair=ascii.read("/Users/dhk/work/cat/NGC_IC/Nair2010/table2.txt")
file=r'/Users/dhk/work/cat/NGC_IC/Steinicke/NI2018.xls'
df = pd.read_excel(file)
file=r'/Users/dhk/work/cat/NGC_IC/ttype.xls'
ttype = pd.read_excel(file)

fm = open('/Users/dhk/Documents/SESE/Dissertation/DKfromTemplate/appendD.tex','w')
cnt = 0

fm.writelines('\chapter{Figure 3-2 for 257 NGC and IC galaxies} \n')
fm.writelines('\label{chap:app-d} \n\n')

for x in range(0,len(ned_tot)):

	if ned_tot['col2'][x] in sha_cat['col1']:
		cnt = cnt + 1
		name = ned_tot['col2'][x]

		sha_match = sha_cat[sha_cat['col1']==name]

		boa = "%.2f" % (ned_tot['col12'][x] / ned_tot['col11'][x])

				######  WCS  read      #######
		w       = WCS('/Users/dhk/work/data/NGC_IC/SHA/SHA_NGC_IC_LONG_mean/'+name+'.fits')
		i,j     = w.all_world2pix(sha_match['col2'][0],sha_match['col3'][0],1)

		semi = ned_tot['col11'][x]                      # read semi-major axis from NED [']
		thumb = int(semi*100)                           # half size of figure size

		if name[0]=='n':
				galnum = name[3:].strip()
				ngc_match = df.loc[(df['N']=='N') & (df['NI']==int(galnum))]
				if len(galnum) ==3:
					galnum='0'+galnum
				elif len(galnum) ==2:
					galnum='00'+galnum
				elif len(galnum) ==1:
					galnum='000'+galnum
				table_name = 'NGC'+galnum
		elif name[0]=='i':
				galnum = name[2:].strip()
				ngc_match = df.loc[(df['N']=='I') & (df['NI']==int(galnum))]
				if len(galnum) == 3:
					galnum = '0' + galnum
				elif len(galnum) == 2:
					galnum = '00' + galnum
				elif len(galnum) == 1:
					galnum = '000' + galnum
				table_name = 'IC' + galnum

		pa = str(ngc_match.iloc[0,21])
		if pa=='nan':
				pa='0'

		#####  xxx.feedme    ########


		fm.writelines('\\begin{figure} \n')
		fm.writelines('\center \n')
		fm.writelines('\includegraphics[width=1.1\\txw]{chapter3_figure/figset2_'+str(cnt)+'.pdf} \\\ \n')
		fm.writelines('\includegraphics[width=1.05\\txw]{galfit/'+name+'.pdf} \n')
		fm.writelines('\caption[Same as Figure~3.2 for '+table_name+']{\small Same as Figure~\\ref{fig:3-2} for '+table_name+'} \n')
		fm.writelines('\label{fig:app-d-'+str(cnt)+'} \n')	# Bad pixel mask
		fm.writelines('\end{figure} \n\n')
		fm.writelines('\clearpage \n\n')

fm.close()
