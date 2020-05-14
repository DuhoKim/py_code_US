from astropy.io import ascii
from astropy.table import Table, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import os
from shutil import copyfile
from pandas import DataFrame, read_csv
import pandas as pd

ned1=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_1.txt")
ned2=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_2.txt")
ned3=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_3.txt")
ned4=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_4.txt")
ned5=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_5.txt")

ned_tot=vstack([ned1,ned2,ned3,ned4,ned5])

s257=ascii.read("/Users/dhk/work/cat/NGC_IC/sha_quarry_batch_257_without_prefix.txt")

rc3=ascii.read("/Users/dhk/work/cat/NGC_IC/RC3/myrc3.dat")
btot=ascii.read("/Users/dhk/work/data/NGC_IC/SHA/galfit/9th_mask/btot_3rd.txt")
#nair=ascii.read("/Users/dhk/work/cat/NGC_IC/Nair2010/table2.txt")
file=r'/Users/dhk/work/cat/NGC_IC/Steinicke/NI2018.xls'
df = pd.read_excel(file)
file=r'/Users/dhk/work/cat/NGC_IC/ttype.xls'
ttype = pd.read_excel(file)

petroR50_V = np.load('/Users/dhk/work/py/pyraf/petroR50_V.npy')
petroR50_V = petroR50_V * 0.6	# pix to ["]

fh=open('/Users/dhk/Documents/SESE/Dissertation/appenC.tex','w')
figs=open('/Users/dhk/Documents/publish/ngcic/fig2.tex','w')

sha2fig=np.zeros(257)
i=0
for x in range(0,len(ned_tot)):
	if ned_tot['col2'][x] in s257['id']:
		boa = "%.2f" % (ned_tot['col12'][x] / ned_tot['col11'][x])
		c = SkyCoord(ned_tot['col4'][x], ned_tot['col5'][x])
	#	dec = ned_tot['col5'][x]
	#	dec = '$'+"{10.6f}".format(c.dec.value)
		
		name = ned_tot['col2'][x]
		if name[0]=='n':
			galnum = name[3:].strip()
			ngc_match = df.loc[(df['N']=='N') & (df['NI']==int(galnum))]
			if len(galnum) == 3:
				galnum='0'+galnum
			elif len(galnum) ==2:
				galnum='00'+galnum
			elif len(galnum) ==1:
				galnum='000'+galnum
			name = 'NGC'+galnum
			table_name = 'NGC '+galnum
		elif name[0]=='i':
			galnum = name[2:].strip()
			ngc_match = df.loc[(df['N']=='I') & (df['NI']==int(galnum))]
			if len(galnum) == 3:
				galnum='0'+galnum
			elif len(galnum)==2:
				galnum='00'+galnum
			elif len(galnum)==1:
				galnum='000'+galnum
			name = 'IC'+galnum
			table_name = 'IC '+galnum

		if len(ngc_match) > 1:
			print(name+'has more than 1 component')

		if len(ngc_match) == 0:
			print(name+'has more no match')

		####  RC3 catalogue #########
		rc3_match = rc3[[j for j,s in enumerate(rc3['name']) if s.strip() == name]]
		if len(rc3_match) != 1:
			print('rc3 match is not one')
		elif rc3_match['T'][0] == '*':
			rc3_match['T'][0] = str(ttype.iloc[i,4])+'*'

		####  Nair+2010 cagalogue ###
		# coord_ned = SkyCoord(ned_tot['col4'][x],ned_tot['col5'][x],frame='icrs')
		# coord_sdss = SkyCoord(nair['ra'],nair['dec'],unit='deg')
		# idx, d2d, d3d = coord_ned.match_to_catalog_sky(coord_sdss)
	
		i=i+1
		sha2fig[np.where(s257['id']==ned_tot['col2'][x])]=i
			
		##### TeX for Table #####
		# Figure number for last column of Table linking to Figures
		# if i % 2:
		# 	nfig=str(i)
		# else:
		# 	nfig=str(i-1)
		nfig='\hyperref[fig:app-d-'+str(i)+']{'+str(i)+'}'

		pa = str(ngc_match.iloc[0,21])
		if pa=='nan':
			pa='*'
		classification = ned_tot['col13'][x].strip()
		classification2 = classification.replace('^','\string^')              
		classification3 = classification2.replace('_','\_')
		hlr = "%4.1f" % petroR50_V[np.where(s257['id']==ned_tot['col2'][x])]
		b2t_idx = np.where(btot['col1']==ned_tot['col2'][x])
		b2t = "%4.2f" % btot['col2'][b2t_idx[0][0]]
		b2t_chi = "%4.1f" % btot['col3'][b2t_idx[0][0]]
		mtot = "--"+"%4.1f" % -btot['col4'][b2t_idx[0][0]]
		mtot_app = "%4.1f" % btot['col6'][b2t_idx[0][0]]
		if c.dec.value > 0:
			dec_str = "{:10.6f}".format(c.dec.value)
		else:
			dec_str = "-"+"{:9.6f}".format(c.dec.value)

		if len(rc3_match) == 1:
			if rc3_match['rBe'][0] == '*':
				rBe = rc3_match['rBe'][0]
			else:
				rBe = "%5.1f" % (float(rc3_match['rBe'][0])*60)
			try:
				if float(rc3_match['T'][0]) > 0:
					T_str = rc3_match['T'][0]
				else:
					T_str = "-" + rc3_match['T'][0]
			except:
				if float(rc3_match['T'][0][:-1]) < 0:
					T_str = "-" + rc3_match['T'][0]
				else:
					T_str = rc3_match['T'][0]
			if rc3_match['eT'][0][0] == '.':
				eT_str = '0'+rc3_match['eT'][0]
			else:
				eT_str = rc3_match['eT'][0]

			fh.writelines(	table_name+' & '+ \
					"{:10.6f}".format(c.ra.value)+' & '+ \
					dec_str+' & '+				\
					str(ned_tot['col8'][x])+' & '+  	\
					str(ngc_match.iloc[0,15])+' & '+	\
					str(ngc_match.iloc[0,16])+' & '+	\
					mtot_app+' & '+	\
					mtot+' & '+ \
					str(ned_tot['col11'][x])+' & '+		\
					str(boa)+' & '+				\
					rBe+' & '+		\
					hlr+' & '+	\
					pa+' & '+				\
					str(ned_tot['col6'][x])+' & '+		\
					b2t+' & '+ 	\
					#b2t_chi+' & '+	\
					T_str+' & '+		\
					eT_str+' & '+		\
					ngc_match.iloc[0,22]+' & '+		\
					str(ttype.iloc[i-1,6])+' & '+		\
					classification3+' & '+nfig+' \\\ \n')
		else:
			fh.writelines(	table_name+' & '+ \
					"{:10.6f}".format(c.ra.value)+' & '+ \
					dec_str+' & '+				\
					str(ned_tot['col8'][x])+' & '+  	\
	                str(ngc_match.iloc[0,15])+' & '+	\
					str(ngc_match.iloc[0,16])+' & '+ \
					mtot_app + ' & ' + \
					mtot + ' & ' + \
					str(ned_tot['col11'][x])+' & '+		\
					str(boa)+' & '+				\
					'* & '+					\
					hlr+' & '+  \
	                pa+' & '+				\
					str(ned_tot['col6'][x])+' & '+		\
                    b2t+' & '+  \
				    #b2t_chi+' & '+  \
					str(ttype.iloc[i-1,4])+'* & '+		\
					'* & '+					\
					ngc_match.iloc[0,22]+' & '+		\
                    str(ttype.iloc[i-1,6])+' & '+           \
					classification3+' & '+nfig+' \\\ \n')

		###### TeX for Figure #######
		# if i % 2:
		# 	prev_name = table_name
		# else:
		figs.write(r'\figsetgrpstart'+'\n')
		figs.write(r'\figsetgrpnum{2.'+str(int(i))+'} \n')
		figs.write(r'\figsetgrptitle{'+table_name+'} \n')
		figs.write(r'\figsetplot{figset2_'+str(int(i))+'.pdf} \n')
		figs.write(r'\figsetgrpnote{Same as Figure~2 for '+table_name+'} \n')
		figs.write(r'\figsetgrpend'+'\n\n')

figs.close()
fh.close()	
np.save('/Users/dhk/work/data/NGC_IC/pdf/sha2fig.npy',sha2fig)
