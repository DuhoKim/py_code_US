#####################################################
# Python script that reads betaV values and T-types
# of 257 NGC/IC galaxies and plot
# written by Duho Kim (06/21/18)        
######################################################
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
from matplotlib import colors
from shutil import copyfile
from pandas import DataFrame, read_csv
import pandas as pd
from operator import itemgetter
import os.path

file=r'/Users/dhk/work/cat/NGC_IC/ttype.xls'
ttype = pd.read_excel(file)

bv0=[	[1.9409089,1.9044532,1.3527486,1.1198820,0.82968015,0.58503551],	# BetaVzero values for SFH3,4,5 and Z=.0001,.0004,.004,.008,.02,.05
	[1.9860456,1.9576220,1.4390702,1.2023316,0.91737698,0.65453906],
	[2.3880801,2.3478914,1.6838646,1.4124115,1.1048444,0.77272439]]

sfh_ttype=[-5,0,5,9]
cols_for_Z=['gray','violet','blue','green','darkorange','red']
text_for_Z=[r'Z=0.0001','Z=0.0004','Z=0.004','Z=0.008','Z=0.02 (Z$_{\odot}$)','Z=0.05']

agn=['Not active','Starburst','LIRG','LINER','Sy 2','Sy 1']
cm=['black','blue','red','darkorange','darkgreen','darkturquoise']

bt=np.load('/Users/dhk/work/py/pyraf/betav_ttype.npy')
bt2=np.load('/Users/dhk/work/py/pyraf/betav_range.npy')
btot=ascii.read("/Users/dhk/work/data/NGC_IC/SHA/galfit/9th_mask/btot_3rd.txt")
sha2fig=np.load('/Users/dhk/work/data/NGC_IC/pdf/sha2fig.npy')
sha_cat=ascii.read('/Users/dhk/work/cat/NGC_IC/sha_ned_match_off_lt_1_SDSS_100psfs_b_over_a_gt_05_from_Tom_2_erase_Tamura_match_segmentation.csv') # Sample NGC/IC numbers

########################################  BetaV vs. T-type #######################################

fig = plt.figure(figsize=(5,5))

a = plt.subplot(111)
cols = ["" for x in range(257)]

for i in range(257):
	#cols[i] = cm[ttype.iloc[int(sha2fig[i]-1),6]]
	if bt[2][i] == 90:
		tt = 12
	elif bt[2][i] == 99:
		tt = 14
	else:
		tt = bt[2][i]
 	e = Ellipse((tt, (bt2[1][i]+bt2[2][i])/2), bt[3][i], (bt2[1][i]-bt2[2][i])/2)
	e.set_alpha(0.2)
	a.add_artist(e)
	e.set_facecolor(cm[ttype.iloc[int(sha2fig[i]-1),6]])

avg_arr = np.zeros(18)
std_arr = np.zeros(18)

for i in range(-5, 13):
	if i < 11:
		idx = np.where( (bt[2] > i-0.5) & (bt[2] < i+0.5) )[0]
	elif i == 11:
		idx = np.where( bt[2]==90 )[0]
	elif i == 12:
		idx = np.where( bt[2]==99 )[0]

	print len(idx)
	arr_1000 = np.zeros(1000)
	for j in range(1000):
		arr_idx = np.zeros(len(idx))
		for k in range(len(idx)):
			arr_idx[k] = np.random.normal((bt2[1][idx[k]]+bt2[2][idx[k]])/2, (bt2[1][idx[k]]-bt2[2][idx[k]])/2)
		arr_1000[j] = np.mean(arr_idx)
	avg_arr[i+5] = arr_1000.mean()
	std_arr[i+5] = arr_1000.std()

print 'E'
idx = np.where(bt[2] < -3)
print np.mean((bt2[1][idx]-bt2[2][idx])/2)
print np.std(bt2[0][idx])
err = np.sqrt(np.mean((bt2[1][idx]-bt2[2][idx])/2)**2 + np.std(bt2[0][idx])**2)
print (2.5 * np.log10(np.mean(bt2[0][idx]) + err) - 2.5 * np.log10(np.mean(bt2[0][idx]) - err)) / 2
print 2.5 * err / np.log(10) / np.mean(bt2[0][idx])
print len(idx[0])
print 'S0'
idx = np.where((bt[2] < 0) & (bt[2] >=-3))
print np.mean((bt2[1][idx]-bt2[2][idx])/2)
print np.std(bt2[0][idx])
err = np.sqrt(np.mean((bt2[1][idx]-bt2[2][idx])/2)**2 + np.std(bt2[0][idx])**2)
print (2.5 * np.log10(np.mean(bt2[0][idx]) + err) - 2.5 * np.log10(np.mean(bt2[0][idx]) - err)) / 2
print 2.5 * err / np.log(10) / np.mean(bt2[0][idx])
print len(idx[0])
print 'Sa'
idx = np.where((bt[2] < 2) & (bt[2] >=0))
print np.mean((bt2[1][idx]-bt2[2][idx])/2)
print np.std(bt2[0][idx])
err = np.sqrt(np.mean((bt2[1][idx]-bt2[2][idx])/2)**2 + np.std(bt2[0][idx])**2)
print (2.5 * np.log10(np.mean(bt2[0][idx]) + err) - 2.5 * np.log10(np.mean(bt2[0][idx]) - err)) / 2
print 2.5 * err / np.log(10) / np.mean(bt2[0][idx])
print len(idx[0])
print 'Sb'
idx = np.where((bt[2] < 4) & (bt[2] >=2))
print np.mean((bt2[1][idx]-bt2[2][idx])/2)
print np.std(bt2[0][idx])
err = np.sqrt(np.mean((bt2[1][idx]-bt2[2][idx])/2)**2 + np.std(bt2[0][idx])**2)
print (2.5 * np.log10(np.mean(bt2[0][idx]) + err) - 2.5 * np.log10(np.mean(bt2[0][idx]) - err)) / 2
print 2.5 * err / np.log(10) / np.mean(bt2[0][idx])
print len(idx[0])
print 'Sbc'
idx = np.where((bt[2] < 5) & (bt[2] >=4))
print np.mean((bt2[1][idx]-bt2[2][idx])/2)
print np.std(bt2[0][idx])
err = np.sqrt(np.mean((bt2[1][idx]-bt2[2][idx])/2)**2 + np.std(bt2[0][idx])**2)
print (2.5 * np.log10(np.mean(bt2[0][idx]) + err) - 2.5 * np.log10(np.mean(bt2[0][idx]) - err)) / 2
print 2.5 * err / np.log(10) / np.mean(bt2[0][idx])
print len(idx[0])
print 'Sc'
idx = np.where((bt[2] < 7) & (bt[2] >=5))
print np.mean((bt2[1][idx]-bt2[2][idx])/2)
print np.std(bt2[0][idx])
err = np.sqrt(np.mean((bt2[1][idx]-bt2[2][idx])/2)**2 + np.std(bt2[0][idx])**2)
print (2.5 * np.log10(np.mean(bt2[0][idx]) + err) - 2.5 * np.log10(np.mean(bt2[0][idx]) - err)) / 2
print 2.5 * err / np.log(10) / np.mean(bt2[0][idx])
print len(idx[0])
print 'Sd'
idx = np.where((bt[2] < 9) & (bt[2] >=7))
print np.mean((bt2[1][idx]-bt2[2][idx])/2)
print np.std(bt2[0][idx])
err = np.sqrt(np.mean((bt2[1][idx]-bt2[2][idx])/2)**2 + np.std(bt2[0][idx])**2)
print (2.5 * np.log10(np.mean(bt2[0][idx]) + err) - 2.5 * np.log10(np.mean(bt2[0][idx]) - err)) / 2
print 2.5 * err / np.log(10) / np.mean(bt2[0][idx])
print len(idx[0])
print 'irr&pec'
idx = np.where(bt[2] >=9)
print np.mean((bt2[1][idx]-bt2[2][idx])/2)
print np.std(bt2[0][idx])
err = np.sqrt(np.mean((bt2[1][idx]-bt2[2][idx])/2)**2 + np.std(bt2[0][idx])**2)
print (2.5 * np.log10(np.mean(bt2[0][idx]) + err) - 2.5 * np.log10(np.mean(bt2[0][idx]) - err)) / 2
print 2.5 * err / np.log(10) / np.mean(bt2[0][idx])
print len(idx[0])



# for i in range(3):		# for SFH3,4,5
# 	for j in range(6):	# for Z=.0001,.0004,.004,.008,.02,.05
# 		plt.plot([sfh_ttype[i],sfh_ttype[i+1]],[bv0[i][j],bv0[i][j]],'-',color=cols_for_Z[j],alpha=0.5)

plt.xlim(-7, 15)
plt.ylim(0.3, 3.0)
plt.xlabel('T-type',size=15)
plt.ylabel(r'$\beta_{V}$',size=15)
custom_ellipse = [	Line2D([0],[0],marker='o',color='w',markerfacecolor=cm[0],alpha=0.3,markersize=15),
			Line2D([0],[0],marker='o',color='w',markerfacecolor=cm[1],alpha=0.3,markersize=15),
			Line2D([0],[0],marker='o',color='w',markerfacecolor=cm[2],alpha=0.3,markersize=15),
			Line2D([0],[0],marker='o',color='w',markerfacecolor=cm[3],alpha=0.3,markersize=15),
			Line2D([0],[0],marker='o',color='w',markerfacecolor=cm[4],alpha=0.3,markersize=15),
			Line2D([0],[0],marker='o',color='w',markerfacecolor=cm[5],alpha=0.3,markersize=15)]
#ells = [Ellipse((-5,2.9),0.3,0.3,facecolor=cm[i]) for i in range(len(agn))]
#black = mpatches.Ellipse(color='black',label=agn[0])
legend1=plt.legend(custom_ellipse,agn,loc=2,bbox_to_anchor=(0.01,0.998), frameon=False, fontsize=12)
custom_lines = [ 	Line2D([0],[0],color=cols_for_Z[0],alpha=0.5),
			Line2D([0],[0],color=cols_for_Z[1],alpha=0.5),
			Line2D([0],[0],color=cols_for_Z[2],alpha=0.5),
			Line2D([0],[0],color=cols_for_Z[3],alpha=0.5),
			Line2D([0],[0],color=cols_for_Z[4],alpha=0.5),
			Line2D([0],[0],color=cols_for_Z[5],alpha=0.5)	]

#plt.scatter(bt[2], bt2[0], color=np.array(cols), alpha=0.2, edgecolor='None')
#plt.errorbar(bt[2], bt2[0], xerr=bt[3]/2, yerr=[bt2[0]-bt2[2], bt2[1]-bt2[0]], fmt='None', alpha=0.2,
#			ecolor=np.array(cols))

#plt.errorbar(3, 2.2, xerr=np.mean(bt[3])/2, yerr=[np.mean(bt2[0]-bt2[2]), np.mean(bt2[1]-bt2[0])])
#plt.legend(custom_lines,text_for_Z,loc=2,fontsize='small')

plt.errorbar(range(-5, 11)+[12, 14], avg_arr, yerr=std_arr)
plt.gca().add_artist(legend1)
plt.tick_params(direction='in',labelsize='large',top=True,right=True)
plt.xticks([-5, 0, 5, 10, 12, 14], [-5, 0, 5, 10, 90, 99])
plt.yticks([0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
fig.savefig('/Users/dhk/Documents/publish/ngcic/betav_ttype.pdf')

########################################  BetaV vs. B/T #######################################

fig = plt.figure(figsize=(5,5))

a = plt.subplot(111)
btot_arr = np.zeros(257)
for i in range(257):
	b2t_idx = np.where(btot['col1'] == sha_cat['col1'][i])
	btot_arr[i] =  btot['col2'][b2t_idx[0][0]]
	#e = Ellipse((btot['col2'][b2t_idx[0][0]],bt[0][i]),btot['col3'][b2t_idx[0][0]],bt[1][i])
	#e.set_alpha(0.2)
	#a.add_artist(e)
plt.scatter(btot_arr,bt[0])
plt.xlim(1,0)
plt.ylim(0,3.5)
plt.xlabel('B/T',size=15)
plt.ylabel(r'$\beta_{V}$',size=15)
plt.tick_params(direction='in',labelsize='large',top=True,right=True)
#plt.xticks([-5,0,5,10])
fig.savefig('/Users/dhk/Documents/publish/ngcic/betav_btot.pdf')


########################################  BetaV_center vs. T-type #######################################

fig2=plt.figure(figsize=(5,5))

fh=open('/Users/dhk/work/data/NGC_IC/cat/betaVcent_SFAGN.txt','w')

b = plt.subplot(111)

type_agn=range(257)

fh.writelines('############################################################\n')
fh.writelines('#  a catalog of 257 NGC/IC galaxies containing             #\n')
fh.writelines('#  1: name,                                                #\n')
fh.writelines('#  2: mean of total pixel-by-pixel values of               #\n')
fh.writelines('#     V/3.6mu flux ratios (betaV),                         #\n')
fh.writelines('#  3: standard deviation of 2                              #\n')
fh.writelines('#  4: mean of central 3x3 pixels values of betaV           #\n')
fh.writelines('#  5: standard deviation of 4                              #\n')
fh.writelines('#  6: SF/AGN types--0: Not Active, 1: Starburst, 2: LIRG,  #\n')
fh.writelines('#     3: LINER, 4: Seyfert2, 5: Seyfert1                   #\n')
fh.writelines('############################################################\n')

for i in range(257):
	e = Ellipse((bt[2][i],bt[4][i]),bt[3][i],bt[5][i])
	e.set_alpha(0.2)
	b.add_artist(e)
	type_agn[i] = ttype.iloc[int(sha2fig[i]-1),6]
	e.set_facecolor(cm[type_agn[i]])
	name = ttype.iloc[int(sha2fig[i]-1),0]
	fh.writelines(name.replace(' ', '') + ' ' +
                    "%5.3f" % (bt[0][i]) + ' ' +
                    "%5.3f" % (bt[1][i]) + ' ' +
                    "%5.3f" % (bt[4][i]) + ' ' +
                    "%5.3f" % (bt[5][i]) + ' ' +
                    str(type_agn[i])+ ' \n')

fh.close()

plt.xlim(-7,12)
plt.ylim(0,3)
plt.xlabel('T-type')
plt.ylabel(r'$\beta_{V}$ in galactic central 3$\times$3 pixels')
custom_ellipse = [      Line2D([0],[0],marker='o',color='w',markerfacecolor=cm[0],alpha=0.2,markersize=15),
       	                Line2D([0],[0],marker='o',color='w',markerfacecolor=cm[1],alpha=0.2,markersize=15),
               	        Line2D([0],[0],marker='o',color='w',markerfacecolor=cm[2],alpha=0.2,markersize=15),
                       	Line2D([0],[0],marker='o',color='w',markerfacecolor=cm[3],alpha=0.2,markersize=15),
                        Line2D([0],[0],marker='o',color='w',markerfacecolor=cm[4],alpha=0.2,markersize=15),
                        Line2D([0],[0],marker='o',color='w',markerfacecolor=cm[5],alpha=0.2,markersize=15)]
plt.legend(custom_ellipse,agn,loc=2)
plt.tick_params(direction='in')
plt.xticks([-5,0,5,10])

fig2.savefig('/Users/dhk/Documents/publish/ngcic/figs/betav_ttype_center.pdf')

########################################  BetaV vs. SF/AGN #######################################

fig3=plt.figure(figsize=(6,8))

yticks=[50,3,1,10,3,2]
for i in range(6):
	###### change the order ########
	# agn=['Not active','Starburst','LIRG','LINER','Sy 2','Sy 1']
	#          0             1         2      3       4      5
	if i == 1:
		ax = plt.axes([0.1, 0.725 - (i + 1) * 0.125, 0.8, 0.125])
	elif i == 2:
		ax = plt.axes([0.1, 0.725 - (i + 3) * 0.125, 0.8, 0.125])
	elif i == 3:
		ax = plt.axes([0.1, 0.725 - (i - 2) * 0.125, 0.8, 0.125])
	elif i == 4 or i == 5:
		ax = plt.axes([0.1, 0.725 - (i - 1) * 0.125, 0.8, 0.125])
	else:
		ax = plt.axes([0.1, 0.725 - i * 0.125, 0.8, 0.125])
	ax.xaxis.set_visible(False)
	ax.set_xlim(0,2)
	ax.set_yticks([0, yticks[i]])
	ax.invert_xaxis()
	ax.tick_params(direction='in',labelsize=15,top=True,right=True)

	agn_arr=np.asarray(type_agn)

	bvs=bt[4][np.where(agn_arr==i)]

	medbv=np.median(bvs[~np.isnan(bvs)])
	ax.hist(bvs[~np.isnan(bvs)],histtype='step',color=cm[i]) # not active
	ax.text(0.02,0.75,agn[i]+' ('+str(len(bvs[~np.isnan(bvs)]))+')',color=cm[i],transform=ax.transAxes,size=15)
	ax.text(0.02,0.45,r'med($\beta_{{V}}$)={0:.2f}'.format(medbv),color=cm[i],transform=ax.transAxes,size=15)
	ax.axvline(x=medbv,color=cm[i],linestyle='dashed')

	if i == 0:
		ax.xaxis.set_visible(True)
		ax.tick_params(direction='in',labelsize=15,top=True,right=True,bottom=False)
		ax.set_xticklabels([])

	if i == 4:
		ax.set_ylabel('N',position=(0,0.9),size=15)

	if i == 2:
		ax.tick_params(direction='in',labelsize=15,right=True,top=False)
		ax.set_xlabel(r'$\beta_{V}$ in central 3x3 pixels',size=15)
		ax.xaxis.set_visible(True)
		ax.set_ylim([0, 2])



fig3.savefig('/Users/dhk/Documents/publish/ngcic/betav_center_agn.pdf')