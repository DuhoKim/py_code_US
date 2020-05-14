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
from astropy.table import Table, vstack
import scipy.stats as stats
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression

model = LinearRegression(fit_intercept=True)


def frange(x, y, jump):
	while x < y:
		yield x
		x += jump

def func(x, a, b):
	return a * x + b

def count_upper(idx, slope, yoff):
	upper = 0
	for idxes in idx:
		if bt[0][idxes] > func(boa[idxes], slope, yoff):
			upper += 1
	return upper

def count_lower(idx, slope, yoff):
	lower = 0
	for idxes in idx:
		if bt[0][idxes] < func(boa[idxes], slope, yoff):
			lower += 1
	return lower



ned1=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_1.txt")
ned2=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_2.txt")
ned3=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_3.txt")
ned4=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_4.txt")
ned5=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_5.txt")

ned_tot=vstack([ned1,ned2,ned3,ned4,ned5])

file=r'/Users/dhk/work/cat/NGC_IC/ttype.xls'
ttype = pd.read_excel(file)

bv0=[	[1.9409089,1.9044532,1.3527486,1.1198820,0.82968015,0.58503551],	# BetaVzero values for SFH3,4,5 and Z=.0001,.0004,.004,.008,.02,.05
	[1.9860456,1.9576220,1.4390702,1.2023316,0.91737698,0.65453906],
	[2.3880801,2.3478914,1.6838646,1.4124115,1.1048444,0.77272439]]

sfh_ttype=[-5,0,5,9]
cols_for_T=['red', 'darkorange', 'orange', 'olive', 'green', 'blue', 'violet', 'darkviolet']
text_for_T=['E', 'S0', 'Sa', 'Sb', 'Sbc', 'Sc', 'Sd', 'Irr&Pec']

agn=['Not active','Starburst','LIRG','LINER','Sy 2','Sy 1']
cm=['black','blue','red','darkorange','darkgreen','darkturquoise']

bt=np.load('/Users/dhk/work/py/pyraf/betav_range.npy')
Ts=np.load('/Users/dhk/work/py/pyraf/betav_ttype.npy')
btot=ascii.read("/Users/dhk/work/data/NGC_IC/SHA/galfit/9th_mask/btot_3rd.txt")
sha2fig=np.load('/Users/dhk/work/data/NGC_IC/pdf/sha2fig.npy')
sha_cat=ascii.read('/Users/dhk/work/cat/NGC_IC/sha_ned_match_off_lt_1_SDSS_100psfs_b_over_a_gt_05_from_Tom_2_erase_Tamura_match_segmentation.csv') # Sample NGC/IC numbers

########################################  BetaV vs. T-type #######################################

zs = np.zeros(257)
boa = np.zeros(257)

cols = ["" for x in range(257)]

for i in range(257):
	ned_idx = np.where(ned_tot['col2'] == sha_cat['col1'][i])
	zs[i] = ned_tot['col8'][ned_idx]
	zerr = ned_tot['col9'][ned_idx]
	boa[i] = ned_tot['col12'][ned_idx] / ned_tot['col11'][ned_idx]

	T = Ts[2][i]

	if T < -3:
		T_idx = 0
	if T < 0 and T >= -3:
		T_idx = 1
	if T < 2 and T >= 0:
		T_idx = 2
	if T < 4 and T >= 2:
		T_idx = 3
	if T < 5 and T >= 4:
		T_idx = 4
	if T < 7 and T >= 5:
		T_idx = 5
	if T < 9 and T >= 7:
		T_idx = 6
	if T >= 9:
		T_idx = 7

#	e = Ellipse((z,bt[0][i]),0.001,bt[1][i])
#	e.set_alpha(0.2)
#	a.add_artist(e)
#	e.set_facecolor(cols_for_T[T_idx])

	cols[i] = cols_for_T[T_idx]



fig = plt.figure(figsize=(10,12))

axes_grid = [[0.1,0.65,0.275,0.275],
			 [0.375, 0.65, 0.275, 0.275],
			 [0.65, 0.65, 0.275, 0.275],
			 [0.1, 0.375, 0.275, 0.275],
			 [0.375, 0.375, 0.275, 0.275],
			 [0.65, 0.375, 0.275, 0.275],
			 [0.1, 0.1, 0.275, 0.275],
			 [0.375, 0.1, 0.275, 0.275]]

# h1 = fig.add_axes([0.1,0.65,0.275,0.275])
# h2 = fig.add_axes([0.375,0.65,0.275,0.275])
# h3 = fig.add_axes([0.65,0.65,0.275,0.275])
# h4 = fig.add_axes([0.1,0.375,0.275,0.275])
# h5 = fig.add_axes([0.375,0.375,0.275,0.275])
# h6 = fig.add_axes([0.65,0.375,0.275,0.275])
# h7 = fig.add_axes([0.1,0.1,0.275,0.275])

custom_ellipse = [	Line2D([0],[0],marker='o',color='w',markerfacecolor=cols_for_T[0],alpha=0.2,markersize=10),
			Line2D([0],[0],marker='o',color='w',markerfacecolor=cols_for_T[1],alpha=0.2,markersize=10),
			Line2D([0],[0],marker='o',color='w',markerfacecolor=cols_for_T[2],alpha=0.2,markersize=10),
			Line2D([0],[0],marker='o',color='w',markerfacecolor=cols_for_T[3],alpha=0.2,markersize=10),
			Line2D([0],[0],marker='o',color='w',markerfacecolor=cols_for_T[4],alpha=0.2,markersize=10),
			Line2D([0],[0],marker='o',color='w',markerfacecolor=cols_for_T[5],alpha=0.2,markersize=10),
			Line2D([0],[0],marker='o',color='w',markerfacecolor=cols_for_T[6],alpha=0.2,markersize=10),
			Line2D([0],[0],marker='o',color='w',markerfacecolor=cols_for_T[7],alpha=0.2,markersize=10)]
#ells = [Ellipse((-5,2.9),0.3,0.3,facecolor=cm[i]) for i in range(len(agn))]
#black = mpatches.Ellipse(color='black',label=agn[0])

#legend1=plt.legend(custom_ellipse,text_for_T,loc=2,bbox_to_anchor=(0.65,0.998), frameon=False, fontsize=12)
cols = np.array(cols)

for i in range(8):
	idx = np.where(cols == cols_for_T[i])
	ax = fig.add_axes(axes_grid[i])
	ax.errorbar(boa[idx], bt[0][idx], yerr=[bt[0][idx]-bt[2][idx], bt[1][idx]-bt[0][idx]], fmt='o', alpha=0.5)
	ax.set_xlim(0.45, 1.05)
	ax.set_ylim(0.0, 3.0)
	ax.tick_params(direction='in', labelsize='20', top=True, right=True, bottom=True, left=True)
	ax.tick_params(which='minor', direction='in', right=True, top=True)
	ax.set_xticks(np.arange(0.5, 1.0, 0.1), minor=True)
	ax.set_xticks([0.5, 1.0])
	ax.set_xticklabels([])
	ax.set_yticks(np.arange(0.0, 3.0, 0.1), minor=True)
	ax.set_yticks([0.0, 1.0, 2.0, 3.0])

	coeff, pvalue = stats.pearsonr(boa[idx],bt[0][idx])

#	ax.text(0.5, 0.1, r'P$_{coeff}$=%.2f' % coeff, size=15)
#	ax.text(0.8, 0.1, r'P$_{pval}$=%.2f' % pvalue, size=15)

	# assign geometric mean values of upper and lower errors as sigma values
	#popt, pcov = curve_fit(func, boa[idx], bt[0][idx], sigma=np.sqrt((bt[1][idx]-bt[0][idx]) * (bt[0][idx]-bt[2][idx])))

	popt, pcov = curve_fit(func, boa[idx], bt[0][idx])
	popt_move, pcov_move = curve_fit(func, boa[idx] - 0.75 , bt[0][idx] - 1.0)

	# model.fit(boa[idx], bt[0][idx])
	#
	# print("Coefficients: \n", regr.coef_)

	# slope, intercept, r_value, p_value, std_err = stats.linregress(boa[idx], bt[0][idx])
	#
	# print("slope: %f intercept: %f r_value: %f p_value: %f std_err: %f" % (slope, intercept, r_value, p_value, std_err))

	perr = np.sqrt(np.diag(pcov_move))

	xdata = np.linspace(0.5, 1.0, 100)

	ax.plot(xdata, func(xdata, *popt))

	sixteen_perc = int( len( idx[0] ) * 0.16 )	# calculate 16% of the total number of galaxy

	if i == 3 or i == 5:
		perr[0] = perr[0] * 1.3			# increase 10% error for the intercept to widen error wedge
		perr[1] = perr[1] * 2  # increase 10% error for the intercept to widen error wedge

	if i == 4 or i == 6:
		perr[0] = perr[0] * 0.6  # decrease 10% error for the intercept to narrow error wedge
		perr[1] = perr[1] * 0.001			# decrease 10% error for the intercept to narrow error wedge

	for yoff in frange(0, 2, 0.1):
		upper_yoff = popt[1] + perr[1] + yoff
		if count_upper(idx[0], popt[0] - perr[0], upper_yoff) <= sixteen_perc:
			break

	for yoff in frange(0, 2, 0.1):
		lower_yoff = popt[1] - perr[1] - yoff
		if count_lower(idx[0], popt[0] + perr[0], lower_yoff) <= sixteen_perc:
			break

	ax.fill_between(xdata, func(xdata, popt[0], popt[1]) + np.sqrt((xdata*perr[0])**2 + perr[1]**2), func(xdata, popt[0], popt[1]) - np.sqrt((xdata*perr[0])**2 + perr[1]**2),
					alpha=0.15, color='blue')

	ax.text(0.55, 2.5, text_for_T[i], size=30)
	# ax.text(0.55, 2, 'perr[0]: %f' % perr[0])
	# ax.text(0.55, 1.5, 'perr[1]: %f' % perr[1])

	if i == 0:
		ax.set_yticks([0.0, 1.0, 2.0, 3.0])
		ax.set_yticklabels([0.0, 1.0, 2.0, 3.0])
	if i == 3:
		ax.set_ylabel(r'$\beta_{V}$', size=30)
		ax.set_yticks([0.0, 1.0, 2.0, 3.0])
		ax.set_yticklabels([0.0, 1.0, 2.0, ''])
	if i == 7:
		ax.set_xlabel('b/a', size=30)
	if i % 3  :
		ax.set_yticklabels([])
	if i == 6:
		ax.set_xticklabels([0.5, 1.0])
		ax.set_yticklabels([0.0, 1.0, 2.0, ''])
	if i == 5 or i == 7:
		#ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
		ax.set_xticklabels([0.5, 1.0])


fig.savefig('/Users/dhk/Documents/publish/ngcic_rev/betav_incl_fit_1sig_err_move_fudge.pdf')
