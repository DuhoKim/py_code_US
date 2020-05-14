from numpy import *
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import Planck15, z_at_value
import astropy.units as u
import numpy as np
import matplotlib.pyplot as pl
import itertools

# AEGIS : All wavelength Extended Groth strip International Survey
# 	U 	B 	G 	V 	R 	Rp 	I 	Ip 	Z 	Zp 	Whitaker+11
# 	-0.22		-0.01		-0.01		0.03		0.05		AEGIS
#	-0.23	-0.07	-0.02	0.12	-0.03	0.03	0.00	-0.34	-0.01	0.11	COSMOS
#	ZP_EAZY = ZP_nominal + offset

#  	g (55,56,57)	r (52,53,54)	i (49,50,51)	z (46,47,48)	J2 (40,41,42)	J3 (37,38,39)	CFHT Megacam
#	4750		6400		7760		9250		11950		12790		NEWFIRM
#	4774		6231		7615		9132						SDSS
#	0.087		-0.181		-0.392		-0.526		0.84		0.96		AB to Vega 	(CFHT)
#	0.085		-0.167		-0.394		-0.530						(SDSS)

# V, z=0 5448, z=0.266 6897, z=0.614 8793, z=1.217 
# magAB = -2.5*log10(flux)+25.0	1.5 arcsec diameter apertures PSF- and astrometrically-matched images

V_ew=5448
opt_ew=[4750,6400,7760,9250,11950,12790]
irac_ew=[3.55,4.493,5.731,7.872]

zs = [i/3.55-1 for i in irac_ew]	# I1 z=0,	I2 z=0.266,	I3 z=0.614,	I4 z=1.217
Vs = [(i+1.0)*V_ew for i in zs]		# Vs = 5448,	6897,		8793,		12078
					# g r		r i		i z		J2 J3

x=np.genfromtxt('/Users/dhk/work/cat/NEWFIRM/AEGIS/aegis-n2.deblend.v5.1.cat')
sp=np.genfromtxt('/Users/dhk/work/cat/NEWFIRM/AEGIS/aegis-n2.deblend.sps/aegis-n2.bc03.del.deblend.v5.1.fout')
ssp=np.genfromtxt('/Users/dhk/work/betav5/data_ascii/beta_V_IRAC_ch1_SSP.txt')
tau1G=np.genfromtxt('/Users/dhk/work/betav5/data_ascii/beta_V_IRAC_ch1_exp1G.txt')
tau10G=np.genfromtxt('/Users/dhk/work/betav5/data_ascii/beta_V_IRAC_ch1_exp10G.txt')
sfh3=np.genfromtxt('/Users/dhk/work/betav4/ascii_sto/SFH3/Zevol/beta_V_IRAC_ch1_sto3_0.txt')
sfh4=np.genfromtxt('/Users/dhk/work/betav4/ascii_sto/SFH4/Zevol/beta_V_IRAC_ch1_sto4_0.txt')
sfh5=np.genfromtxt('/Users/dhk/work/betav4/ascii_sto/SFH5/Zevol/beta_V_IRAC_ch1_sto5_0.txt')

Zsol=min(ssp[:,0],key=lambda x:abs(x-0))
ind=(ssp[:,0]==Zsol) & (ssp[:,1] > 3.6)
zz=[z_at_value(cosmo.age,10**(i-3.0)*u.Gyr+0.545*u.Gyr) for i in ssp[ind,1]]

ind1 = (0.0   < x[:,73]) & (x[:,73] < 0.100) & (x[:,16] > 0.0) & (x[:,55] > 0.0) & (x[:,52] > 0.0)	# 0.0     z < 0.1
ind2 = (0.166 < x[:,73]) & (x[:,73] < 0.366) & (x[:,13] > 0.0) & (x[:,52] > 0.0) & (x[:,49] > 0.0)	# 0.166 < z < 0.266
ind3 = (0.514 < x[:,73]) & (x[:,73] < 0.714) & (x[:,10] > 0.0) & (x[:,49] > 0.0) & (x[:,46] > 0.0)	# 0.614 < z < 0.714
ind4 = (1.117 < x[:,73]) & (x[:,73] < 1.317) & (x[:,7] > 0.0)  & (x[:,40] > 0.0) & (x[:,37] > 0.0)	# 1.117 < z < 1.317

ind5=ind1 | ind2 | ind3 | ind4

rV1=x[ind1,55]+(x[ind1,52]-x[ind1,55])*(Vs[0]-opt_ew[0])/(opt_ew[1]-opt_ew[0])
rV2=x[ind2,52]+(x[ind2,49]-x[ind2,52])*(Vs[1]-opt_ew[1])/(opt_ew[2]-opt_ew[1])
rV3=x[ind3,49]+(x[ind3,46]-x[ind3,49])*(Vs[2]-opt_ew[2])/(opt_ew[3]-opt_ew[2])
rV4=x[ind4,40]+(x[ind4,37]-x[ind4,40])*(Vs[3]-opt_ew[4])/(opt_ew[5]-opt_ew[4])

rVs=np.concatenate([rV1,rV2,rV3,rV4])
rLs=np.concatenate([x[ind1,16],x[ind2,13],x[ind3,10],x[ind4,7]])
beta=rVs/rLs

l1=pl.plot(zz,ssp[ind,2])
l2=pl.plot(zz,tau1G[ind,2])
l3=pl.plot(zz,tau10G[ind,2])
l4=pl.plot(zz,sfh3[ind,2])
l5=pl.plot(zz,sfh4[ind,2])
l6=pl.plot(zz,sfh5[ind,2])
o1=pl.scatter(x[ind1,73],rV1*(10**(sp[ind1,5]/2.5))/x[ind1,16],color='r',s=10,marker='o',alpha=.2)
o2=pl.scatter(x[ind2,73],rV2*(10**(sp[ind2,5]/2.5))/x[ind2,13],color='b',s=10,marker='o',alpha=.2)
o3=pl.scatter(x[ind3,73],rV3*(10**(sp[ind3,5]/2.5))/x[ind3,10],color='g',s=10,marker='o',alpha=.2)
o4=pl.scatter(x[ind4,73],rV4*(10**(sp[ind4,5]/2.5))/x[ind4,7],color='c',s=10,marker='o',alpha=.2)
o5=pl.scatter(x[ind5,73],rVs/rLs,color='grey',s=5,marker='o',alpha=.1,label='A_V uncorrected')
legend1=pl.legend([l1,l2,l3,l4,l5,l6],["SSP, V/I1","tau = 1Gyr, V/I1","tau = 10Gyr, V/I1","SFH3, V/I1","SFH4, V/I1","SFH5, V/I1"],loc=2)
pl.gca().add_artist(legend1)
pl.legend([o1,o2,o3,o4,o5],["0.00 < z < 0.10, (g r)/I1","0.17< z < 0.37, (r i)/I2","0.51< z < 0.71, (i z)/I3","1.12 < z < 1.32, (J2J3)/I4"],loc=1)
pl.ylim(0,10)
pl.xlabel('z')
pl.ylabel('BetaV')
#pl.gca().add_artist(legend1)
#pl.gca().add_artist(legend2)
pl.show()               
