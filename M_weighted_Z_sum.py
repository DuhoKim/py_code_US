#####################################################
# Python script that calculate mass weighted metallicity
# based on the star-formation history from Behroozi+03
# and the metallicity set from Kim+17
# written by Duho Kim (9/12/18)
######################################################
# from pyraf import iraf
from astropy.io import fits
from astropy.io import ascii
import numpy as np
from astropy.cosmology import Planck15 as Cosmo
from operator import itemgetter
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

beh = ascii.read("/Users/dhk/work/data/Behroozi13/release-sfh_z0_z8_052913/sfh/sfh_release.dat")

Zs = np.array([0.05, 0.02, 0.008, 0.004, 0.0004, 0.0001])
Zs_log = np.log10(Zs / 0.02)
Zevol_coeff = [-0.18, -0.30, -0.58]

cols_for_SFH = ['red', 'green', 'blue']
linestyles_for_SFHs = ['-', '--', ':']

m15_init_idx = min(enumerate(np.abs(beh['col2']-15)),key=itemgetter(1))[0]		# find halo mass 10^15 Msun for SFH3
m13_init_idx = min(enumerate(np.abs(beh['col2']-13)),key=itemgetter(1))[0]		# find halo mass 10^13 Msun for SFH4
m11_init_idx = min(enumerate(np.abs(beh['col2']-11.4)),key=itemgetter(1))[0]		# find halo mass 10^11.4 Msun for SFH5

m15_idx = [beh['col2'] == beh['col2'][m15_init_idx]]
m13_idx = [beh['col2'] == beh['col2'][m13_init_idx]]
m11_idx = [beh['col2'] == beh['col2'][m11_init_idx]]

nz = np.sum(m15_idx)
zs = beh['col1'][m15_idx] - 1
dts = np.zeros(nz-1)
lbt = Cosmo.lookback_time(beh['col1'][m15_idx]-1)
for i in range(0, nz-1):
    dts[i] = lbt[i+1].value - lbt[i].value

Zevol_record = np.zeros((len(Zevol_coeff), nz, len(Zs)))
################ M weighted Z calculation  #######################
for k in range(0, 3):                       # for each SFH
    if k == 0:
        sfr = 10 ** beh['col3'][m15_idx]
    elif k ==1:
        sfr = 10 ** beh['col3'][m13_idx]
    else:
        sfr = 10 ** beh['col3'][m11_idx]
    for i in range(0, nz-1):              # for each z
        nomi = np.zeros(len(Zs))            # Z * dz * SFR
        denomi = 0                          # dz * SFR
        for j in range(i, nz-1):
            Z = Zs_log + zs[j] * Zevol_coeff[k]
            for l in range(0, 6):
                if Z[l] < np.min(Zs_log):
                    Z[l] = np.min(Zs_log)
            nomi = nomi + Z  * sfr[j] * dts[j]
            denomi = denomi + sfr[j] * dts[j]
        Zevol_record[k, i, :] = nomi / denomi

fig = plt.figure(figsize=(8, 6))

h1 = fig.add_axes([ 0.15,    0.15,   0.8,  0.7])

for i in range(0, 3):           # for each SFHs
    for j in range(0, 3):       # for each Zs
        h1.plot(zs, Zevol_record[i, :, j], color=cols_for_SFH[i], linestyle=linestyles_for_SFHs[i], linewidth=3.0)

h1.set_ylabel(r'<log(Z/Z$_{\odot})>_{\mathrm{M}}$', fontsize=30)
h1.set_xlim([0, 7])
h1.set_ylim(-2.5,0)
h1.tick_params(direction='in',top=True,right=True,labelsize=20)
h1.set_xlabel('Redshift', fontsize=30)

custom_lines1 = [       Line2D([0],[0],color=cols_for_SFH[0], linestyle='-', linewidth=3),
                        Line2D([0],[0],color=cols_for_SFH[1], linestyle='--', linewidth=3),
                        Line2D([0],[0],color=cols_for_SFH[2], linestyle=':', linewidth=3)  ]
custom_lines2 = [       Line2D([0],[0],color=cols_for_SFH[0], linestyle='-', linewidth=3),
                        Line2D([0],[0],color=cols_for_SFH[1], linestyle='--', linewidth=3),
                        Line2D([0],[0],color=cols_for_SFH[2], linestyle=':', linewidth=3)  ]
custom_lines3 = [       Line2D([0],[0],color=cols_for_SFH[0], linestyle='-', linewidth=3),
                        Line2D([0],[0],color=cols_for_SFH[1], linestyle='--', linewidth=3),
                        Line2D([0],[0],color=cols_for_SFH[2], linestyle=':', linewidth=3)   ]
# custom_lines4 = [       Line2D([0],[0],color='red',  linewidth=3, linestyle='-'),
#                         Line2D([0],[0],color='green',  linewidth=3, linestyle='--'),
#                         Line2D([0],[0],color='blue',  linewidth=3, linestyle=':')   ]
# custom_lines5 = [       Line2D([0],[0],color='silver', linestyle='-', alpha=0.5, linewidth=3),
#                         Line2D([0],[0],color='silver', linestyle='--', alpha=0.5, linewidth=3),
#                         Line2D([0],[0],color='silver', linestyle=':', alpha=0.5, linewidth=3)   ]

legend_texts1 = [r'<Zevol3(Z0=2.5, 1.0, 0.4$\,$Z$\odot$)>$_{\mathrm{M}}$', r'<Zevol4(Z0=2.5, 1.0, 0.4$\,$Z$\odot$)>$_{\mathrm{M}}$', r'<Zevol5(Z0=2.5, 1.0, 0.4$\,$Z$\odot$)>$_{\mathrm{M}}$']
#legend_texts2 = [r'<Zevol3(Z0=Z$\odot$)>$_{\mathrm{M}}$', r'<Zevol4(Z0=Z$\odot$)>$_{\mathrm{M}}$', r'<Zevol5(Z0=Z$\odot$)>$_{\mathrm{M}}$']
#legend_texts3 = [r'<Zevol3(Z0=0.4$\,$Z$\odot$)>$_{\mathrm{M}}$', r'<Zevol4(Z0=0.4$\,$Z$\odot$)>$_{\mathrm{M}}$', r'<Zevol5(Z0=0.4$\,$Z$\odot$)>$_{\mathrm{M}}$']
# legend_texts4 = [       r'Zevol3 = -0.18 z + log Z0',
#                         r'Zevol4 = -0.30 z + log Z0',
#                         r'Zevol5 = -0.58 z + log Z0'  ]
# legend_texts5 = [       r'SFH3 (log(M/M$\odot$)=11.4',
#                         r'SFH4 (log(M/M$\odot$)=10.95',
#                         r'SFH5 (log(M/M$\odot$)=9.5'  ]

# h1.text(0.3, 0.9,'Behroozi+13', fontsize=30)
h1.text(0.5,-2.3,'K17', fontsize=30)

legend1 = h1.legend(custom_lines1,legend_texts1,loc='lower left',bbox_to_anchor=(0.55,0.8),frameon=False)
#legend2 = h1.legend(custom_lines2,legend_texts2,loc='lower left',bbox_to_anchor=(0.33,1),frameon=False)
#legend3 = h1.legend(custom_lines3,legend_texts3,loc='lower left',bbox_to_anchor=(0.66,1),frameon=False)
# legend4 = h1.legend(custom_lines4,legend_texts4,loc='lower left',bbox_to_anchor=(0.45,0.75),frameon=False, fontsize=15)
# legend5 = h1.legend(custom_lines5,legend_texts5,loc='lower left',bbox_to_anchor=(0.47,0.75),frameon=False, fontsize=15)

plt.gca().add_artist(legend1)
#plt.gca().add_artist(legend2)
#plt.gca().add_artist(legend3)
# plt.gca().add_artist(legend4)
# plt.gca().add_artist(legend5)

# zevol008 = Zs_log[0] + Zevol_coeff[2] * zs
# zevol008[np.where(zevol008 < Zs_log[5])] = Zs_log[5]
# h1.plot(zs,Zs_log[0] + Zevol_coeff[0] * zs, color='red',  linewidth=3)
# h1.plot(zs,Zs_log[0] + Zevol_coeff[1] * zs, color='green', linewidth=3, linestyle='--')
# h1.plot(zs,zevol008, color='blue',  linewidth=3, linestyle=':')
#
# zevol008 = Zs_log[1] + Zevol_coeff[2] * zs
# zevol008[np.where(zevol008 < Zs_log[5])] = Zs_log[5]
# h1.plot(zs,Zs_log[1] + Zevol_coeff[0] * zs, color='red', linewidth=3)
# h1.plot(zs,Zs_log[1] + Zevol_coeff[1] * zs, color='green',  linewidth=3, linestyle='--')
# h1.plot(zs,zevol008, color='blue', linewidth=3, linestyle=':')
#
# zevol008 = Zs_log[2] + Zevol_coeff[2] * zs
# zevol008[np.where(zevol008 < Zs_log[5])] = Zs_log[5]
# h1.plot(zs,Zs_log[2] + Zevol_coeff[0] * zs, color='red',  linewidth=3)
# h1.plot(zs,Zs_log[2] + Zevol_coeff[1] * zs, color='green', linewidth=3, linestyle='--')
# h1.plot(zs,zevol008, color='blue',  linewidth=3, linestyle=':')

# h2 = h1.twinx()

# h1.plot(zs,1e1 ** (beh['col3'][m15_idx]+9.0-11.4), color='red',  linewidth=3)
# h1.plot(zs,1e1 ** (beh['col3'][m13_idx]+9.0-10.95), color='green', linestyle='--',   linewidth=3)
# h1.plot(zs,1e1 ** (beh['col3'][m11_idx]+9.0-9.5), color='blue', linestyle=':', linewidth=3)

# h1.set_ylabel('sSFR (arbitrary)', fontsize=30)
# h2.tick_params(labelsize='large', direction='in')


fig.savefig('/Users/dhk/Documents/publish/ngcic/Zevol_sum.pdf')
plt.close(fig)