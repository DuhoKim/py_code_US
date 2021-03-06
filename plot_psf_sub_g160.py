from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

box_size=15

hgrp=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/g.r.pp.160.fits')
hgrps=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/g.r.pp.160.sub.1.fits')
hgpsf=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/g.r.psf.1.out.160.fits')
hipsf=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/i1.psf.1.out.160.fits')

grp=hgrp[0].data
grps=hgrps[0].data
gpsf=hgpsf[0].data
ipsf=hipsf[0].data

cat2=ascii.read("/Users/dhk/work/data/NGC_IC/pilot4/g.r.p.cat")
cat=cat2[(cat2['FLAGS'] == 0) & (cat2['FWHM_IMAGE'] < 6) & (cat2['FWHM_IMAGE'] > 1.5)]
cat.sort(['MAG_AUTO'])

submean=[]
substd=[]
subcol=[]

fig = plt.figure(figsize=(8,8))

h=fig.add_axes([0,0.75,0.25,0.25])
h.text(-0.4,600,'30x30 pixels')
h.set_xlabel('pixel value')
h.text(-0.4,800,'g_r_p_D160-I1_PSF_D160')
h.set_xlim(-0.5,0.5)

s=fig.add_axes([0.3,0.75,0.25,0.25])
s.set_xlabel('mean')
s.set_ylabel('std')
s.set_xlim(-0.1,0.5)
s.set_ylim(0,4)

g=fig.add_axes([0.6,0.75,0.2,0.2])
g.imshow(gpsf,cmap='gray_r',clim=[0,np.mean(gpsf)])
g.text(5,-5,'g PSF w/o mask out')

i2=fig.add_axes([0.8,0.75,0.2,0.2])
i2.imshow(ipsf,cmap='gray_r',clim=[0,np.mean(ipsf)])
i2.text(5,-5,'I1 PSF w/o mask out')


for i in range(0,30):
        pos1 = [ (i%5)*0.2, 0.55 - (i/5)*0.1, 0.1, 0.1]
        pos2 = [ 0.1+(i%5)*0.2, 0.55 - (i/5)*0.1, 0.1, 0.1]
	a1 = fig.add_axes(pos1)
	a2 = fig.add_axes(pos2)
	a1.axes.get_xaxis().set_visible(False)
	a2.axes.get_xaxis().set_visible(False)
	a1.axes.get_yaxis().set_visible(False)
	a2.axes.get_yaxis().set_visible(False)
	thumb=grp[int(cat['Y_IMAGE'][i])-box_size:int(cat['Y_IMAGE'][i])+box_size,int(cat['X_IMAGE'][i])-box_size:int(cat['X_IMAGE'][i])+box_size]
	thumbs=grps[int(cat['Y_IMAGE'][i])-box_size:int(cat['Y_IMAGE'][i])+box_size,int(cat['X_IMAGE'][i])-box_size:int(cat['X_IMAGE'][i])+box_size]
	a1.imshow(thumb, cmap='gray_r', clim=[0,np.mean(thumb)])
	a2.imshow(thumbs, cmap='gray_r', clim=[0,np.mean(thumb)])
	y,binEdges = np.histogram(thumbs,bins=np.arange(-1,1,0.1))
	bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
	base_line, =	h.plot(bincenters,y,'-') 
	color_cycle = base_line.get_color()
	submean.append(np.mean(thumbs))
	substd.append(np.std(thumbs))
	a1.text(3,5,'%d' % i,color=color_cycle)
	a2.text(3,5,'%4.2f' % cat['FWHM_IMAGE'][i],color='red')
	subcol.append(color_cycle)
	


s.scatter(submean,substd)
for i in range(0,30):
	s.annotate(i,(submean[i],substd[i]),color=subcol[i])


fig.savefig("/Users/dhk/work/data/NGC_IC/pilot4/plot_grp160_mean.pdf")
			
