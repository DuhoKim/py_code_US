from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

box_size=15

hgrp=fits.open('/Users/dhk/work/data/NGC_IC/pilot3/r.r.fits')
hgrps=fits.open('/Users/dhk/work/data/NGC_IC/pilot3/r.r.pp.def.sub.1.fits')

grp=hgrp[0].data
grps=hgrps[0].data

cat2=ascii.read("/Users/dhk/work/data/NGC_IC/pilot3/test.cat")
cat=cat2[(cat2['X_IMAGE'] > 889) & (cat2['X_IMAGE'] < 1485) & (cat2['Y_IMAGE'] > 1390) & (cat2['Y_IMAGE'] < 1983)]
cat.sort(['MAG_AUTO'])

submean=[]
substd=[]
subcol=[]

fig = plt.figure(figsize=(8,8))

h=fig.add_axes([0,0.75,0.25,0.25])
h.text(-0.9,600,'30x30 pixels')
h.set_xlabel('pixel value')
h.text(-0.9,800,'r_r_p25-I1_PSF25')

s=fig.add_axes([0.3,0.75,0.25,0.25])
s.set_xlabel('mean')
s.set_ylabel('std')
for i in range(0,15):
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
	a1.imshow(thumb, cmap='gray_r', clim=[0,1.0])
	a2.imshow(thumbs, cmap='gray_r', clim=[0,1.0])
	y,binEdges = np.histogram(thumbs,bins=np.arange(-1,1,0.1))
	bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
	base_line, =	h.plot(bincenters,y,'-') 
	color_cycle = base_line.get_color()
	submean.append(np.mean(thumbs))
	substd.append(np.std(thumbs))
	a1.text(3,5,'%d' % i,color=color_cycle)
	subcol.append(color_cycle)
	


s.scatter(submean,substd)
for i in range(0,15):
	s.annotate(i,(submean[i],substd[i]),color=subcol[i])


fig.savefig("/Users/dhk/work/data/NGC_IC/pilot3/plot_rrp25.pdf")
			
