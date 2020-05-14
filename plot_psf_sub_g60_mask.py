from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

box_size=30
n=7      # a threshold throttle (with satisfactory values between 5 and 10) Jarrett et al. 2003
scl=(-0.05,0.05)

hgrp=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/gr.pp.60.fits')
hgrps=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/gr.pp.sub.60.1.fits')
hgpsf=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/gr.psf.160.out.5.60D.fits')
hipsf=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/i2.psf.160.out.6.60D.fits')

grp=hgrp[0].data
grps=hgrps[0].data
gpsf=hgpsf[0].data
ipsf=hipsf[0].data

cat2=ascii.read("/Users/dhk/work/data/NGC_IC/pilot4/g.r.p.cat")
cat=cat2[(cat2['FLAGS'] == 0) & (cat2['FWHM_IMAGE'] < 3.5) & (cat2['FWHM_IMAGE'] > 3)]
cat.sort(['MAG_AUTO'])

submean=[]
substd=[]
subcol=[]

fig = plt.figure(figsize=(8,8))

h=fig.add_axes([0,0.75,0.25,0.25])
h.set_xlabel('pixel value')
h.set_xlim(-0.3,0.3)

s=fig.add_axes([0.3,0.75,0.25,0.25])
s.set_xlabel('mean')
s.set_ylabel('std')
s.set_xlim(-0.01,0.05)
s.set_ylim(0,1)
s.text(0,0.9,'60x60 pixels')
s.text(0,0.8,'g_r_p_D60-I1_PSF_D60')

g=fig.add_axes([0.6,0.75,0.2,0.2])
g.axes.get_xaxis().set_visible(False)
g.axes.get_yaxis().set_visible(False)
g.imshow(gpsf,cmap='gray_r',clim=scl)
g.text(5,-5,'g PSF w/ mask out')

i2=fig.add_axes([0.8,0.75,0.2,0.2])
i2.axes.get_xaxis().set_visible(False)
i2.axes.get_yaxis().set_visible(False)
i2.imshow(ipsf,cmap='gray_r',clim=scl)
i2.text(5,-5,'I1 PSF w/ mask out')

gx=fig.add_axes([0.05,0.45,0.2,0.2])
gx.plot(range(0,60),np.log10(np.abs(grp[int(cat['Y_IMAGE'][0])-box_size:int(cat['Y_IMAGE'][0])+box_size,int(cat['X_IMAGE'][0])])),'-')
gx.set_xlabel('No 0 horizontal')
gx.set_ylim(-4,3)
gx.grid(True)
gx.set_xticks(range(0,61,15))

gy=fig.add_axes([0.3,0.45,0.2,0.2])
gy.grid(True)
gy.plot(range(0,60),np.log10(np.abs(grp[int(cat['Y_IMAGE'][0]),int(cat['X_IMAGE'][0])-box_size:int(cat['X_IMAGE'][0])+box_size])),'-')
gy.set_xlabel('No 0 vertical')
gy.set_ylim(-4,3)
gy.set_xticks(range(0,61,15))

i2x=fig.add_axes([0.55,0.45,0.2,0.2])
i2x.plot(range(0,60),np.log10(np.abs(grp[int(cat['Y_IMAGE'][1])-box_size:int(cat['Y_IMAGE'][1])+box_size,int(cat['X_IMAGE'][1])])),'-')
i2x.set_xlabel('No 1 horizontal')
i2x.set_ylim(-4,3)
i2x.grid(True)
i2x.set_xticks(range(0,61,15))

i2y=fig.add_axes([0.8,0.45,0.2,0.2])
i2y.grid(True)
i2y.plot(range(0,60),np.log10(np.abs(grp[int(cat['Y_IMAGE'][1]),int(cat['X_IMAGE'][1])-box_size:int(cat['X_IMAGE'][1])+box_size])),'-')
i2y.set_xlabel('No 1 vertical')
i2y.set_ylim(-4,3)
i2y.set_xticks(range(0,61,15))

for i in range(0,15):
        pos1 = [ (i%5)*0.2, 0.2 - (i/5)*0.1, 0.1, 0.1]
        pos2 = [ 0.1+(i%5)*0.2, 0.2 - (i/5)*0.1, 0.1, 0.1]
	a1 = fig.add_axes(pos1)
	a2 = fig.add_axes(pos2)
	a1.axes.get_xaxis().set_visible(False)
	a2.axes.get_xaxis().set_visible(False)
	a1.axes.get_yaxis().set_visible(False)
	a2.axes.get_yaxis().set_visible(False)
	thumb=grp[int(cat['Y_IMAGE'][i])-box_size:int(cat['Y_IMAGE'][i])+box_size,int(cat['X_IMAGE'][i])-box_size:int(cat['X_IMAGE'][i])+box_size]
	thumbs=grps[int(cat['Y_IMAGE'][i])-box_size:int(cat['Y_IMAGE'][i])+box_size,int(cat['X_IMAGE'][i])-box_size:int(cat['X_IMAGE'][i])+box_size]
	a1.imshow(thumb, cmap='gray_r',clim=scl)
	a2.imshow(thumbs, cmap='gray_r',clim=scl)
	y,binEdges = np.histogram(thumbs,bins=np.arange(-1,1,0.05))
	bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
	base_line, =	h.plot(bincenters,y,'-') 
	color_cycle = base_line.get_color()
#	h.hist(thumbs,bins=np.arange(-1,1,0.1),histtype='bar')
	submean.append(np.mean(thumbs))
	substd.append(np.std(thumbs))
	a1.text(3,10,'%d' % i,color=color_cycle)
	a2.text(3,10,'%4.2f' % cat['FWHM_IMAGE'][i],color='red')
	subcol.append(color_cycle)
	


s.scatter(submean,substd)
for i in range(0,15):
	s.annotate(i,(submean[i],substd[i]),color=subcol[i])


fig.savefig("/Users/dhk/work/data/NGC_IC/pilot4/plot_grp60_order.pdf")
			
