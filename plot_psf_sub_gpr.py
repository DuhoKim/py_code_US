from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

box_size=15

hgpr=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/g.p.r.fits')
hgprs=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/g.p.r.sub.1.fits')

gpr=hgpr[0].data
gprs=hgprs[0].data

cat2=ascii.read("/Users/dhk/work/data/NGC_IC/pilot4/g.p2.r.cat")
cat=cat2[(cat2['FLAGS'] == 0) & (cat2['FWHM_IMAGE'] < 6) & (cat2['FWHM_IMAGE'] > 3)]
cat.sort(['MAG_AUTO'])

submean=[]
substd=[]
subcol=[]

fig = plt.figure(figsize=(8,8))

h=fig.add_axes([0,0.75,0.25,0.25])
h.text(-0.9,600,'30x30 pixels')
h.set_xlabel('pixel value')
h.text(-0.9,800,'g_p_r-I1_PSF')

s=fig.add_axes([0.3,0.75,0.25,0.25])
s.set_xlabel('mean')
s.set_ylabel('std')
s.set_xlim(-0.05,0.05)
s.set_ylim(0,0.2)
for i in range(0,30):
        pos1 = [ (i%5)*0.2, 0.55 - (i/5)*0.1, 0.1, 0.1]
        pos2 = [ 0.1+(i%5)*0.2, 0.55 - (i/5)*0.1, 0.1, 0.1]
	a1 = fig.add_axes(pos1)
	a2 = fig.add_axes(pos2)
	a1.axes.get_xaxis().set_visible(False)
	a2.axes.get_xaxis().set_visible(False)
	a1.axes.get_yaxis().set_visible(False)
	a2.axes.get_yaxis().set_visible(False)
	thumb=gpr[int(cat['Y_IMAGE'][i])-box_size:int(cat['Y_IMAGE'][i])+box_size,int(cat['X_IMAGE'][i])-box_size:int(cat['X_IMAGE'][i])+box_size]
	thumbs=gprs[int(cat['Y_IMAGE'][i])-box_size:int(cat['Y_IMAGE'][i])+box_size,int(cat['X_IMAGE'][i])-box_size:int(cat['X_IMAGE'][i])+box_size]
	a1.imshow(thumb, cmap='gray_r', clim=[0,0.5])
	a2.imshow(thumbs, cmap='gray_r', clim=[0,0.5])
	y,binEdges = np.histogram(thumbs,bins=np.arange(-0.5,0.5,0.05))
	bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
	base_line, =	h.plot(bincenters,y,'-') 
	color_cycle = base_line.get_color()
	submean.append(np.mean(thumbs))
	substd.append(np.std(thumbs))
	a1.text(3,5,'%d' % i,color=color_cycle)
	subcol.append(color_cycle)
	


s.scatter(submean,substd)
for i in range(0,30):
	s.annotate(i,(submean[i],substd[i]),color=subcol[i])


fig.savefig("/Users/dhk/work/data/NGC_IC/pilot4/plot_gpr.pdf")
			
