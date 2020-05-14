from astropy.io import fits
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

box_size=15

hi1=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/i1.fits')
hi1s1=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/i1.sub.1.fits')
hi1s=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/i1.sub.2.fits')
hi1s3=fits.open('/Users/dhk/work/data/NGC_IC/pilot4/i1.sub.3.fits')
i1=hi1[0].data
i1s=hi1s[0].data

cat=ascii.read("/Users/dhk/work/data/NGC_IC/pilot4/i1.cat")
cat.sort(['MAG_AUTO'])

submean=[]
substd=[]
subcol=[]

fig = plt.figure(figsize=(8,8))

h=fig.add_axes([0,0.75,0.25,0.25])
h.text(-0.9,500,'30x30 pixels')
h.set_xlabel('pixel value')
h.text(-0.9,550,'I1-I1_PSF_32x32')

s=fig.add_axes([0.3,0.75,0.25,0.25])
s.set_xlabel('mean')
s.set_ylabel('std')
s.set_xlim(-0.15,0.15)
s.set_ylim(0,1)
for i in range(3,30):
        pos1 = [ (i%5)*0.2, 0.55 - (i/5)*0.1, 0.1, 0.1]
        pos2 = [ 0.1+(i%5)*0.2, 0.55 - (i/5)*0.1, 0.1, 0.1]
	a1 = fig.add_axes(pos1)
	a2 = fig.add_axes(pos2)
	a1.axes.get_xaxis().set_visible(False)
	a2.axes.get_xaxis().set_visible(False)
	a1.axes.get_yaxis().set_visible(False)
	a2.axes.get_yaxis().set_visible(False)
	thumb=i1[int(cat['Y_IMAGE'][i])-box_size:int(cat['Y_IMAGE'][i])+box_size,int(cat['X_IMAGE'][i])-box_size:int(cat['X_IMAGE'][i])+box_size]
	thumbs=i1s[int(cat['Y_IMAGE'][i])-box_size:int(cat['Y_IMAGE'][i])+box_size,int(cat['X_IMAGE'][i])-box_size:int(cat['X_IMAGE'][i])+box_size]
	a1.imshow(thumb, cmap='gray_r', clim=[0.0,0.1])
	a2.imshow(thumbs, cmap='gray_r', clim=[0.0,0.1])
	y,binEdges = np.histogram(thumbs,bins=np.arange(-1,1,0.1))
	bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
	base_line, =	h.plot(bincenters,y,'-') 
	color_cycle = base_line.get_color()
#	h.hist(thumbs,bins=np.arange(-1,1,0.1),histtype='bar')
	submean.append(np.mean(thumbs))
	substd.append(np.std(thumbs))
	a1.text(3,5,'%d' % i,color=color_cycle)
	subcol.append(color_cycle)
	


s.scatter(submean,substd)
for i in range(3,30):
	s.annotate(i,(submean[i-3],substd[i-3]),color=subcol[i-3])


fig.savefig("/Users/dhk/work/data/NGC_IC/pilot4/plot_i1_32d.pdf")
			
