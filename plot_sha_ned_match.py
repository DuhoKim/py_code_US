import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from astropy.io import ascii
from matplotlib.ticker import NullFormatter

match=ascii.read("sha_ned_match.txt")
x = match['off']
y = match['a']

nullfmt = NullFormatter()

left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]
rect_inset = [0.5, 0.5, 0.2, 0.2]

plt.figure(1, figsize=(8,8))

axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)
axInset = plt.axes(rect_inset)

axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

axScatter.scatter(x,y,s=1)

axInset.hist(match['fov'],[0,5,10,15,20,25,30,35,40,45,50,55,60,65])
axInset.set_xlabel('s_fov [\']')

axScatter.set_ylabel('Major axis 2a [\']')
axScatter.set_xlabel('Offset [\']')
axHistx.set_title('NGC/IC w/ IRAC1 mosaic data & Major axis value from NED')

xbinwidth = 0.1
ybinwidth = 0.1
#xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
#lim = (int(xymax/binwidth) + 1) * binwidth

axScatter.set_xlim((0, 5))
axScatter.set_ylim((0, 10))

xbins = np.arange(0, 5 + xbinwidth, xbinwidth)
ybins = np.arange(0, 10 + ybinwidth, ybinwidth)
axHistx.hist(x, bins=xbins)
#yy=mlab.normpdf(xbins,np.mean(x),np.std(x))
#axHistx.plot(xbins,yy*5000,'k--',linewidth=1.5)
axHistx.axvline(x=1.0,color='r')
#axHistx.axvline(x=np.mean(x),color='k',linewidth=1.0)
#axHistx.plot([np.mean(x)-np.std(x),np.mean(x)+np.std(x)],[55,55],'k--',linewidth=1.0)

#axScatter.axvline(x=0)
#axScatter.axhline(y=12)
#axScatter.text(-35,38,'%d' % (len(x[np.logical_and(x<0,y>12)])))
#axScatter.text(-35,2,'%d' % (len(x[np.logical_and(x<0,y<12)])))
#axScatter.text(50,38,'%d' % (len(x[np.logical_and(x>0,y>12)])))
#axScatter.text(50,2,'%d' % (len(x[np.logical_and(x>0,y<12)])))

#axHistx.text(40,150,'mean: %3.1f' % (np.mean(x)))
#axHistx.text(40,120,'stddev: %4.1f' % (np.std(x)))
#axHistx.text(40,90,'min: %d' % (min(x)))
#axHistx.text(40,60,'max: %d' % (max(x)))
axHistx.text(0.5,100,'x<1 : %d' % (len(x[x<1.0])))
#axHistx.text(-10,20,'%d' % (len(x[x<0])),color='r')
#axHistx.text(10,20,'%d' % (len(x[x>0])),color='r')

axHisty.hist(y, bins=ybins, orientation='horizontal')

axHistx.set_xlim(axScatter.get_xlim())
axHisty.set_ylim(axScatter.get_ylim())

plt.savefig('a_vs_off.png')

plt.show()
