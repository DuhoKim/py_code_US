#####################################################
# Python script of generating DS9 region file
# from IRAF output mag.1 file
# written by Duho Kim (1/18/18)        
######################################################
import numpy as np
from astropy.io import ascii
from astropy.table import Table
from regions import DS9Parser, write_ds9

work_dir = '/Users/dhk/work/data/NGC_IC/pilot_script/'

mag1 = ascii.read(work_dir+'ic1065.mag.1')


reg_string = 'galactic\ncircle(42,43,3) # color=green'

parser = DS9Parser(reg_string)
parser.run()


reg = Table(names=('shape','x','y','rad'),dtype=('S6','f8','f8','f8'))

for i in range(len(mag1)):
	write_ds9('image\ncircle('+mag1['XCENTER'][i]+','+mag1['YCENTER'][i]+')')
	reg.add_row( ('circle',mag1['XCENTER'][i],mag1['YCENTER'][i],np.sqrt(mag1['AREA'][i]/np.pi)) )

reg.write(work_dir+'ic1065.reg',format='ascii.commented_header',comment='# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1')


