#####################################################
# Python script of reading output file generated during
# building PSF for given FITS images 
# written by Duho Kim (11/30/17)        
######################################################
from astropy.io import ascii
from astropy.table import Table
import os
import numpy as np

data_dir = '/Users/dhk/work/data/NGC_IC/SDSS/'
cat_dir = '/Users/dhk/work/cat/NGC_IC/'
work_dir = data_dir+'g_bench1'

os.chdir(work_dir)

sha_cat=ascii.read(cat_dir+'sha_quarry_batch_257_without_prefix.txt')

out = Table(names=('ID','NUMSTARS','FUNC1','NORM1','MAG1','FUNC2','NORM2','MAG2','FUNC3','NORM3','MAG3'),dtype=('S7','i4','S8','f8','f8','S8','f8','f8','S8','f8','f8'))

#for x in range(0,len(sha_cat)):
for x in range(0,231):
	galid = sha_cat['id'][x]+'-g'
	if os.path.isfile(galid+'.psf1.out'):
		out1 = open(galid+'.psf1.out','r').read()
		numstars = int(out1.split('PSF stars',1)[0].split()[-1])
		try:
			func1 = out1.split('Best fitting function is ',1)[1].split()[0]
			norm1 = float(out1.split(func1+' norm scatter = ',1)[1].split()[0])
			mag1 = float(out1.split('Psfmag: ',1)[1].split()[0])
		except:
			out.add_row([galid,numstars,'',np.nan,np.nan,'',np.nan,np.nan,'',np.nan,np.nan])
			continue
		if os.path.isfile(galid+'.psf2.out'):
			out2 = open(galid+'.psf2.out','r').read()
			try:
	                	func2 = out2.split('Best fitting function is ',1)[1].split()[0]
        	        	norm2 = float(out2.split(func1+' norm scatter = ',1)[1].split()[0])
                		mag2 = float(out2.split('Psfmag: ',1)[1].split()[0])
			except:
				out.add_row([galid,numstars,func1,norm1,mag1,'',np.nan,np.nan,'',np.nan,np.nan])
				continue
	                if os.path.isfile(galid+'.psf3.out'):
				out3 = open(galid+'.psf3.out','r').read()
				try:
	        	                func3 = out3.split('Best fitting function is ',1)[1].split()[0]
        	        	        norm3 = float(out3.split(func1+' norm scatter = ',1)[1].split()[0])
                	        	mag3 = float(out3.split('Psfmag: ',1)[1].split()[0])
					out.add_row([galid,numstars,func1,norm1,mag1,func2,norm2,mag2,func3,norm3,mag3])		
				except:
					out.add_row([galid,numstars,func1,norm1,mag1,func2,norm2,mag2,'',np.nan,np.nan])
					continue
			else:
				out.add_row([galid,numstars,func1,norm1,mag1,func2,norm2,mag2,'',np.nan,np.nan])
		else:
			out.add_row([galid,numstars,func1,norm1,mag1,'',np.nan,np.nan,'',np.nan,np.nan])
	else:
		out.add_row([galid,np.nan,'',np.nan,np.nan,'',np.nan,np.nan,'',np.nan,np.nan])

out.write('numstars_fit_func_norm_mag.txt',format='ascii.rst',overwrite=True)
