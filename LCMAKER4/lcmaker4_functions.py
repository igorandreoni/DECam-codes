##########################
####### FUNCTIONS ########
##########################
import pdb
from astropy import wcs
from astropy.wcs import WCS
import sewpy   # SExtractor wrapper
import pyfits
from astropy.io import fits
import random as random
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import scipy.signal

# Convert RA, DEC to x,y coordinates
def wcstoxy(image,ra,dec):

	w=WCS(image)
	x, y = w.all_world2pix(ra, dec, 0)
	#add 1px to be consistent with ds9:
	x1=x
	y1=y
	return (x1,y1)

# Convert x,y to RA,DEC coordinates
def xytowcs(image,x,y):

	w=WCS(image)
	ra, dec = w.all_pix2world(x, y, 0)

	return (ra,dec)

#################################
###################################
# Convert flagmap into weightmap
###################################

def convertflag(outflagname,hdulistmask,weightimagename):
	hdulistmask = fits.open(outflagname)
	maskarr=hdulistmask[0].data
	for j in range(len(maskarr)):
		for k in range(len(maskarr[1])):
			if maskarr[j,k] <= 64:
				maskarr[j,k] = 1
			else:
				maskarr[j,k] = 0
	hdulistmask.writeto(weightimagename)
	hdulistmask.close()	
	return

#################################
# Add mock stars
#################################


# First define a function to create a PSF
def makePSF(size, fwhm=1.0, center=None):
# makes a Gaussian PSF with fwhm
    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]
    if center is None:
        x0 = y0 = size // 2.0
    else:
        x0 = center[0]
        y0 = center[1]
    return np.exp(-4.0*np.log(2)*((x-x0)**2 +(y-y0)**2)/fwhm**2)

def addstar(img,image,psffolder,imalignfolder,ramock0, decmock0):
	hdulistimg = pyfits.open(img)
	hduimg=hdulistimg[0].header # Primary header
	dataimg=hdulistimg[0].data

	# Create an image with mock stars, that will be added to the input image.
	hdu = fits.PrimaryHDU()
	npsf = 22

	ny = len(dataimg[0]) -4 + npsf 
	nx = len(dataimg)-4 + npsf 

	# create seeing PSF
	hdu.data = makePSF(19, fwhm=3.90)
	hdu.writeto(psffolder+'/psf.fits', clobber=True)
	psf = hdu.data

	# define CCD frame size
	hdu.data = np.zeros((nx,ny))

	# add stars around the frame
	nstars = 3
	alpha = 1.0
	zeropoint=25

	# Defined position of the stars:  
	#Give default values if the image is the first one: 
	xstar,ystar=wcstoxy(img,ramock0, decmock0)

#Open the file listing the mock stars
	#mocklist='/lustre/projects/p025_swin/pipes/DWF_PIPE/MARY_WORK/4hr_151218_m18oldall/ccd4/images_resampled/mocklist.cat'
	#m = open(mocklist, 'w')
	for i in range(0,nstars):
# power law distribution of sflux with slope -alpha
   		#sflux = 5.0*np.abs(np.random.rand())**(-1.0/alpha)
   		sflux = 15000
   		
    	#xstar=nx*np.random.rand()
    	#ystar=ny*np.random.rand()
    	hdu.data[int(xstar),int(ystar)] = sflux
    	#outputline=str(i) + ' ' + str(ystar) + ' ' + str(xstar) + ' ' + str(sflux) + ' ' + str(zeropoint-2.5*np.log10(sflux)) + ' \n'
    	#m.write(outputline)
    	xstar=xstar+150
    	ystar=ystar+150
	#m.close()
# convolve model with seeing
	hdu.data = scipy.signal.convolve2d(hdu.data, psf, mode='valid')

	# write out the model (i.e. noise free) as a fits file
	modelname=psffolder+'/'+image+'_model.fits'
	hdu.writeto(modelname, clobber=True)

# Now add the model to the input image.
	mockstars = hdu.data 
	mock=dataimg+mockstars
	header = hduimg
# write out the data file with a coordinate frame
###pyfits.writeto('mock.fits', hdu.data, header, clobber=True)
	outname=imalignfolder + '/' + image + '_mock.fits'
	pyfits.writeto(outname, mock, header, clobber=True)
	return outname


#################################
# SExtractor (wrapped with sewpy)
#################################

def sewfunction(insex,incatalog,seeing,zeropoint,DETECT_MINAREA,THRESH_TYPE,DETECT_THRESH,ANALYSIS_THRESH,sexconfigfile,prefixname):
	sew = sewpy.SEW(params=["VECTOR_ASSOC(3)",\
	"FLUX_AUTO", \
	"FLUXERR_AUTO",   \
	"MAG_AUTO", \
	"MAGERR_AUTO",\
	"FLUX_ISO", \
	"FLUX_ISOCOR",\
	"MAG_ISOCOR",\
	"MAG_ISOCOR", \
	"X_IMAGE",\
	"Y_IMAGE",\
	"X_WORLD",\
	"Y_WORLD", \
	"BACKGROUND"],\
	config={\
#	"ASSOC_NAME": incatalog,\
	"ASSOC_DATA": 0,\
#	"ASSOC_PARAMS": '10,11',\
	"ASSOC_RADIUS": 4,\
#"PARAMETERS_NAME":sexparfile, \
#	"CATALOG_NAME":sexcatalogname, \
	"SEEING_FWHM":str(seeing), \
	"MAG_ZEROPOINT": str(zeropoint), \
	"DETECT_MINAREA": DETECT_MINAREA, \
	"THRESH_TYPE": THRESH_TYPE, \
	"DETECT_THRESH": DETECT_THRESH, \
	"ANALYSIS_THRESH": ANALYSIS_THRESH},\
	configfilepath=sexconfigfile,sexpath='/home/tpritcha/bin/sex'\
# prefix="new"\, assoc_cat=outbase, assoc_xname="X_IMAGE", assoc_yname="Y_IMAGE"
	) 
	
	out = sew(insex, prefix=prefixname, assoc_cat=incatalog, assoc_xname="X_IMAGE", assoc_yname="Y_IMAGE")

	return out