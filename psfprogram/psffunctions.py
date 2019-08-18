import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
import scipy.signal
from PythonPhot import getpsf
from PythonPhot import aper
from PythonPhot import pkfit
from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u
import pdb 
import pyfits
import sys


########################################
# Convert RA, DEC to x,y coordinates
def wcstoxy(imagename,ra,dec):

	w=WCS(imagename)
	x, y = w.all_world2pix(ra, dec, 1)
	#add 1px to be consistent with ds9:
	x1=x
	y1=y
	return (x1,y1)

########################################
# Convert x,y to RA,DEC coordinates

def xytowcs(imagename,x,y):
	
	w=WCS(imagename)
	ra, dec = w.all_pix2world(x, y, 0)

	return (ra,dec)


#################################################
#Make a list of stars tho extract the PSF from.
def makelist(image, imagename,imagenameshort, RAcand, DECcand,pathstars='./starstore',magrange=[16.5,19.0],radiussearch='3',fwhmnominal=4.5):
	#find the dimensions of the image:
	xdimension=int(len(image[0]))
	ydimension=int(len(image))

	#Query the GSC catalog to find stars usable as reference:
	Vizier.ROW_LIMIT=999999999  #Remove the limit of 50 rows
	gsc = Vizier.query_region(coord.SkyCoord(ra=RAcand, dec=DECcand,unit=(u.deg, u.deg)),radius=radiussearch+'m',catalog="GSC2.3")
	####if imagenameshort=='c4d_151218_081209_ooi_g_3':
		#pdb.set_trace()
	# Has any source been found during the query? True if yes, False if not
	presentsource = (gsc!=[])

	# Verify that at least one of the queried sources falls within the CCD and within the magnitude range and is not coincident witht the candidate:
	z=0
	intheCCD=False
	inthemagrange=False
	notcoincidentwiththecand=False
	
	if presentsource:
		while (z < len(gsc[0]["RAJ2000"])  and  intheCCD == False  and  inthemagrange == False):
			Xstar, Ystar=wcstoxy(imagename+".fits",gsc[0]["RAJ2000"][z],gsc[0]["DEJ2000"][z])  
			#Fill possible masked value with 0.0
			gsc['I/305/out']["Fmag"].fill_value=0.0
			
			#Check location and magnitude to fall within acceptable values:
			checklocation=(Xstar > 2*fwhmnominal and Xstar < (xdimension-2*fwhmnominal) and  Ystar > 2*fwhmnominal and Ystar < (ydimension-2*fwhmnominal))  #Frame of 2*fwhm of acceptability for the location.
			checkmagnitude=((gsc[0]["Fmag"][z]) >= magrange[0] and gsc[0]["Fmag"][z] <= magrange[1] and  gsc[0]["Class"][z]==0)
			#Distance between the source and the candidate must be more than 2 arcsec
			dist=np.sqrt((np.mean(gsc[0]["RAJ2000"][z])-RAcand)**2 + (np.mean(gsc[0]["DEJ2000"][z]) - DECcand)**2)
			checkcoincidentwiththecand=dist >= 2.0/3600.0 
			
			#Change the values of the variables:
			if (checklocation == True):
				intheCCD=True	
			if (checklocation == True):
				inthemagrange=True
			if (checkcoincidentwiththecand == True):
				notcoincidentwiththecand=True
			
			z=z+1
			
	#Has any source good for the PSF extraction been found?
	if (presentsource and intheCCD and inthemagrange and notcoincidentwiththecand):
		OKgo=True
		print ">>> Your query returned at least one star suitable for the PSF extraction."
	else:
		OKgo=False
		print "!Watch out! Your query returned no star suitable for the PSF extraction. Try again with a larger radius or a different range of magnitudes. "
		print " "
						
			
	#If there is at least one source good for PSF extraction is found during the query, append the coordinates and the magnitudes to the lists
	if OKgo:		
		#Fill possible masked value with 0.0
		gsc['I/305/out']["Fmag"].fill_value=0.0

		#Open the file, that will be in the same directory of the images and will end with _stars.txt
		fs = open(pathstars+'/' + imagenameshort +"_stars.txt", 'w')
		for c in range(0,len(gsc[0]["Fmag"])):
			
			#Check that the source lies within the specified range of magnitude and that it is catalogued as "star".
			if ((gsc[0]["Fmag"][c]) >= magrange[0] and gsc[0]["Fmag"][c] <= magrange[1] and  gsc[0]["Class"][c]==0):
				
				#Check that the source actually falls INSIDE the CCD:
				#First, convert the RA, DEC into X, Y
				Xstar, Ystar=wcstoxy(imagename+".fits",gsc[0]["RAJ2000"][c],gsc[0]["DEJ2000"][c]) 
				#...then check that the star is within the limits of the CCD, assuming a frame of 20px from the borders:
				checklocation=(Xstar > 20 and Xstar < (xdimension-20) and  Ystar > 20 and Ystar < (ydimension-20))
				if checklocation: 
					fs.write(str(gsc[0]["RAJ2000"][c]) + " " + str(gsc[0]["DEJ2000"][c]) + "\n") 
		fs.close()
					



########################################
########################################
#Gets the PSF of the image.

def getPSFimage(image,imagename,starslist, Xcand,Ycand,pathpsf='./psfstore',fwhmnominal=5):
	
	#Convert the RA-DEC of the reference stars into XY coordinates:
	xpos,ypos=wcstoxy(imagename+".fits",starslist["RAstar"],starslist["DECstar"])
	
	#Check that all the stars are within the image. Report the number of removed stars in the log.
	#This operation may or may not be reduntant, but it is good to double-check anyway.
	#Dimension of the image
	xdimension=int(len(image[0]))
	ydimension=int(len(image))
	
	for Xstar,Ystar,index in zip(xpos,ypos,range(len(xpos))):
			#Check location to fall within acceptable values:
			checklocation=(Xstar > 2*fwhmnominal and Xstar < (xdimension-2*fwhmnominal) and  Ystar > 2*fwhmnominal and Ystar < (ydimension-2*fwhmnominal))  #Frame of 2*fwhm of acceptability for the location.
			if checklocation == False:
				#Remove the star if it is outside the image:
				np.delete(xpos, index)
				np.delete(ypos, index)
			
	mag,magerr,flux,fluxerr,sky,skyerr,badflag,outstr = aper.aper(image,xpos,ypos,phpadu=1,apr=5,zeropoint=25,skyrad=[40,50],badpix=[-12000,60000],exact=True, verbose=False)

	# use the stars at those coords to generate a PSF model (gauss)
	imagenameshort=imagename.split("/")[-1]
	gauss,psf,psfmag = getpsf.getpsf(image,xpos,ypos,mag,sky,1,1,np.arange(len(xpos)), 15,3, pathpsf + '/' + imagenameshort+'_psf_residuals.fits', verbose=False)
	
	#print 'gauss= ', gauss
	return gauss, psf, psfmag


###To query the Gaia catalog:
###info = Vizier.query_region(coord.SkyCoord(ra=50.0000, dec=-64.000,unit=(u.deg, u.deg)),radius='1m', catalog="I/337/gaia")
###info[0]["__Gmag_"]
