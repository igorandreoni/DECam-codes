import parameters as p
import psffunctions as pf
import lcmaker4_functions as fu 
import pdb
import sys
import os


#This program creates light curves of point sources
#Over an arbitrary number of images using image subtraction. 


#Debug
#pdb.set_trace()


# Initial statements to summarise the chosen parameters
print ' '

#In case the input is not valid:
if p.imaligned != 'NO' and p.imaligned != 'YES':
	print 'ERROR: Please enter correctly if the images are aligned or not (YES/NO) in the parameters file!'
	print ' '
	exit()

if p.imaligned == 'YES':
	print 'Images are already aligned and stored at ' + p.imalignedfolder
	print ' '
if p.imaligned == 'NO':
	print 'The images are not already aligned. They are located at ' + p.imfolder
	print 'After alignment they will be stored at ' + p.imalignedfolder

print 'You decided to ' + p.imsubchoice + ' the subtracted images.'
print ''
print 'The template images is ' + p.tmpl 
print ''


##########################
# Import useful packages
#########################
import numpy as np
import matplotlib.pyplot as plt
import subprocess
from subprocess import call
from astropy.io import fits
from astropy.table import Table
from PythonPhot import getpsf
from PythonPhot import aper
from PythonPhot import pkfit
import pyfits
import math
from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u

###########################
# Useful definitions
###########################
hotpants='/lustre/projects/p025_swin/iandreoni/hotpants/hotpants-master/hotpants'

swarp='/lustre/projects/p025_swin/pipes/arest/photpipe/Cfiles/bin/linux/swarp'
swarpconfigfile='/lustre/projects/p025_swin/iandreoni/LCMAKER4/config.swarp'

temp=p.tmpl
temproot=temp.split('/')[-1] #[0:-5]  this is for files finishing with .fits

############################################################################################
#First, create a fake base catalog, with just the RA,DEC coords of the listed candidate(s):
############################################################################################

candlist=Table.read(p.candlist,format='ascii',names=['ID','ccdnum','candnum','candRA','candDEC'])

#Assign to two variables the arrays with RA and DEC of candidates
IDcand=str(candlist["ID"][0])
RAcand=np.array(candlist["candRA"])
DECcand=np.array(candlist["candDEC"])

############################################################################################

#Open the log file
logfile = open(p.logfilename, 'a')

#Open the output file:
#Open the output file that will include the results:
outputfile = open(p.outputfilename + '_' + str(IDcand) + '.txt', 'w')
outputfile.write('# Magnitudes calibrated using the %s	catalog. %s'%(p.catmag, '\n' ))
outputfile.write('#%s %s %s %s %s %s %s %s %s %s %s'%('image','MJD','X', 'Y', 'MAG', 'EMAG', 'FLUX', 'EFLUX', 'ULMAG', 'ULFLUX', '\n' ))


################################
# Calibrate the template image
################################

imagename=temp
#Read the image as an array
hdulist = pyfits.open(imagename+'.fits')
image=hdulist[0].data
#Read the header
headerimage=hdulist[0].header
imagenameshort=imagename.split("/")[-1]


#Compute the dimension of the image
xdimension=int(len(image[0]))
ydimension=int(len(image))

#Compute the XY coordinates of the candidate in image
Xcand,Ycand=pf.wcstoxy(imagename+'.fits',RAcand,DECcand)

#If the candidate falls outside the image, then stop the program
outside=(Xcand <= p.fwhmnominal or Xcand >= (xdimension-p.fwhmnominal) or  Ycand <= p.fwhmnominal or Ycand >= (ydimension-p.fwhmnominal))
if outside:
	print 'The candidate is outside the template image! Stopping the program. '
	#Write the outcome in the log
	logfile.write('The candidate is outside the template image! Stopping the program. \n' )
	
	pdb.set_trace()
	#sys.exit("The candidate is outside the template image! Stopping the program. ")


#Extract the information about the date:
#########################MJD=headerimage["MJD-OBS"]
#CHANGE!!! Put a real MJD
MJD=0

#Prepare a list of stars to compute the PSF of the image, if the list is not given
if (p.myownlistofstars == False):
	
	#Create the list of stars for the template
	#depending on the choice of catalog (GSC, SDSS-DR12)
	if (p.whichcatforstars=='SDSS-DR12'):
		pf.makelistsdss(image, imagename,imagenameshort, RAcand, DECcand,pathstars=p.pathstars, magrange=p.magrange,radiussearch=str(p.radiussearch),fwhmnominal=p.fwhmnominal)
		
	if (p.whichcatforstars=='GSC'):
		pf.makelist(image, imagename,imagenameshort, RAcand, DECcand,pathstars=p.pathstars, magrange=p.magrange,radiussearch=str(p.radiussearch),fwhmnominal=p.fwhmnominal)
		
	#Safety check
	if (p.whichcatforstars != 'GSC' and p.whichcatforstars != 'SDSS-DR12'):
		print 'Select a valid catalog for the search for stars for the PSF extraction!'
		print 'In the parameters.py file: whichcatforstars is SDSS-DR12 or GSC'
		pdb.set_trace()
	
	listofstarsname=p.pathstars+'/' + imagenameshort +"_stars.txt"
if (p.myownlistofstars == True):
	listofstarsname=p.listofstarsname

#Skip the candidate if the star list is not produced. Add the info to a log file.
checkstars=os.popen('more '+ listofstarsname).read().splitlines()
if (checkstars == []):
	
	print 'The template has not enough good stars. Stopping the program. '
	#Write the outcome in the log
	logfile.write('The candidate is outside the template image! Stopping the program. \n' )
	
	
	#Remove the new extension temporarily created
	#####os.system("/bin/rm " + imagename + ".fits") 
	
	pdb.set_trace()
	#Remove the new extension temporarily created
	#os.system("/bin/rm " + imagename + ".fits") 
	

#Read the list of stars to use for PSF extraction
starslist = Table.read(listofstarsname, format='ascii') 
starslist.rename_column('col1','RAstar')
starslist.rename_column('col2','DECstar')

#Check again that enough reference stars are present
#Create a new (empty) list of stars
starslistnew = Table([[], [], [], []], names=('RAstar', 'DECstar','Xstar','Ystar'))
nremovedstars=0	
for rastar, decstar in zip(starslist["RAstar"], starslist["DECstar"]):
	#Find x,y coordinates
	xstar,ystar=pf.wcstoxy(imagename+'.fits',rastar,decstar)
		
	#Check location to fall within acceptable values:
	checklocation=(xstar > 2*p.fwhmnominal and xstar < (xdimension-2*p.fwhmnominal) and  ystar > 2*p.fwhmnominal and ystar < (ydimension-2*p.fwhmnominal))  #Frame of 2*fwhm of acceptabil
	if checklocation == False:
		nremovedstars=nremovedstars+1
		continue
	else:
		#Add the star (that is within the image borders) to the new list of stars
		starslistnew.add_row([rastar,decstar,xstar,ystar])
	
#Skip the image in case there are no enough good stars:
if (len(starslistnew["RAstar"]) ==0):
	print 'Not enough reference stars in the template, stopping the program.'
	
	#Write the outcome in the log
	logfile.write('Not enough reference stars in the template, stopping the program. \n' )

	#Remove the PSF model
	#os.system("/bin/rm " + pathpsf + "/*.fits") 

	#Remove the new extension temporarily created
	#os.system("/bin/rm " + imagename + ".fits") 
	pdb.set_trace()

#''''''''''''''''''''''''''''''''''''''''''''''

#Acquire information about the PSF of the image.
#The coordinates of the candidate will be appended at the end of the list.
gauss,psfresiduals,psfmag=pf.getPSFimage(image,imagename,starslistnew,Xcand,Ycand,pathpsf=p.pathpsf,zeropoint=p.zeropointvalue,skyRmin=p.skyRmin,skyRmax=p.skyRmax)

#Compute the FWHM of the template
#gauss[3]=sigma(X)
#gauss[4]=sigma(Y)
tempsigma=np.mean([gauss[3],gauss[4]])
fwhmtemplate=2.355*tempsigma


#Aperture photometry to get mags and sky values for specified coordinates
mag,magerr,flux,fluxerr,sky,skyerr,badflag,outstr = \
aper.aper(image,Xcand,Ycand,phpadu=p.phpadu,apr=fwhmtemplate,zeropoint=p.zeropointvalue, \
skyrad=[p.skyRmin,p.skyRmax],badpix=[-12000,60000],exact=True, verbose=False)


#Check that the PSF was successfully generated. If not, skip the image and report in the log.
lsPSFresult=os.popen("ls " + p.pathpsf+'/' +imagenameshort + '_psf_residuals.fits').read().splitlines()
if (lsPSFresult==[]):
	print 'PSF not successfully generated for the template. Stopping the program.'
	#Write the outcome in the log
	logfile.write('PSF not successfully generated for the template. Stopping the program. \n' )
	
	#Remove the new extension temporarily created
	#os.system("/bin/rm " + imagename + ".fits") 				
	pdb.set_trace()

hpsf = pyfits.getheader(p.pathpsf+'/' +imagenameshort + '_psf_residuals.fits')


#Perform the PSF fitting
if (p.recenter == "YES"): 
	#Load the pkfit class (here one can add mask and noise image)
	pk = pkfit.pkfit_class(image,gauss,psfresiduals,p.rnois,p.phpadu) #Here you can add noiseim, maskim (noise image, mask image)
	errmag,chi,sharp,niter,scale = pk.pkfit(1,Xcand,Ycand,sky,5) 
else:
	from PythonPhot import pkfit_norecenter
	#Load the pkfit class (here one can add mask and noise image)
	pk = pkfit_norecenter.pkfit_class(image,gauss,psfresiduals,p.rnois,p.phpadu) #Here you can add noiseim, maskim (noise image, mask image)
	errmag,chi,sharp,niter,scale = pk.pkfit_norecenter(1,Xcand,Ycand,sky,5) 

#Alternative versions: 
#pk.pkfit_norecenter
#pk.pkfit_noise
#pk.pkfit_norecent_noise 
#pk.pkfit_norecent_noise_smp

#Compute flux and magnitude
fluxtemplate = scale*10**(0.4*(p.zeropointvalue-hpsf['PSFMAG']))
dfluxtemplate = errmag*10**(0.4*(p.zeropointvalue-hpsf['PSFMAG']))
magvaluetemplate=p.zeropointvalue-2.5*np.log10(fluxtemplate)
#Estimate a conservative measure for the error
magerrortemplate=np.absolute(p.zeropointvalue-2.5*np.log10(fluxtemplate-dfluxtemplate) - magvaluetemplate)

#If the error of the magnitude is not a number, assign the 99 value.
if math.isnan(magerrortemplate):
	magerrortemplate=99


#Zeropoint calibration: Compute the offset to apply to the magnitudes
#using the USNO-B1 or Gaia catalogue.
#First initialise an arrey for the offsets
offsetarray=[]
		
#Initialise the number of stars deleted from the list because outside the image:

for rastar, decstar,xstar,ystar in zip(starslistnew["RAstar"], starslistnew["DECstar"],starslistnew["Xstar"],starslistnew["Ystar"]):
				
	if (p.catmag == 'USNO-B1'): 
		#Query the USNO-B1 catalogue
		usno = Vizier.query_region(coord.SkyCoord(ra=rastar, dec=decstar,unit=(u.deg, u.deg)),radius='3s',catalog="USNO-B1")
		if (usno != []):	
			#Check if the values are masked or not
			#If both are unmasked take the mean value
							
			if ( any(usno[0]["B2mag"]) == True  and any(usno[0]["B1mag"]) == True):
				USNOinfoB=np.mean([np.mean(usno[0]["B2mag"]),np.mean(usno[0]["B1mag"])])
			elif (any(usno[0]["B2mag"]) == True  and any(usno[0]["B1mag"]) == False):
				USNOinfoB=np.mean(usno[0]["B2mag"])
			elif (any(usno[0]["B2mag"]) == False  and any(usno[0]["B1mag"]) == True):
				USNOinfoB=np.mean(usno[0]["B1mag"])
			#Store the result in the variable good for every chosen catalog.	
			magcatalog=USNOinfoB
		#In case there in no catalog information, report the lack of info in the log and move on
		else:
			############Complete the log!!!
			continue

	if (p.catmag == 'GAIA'): 
	#Query the Gaia catalogue
		gaia = Vizier.query_region(coord.SkyCoord(ra=rastar, dec=decstar,unit=(u.deg, u.deg)),radius='1s', catalog="I/337/gaia")

		if (gaia != []):        
	        #Check if the values are masked or not
	        #If both are unmasked take the mean value		       
			if ( any(gaia[0]["__Gmag_"]) == True):
				gaiainfoG=np.mean(gaia[0]["__Gmag_"])
				#Store the result in the variable good for every chosen catalog.
				magcatalog=gaiainfoG
		else:
			print "No Gaia information for star at " + str(rastar) + ' ' + str(decstar)
			continue

	if (p.catmag == 'SDSS-DR12'): 
	#Query the Gaia catalogue
		sdss = Vizier.query_region(coord.SkyCoord(ra=rastar, dec=decstar,unit=(u.deg, u.deg)),radius='1s', catalog="V/147/sdss12")

		if (sdss != []):        
	        #Check if the values are masked or not
	        #If both are unmasked take the mean value		       
			if ( any(sdss[0]["gmag"]) == True):
				sdssinfoG=np.mean(sdss[0]["gmag"])
				#Store the result in the variable good for every chosen catalog.
				magcatalog=sdssinfoG
		else:
			print "No SDSS-DR12 information for star at " + str(rastar) + ' ' + str(decstar)
			continue


	#Perform the aperture and PSF photometry to get a measured value
		
	#Aperture photometry to get mags and sky values for specified coordinates.
	#Aperture photometry
	magstarap,magerrstarap,fluxstarap,fluxerrstarap,skystar,skyerrstar,badflagstar,outstrstar = \
	aper.aper(image,xstar,ystar,phpadu=p.phpadu,apr=fwhmtemplate,zeropoint=p.zeropointvalue, \
	skyrad=[p.skyRmin,p.skyRmax],badpix=[-12000,30000],exact=True, verbose=False)
		
	#PSF photometry
	pk = pkfit.pkfit_class(image,gauss,psfresiduals,p.rnois,p.phpadu) #Here you can add noiseim, maskim (noise image, mask image)
	errmag,chi,sharp,niter,scale = pk.pkfit(1,xstar,ystar,skystar,5) 
		
	#Compute flux and magnitude
	fluxstar = scale*10**(0.4*(p.zeropointvalue-hpsf['PSFMAG']))
	dfluxstar = errmag*10**(0.4*(p.zeropointvalue-hpsf['PSFMAG']))
	magvaluestar=p.zeropointvalue-2.5*np.log10(fluxstar)
	
	#Compute the offset and append to the main aray
	offsetarray.append(magcatalog-magvaluestar)


	
#Compute the offset to apply for the zeropoint correction of the image
offset=np.median(offsetarray)

########Very important: Zero point of the template, to be applied in all the images:
ZEROPOINT_TEMP=p.zeropointvalue+offset

magvaluecalibrated=magvaluetemplate+offset


#Is the magnitude value of the template an upper limit?
#Use the magerr<0.198 method for 5sigma upper limit
ulmagtemplate=0
ulfluxtemplate=0
ultemplate=False
#if (magerrortemplate >= 0.198):
if (magerrortemplate >= 0.17):
	#Assign new non-zero value to the upper limit
	ulmagtemplate=ZEROPOINT_TEMP-2.5*np.log10(5*skyerr*np.sqrt(3.14)*2*fwhmtemplate) 
	#Bring flux and flux error for the teplate to zero
	fluxtemplate=0 
	dfluxtemplate=0
	magvaluecalibrated=0
	magerrortemplate=0
	ulfluxtemplate=10**(0.4*(ZEROPOINT_TEMP-ulmagtemplate))
	ultemplate=True

#```````````````````````````````````````````````````````````````````````````````
#Print the result taking into account possible upper limits

if ultemplate:
	print('PSF fit for image %s gives an upper limit of mag %.3f'%(imagenameshort+'.fits',ulmagtemplate))
else:
	print('PSF fit for image %s gives mag %.3f +/- %.4f'%(imagenameshort+'.fits',magvaluecalibrated,magerrortemplate))
	print('PSF fit for image %s gives flux %.2f +/- %.3f'%(imagenameshort+'.fits',fluxtemplate,dfluxtemplate))
#print('PSF fit to coords %.2f,%.2f gives mag %s +/- %s'%(magvalue,errmag))
#print('PSF fit to coords %.2f,%.2f gives flux %s +/- %s'%(flux,dflux))

outputfile.write('%s %.8f %.2f %.2f %.3f %.4f %.2f %.3f %.3f %.2f %s'%(imagename+'.fits',MJD,Xcand, Ycand, magvaluecalibrated, magerrortemplate, fluxtemplate, dfluxtemplate, ulmagtemplate, ulfluxtemplate,'\n' ))

#Write the outcome in the log
logfile.write('%s	OK	%s'%(imagenameshort+'.fits', '\n' ))

#Report IF any reference star was deleted from the list because outside the image
totnumberstars=len(starslist["RAstar"])
if nremovedstars != 0:
	logfile.write('%s / %s	stars were deleted from the initial list.	%s'%(str(nremovedstars),str(totnumberstars), '\n' )) 


########################################




# Read the list of images
f = open(p.imlist, 'r')
	#imgtable= f.read()
imlist2=[]
for line in f:
	singleimg=line.strip()
	imlist2.append(singleimg)

f.close()

###############################
#Execute the full swarp+hotpants
#process for every image. In this
#Way even dithered images can be 
#accounted for.
###############################

for imagename in imlist2:

	#CHANGE!!
	MJD=MJD+1
	#Extract the information about the date for multiextension files:
	if p.multiextension:
		if p.fz:
			imagenametot=imagename + ".fits.fz"
		else:
			imagenametot=imagename + ".fits"
		hdulisttot = pyfits.open(imagenametot)
		headerimagetot=hdulisttot[0].header
		#CHANGE BACK FOR DECAM
		#########MJD=headerimagetot["MJD-OBS"]	
		#CHANGE!!!
		###########MJD=headerimagetot["MJD"]
	
#Case of compressed images
	CCDcand=candlist["ccdnum"][0]  ##########CHANGE: Number to go from 0 to a cycle for, for all the candidates in the list.
	if p.multiextension:
	
		if p.fz:
			#copy the extension from the compressed *fz in a separate FITS file, that will be delated at the end.
			os.system("/usr/local/x86_64/gnu/astrometry.net-0.50/bin/imcopy " + imagename + ".fits.fz["+str(CCDcand)+"] " + imagename + "_" + str(CCDcand) +".fits" )
			#Rename the variable with the name of the image
			imagename=imagename + "_" + str(CCDcand)
		else:			
			#copy the extension in a separate FITS file, that will be delated at the end.
			os.system("/usr/local/x86_64/gnu/astrometry.net-0.50/bin/imcopy " + imagename + ".fits["+str(CCDcand)+"] " + imagename + "_" + str(CCDcand) +".fits" )
			#Rename the variable with the name of the image
			imagename=imagename + "_" + str(CCDcand)



####################################################
# Does the candidate lye within the image?
####################################################


	#Read the image as an array
	hdulist = pyfits.open(imagename+'.fits')
	image=hdulist[0].data
	#Read the header
	headerimage=hdulist[0].header
	imagenameshort=imagename.split("/")[-1]
	
	#Extract the information about the date (if the file is NOT multiextension)
	
	
	#CHANGE!!!!!!!!!!
	if (p.multiextension == False):
	#########CHANGE BACK FOR DECAM
	##########	MJD=headerimage["MJD-OBS"]
	##########MJD=headerimage["MJD"]
		MJD=MJD+1

	#Compute the dimension of the image
	xdimension=int(len(image[0]))
	ydimension=int(len(image))

	#Compute the XY coordinates of the candidate in image
	Xcand,Ycand=pf.wcstoxy(imagename+'.fits',RAcand,DECcand)

	#If the candidate falls outside the image, then stop the program
	outside=(Xcand <= p.fwhmnominal or Xcand >= (xdimension-p.fwhmnominal) or  Ycand <= p.fwhmnominal or Ycand >= (ydimension-p.fwhmnominal))
	if outside:
		print 'Skipping image ' + imagename +'.fits, candidate outside the image'
		#Write the outcome in the log
		logfile.write('%s	Skipped: candidate outside the image.	%s'%(imagenameshort+'.fits', '\n' ))

		#Remove the new extension temporarily created, if multiple extensions are present:
		if p.multiextension:
			os.system("/bin/rm " + imagename + ".fits") 
		
		continue

#####################################
# Check list of stars
#####################################
	 
	#Check again that enough reference stars are present in each original image
	#Create a new (empty) list of stars
	starslistnew = Table([[], [], [], []], names=('RAstar', 'DECstar','Xstar','Ystar'))
	nremovedstars=0	
	for rastar, decstar in zip(starslist["RAstar"], starslist["DECstar"]):
		#Find x,y coordinates
		xstar,ystar=pf.wcstoxy(imagename+'.fits',rastar,decstar)
		
		#Check location to fall within acceptable values:
		checklocation=(xstar > 2*p.fwhmnominal and xstar < (xdimension-2*p.fwhmnominal) and  ystar > 2*p.fwhmnominal and ystar < (ydimension-2*p.fwhmnominal))  #Frame of 2*fwhm of acceptabil
		if checklocation == False:
			nremovedstars=nremovedstars+1
			continue
		else:
			#Add the star (that is within the image borders) to the new list of stars
			starslistnew.add_row([rastar,decstar,xstar,ystar])
	
#Skip the image in case there are no enough good stars:
	if (len(starslistnew["RAstar"]) ==0):
		print 'Skipping image ' + imagename +'.fits, not enough stars'

		#Write the outcome in the log
		logfile.write('%s	Skipped: not enough stars.	%s'%(imagenameshort+'.fits', '\n' ))
		
		if p.multiextension:
			#Remove the new extension temporarily created
			os.system("/bin/rm " + imagename + ".fits") 
		continue


#################################
#################################
# Run Swarp
#################################


# Go directly to the subtraction if the images are already aligned with the template
	if p.imaligned == 'YES':
		imaligned=p.imalignedfolder + '/' + imagenameshort + '.fits'

#Swarp (resample) the images with the template in case they are not aligned
	if p.imaligned == 'NO':
		sci=imagename + '.fits'
		
		swarpcommand= swarp + ' ' + sci + ' ' + temp+'.fits' + ' -c ' + swarpconfigfile + ' -IMAGE_SIZE ' + p.size +\
		' -RESAMPLE_DIR ' + p.imalignedfolder  + ' -CENTER_TYPE MOST '
		
		
		###' -CENTER_TYPE MANUAL -CENTER ' + " ' " +str(RAcand[0]) + "," +str(DECcand[0]) + " ' "   #doesn't align well images if candidate not at the centre.
		#####CHANGE important...RAcand is an array, you may want to insert the index
		 
		if p.doswarp == 'YES':
			subprocess.call(swarpcommand, shell=True)
# Name the resampled images
		imaligned=p.imalignedfolder + '/' + imagenameshort + '.resamp.fits'
		tempaligned=p.imalignedfolder + '/' + temproot + '.resamp.fits'         
# Now I have an image(imaligned).
#Remove the weightmaps if wresampchoice='DELETE'
		if p.wresampchoice=='DELETE':
			subprocess.call('rm ' + p.imalignedfolder +'/'+ imagenameshort + '.resamp.weight.fits', shell=True)

	
##########################
# Run HOTPANTS.
##########################
#		
# Now considering default upper and lower ADU limits
	ilvalue='-1000'
	tlvalue='-1000'
# Name of subtracted image
	imagesubname=p.imsubfolder + '/' + imagenameshort + '_sub.fits'
	
# Name of output fits file with kernel informatiom
	############outkernelname=p.imsubfolder + '/' + imagenameshort + '_ker.fits'
	outflagname=p.imsubfolder + '/' + imagenameshort + '_flag.fits'
	#Force the convolution of the template and the normalisation to the template gain value
	#  -oki    print extensive kernel info to output image header (0)
	#  -n t   normalize to (t)emplate, (i)mage, or (u)nconvolved (t)
	
	hotpantscommand=hotpants + ' -inim ' + imaligned + ' -tmplim ' + tempaligned + ' -outim ' + imagesubname  + ' -il ' + ilvalue + ' -tl ' + tlvalue + ' -oki 1 ' + ' -omi ' + outflagname + ' -n t '  
	if p.dohotpants == 'YES':
		subprocess.call(hotpantscommand, shell=True)
	if p.imsubchoice=='DELETE':
		subprocess.call('rm ' + p.imalignedfolder +'/'+ imagenameshort + 'resamp.fits', shell=True)

##################################
# Convert flagmap into weightmap
##################################
	if p.usedataqualitymap == 'YES':
		#Name the weightmap
		weightimagename=p.imsubfolder + '/' + imagenameshort + '_wmap.fits'
		#Call the function to convert a fagmap into a boolean weightmap
		fu.convertflag(outflagname,hdulistmask,weightimagename)
		
	
##################################	
##################################
# Perform FORCED PSF PHOTOMETRY
##################################	
##################################

	#Recap: 
	# - We calibrated the template image and obtained a zero point: ZEROPOINT_TEMP
	# - The image subtraction was performed normalising the values to the template.
	#   Thus we just need to use the new ZEROPOINT_TEMP for the future photometry.
	# - The name of the original image is: imagename  and we want to extract the PSF from it.
	# - To extract the PSF we can use the starslist of the template, but checking again for them to be inside the image.
	# - We want to performed forced photometry on the subtracted image (imagesubname)
	# - The location of the transient is RAcand, DECcand and XY values must be computed in the sub image. 
	

############################################################
# Check that the candidate lyes within the subtracted image
############################################################



	#Read the image as an array
	hdulistsub = pyfits.open(imagesubname)
	imagesub=hdulistsub[0].data
	#Read the header
	headerimagesub=hdulistsub[0].header


	#Compute the dimension of the image
	xdimensionsub=int(len(imagesub[0]))
	ydimensionsub=int(len(imagesub))

	#Compute the XY coordinates of the candidate in image
	Xcandsub,Ycandsub=pf.wcstoxy(imagesubname,RAcand,DECcand)

	#If the candidate falls outside the image, then stop the program
	outside=(Xcandsub <= p.fwhmnominal or Xcandsub >= (xdimensionsub-p.fwhmnominal) or  Ycandsub <= p.fwhmnominal or Ycandsub >= (ydimensionsub-p.fwhmnominal))
	if outside:
		print 'Skipping image ' + imagename +'.fits, candidate outside the SUBTRACTED image'
		#Write the outcome in the log
		logfile.write('%s	Skipped: candidate outside the SUBTRACTED image.	%s'%(imagenameshort+'.fits', '\n' ))
		
		#Remove the subtracted image
		os.system("/bin/rm " + imagesubname)
		
		#Remove the new extension temporarily created, if multiple extensions are present:
		if p.multiextension:
			os.system("/bin/rm " + imagename + ".fits") 
		
		continue

	 
###################################	 
# Extract PSF
###################################

	#To extract the PSF in the individual images we can take advantage of the same list
	#of stars used to calibrate the template, that have already been checked before 
	#running Swarp.
	
	#Acquire information about the PSF of the image.
	#The coordinates of the candidate will be appended at the end of the list.
	gauss,psfresiduals,psfmag=pf.getPSFimage(image,imagename,starslistnew,Xcand,Ycand,pathpsf=p.pathpsf, zeropoint=ZEROPOINT_TEMP,skyRmin=p.skyRmin,skyRmax=p.skyRmax)

	#Compute the FWHM of the image
	#gauss[3]=sigma(X)
	#gauss[4]=sigma(Y)
	imsigma=np.mean([gauss[3],gauss[4]])
	fwhmimage=2.355*imsigma


	#Aperture photometry to get mags and sky values for specified coordinates
	mag,magerr,flux,fluxerr,sky,skyerr,badflag,outstr = \
	aper.aper(image,Xcand,Ycand,phpadu=p.phpadu,apr=fwhmimage,zeropoint=ZEROPOINT_TEMP, \
	skyrad=[p.skyRmin,p.skyRmax],badpix=[-12000,60000],exact=True, verbose=False)


	#Check that the PSF was successfully generated. If not, skip the image and report in the log.
	lsPSFresult=os.popen("ls " + p.pathpsf+'/' +imagenameshort + '_psf_residuals.fits').read().splitlines()
	if (lsPSFresult==[]):
		#Write the outcome in the log
		logfile.write('%s	Skipped: no PSF generated.	%s'%(imagenameshort+'.fits', '\n' ))
		
		if p.multiextension:
			#Remove the new extension temporarily created
			os.system("/bin/rm " + imagename + ".fits") 				
		continue

	hpsf = pyfits.getheader(p.pathpsf+'/' +imagenameshort + '_psf_residuals.fits')

##########################################################
# Forced (or not) PSF photometry on the subtracted image
##########################################################

	#Aperture photometry to get mags and sky values for specified coordinates
	mag,magerr,flux,fluxerr,sky,skyerr,badflag,outstr = \
	aper.aper(imagesub,Xcandsub,Ycandsub,phpadu=p.phpadu,apr=fwhmimage,zeropoint=ZEROPOINT_TEMP, \
	skyrad=[p.skyRmin,p.skyRmax],badpix=[-12000,60000],exact=True, verbose=False)



	#Perform the PSF fitting
	if (p.recenter == "YES"): 
		#Load the pkfit class (here one can add mask and noise image)
		pk = pkfit.pkfit_class(imagesub,gauss,psfresiduals,p.rnois,p.phpadu) #Here you can add noiseim, maskim (noise image, mask image)
		errmag,chi,sharp,niter,scale = pk.pkfit(1,Xcandsub,Ycandsub,sky,5) 
	else:
		from PythonPhot import pkfit_norecenter
		#Load the pkfit class (here one can add mask and noise image)
		pk = pkfit_norecenter.pkfit_class(imagesub,gauss,psfresiduals,p.rnois,p.phpadu) #Here you can add noiseim, maskim (noise image, mask image)
		errmag,chi,sharp,niter,scale = pk.pkfit_norecenter(1,Xcandsub,Ycandsub,sky,5) 

	#Alternative versions: 
	#pk.pkfit_norecenter
	#pk.pkfit_noise
	#pk.pkfit_norecent_noise 
	#pk.pkfit_norecent_noise_smp

	#Compute flux (adding the flux of the template) and magnitude
	##########flux = fluxtemplate+scale*10**(0.4*(ZEROPOINT_TEMP-hpsf['PSFMAG']))
	flux = scale*10**(0.4*(ZEROPOINT_TEMP-hpsf['PSFMAG']))
	dflux = errmag*10**(0.4*(ZEROPOINT_TEMP-hpsf['PSFMAG']))
	#Correct the dflux for the flux error of the template (which is zero in case of upper limit in the template)
	
	#########dflux= np.sqrt( dflux**2 + dfluxtemplate**2 )
	
	magvalue=ZEROPOINT_TEMP-2.5*np.log10(flux)
	#Estimate a conservative measure for the error
	magerror=np.absolute(ZEROPOINT_TEMP-2.5*np.log10(flux-dflux) - magvalue)
	
	
	#If the error of the magnitude is not a number, assign the 0 value.
	if math.isnan(magerror):
		magerror=0
		
	#Is the magnitude value of the candidate an upper limit?
	#Use the magerr<0.198 method for 5sigma upper limit
	ulmag=0
	ulflux=0
	ul=False
	#Modified from magerror >= 0.198 to 1.7
	if (magerror >= 0.198):
		#Assign new non-zero value to the upper limit
		ulmag=ZEROPOINT_TEMP-2.5*np.log10(5*skyerr*np.sqrt(3.14)*2*fwhmimage) 
		#Bring flux and flux error to zero
		flux=0 
		dflux=0
		magvalue=0
		magerror=0
		ulflux=	10**(0.4*(ZEROPOINT_TEMP-ulmag))
		ul=True


		
#``````````````````````````````````````````````````````````````````````````````````````````````
	#Print the results taking into account the possibility of an upper limit
	if ul:
		print('PSF fit for image %s gives an upper limit at mag %.3f'%(imagenameshort+'.fits',ulmag))
	else:
		print('PSF fit for image %s gives mag %.3f +/- %.4f'%(imagenameshort+'.fits',magvalue,magerror))
		print('PSF fit for image %s gives flux %.2f +/- %.3f'%(imagenameshort+'.fits',flux,dflux))
	
	
	#print('PSF fit to coords %.2f,%.2f gives mag %s +/- %s'%(magvalue,errmag))
	#print('PSF fit to coords %.2f,%.2f gives flux %s +/- %s'%(flux,dflux))

	outputfile.write('%s %.8f %.2f %.2f %.3f %.4f %.2f %.3f %.3f %.2f %s'%(imagenameshort+'.fits',MJD,Xcand, Ycand, magvalue, magerror, flux, dflux, ulmag, ulflux, '\n' ))

	#Write the outcome in the log
	logfile.write('%s	OK	%s'%(imagenameshort+'.fits', '\n' ))

	#Report IF any reference star was deleted from the list because outside the image
	totnumberstars=len(starslist["RAstar"])
	if nremovedstars != 0:
		logfile.write('%s / %s	stars were deleted from the initial list.	%s'%(str(nremovedstars),str(totnumberstars), '\n' )) 


	#Remove the PSF model
	os.system("/bin/rm " + p.pathpsf+'/' +imagenameshort + '_psf_residuals.fits') 
	if p.multiextension:
	
	
		#Remove the new extension temporarily created
		os.system("/bin/rm " + imagename + ".fits") 

	
outputfile.close()
