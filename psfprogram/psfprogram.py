#This program is based on the PythonPhot package,
#that perform DAOPHOT-type PSF photometry.
#More info here: https://github.com/djones1040/PythonPhot/tree/master/PythonPhot

#This version relates to NOAO-portal images, that come as multi-extension files.
#This version is accounts for multiple targets at the same time, too.

#This program first creates (or reads) a list of stars suitable for the PSF photometry. 
#If a list is not provided, it queries the Vizier database within a user-defined radius. 
#Then PSF photometry is performed at the coordiantes of the candidate. 

#This program goes with the "psffunctions.py" auxiliary program, that contains
#a set of useful functions.

#Preparation:
#  1. Prepare a list of images named "list_images.txt". 
#     All the images must be listed here with full path (unless they are present in the
#     directory where the program is called from) and WITHOUT the .fits extension.
#
#  2. Edit the "parameters to set manually" section below. They are individually explained.
#
#  3. Call the program from Python or simply as: python psfprogram.py 
#


from PythonPhot import getpsf
from PythonPhot import aper
from PythonPhot import pkfit
import pyfits
import numpy as np
from astroquery.vizier import Vizier
from astropy.table import Table
import astropy.coordinates as coord
import astropy.units as u
import pdb
import os
import math

#For the zeropoint calibration:
from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PARAMETERS TO SET MANUALLY:


#Name of the log file
logfilename="log_psfprogram_SN19941.log"


#List of candidates: ccdNumber (I), RA (D), DEG (D) 
listcandname='list_candidates_SN19941.txt'

#Name of the list with input images
#NOTE: the images mast include full path and NO extention.
#For example, /my/path/image.fits becomes:  /my/path/image 
listname='list_images_SN19941_z.txt' 

#Are the images multi-extension images? (True/False). 
#If True the CCD number will correspond to the extension number.
multiextension=False

#Are the images compressed as *.fz files? (True/False)
fz=False

#Name of the output file
outputfilename='results'  # The output will become: outputfilename_ID.txt

#Path where to store the PSF models and the lists of stars:
pathpsf='./psfstore'
pathstars='./starstore'

#---------- STARS ----------
#Do you want to give a defined list of stars for the PSF fitting? (True/False)
myownlistofstars=True

#Path to the list of stars
listofstarsname="./list_stars_SN19941.txt"
#Note: the list must be in the form "RA	DEC" with values in degrees.

#Do you want the first list of stars to be kept the same? (True/False)
keepfirstlist=False
#Note: This functionality allows the light curves to be
#generated faster and more self-consistent.  

#Radius (in arcmin) for the query of stars suitable for the PSF extraction
radiussearch=1.5  #arcmin


#Do you want to calibrate to DECam magnitudes? (True/False)
#If True, then no calibration using catalogs will be performed.
calibDECam=False


#Which catalog to extract the magnitude of the stars from? ('USNO-B1'/'GAIA')
#Note: this will not be active if calibDECam=True
catmag='USNO-B1'

#Range of F filter magnitudes (red) good to select stars
#suitable for the PSF extraction (basically, not too bright, not too faint)
magrange=[16.5,19.5]

#-------------------------------
#Recentering during the PSF fitting? ("YES"/"NO")
recenter="NO"

#filter 'g'/'r' ('u','i','z','Y' only if calibDECam=True)
filt='r'

#Zero point:
zeropointvalue=25.0

#Rough value for the FWHM (in pixel)
fwhmnominal=3.0  #pixel, mainly used to establish the distance of the source from the border.

#Maximum distance for coincident sources (radius, in arcsec)
radmatch=1.2

#Photons per ADU:
phpadu=120.4

#Readout noise per pixel (scalar)
rnois=5.

# END OF PARAMETERS TO SET MANUALLY:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ACTION!

#Import the functions
import psffunctions as pf


#Read the list of candidates
listcand = Table.read(listcandname, format='ascii',delimiter=' ')
listcand.rename_column('col1','ID')
listcand.rename_column('col2','CCD')
listcand.rename_column('col3','RA')
listcand.rename_column('col4','DEC')

#Sort the table with increasing CCD number (useful for when two candidates are on the same CCD)
#listcand.sort('CCD')


#Read the list of images
with open(listname) as f:
    imagenamelist = f.read().splitlines()

#Convert the search radius to STRING
radiussearch=str(radiussearch)

#For cycle over the candidates
for IDcand, CCDcand, RAcand, DECcand in zip (listcand['ID'],listcand['CCD'],listcand['RA'],listcand['DEC']):
	
	#Open the log file
	logfile = open(logfilename, 'a')

	#Open the output file that will include the results:
	outputfile = open(outputfilename + '_' + str(IDcand) + '.txt', 'w')
	###outputfile.write('# Magnitudes calibrated using the %s	catalog. %s'%(catmag, '\n' ))
	outputfile.write('#%s	%s	%s	%s	%s	%s	%s	%s	%s	%s'%('image','MJD','X', 'Y', 'MAG', 'EMAG', 'FLUX', 'EFLUX', 'ZPT', '\n' ))
	
	#Write the name of the candidate
	logfile.write('candidate %s	%s'%(IDcand, '\n' ))
	
	
	#Minimise the queries. Create a dictionary with RA, DEC, and catalogue info for each 
	#reference star.		
	tstars = Table([[], [], []], names=('RA', 'DEC', 'catalog'))
	
	##############CHANGE
	MJD0=0
	#For cycle over the images
	for imagename in imagenamelist:

		#```````````````````````````````````````````
		#Extract the information about the date:
	
		#Read the image as an array
		if (fz == True):
			hdulist0 = pyfits.open(imagename+'.fits.fz')
		else:
			hdulist0 = pyfits.open(imagename+'.fits')
		#Read the header
		headerimage0=hdulist0[0].header
		#CHANGE
		MJD=headerimage0["MJD-OBS"]
		###########MJD=MJD0+1   ########headerimage["MJD-OBS"]
		#########i=i+1	
		
		
		#```````````````````````````````````````````
	
		if multiextension:
			#copy the extension in a separate FITS file, that will be delated at the end.
			if (fz == True):
				#CHANGE!!!!!!!!!!!
				###os.system("/usr/local/x86_64/gnu/astrometry.net-0.50/bin/imcopy " + imagename + ".fits.fz["+str(CCDcand)+"] " + imagename + "_" + str(CCDcand) +".fits" )
				os.system("/usr/local/x86_64/gnu/astrometry.net-0.50/bin/imcopy " + imagename + ".fits.fz["+str(CCDcand)+"] ./tempimage_"+ str(CCDcand) +".fits" )
				imagename="./tempimage"
			else:
				os.system("/usr/local/x86_64/gnu/astrometry.net-0.50/bin/imcopy " + imagename + ".fits["+str(CCDcand)+"] " + imagename + "_" + str(CCDcand) +".fits" )
		
			#Rename the variable with the name of the image
			imagename=imagename + "_" + str(CCDcand)
		
		#Read the image as an array
		hdulist = pyfits.open(imagename+'.fits')
		
		image=hdulist[0].data	
		#HACK: the values are too low. Multiply times 10**12
		image=image*10**12
		
		#Read the header
		headerimage=hdulist[0].header

		#Extract the name of the image without full path
		imagenameshort=imagename.split("/")[-1]


		#Compute the dimension of the image
		xdimension=int(len(image[0]))
		ydimension=int(len(image))
		
		#Compute the XY coordinates of the candidate in image
		Xcand,Ycand=pf.wcstoxy(imagename+".fits",RAcand,DECcand)
		
		#If the candidate falls outside the image, then skip the image
		outside=(Xcand <= fwhmnominal or Xcand >= (xdimension-fwhmnominal) or  Ycand <= fwhmnominal or Ycand >= (ydimension-fwhmnominal))
		if outside:
			print 'Skipping image ' + imagename +'.fits, candidate outside the image'
			#Write the outcome in the log
			logfile.write('%s	Skipped: candidate outside the image.	%s'%(imagenameshort+'.fits', '\n' ))
	
			#Remove the new extension temporarily created
			if multiextension:
				os.system("/bin/rm " + imagename + ".fits") 
			
			continue
		
		
		#Prepare a list of stars to compute the PSF of the image, if the list is not given
		if (myownlistofstars == False):
			
			#Create the list of stars just for the first image if the "keep the first list of stars" option is activated.
			if (keepfirstlist and (imagename == imagenamelist[0]+ "_" + str(CCDcand))):  
				pf.makelist(image, imagename,imagenameshort, RAcand, DECcand,pathstars=pathstars, magrange=magrange,radiussearch=radiussearch,fwhmnominal=fwhmnominal)
				listofstarsname=pathstars+'/' + imagenameshort +"_stars.txt"
				
			#Create every time a new list of stars if the "keep the first list of stars" option is NOT activated.
			elif (keepfirstlist == False):
				pf.makelist(image, imagename,imagenameshort, RAcand, DECcand,pathstars=pathstars, magrange=magrange,radiussearch=radiussearch,fwhmnominal=fwhmnominal)
				listofstarsname=pathstars+'/' + imagenameshort +"_stars.txt"
		
		
		#Skip the candidate if the star list is not produced. Add the info to a log file.
		checkstars=os.popen('more '+listofstarsname).read().splitlines()
		if (checkstars == []):
			
			print 'Skipping image ' + imagename +'.fits, not enough stars'
			
			#Write the outcome in the log
			logfile.write('%s	Skipped: not enough stars.	%s'%(imagenameshort+'.fits', '\n' ))
		
			#Remove the PSF model
			os.system("/bin/rm " + pathpsf + "/*.fits") 

			#Remove the new extension temporarily created
			if multiextension:
				#CHANGE!!!!!!!!!
				####os.system("/bin/rm " + imagename + ".fits") 
				os.system("/bin/rm ./tempimage_"+ str(CCDcand) +".fits") 
			continue
		
		#Read the list of stars to use for PSF extraction
		starslist = Table.read(listofstarsname, format='ascii') 
		starslist.rename_column('col1','RAstar')
		starslist.rename_column('col2','DECstar')
		
		#Check again that enough reference stars are present
		#Create a new (empty) list of stars
		starslistnew=t = Table([[], [], [], []], names=('RAstar', 'DECstar','Xstar','Ystar'))
		nremovedstars=0	
		for rastar, decstar in zip(starslist["RAstar"], starslist["DECstar"]):
			#Find x,y coordinates
			xstar,ystar=pf.wcstoxy(imagename + '.fits',rastar,decstar)
				
			#Check location to fall within acceptable values:
			checklocation=(xstar > 2*fwhmnominal and xstar < (xdimension-2*fwhmnominal) and  ystar > 2*fwhmnominal and ystar < (ydimension-2*fwhmnominal))  #Frame of 2*fwhm of acceptability for the location.
			if checklocation == False:
				nremovedstars=nremovedstars+1
				continue
			else:
				#Add the star (that is within the image borders) to the new list of stars
				starslistnew.add_row([rastar,decstar,xstar,ystar])
			
		#Skip the image in case there are not enough good stars:
		if (len(starslistnew["RAstar"]) ==0):
			print 'Skipping image ' + imagename +'.fits, not enough reference stars.'
			
			#Write the outcome in the log
			logfile.write('%s	Skipped: not enough reference stars.	%s'%(imagenameshort+'.fits', '\n' ))
		
			#Remove the PSF model
			os.system("/bin/rm " + pathpsf + "/*.fits") 

			#Remove the new extension temporarily created
			if multiextension:
				os.system("/bin/rm " + imagename + ".fits") 
			continue
		
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		
		#Acquire information about the PSF of the image.
		#The coordinates of the candidate will be appended at the end of the list.
		aux=True
		aux2=0
		while aux==True:
			aux2=aux2+1	
			try:		
				gauss,psfresiduals,psfmag=pf.getPSFimage(image,imagename,starslistnew,Xcand,Ycand,pathpsf=pathpsf)
				break
			
			except ValueError:
				logfile.write('%s       Skipped: no gaussian computed.      %s'%(imagenameshort+'.fits', '\n'))
                        	#Remove the new extension temporarily created
				if multiextension:
					os.system("/bin/rm " + imagename + ".fits")
				aux=False
			if aux2==2:
				#Break after 2 attempts
				aux=False
		#Skip image if aux==False
		if aux==False:
			continue

		#Compute the real FWHM of the image
		#gauss[3]=sigma(X)
		#gauss[4]=sigma(Y)
		sigmaimage=np.mean([gauss[3],gauss[4]])
		fwhmimage=2.355*sigmaimage
		
		#Check that the FWHM is physical
		if fwhmimage < 0.3:
			logfile.write('%s       Skipped: bad FWHM computed.      %s'%(imagenameshort+'.fits', '\n'))
			#Remove the new extension temporarily created
			if multiextension:
                        	os.system("/bin/rm " + imagename + ".fits")
                        continue
		#Define the sky radius
		skyRmin = 3*fwhmimage
		skyRmax = 5*fwhmimage
		
		#Aperture photometry to get mags and sky values for specified coordinates
		mag,magerr,flux,fluxerr,sky,skyerr,badflag,outstr = \
		aper.aper(image,Xcand,Ycand,phpadu=phpadu,apr=fwhmimage,zeropoint=zeropointvalue, \
				
		skyrad=[skyRmin,skyRmax],badpix=[-12000,60000],exact=True, verbose=False)
		
		
		#Check that the PSF was successfully generated. If not, skip the image and report in the log.
		lsPSFresult=os.popen("ls " + pathpsf+'/' +imagenameshort + '_psf_residuals.fits').read().splitlines()
		if (lsPSFresult==[]):
			#Write the outcome in the log
			logfile.write('%s	Skipped: no PSF generated.	%s'%(imagenameshort+'.fits', '\n' ))
			
			#Remove the new extension temporarily created
			if multiextension:
				os.system("/bin/rm " + imagename + ".fits") 				
			continue
		
		hpsf = pyfits.getheader(pathpsf+'/' +imagenameshort + '_psf_residuals.fits')
		
	
		#Perform the PSF fitting
		if (recenter == "YES"): 
			#Load the pkfit class (here one can add mask and noise image)
			pk = pkfit.pkfit_class(image,gauss,psfresiduals,rnois,phpadu) #Here you can add noiseim, maskim (noise image, mask image)
			errmag,chi,sharp,niter,scale = pk.pkfit(1,Xcand,Ycand,sky,5) 
		else:
			from PythonPhot import pkfit_norecenter
			#Load the pkfit class (here one can add mask and noise image)
			pk = pkfit_norecenter.pkfit_class(image,gauss,psfresiduals,rnois,phpadu) #Here you can add noiseim, maskim (noise image, mask image)
			errmag,chi,sharp,niter,scale = pk.pkfit_norecenter(1,Xcand,Ycand,sky,5) 
		
		#Alternative versions: 
		#pk.pkfit_norecenter
		#pk.pkfit_noise
		#pk.pkfit_norecent_noise 
		#pk.pkfit_norecent_noise_smp
		
		#Compute flux and magnitude
		flux = scale*10**(0.4*(zeropointvalue-hpsf['PSFMAG']))
		dflux = errmag*10**(0.4*(zeropointvalue-hpsf['PSFMAG']))
		magvalue=zeropointvalue-2.5*np.log10(flux)
		#Estimate a conservative measure for the error
		magerror=np.absolute(zeropointvalue-2.5*np.log10(flux-dflux) - magvalue)
		
		#If the error of the magnitude is not a number, assign the 0 value.
		if math.isnan(magerror):
			magerror=0




#``````````````````````````````````````````
		#Calibration to DECam magnitudes.
		#This step is optional but, if chosen (calibDECam=True), no zeropoint
		#calibration with the USNO/GAIA catalogs will be performed.
		
		if (calibDECam==True):
			tabDECam=Table.read('/home/iandreon/psfprogram/standard_star_photometry_DECam.txt', format='ascii') # to be moved outside the for cycle
			
			#Get CCDNUM and airmass values from the headers (indivisual CCD and full image)
			CCDNUMimage=headerimage["CCDNUM"]
			airmass=headerimage0["AIRMASS"] 
			exptime=headerimage0["EXPTIME"]
						
			#Get the coefficients for the standard star photometry
			checktable=0
			for ccdnumtable,filtable,atable,aerrtable,btable,berrtable,ktable,kerrtable in zip(tabDECam["CCDNUM"],tabDECam["filter"],tabDECam["a"],tabDECam["aerr"],tabDECam["b"],tabDECam["berr"],tabDECam["k"],tabDECam["kerr"]):
				if ( filt == filtable and CCDNUMimage == ccdnumtable):
					checktable=checktable+1
					aimage=atable
					aerrimage=aerrtable
					bimage=btable
					berrimage=berrtable
					kimage=ktable
					kerrimage=kerrtable
			#Verify that everything is OK
			if (checktable != 1):
				print '>>> There is a problem in associating the CCD to the standard star photometric calibration!'
				pdb.set_trace()
			
			#Calculate the magnitude with the formula
			#Main formula: g = g_instr - a_g - b_g*( (g-r) - (g-r)0 ) - k_g * X
			#http://www.ctio.noao.edu/noao/content/Mean-Photometric-Standard-Star-Module-PSM-Solutions-mean-zeropoints-color-terms-extinctions#g-band%20fit
			#No color term for now.
			magvaluecalibrated=magvalue-25.-aimage-kimage*airmass + 2.5*np.log10(exptime)
			
			#zeropoint for DECam calib (check if it is correct)
			zptDECam=-aimage-kimage*airmass + 2.5*np.log10(exptime)
 			
			#Correct the errorbar:
			magerror=magerror+aerrimage+(kerrimage*airmass)
			
			print('PSF fit for image %s gives mag %.3f +/- %.4f'%(imagenameshort+'.fits',magvaluecalibrated,magerror))
			print('PSF fit for image %s gives flux %.2f +/- %.3f'%(imagenameshort+'.fits',flux,dflux))
			#print('PSF fit to coords %.2f,%.2f gives mag %s +/- %s'%(magvalue,errmag))
			#print('PSF fit to coords %.2f,%.2f gives flux %s +/- %s'%(flux,dflux))
	
			outputfile.write('%s	%.8f	%.2f	%.2f	%.3f	%.4f	%.2f	%.3f	%.2f	%s'%(imagenameshort+'.fits',MJD,Xcand, Ycand, magvaluecalibrated, magerror, flux, dflux, zptDECam,'\n' ))
		
			#Write the outcome in the log
			logfile.write('%s	OK	%s'%(imagenameshort+'.fits', '\n' ))
		
			#Report IF any reference star was deleted from the list because outside the image
			totnumberstars=len(starslist["RAstar"])
			if nremovedstars != 0:
				logfile.write('%s / %s	stars were deleted from the initial list.	%s'%(str(nremovedstars),str(totnumberstars), '\n' )) 
		
			if multiextension:
				#Remove the PSF model
				os.system("/bin/rm " + pathpsf + "/*.fits") 
		
				#Remove the new extension temporarily created
				#CHANGE!!!!!!!!!
				####os.system("/bin/rm " + imagename + ".fits") 
				os.system("/bin/rm ./tempimage_"+ str(CCDcand) +".fits") 

			
			continue








		
		
#```````````````````````````````````````````
		#Zeropoint calibration: Compute the offset to apply to the magnitudes
		#using the USNO-B1 or GAIA catalogue.
		#First initialise an arrey for the offsets
		offsetarray=[]
				
		#Initialise the number of stars deleted from the list because outside the image:

		for rastar, decstar,xstar,ystar in zip(starslistnew["RAstar"], starslistnew["DECstar"],starslistnew["Xstar"],starslistnew["Ystar"]):

			
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#Catalogs query	
			#If the table with catalogued values is not empty, check if the star was already queried
			if ( len(tstars['RA']) != 0 ):
				match=False #Initialisation
				for ratable,dectable,cattable in zip(tstars['RA'],tstars['DEC'], tstars["catalog"]):
					dist=np.sqrt((rastar-ratable)**2 + (decstar-dectable)**2)
					if (dist <= radmatch/3600.0):  #Radius converted to degrees
						match=True
						magcatalog=cattable
						break
						
			#If the star was not queried, yet, then perform the query	 
			if ( len(tstars['RA']) == 0 ) or ( match==False ): 
				if (catmag == 'USNO-B1'): 
					#Query the USNO-B1 catalogue
					usno = Vizier.query_region(coord.SkyCoord(ra=rastar, dec=decstar,unit=(u.deg, u.deg)),radius='3s',catalog="USNO-B1")
					if (usno != []   and   filt == 'g'):	
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
						
					if (usno != []   and   filt == 'r'):
						if ( any(usno[0]["R2mag"]) == True  and any(usno[0]["R1mag"]) == True):
							USNOinfoR=np.mean([np.mean(usno[0]["R2mag"]),np.mean(usno[0]["R1mag"])])
						elif (any(usno[0]["R2mag"]) == True  and any(usno[0]["R1mag"]) == False):
							USNOinfoR=np.mean(usno[0]["R2mag"])
						elif (any(usno[0]["R2mag"]) == False  and any(usno[0]["R1mag"]) == True):
							USNOinfoR=np.mean(usno[0]["R1mag"])
						#Store the result in the variable good for every chosen catalog.	
						magcatalog=USNOinfoR
						
					#In case there in no catalog information, report the lack of info in the log and move on
					if (usno == [] ):
						############Complete the log!!!
						print 'No catalog information!!!'
						#pdb.set_trace()
						continue
				
				if (catmag == 'GAIA'): 
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
		               
				
				#Add the queried results to the table of queried stars
				tstars.add_row([rastar,decstar,magcatalog])	
				
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
			
			#Perform the aperture and PSF photometry to get a measured value
				
			#Aperture photometry to get mags and sky values for specified coordinates.
			#Aperture photometry
			magstarap,magerrstarap,fluxstarap,fluxerrstarap,skystar,skyerrstar,badflagstar,outstrstar = \
			aper.aper(image,xstar,ystar,phpadu=phpadu,apr=fwhmimage,zeropoint=zeropointvalue, \
			skyrad=[skyRmin,skyRmax],badpix=[-12000,30000],exact=True, verbose=False)
				
			#PSF photometry
			pk = pkfit.pkfit_class(image,gauss,psfresiduals,rnois,phpadu) #Here you can add noiseim, maskim (noise image, mask image)
			errmag,chi,sharp,niter,scale = pk.pkfit(1,xstar,ystar,skystar,5) 
				
			#Compute flux and magnitude
			fluxstar = scale*10**(0.4*(zeropointvalue-hpsf['PSFMAG']))
			dfluxstar = errmag*10**(0.4*(zeropointvalue-hpsf['PSFMAG']))
			magvaluestar=zeropointvalue-2.5*np.log10(fluxstar)
			
			#Compute the offset and append to the main aray
			offsetarray.append(magcatalog-magvaluestar)

			
		#Compute the offset to apply for the zeropoint correction of the image
		offset=np.median(offsetarray)
		
		
		magvaluecalibrated=magvalue+offset
		
#``````````````````````````````````````````````````````````````````````````````````````````````
	
		print('PSF fit for image %s gives mag %.3f +/- %.4f'%(imagenameshort+'.fits',magvaluecalibrated,magerror))
		print('PSF fit for image %s gives flux %.2f +/- %.3f'%(imagenameshort+'.fits',flux,dflux))
		#print('PSF fit to coords %.2f,%.2f gives mag %s +/- %s'%(magvalue,errmag))
		#print('PSF fit to coords %.2f,%.2f gives flux %s +/- %s'%(flux,dflux))
	
		outputfile.write('%s	%.8f	%.2f	%.2f	%.3f	%.4f	%.2f	%.3f	%.2f	%s'%(imagenameshort+'.fits',MJD,Xcand, Ycand, magvaluecalibrated, magerror, flux, dflux, zeropointvalue+offset, '\n' ))
		
		#Write the outcome in the log
		logfile.write('%s	OK	%s'%(imagenameshort+'.fits', '\n' ))
		
		#Report IF any reference star was deleted from the list because outside the image
		totnumberstars=len(starslist["RAstar"])
		if nremovedstars != 0:
			logfile.write('%s / %s	stars were deleted from the initial list.	%s'%(str(nremovedstars),str(totnumberstars), '\n' )) 
		
		if multiextension:
			#Remove the PSF model
			os.system("/bin/rm " + pathpsf + "/*.fits") 
		
			#Remove the new extension temporarily created
			os.system("/bin/rm " + imagename+".fits") 
		
	logfile.close()
		
	outputfile.close()

