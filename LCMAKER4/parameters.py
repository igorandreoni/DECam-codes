# Set the useful parameters to generate differential photometry lightcurves.

#are the images already aligned with the template, and of the same size? ('YES'/'NO')
imaligned='NO'

#Do you want to run Swarp? ('YES'/'NO')
#This is a shortcut to avoid re-aligning without changing the images names.
doswarp='YES'

#Folder for the NOT aligned images (discard this entry if the images are already aligned) 
imfolder='/lustre/projects/p025_swin/iandreoni/LCMAKER4/DWF17a3962systembis17/IM'

#Folder for the aligned images (template icluded)
imalignedfolder='/lustre/projects/p025_swin/iandreoni/LCMAKER4/DWF17a3962systembis17/IMALIGN'

#List of images on which perform the differential photometry (do NOT include the template)
#These images will be added *resamp.fits automatically in case they are not already aligned.
imlist='/lustre/projects/p025_swin/iandreoni/LCMAKER4/DWF17a3962systembis17/list_images.txt'

#Template image (must finish without .fits)
tmpl='/lustre/projects/p025_swin/pipes/DWF_PIPE/MARY_WORK/Antlia_170202_msystembis17_1/ccd39/images_resampled/temp_39.resamp'

# Size (a little less) of the non-aligned images?
size='3900,1900'

#Folder for the subtracted images
imsubfolder='/lustre/projects/p025_swin/iandreoni/LCMAKER4/DWF17a3962systembis17/SUB'

#Folder for the catalogs
psffolder='/lustre/projects/p025_swin/iandreoni/LCMAKER4/DWF17a3962systembis17/PSF'

#Name of the list of candidates (that must include: ccdnum,candnum,ra,dec)
candlist='/lustre/projects/p025_swin/iandreoni/LCMAKER4/DWF17a3962systembis17/cand.txt'

#Do you want to run HOTPANTS? ('YES'/'NO')
dohotpants='YES'

#Do you want to create and use data quality maps? (YES/NO)
usedataqualitymap='NO' 

# Do you want to keep or delete the subtracted images? ('KEEP'/'DELETE')
imsubchoice='KEEP'

# Do you want to keep or delete the weightmaps created during the resampling? ('KEEP'/'DELETE')
wresampchoice='DELETE'


#Are the images multi-extension images? (True/False). 
#If True the CCD number will correspond to the extension number.
multiextension=False

#Are the images compressed in *fits.fz format?
fz=False

#Name of the log file
logfilename="DWF17a3962systembis17/log_DWF17a3962systembis17.log"

#Name of the output file
outputfilename='results'  # The output will become: outputfilename_ID.txt


#Path where to store the PSF models and the lists of stars:
pathpsf='./psfstore'
pathstars='./starstore'

#---------- STARS ----------
#Do you want to give a defined list of stars for the PSF fitting? (True/False)
myownlistofstars=False

#If automatic search for stars, do you want to select them using the GSC or SDSS-DR12 catalog? ('GSC'/'SDSS-DR12')
whichcatforstars='GSC'

#Path to the list of stars
listofstarsname="./DWF17a3962systembis17/list_stars.txt"
#Note: the list must be in the form "RA	DEC" with values in degrees.

#Do you want the first list of stars to be kept the same? (True/False)
keepfirstlist=True
#Note: This functionality allows the light curves to be
#generated faster and more self-consistent.  

#Radius (in arcmin) for the query of stars suitable for the PSF extraction
radiussearch=4  #arcmin

#Which catalog to extract the magnitude of the stars from? ('USNO-B1'/'GAIA','SDSS-DR12')
catmag='USNO-B1'

#Range of F filter magnitudes (red) good to select stars
#suitable for the PSF extraction (basically, not too bright, not too faint)
magrange=[17,20]
#-------------------------------
#Recentering during the PSF fitting? ("YES"/"NO")
recenter="NO"

#Zero point:
zeropointvalue=25.0

#Sky annulus (pix)
skyRmin=15
skyRmax=25

#Rough value for the FWHM (in pixel)
fwhmnominal=4.8  #pixel

#Maximum distance for coincident sources (radius, in arcsec)
radmatch=1.0

#Photons per ADU:
phpadu=2.4

#Readout noise per pixel (scalar)
rnois=5.9 #Rough re



#Name of the output file
outputfilename='results'  # The output will become: outputfilename_ID.txt

#Path where to store the PSF models and the lists of stars:
pathpsf='./psfstore'
pathstars='./starstore'




# SExtractor threshold and parameters

DETECT_MINAREA = '8'             # minimum number of pixels above threshold
THRESH_TYPE='RELATIVE'           # threshold type: RELATIVE (in sigmas)
                                 # or ABSOLUTE (in ADUs)
DETECT_THRESH = '1.1'           # <sigmas> or <threshold>,<ZP> in mag.arcsec-2
ANALYSIS_THRESH = '1.1' 	 # <sigmas> or <threshold>,<ZP> in mag.arcsec-2

seeing=1.0 				#FWHM in arcsec
zeropoint=25
