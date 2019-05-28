
'''
Run SExtractor on all the original images to identify saturating sources.

Find the XY position of the possible crosstalk and convert to WCS coordinates.

create a file including the world coords of possible regions affected by crosstalk.
'''

import numpy as np
import subprocess
from subprocess import call
import os
from astropy.io import ascii
from astropy.io import fits
from astropy import wcs
from astropy.table import Table
import pdb

def flag_crosstalk(path_images, image_name, callsex, path_config, set_saturation = True, saturation_adu = 5000, scale = 0.2637):

    # Initialise the lists:
    listRA = []
    listDEC = []
    listRADII = []

    # Initialise a list of saturating sources (different catalog)
    listSATra = []
    listSATdec = []
		
    # Run SExtractor to find the sources flagged as saturating.
    catname=crosstalk_filename=f"{path_images}/{image_name.replace('.fits','_xtalk_temp.cat')}"
	
    sexcommand= f"{callsex} {path_images}/{image_name} -c {path_config}/config_crosstalk.sex -CATALOG_NAME {catname} -PARAMETERS_NAME {path_config}/sex_crosstalk.param"

    # Fix the saturation level in case the user does not want to use the one given in the header:
    if set_saturation:
        sexcommand=sexcommand + f" -SATUR_LEVEL {saturation_adu} -SATUR_KEY  NONE "
    else:
        sexcommand=sexcommand + " -SATUR_KEY SATURATION "

    # Call the SExtractor command:
    os.system(sexcommand)
		
    # Select the saturating sources (flag = 4,5,6,7)		
    saturation_tbl = ascii.read(catname, format='sextractor') 
    saturation_tbl = saturation_tbl[saturation_tbl['FLAGS'] >= 4]
    saturation_tbl = saturation_tbl[saturation_tbl['FLAGS'] <= 7]    
	
    #Return if there are no saturating sources
    if len(saturation_tbl) == 0:
        print('No saturating sources. The crosstalk flagging file will be empty.')
        return

    # Read in the image (fits) and get exact size information
    image_hdul = fits.open(f"{path_images}/{image_name}")
    image_data = image_hdul[0].data   
    NAXIS2, NAXIS1 = len(image_data[:,0]), len(image_data[0,:])
    
    # Get the WCS solution
    w = wcs.WCS(image_hdul[0].header)
    
    # Initialize new table for the crosstalk
    crosstalk_tbl = Table([[],[],[],[],[]], names = ('RA_sat','Dec_sat', 'RA_xtalk', 'Dec_xtalk', 'radius')) 

    # For each saturating star, compute the location of the crosstalk
    for Xsat, Ysat, RAsat, DECsat, kr, a_image in zip(saturation_tbl['X_IMAGE'], saturation_tbl['Y_IMAGE'], saturation_tbl['X_WORLD'], saturation_tbl['Y_WORLD'], saturation_tbl['KRON_RADIUS'], saturation_tbl['A_IMAGE']):
        xcross = (NAXIS1 / 2.) + ((NAXIS1 / 2.) - Xsat)
        ycross = Ysat
        world_coords_xtalk = w.wcs_pix2world(np.array([[xcross, ycross]]), 1)
        RAcross, DECcross = world_coords_xtalk[0][0], world_coords_xtalk[0][1]
        newradius = kr * a_image
        # Add a row to the output table, converting the radius from px to arcsec
        crosstalk_tbl.add_row([RAsat, DECsat, RAcross, DECcross, newradius*scale])
        print(Xsat, Ysat, xcross,ycross)

    # File that will include RA, Dec, RADIUS of the regions to flag as possible crosstalk, along with bright sources.
    crosstalk_filename=f"{path_images}/{image_name.replace('.fits','_crosstalk.csv')}"
    crosstalk_tbl.write(crosstalk_filename, format='csv')

    return crosstalk_tbl


if __name__ == '__main__':
    path_images = 'test'
    image_name = '4hr.g.ut151219.503748_23.fits'
    callsex = '/Users/igor/Software/ureka/Ureka/bin/sex'
    path_config = '.'
    crosstalk = flag_crosstalk(path_images, image_name, callsex, path_config, set_saturation = True, saturation_adu = 5000, scale = 0.2637)
