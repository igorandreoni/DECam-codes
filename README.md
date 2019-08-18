# DECam python codes

## match_db_decam 

usage: match_db_decam.py [-h] [--c RADEC RADEC] [--r MATCH_RADIUS] [--n NAME] <br>

Match coordinates or retrieve info of candidates stored on the DECam database. <br>

optional arguments:
<ul>
   -h, --help        show this help message and exit <br>
   --c RADEC RADEC   Coordinates of the object you want to find. RA Dec (deg). <br>                    
   --r MATCH_RADIUS  Match radius (arcsec) <br>
   --n NAME          Name of the candidate that you want to find, for example DGabcd <br>
</ul>

  
<p> Examples: </p> 
    <p> match_db_decam.py --c 13.23422 -67.32677 </p>
    <p> match_db_decam.py --c 13.23422 -67.32677 --r 5 </p>
    <p> match_db_decam.py --n DG19bexl </p>


## LCMAKER4: Forced PSF photometry with image subtraction
The main code is lcmaker4.py <br>  
More documentation will be provided. <br>
The image subtraction is performed with HOTPANTS and the PSF photometry is based on PythonPhot: https://github.com/djones1040/PythonPhot 

## psfprogram: PSF photometry 
The main code is psfprogram.py <br>  
More documentation will be provided. <br>
PSF photometry is based on PythonPhot: https://github.com/djones1040/PythonPhot 

## Crosstalk
The main code is crosstalk.py <br>  
The code identifies saturating sources and flags those regions of the CCD where crosstalk ghost images are likely to be found.  It also indicates the radius (in arcsec) of the region that should be flagged as possible crosstalk. 
