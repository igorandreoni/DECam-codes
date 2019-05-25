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
The main code is lcmaker4.py.  More documentation will be provided. <br>
Based on the PythonPhot: https://github.com/djones1040/PythonPhot 


