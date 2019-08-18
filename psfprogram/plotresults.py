#This program plots the light curves resulting from the psfprogram.py

import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import pdb
import os


#Function fo plot the results
def plotlc(resultsname, saveplot=False):

	tbl=ascii.read(resultsname)

	#########################################################
	#Go through the whole table and divide detections
	#from upper limits (5sigma, emag>0.198).
	#The upper limits are defined as the magnitude value minus its error.

	#Initialise the arrays
	magarr=[]
	emagarr=[]
	ularr=[]

	#Also initialise arrays for the plot axis limits
	ymin0=[]
	ymax0=[]
	
	for mag,emag in zip(tbl["MAG"],tbl["EMAG"]):
		detection=(emag <= 0.198 and emag !=0)
	
		if detection:
			magarr.append(mag)
			emagarr.append(emag)
			ularr.append(0)

		else:
			magarr.append(0)
			emagarr.append(0)
			ularr.append(mag-emag)
	
		ymin0.append(mag+emag)
		ymax0.append(mag-emag)

	#Limits for the y axis
	
	ymin=np.amax(ymin0)
	ymax=np.amin(ymax0)

	########################################################

	#Prepare the time axis.
	#If all the observations happen on the same night,
	#then minutes timescale. Otherwise MJD.

	time=[]

	if ((tbl["MJD"][-1]-tbl["MJD"][0]) < 1):
		MJD0=tbl["MJD"][0]
		for MJD in tbl["MJD"]:
			time.append(  (MJD-MJD0)*24*60  )
			timename="Minutes since MJD="+str(MJD0)	
	else:
		MJD0=int(tbl["MJD"][0])-1  
		for MJD in tbl["MJD"]:
			time.append(  (MJD-MJD0)  )
			timename="MJD-"+str(MJD0)	

	#Set the limits of the plot axis
	xmin=time[0]-1
	xmax=time[-1]+1

##########################################################



	#Plot
	plt.figure(figsize=(15,10))
	plt.errorbar(time, magarr, yerr=emagarr, fmt="bo")
	#plt.plot(time, magarr, lw=1,color='black')
	plt.plot(time, magarr, 'bo', color='black')
	plt.plot(time, ularr, "rv")
	plt.xlabel(timename, fontsize=22)
	plt.ylabel('g mag', fontsize=22)
	plt.axis([xmin,xmax,ymin,ymax])
	plt.tick_params(labelsize=18)
	titlename0=resultsname.replace(".txt", "")
	titlename1=titlename0.replace("results_", "")
	#########plt.title(titlename1)

	if saveplot:
		return plt.savefig("png/" + resultsname.replace(".txt", "_lc.png"))
	
	else: 
		return plt.show()



resultnames=os.popen('ls diff_results_SN27920*nr_r_60d*').read().splitlines()

for resultname in resultnames:
	plotlc(resultname,saveplot=False)

