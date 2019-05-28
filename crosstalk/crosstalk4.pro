pro crosstalk4,ccdnumber,callsex,callimcopy,path_mary=path_mary,path_workspace=path_workspace,path_images=path_images,path_resamp=path_resamp,path_diff=path_diff,path_cat=path_cat,path_original=path_original,scale=scale,zpt=zpt,sequencenumber,field,date,useoldtemplate=useoldtemplate,codetemp=codetemp,codesci=codesci,datetemp=datetemp,filter=filter,satlevelchoice=satlevelchoice,fixsatlevel=fixsatlevel,pipesetup=pipesetup,usemarytemplate=usemarytemplate

; "imcopy" the original images to single files.
; This operation is necessary to have correct coordinates as output.
;
;Run SExtractor on all the original images to identify saturating sources.
;
;Find the XY position of the possible crosstalk and convert to WCS coordinates.
;
;create a file including the world coords of possible regions affected by crosstalk.
;
;Delete the imcopy-ed image to save space on disk. 

;;; Initialise the lists:
listCCD=[]   ; Useful to wrap up the results of all CCDs, but not necessary.
listRA=[]
listDEC=[]
listRADII=[]

;;; Initialise a list of saturating sources (different catalog)
listSATra=[]
listSATdec=[]

;;; Open the file that will include CCD#, RA, DEC, RADIUS of the regions to flag as possible crosstalk.

crosstalkname=path_cat + '/crosstalk_' + ccdnumber + '.cat'
OpenW, lun, /get_lun, crosstalkname
;;;;;;;;;;;;OpenW, lun2, /get_lun, saturatename  ;;;No list of saturating sources created here.


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; First attach to the (empty) list the position extracted from templates.
;;; This information can be useful to understand real time analysis anomalies.
;;; In case an old template is used, this operation will NOT be performed.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (useoldtemplate eq 'NO' and usemarytemplate eq 'NO') then begin

	for j=0, n_elements(codetemp)-1 do begin
       		
;;; imcopy the images (temporarily)
		     
		;;;Images downloaded from NOAO portal:     
     		if (pipesetup eq 'NOAO') then img=strcompress(path_original + '/' + 'c4d_' + datetemp + '_' + codetemp[j] + '_ooi_' + filter + '_v1.fits[' + ccdnumber + ']',/remove_all)
		;;;Images processed in real time:     
     		if (pipesetup eq 'RT') then img=strcompress(path_workspace + '/' + 'ut' + datetemp + '/' + ccdnumber + '/' + field + '.' + filter + '.ut' + datetemp + '.' + codetemp[j]+'_'+ccdnumber+'.fits', /remove_all)     
		
		imcopyfile=path_images + '/temporaryimg.fits'
		imcopycommand=callimcopy+' ' + img + ' ' + imcopyfile
		spawn, imcopycommand, imcoupyout
		
;;; Run SExtractor to find the sources flagged as saturating.
		catname=path_images + '/temporaryimg.cat'
		sexcommand= callsex +' '+ imcopyfile + ' -c ' + path_mary + '/config_crosstalk.sex -CATALOG_NAME ' + catname + ' -PARAMETERS_NAME '+path_mary + '/sex_crosstalk.param'
;;; Fix the saturation level in case the user does not want to use the one given in the header:
		if (satlevelchoice eq 'YES') then begin
			fixsatlevel1=string(fixsatlevel)
			sexcommand=sexcommand + ' -SATUR_LEVEL ' + fixsatlevel1 + ' -SATUR_KEY  AAAA '
		endif
;;; Spawn the SExtractor command:	
		spawn, sexcommand, sexout
		
		
;;; Select the saturating sources (flag = 4,5,6,7)		
		readcol, format='(F,F,D,D,I,F,F)', catname, Xsat,Ysat,RAsat,DECsat,FLAGS,KRON_RADIUS,A_IMAGE, /SILENT      
		
		satindex=where(flags ge 4 AND flags le 7)
;;;If there is any saturating source
		if (satindex[0] ge 0) then begin
		
;;; Read the image and extract astrometric information
			imgarr=readfits(imcopyfile, hdr, 1)
			naxis1=n_elements(imgarr[*,1])
			EXTAST, hdr, astr
			for z=0, n_elements(satindex)-1 do begin
				listCCD=[listCCD, ccdnumber]
;;;Compute x,y coordinates of the crosstalk
				xcross=(naxis1 / 2) + ((naxis1 / 2) - Xsat(satindex[z]))
				ycross=Ysat(satindex[z])
;;;Compute WCS coordinates of the crosstalk
				XY2AD, xcross,ycross,astr,RAcross,DECcross
;;;Update lists of coordinates and radii if the exclusion regions
				listRA=[listRA, RAcross]
				listDEC=[listDEC, DECcross]
				newradius=KRON_RADIUS(satindex[z])*A_IMAGE(satindex[z])
				listRADII=[listRADII,newradius]
;;; Update list of saturating sources
				newSATra=RAsat(satindex[z])
				newSATdec=DECsat(satindex[z])
				listSATra=[listSATra,newSATra]
				listSATdec=[listSATdec,newSATdec]
			endfor		
		
		
		endif 
			
;;; Remove the temporary image and the temporary catalog.
		spawn,'rm ' + imcopyfile, outrmimg
		spawn,'rm ' + catname, outrmcat
	
    	endfor


endif




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Find and list possible crosstalk in each "science" image.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




for j=0, n_elements(codesci)-1 do begin
       		
;;; imcopy the images (temporarily)
	
	;;;Images downloaded from NOAO portal:     
     	if (pipesetup eq 'NOAO') then img=strcompress(path_original + '/' + 'c4d_' + date + '_' + codesci[j] + '_ooi_' + filter + '_v1.fits[' + ccdnumber + ']',/remove_all)
	;;;Images processed in real time:     
     	if (pipesetup eq 'RT') then img=strcompress(path_workspace + '/' + 'ut' + date + '/' + ccdnumber + '/' + field + '.' + filter + '.ut' + date + '.' + codesci[j]+'_'+ccdnumber+'.fits', /remove_all)   

	;;;imcopyfile='$PBS_JOBFS/temporaryimg_'+field+'_'+date+'_m'+sequencenumber+'_'+ccdnumber+'.fits'
	imcopyfile=strcompress(path_images+'/temporaryimg_'+field+'_'+date+'_m'+sequencenumber+'_'+ccdnumber+'.fits', /remove_all)
	
	imcopycommand=callimcopy +' '+ img + ' ' + imcopyfile
	spawn, imcopycommand, imcoupyout
		
;;; Run SExtractor to find the sources flagged as saturating.
	;;;catname='$PBS_JOBFS/temporaryimg_'+field+'_'+date+'_m'+sequencenumber+'_'+ccdnumber+'.cat'
	catname=strcompress(path_images+'/temporaryimg_'+field+'_'+date+'_m'+sequencenumber+'_'+ccdnumber+'.cat', /remove_all)
	
	sexcommand= callsex +' '+ imcopyfile + ' -c ' + path_mary + '/config_crosstalk.sex -CATALOG_NAME ' + catname + ' -PARAMETERS_NAME '+path_mary + '/sex_crosstalk.param'
;;; Fix the saturation level in case the user does not want to use the one given in the header:
	if (satlevelchoice eq 'YES') then begin
		fixsatlevel1=string(fixsatlevel)
		sexcommand=sexcommand + ' -SATUR_LEVEL ' + fixsatlevel1 + ' -SATUR_KEY  AAAAA ' 
		
	endif
;;; Spawn the SExtractor command:
	spawn, sexcommand, sexout
		
;;; Select the saturating sources (flag = 4,5,6,7)		
	readcol, format='(F,F,D,D,I,F,F)', catname, Xsat,Ysat,RAsat,DECsat,FLAGS,KRON_RADIUS,A_IMAGE, /SILENT        
		
	satindex=where(flags ge 4 AND flags le 7)
;;;If there is any saturating source
	if (satindex[0] ge 0) then begin
		
;;; Read the image and extract astrometric information
		imgarr=readfits(imcopyfile, hdr, 1)
		naxis1=n_elements(imgarr[*,1])
		EXTAST, hdr, astr
		for z=0, n_elements(satindex)-1 do begin
			listCCD=[listCCD, ccdnumber]
;;;Compute x,y coordinates of the crosstalk
			xcross=(naxis1 / 2) + ((naxis1 / 2) - Xsat(satindex[z]))
			ycross=Ysat(satindex[z])
;;;Compute WCS coordinates of the crosstalk
			XY2AD, xcross,ycross,astr,RAcross,DECcross
;;;Update lists of coordinates and radii if the exclusion regions
			listRA=[listRA, RAcross]
			listDEC=[listDEC, DECcross]
			newradius=KRON_RADIUS(satindex[z])*A_IMAGE(satindex[z])
			listRADII=[listRADII,newradius]
			
;;; Update list of saturating sources
				newSATra=RAsat(satindex[z])
				newSATdec=DECsat(satindex[z])
				listSATra=[listSATra,newSATra]
				listSATdec=[listSATdec,newSATdec]
			
		endfor		
		
		
	endif 
	
	
	spawn,'rm ' + imcopyfile, outrmimg
	spawn,'rm ' + catname, outrmcat	
		
endfor


;;;convert the radii from PIXELS to ARCSEC
if (satindex[0] ge 0) then listRADIIarcsec=listRADII*scale


;;; If the lists are empty (no possible crosstalk) then print an empty space in the file.
if (listCCD eq !NULL) then printf, lun, ' '


;;; Print the lists:

if (listCCD ne !NULL  and  satindex[0] ge 0) then begin
	for t=0,n_elements(listCCD)-1 do begin
		printf, format='(A,D,D,D)', lun,listCCD[t], listRA[t], listDEC[t], listRADIIarcsec[t]
		;;;;;;;printf, format='(A,D,D,D)', lun2,listCCD[t], listSATra[t], listSATdec[t], listRADIIarcsec[t]
	endfor 
endif


close, lun
free_lun, lun

;;;;;;close, lun2
;;;;;;free_lun, lun2


;################################
; CASE OF OTHER MARY AS TEMPLATE
;################################

;;;;;;;;if (usemarytemplate eq 'YES') then begin

;;;;;;;;;;endif
print, '>>>>>>>>>>>>>>>>MADE IT TO THE END OF CROSSTALK'

end


