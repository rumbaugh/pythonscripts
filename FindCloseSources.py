import numpy as np
import math as m
import sys

minsrchrad = 5

def SphDist(RAi,Deci,Raf,Decf):
	#Distance in arcminutes
	dist = 2.0*m.asin(m.sqrt((m.sin(0.5*m.pi*(Decf-Deci)/180))**2+m.cos(m.pi*Decf/180.0)*m.cos(m.pi*Deci/180.0)*(m.sin(0.5*m.pi*(Raf-RAi)/180))**2))
	dist *= 180.0*60/m.pi
	return dist

def WindingNum(RA,Dec,boundRAs,boundDecs):
	#Computes winding number by counting number of quadrants moved through counterclockwise
	NumCCW = 0
	for i in range(0,len(boundRAs)):
		if boundRAs[i] > RA:
			if boundDecs[i] > Dec:
				CQ = 1
			else: 
				CQ = 4
		else:
			if boundDecs[i] > Dec:
				CQ = 2
			else:
				CQ = 3
		if i != 0:
			DelQ = CQ-PQ
			if m.fabs((m.fabs(DelQ)-2)) < 0.1:
				testDec = (boundDecs[i]-boundDecs[i-1])*(RA-boundRAs[i-1])/(boundRAs[i]-boundRAs[i-1])+boundDecs[i-1]
				if testDec < Dec:
					if ((PQ == 2) or (PQ == 3)):
						NumCCW += 2
					else: 
						NumCCW -= 2
				else: 
					if ((PQ == 4) or (PQ == 1)):
						NumCCW += 2
					else:
						NumCCW -= 2
			else:
				if (m.fabs(DelQ)) > 2.9: DelQ /= -3
				NumCCW += DelQ
		PQ = CQ
	return NumCCW

def FindCloseSources(ra,dec,tol,ra_opt,dec_opt,usemin):
	#finds sources from the ra_opt and dec_opt catalogs close to the single
	# RA and Dec given by the variables ra and dec
	if ((tol < minsrchrad) & (usemin != 0)): tol = minsrchrad
	tolx = tol/m.cos(dec*m.pi/180.0)
	ra_box_temp_ind1 = np.where((ra_opt >= ra-tolx/3600.0) & (ra_opt <= ra+tolx/3600.0))
	if len(ra_box_temp_ind1) > 0: 
		ra_box_temp_ind1 = ra_box_temp_ind1[0]
		ra_box_temp1 = ra_opt[ra_box_temp_ind1]
		dec_box_temp1 = dec_opt[ra_box_temp_ind1]
		dec_box_temp_ind2 = np.where((dec_box_temp1 >= dec-tol/3600.0) & (dec_box_temp1 <= dec+tol/3600.0))
		if len(dec_box_temp_ind2) > 0:
			dec_box_temp_ind2 = dec_box_temp_ind2[0]
			dec_box_temp3 = dec_box_temp1[dec_box_temp_ind2]
			ra_box_temp3 = ra_box_temp1[dec_box_temp_ind2]
			temp_dist_ar = np.zeros(len(dec_box_temp3))
			for i in range(0L,len(temp_dist_ar)):
				temp_dist_ar[i] = SphDist(ra_box_temp3[i],dec_box_temp3[i],ra,dec)
			inside_tol_ind = np.where(temp_dist_ar*60.0 <= tol)
			if len(inside_tol_ind) > 0: 
				inside_tol_ind = inside_tol_ind[0]
				if len(inside_tol_ind) > 0:
					return ra_box_temp_ind1[[dec_box_temp_ind2[inside_tol_ind]]]
				else:
					return np.zeros(0)
			else:
				return np.zeros(0)
		else:
			return np.zeros(0)
	else:
		return np.zeros(0)

def mkinArr(imgrid,by,bx):
	#When loading from the boundary RA and Dec (actually image X and Y), 
	#but whatever, the second column(Y) is by. This is because python
	# reverses rows and columns or something like that. The output g can
	#be applied right to a table made with get_piximgvals used on an image
	#crate. The indices should be right. Don't forget that the bx and by
	#need to be close(1st and last index should be the same)
	try:
		MIAstep
	except NameError:
		MIAstep = 5
	try:
		ymaxmin
	except NameError:
		ymaxmin = 15
    	ishp = imgrid.shape
    	inArr = np.zeros((ishp[0],ishp[1]))
    	maxy = ishp[0]
    	miny = 0
    	for x in range(0,ishp[1]):
       	    y = maxy + MIAstep
       	    if ((y < ymaxmin) or (x < 2)): y = ishp[0]
       	    if y > ishp[0]: y = ishp[0]
       	    ccw = 0
       	    while ((ccw == 0) & (y > 0)):
       	    	  y -= 1
       	     	  ccw = WindingNum(x,y,bx,by)
            maxy,ccw = y,0
       	    y = miny - MIAstep
       	    if ((y < 0) or (x < 2)): y = -1
       	    while ((ccw == 0) & (y < maxy)):
       	     	  y += 1
       	     	  ccw = WindingNum(x,y,bx,by)
            miny = y
       	    if maxy > 1: inArr[miny:maxy+1,x] = 1
    	g = np.where(inArr != 0)
	return g
