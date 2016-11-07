import os
import numpy as np
import time
import gc

def FindPeaks(imgvals):
    imshp = imgvals.shape
    peaks = np.zeros(0)
    xs,ys = np.zeros(0),np.zeros(0)
    g = np.where(imgvals > (0.66*imgvals.max()+0.34*imgvals.mean()))
    for i in range(0,len(g[0][:])):
    	x = g[0][i]
	y = g[1][i]
	t = imgvals[x][y]
    	sxmax = x + 20
    	if sxmax >= imshp[0]: imshp[0]-1
    	sxmin = x - 20
   	if sxmin < 0: sxmin = 0
        symax = y + 20
        if symax >= imshp[1]: symax = imshp[1]-1
        symin = y-20
        if symin < 0: symin=0
        tt = imgvals[sxmin:sxmax,symin:symax]
        if t >= tt.max(): 
            peaks,xs,ys = np.append(peaks,t),np.append(xs,x),np.append(ys,y)
    return peaks,xs,ys
  
def FindPeaksAll(imgvals):
    imshp = imgvals.shape
    peaks = np.zeros(0)
    xs,ys = np.zeros(0),np.zeros(0)
    for x in range(10,imshp[0]-10):
        for y in range(10,imshp[1]-10):
            t = imgvals[x][y]
            sxmax = x + 20
            if sxmax >= imshp[0]: imshp[0]-1
            sxmin = x - 20
            if sxmin < 0: sxmin = 0
            symax = y + 20
            if symax >= imshp[1]: symax = imshp[1]-1
            symin = y-20
            if symin < 0: symin=0
            tt = imgvals[sxmin:sxmax,symin:symax]
            if t >= tt.max(): 
                peaks,xs,ys = np.append(peaks,t),np.append(xs,x),np.append(ys,y)
    return peaks,xs,ys
  
