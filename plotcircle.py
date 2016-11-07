import matplotlib
import matplotlib.pylab as pylab
import numpy as np
import math as m
import sys
execfile("/home/rumbaugh/FindCloseSources.py")

#Given a length, a start point, and the latitude(declination) of an end point,
#circFindLamda will calculate the RA of the end point. radius in degrees. 
#answer also in degrees. Keep in mind the answer is degenerate
def circFindLamda(startRA,startDec,endDec,rad):
    endRA =startRA*m.pi/180 + 2.0*m.asin(m.sqrt(((m.sin(m.pi*rad/360.0))**2-(m.sin((endDec-startDec)*m.pi/360.0))**2)/(m.cos(startDec*m.pi/180)*m.cos(endDec*m.pi/180))))
    endRA *= 180.0/m.pi
    return endRA

def plotcircle(RAcen=0.0,DECcen=0.0,rad=0.01,lstyle='dashed',pcolor='k',npoints=200):
    #rad should be in degrees
    if m.floor(npoints/2.0) != npoints/2: sys.exit("npoints must be even for plotcircle")
    xpoints,ypoints=np.zeros(npoints),np.zeros(npoints)
    theta = m.pi/2.0
    thstep = 2.0*m.pi/npoints
    xmax = circFindLamda(RAcen,DECcen,DECcen,rad)
    for ipc in range(0,npoints/2.0):
        ypoints[ipc] = 0.99999*rad*m.sin(theta) + DECcen
        ypoints[npoints-1-ipc] = ypoints[ipc]
        eRA = circFindLamda(RAcen,DECcen,ypoints[ipc],rad)
        xpoints[ipc] = eRA
        xpoints[npoints-1-ipc] = 2*RAcen-eRA
        theta -= thstep
    pylab.plot(xpoints,ypoints,color=pcolor,linestyle=lstyle)
