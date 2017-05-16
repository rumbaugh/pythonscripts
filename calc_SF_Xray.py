import numpy as np
import scipy.signal as signal

def make_tophat(R,xsize,ysize):
    xcen,ycen=0.5*(xsize-1),0.5*(ysize-1)
    tophat_arr=np.zeros((ysize,xsize))
    xcoord,ycoord=np.arange(xsize)*np.ones(ysize).reshape((ysize,1)),np.ones(xsize)*np.arange(ysize).reshape((ysize,1))
    tophat_arr[np.sqrt((xcoord-xcen)**2+(ycoord-ycen)**2)<=R]=1
    return tophat_arr

def calc_SF_Xray(InArray,smooth='tophat',smooth_param=[5]):
    xsize,ysize=np.shape(InArray)[1],np.shape(InArray)[0]
    if smooth=='tophat':
        R=smooth_param[0]
        smooth_arr=make_tophat(R,R,R)/(np.pi*R**2)
    else:
        print "smooth currently only accepts 'tophat'"
        return
    SFarr=signal.convolve2d(InArray,smooth_arr,mode='same')
    return SFarr
    
