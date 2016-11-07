import numpy as np
import math as m
import sys

def lum2lum(lum=1.0,kt=3.0,band='0.5-2.0'):
    #This is a conversion to python of Dale's lum2lum.pro script
    #which calculated bolometric luminosity for a raymond-smith
    #profile based on temperature and luminosity in a given band
    #Right now, this only works with the 0.5-2.0 band
    if band == '0.5-2.0':
        kfile = '/home/rumbaugh/kcor.0.5-2.0'
        lfile = '/home/rumbaugh/lum2lum.0.5-2.0_bolometric.dat'
    else:
        sys.exit("Input band must be 0.5-2.0 for now")

    cr = read_file(lfile)
    lcor = get_colvals(cr,'col1')
    dkt = 0.25
    kt = (np.arange(80)/4.0)+0.25
    conv2bol = interp1(kt,lcor,[t])
 
    return, lum*conv2bol
