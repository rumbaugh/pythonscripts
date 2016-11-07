#!/usr/bin/env python

import sys,pyfits as pf

input = sys.argv[1]
x = eval(sys.argv[2])
y = eval(sys.argv[3])
box = eval(sys.argv[4])

cut = pf.getdata(input)[y-box-1:y+box,x-box-1:x+box]

fits = pf.PrimaryHDU(cut)
h = fits.header
d = fits.data

h2 = pf.getheader(input)
try:
    test = h2['cd1_1']
    h.update('wcsaxes',2)
    h.update('ctype1','RA---TAN')
    h.update('ctype2','DEC--TAN')
    h.update('cd1_1',h2['cd1_1'])
    h.update('cd1_2',h2['cd1_2'])
    h.update('cd2_1',h2['cd2_1'])
    h.update('cd2_2',h2['cd2_2'])
    h.update('crval1',h2['crval1'])
    h.update('crval2',h2['crval2'])
    h.update('crpix1',h2['crpix1']-x+box+1)
    h.update('crpix2',h2['crpix2']-y+box+1)
    h.update('ltv1',h2['ltv1'])
    h.update('ltv2',h2['ltv2'])
    h.update('ltm1_1',h2['ltm1_1'])
    h.update('ltm2_2',h2['ltm2_2'])
except:
    pass

try:
    #fits.writeto('%s.fits'%sys.argv[5],clobber=True)
    pf.writeto('%s.fits'%sys.argv[5],d,h,clobber=True)
except:
    #fits.writeto('stamp.fits',clobber=True)
    pf.writeto('stamp.fits',d,h,clobber=True)
