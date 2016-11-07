"""
Script to make a 3-panel plot showing 3 postage stamp cutouts of a lens system
in different observed bands.  The cutouts are defined by an image size in
arcseconds and image centers, given as RA and Dec in decimal degrees.
"""

import numpy as n
import imfuncs as imf
import pyfits as pf
from matplotlib import pyplot as plt

try:
   infiles
except NameError:
   infiles = None

""" Set up the object names, coordinates, etc. """
try:
   lensroot
except NameError:
   lensroot = '1608'
try:
   inst
except NameError:
   inst = ['acs_c13','acs_c13','nic1_arch']
try:
   band
except NameError:
   band = [None,'f814w','f160w']
try:
   ralist
except NameError:
   ralist  = [242.308150, 242.308150, 242.308150]
try:
   declist
except NameError:
   declist = [65.5411111, 65.5411111, 65.5411111]
try:
   imsize
except NameError:
   imsize  = 5.0 # Image size in arcsec
""" 
Also put in upper limit of display range to make the images look nice.
The default display range is determined by finding the clipped mean, mu_clip,
and then setting min = mu_clip - 1 sigma_clip and max = mu_clip + 10 sigma_clip.
The values below override the 10x factor for the max value.
"""
try:
   sighi
except NameError:
   sighi = [15., 50., 10.]

try:
   outfile
except NameError:
   outfile = None

""" Set up plotting basics """
plt.figure(figsize=(9,3))
xsize = 0.3 # 3 images per row
ysize = 1   # 3 images per row
ycount = 1

""" Loop through the images """

for i in range(len(band)):

   """ If there is no image for this particular band, then skip this step """
   if band[i] is None:
      continue

   """ Open the image """
   if infiles != None:
      infile = infiles[i]
   else:
      infile = '%s_%s_%s.fits' % (lensroot,inst[i],band[i])
   f = imf.Image(infile)

   """ Define the image location """
   xstart = 0.0025 + 0.3 * i
   ystart = 0.0025
   imax = plt.axes([xstart,ystart,xsize,ysize])

   """ Make the png images """
   plt.axis('off')
   f.display(cmap='gaia',sighigh=sighi[i],subimdef='radec',
             subimcent=[ralist[i],declist[i]],subimsize=[imsize,imsize])


   """ Label the plots """
   plt.text(0.5,0.85,band[i].upper(),horizontalalignment='center',color='w',
            transform=imax.transAxes,fontsize=10)

   del f

plt.show()
if outfile != None: plt.savefig(outfile)
