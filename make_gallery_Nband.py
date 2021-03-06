"""
Script to make a N-panel plot showing N postage stamp cutouts of a lens system
in different observed bands.  The cutouts are defined by an image size in
arcseconds and image centers, given as RA and Dec in decimal degrees,
except for the NIRC2 data in which case the centers are defined in terms
of (x,y).
"""

import numpy as n
import imfuncs as imf
import pyfits as pf
from matplotlib import pyplot as plt

""" Set up some shared basic paramters """
nband = 4
imsize = 5.0
panelsize = 2.0 # Size of an individual panel, in inches

""" Set up the lens-based paramters"""
lensinfo = n.array([
#   ('SDSS_J0252+0039','J0252',xx,xx,20.)
   ('B0128+437','B0128',
    ('wfpc2','wfpc2','nic2','nirc2_n'),
    (None,'F814W','F160W','Kp'),
    (22.8059167,22.806102,22.806374,300.),
    (43.9702833,43.970116,43.970403,300.),
    (20.,10.,40.,20.)),
   ('HE0435-1223','HE0435',
    ('acs','acs','nic2','nirc2_n'),
    ('f555w','f814w','F160W','Kp'),
    (69.5619625,69.5619625,69.561607,300.),
    (-12.2874889,-12.2874889,-12.287333,300.),
    (20.,20.,20.,12.)),
   ('B0445+123','B0445',
    ('acs','acs','nic2','nirc2_n'),
    ('f555w','f814w','F160W','Kp'),
    (72.0918708,72.0918708,72.0918708,300.),
    (12.4655528,12.4655528,12.4655528,300.),
    (20.,40.,40.,20.)),
   ('B0631+519','B0631',
    ('acs','acs','nic2','nirc2_n'),
    (None,'f814w','F160W','Kp'),
    (98.8013083,98.8013083,98.8013083,300.),
    (51.9505008,51.9505008,51.9505008,300.),
    (20.,40.,40.,20.)),
   ('B0712+472','B0712',
    ('wfpc2','wfpc2','nic1','nirc2_n'),
    (None,'F814W',None,'Kp'),
    (109.014904,109.014904,109.014904,300.),
    (47.1472639,47.1472639,47.1472639,300.),
    (20.,10.,40.,20.)),
   ('MG0751+2715','MG0751',
    ('wfpc2','wfpc2','nic2','nirc2_n'),
    (None,'F814W','F160W','Kp'),
    (117.922750,117.922750,117.922750,300.),
    (27.275375,27.275375,27.275375,300.),
    (20.,10.,40.,20.)),
   ('SDSS_J0924+0219','J0924',
    ('acs','acs','nic2','nirc2_n'),
    ('f555w','f814w',None,'Kp'),
    (141.23246,141.23246,141.23246,300.),
    (2.32349,2.32349,2.32349,300.),
    (50.,50.,50.,20.))
#   ('RX_J1131-1231','J1131',172.9644583,-12.5323444,20.)
   ],
   dtype=[('fullname','S20'), ('root','S20'), ('inst','S20',nband),
          ('band','S5',nband), ('ra',float,nband), ('dec',float,nband), 
          ('sighi',float,nband)]
   )

lensinfo = lensinfo.view(n.recarray)


""" Set things up for 3 panels per lens, for now """
plt.figure(figsize=(nband*panelsize,lensinfo.size*panelsize))
xsize = 1.0/nband - 0.005
ysize = 1.0/lensinfo.size - 0.005

ycount = 1

print lensinfo.size
for i in range(lensinfo.size):

   """ Loop through the images """
   
   print ''
   for j in range(nband):

      """ Define the image location """
      xstart = 0.0025 + 1.0 * (j % nband) / nband
      ystart = 1.0025 - 1.0 * ycount / lensinfo.size
      if (j % nband) == nband-1:
         ycount += 1
   
      """ If there is no image for this particular band, then skip this step """
      if lensinfo.band[i][j] is None or lensinfo.band[i][j] == 'None':
         continue
   
      """ Open the image """
      indir = 'Data/%s/%s/Cutouts/%s' % (lensinfo.inst[i][j].upper(),
                                         lensinfo.band[i][j],
                                         lensinfo.fullname[i])
      infile = '%s/%s_%s_%s_6x6.fits' % (indir,lensinfo.root[i],
                                         lensinfo.inst[i][j],
                                         lensinfo.band[i][j])
      f = imf.Image(infile)
   
      """ Make the png images """
      if lensinfo.inst[i][j] == 'nirc2_n':
         imdef = 'xy'
         imsizen = imsize * 100.
         band = lensinfo.band[i][j]
      else:
         imdef = 'radec'
         imsizen = imsize
         band = lensinfo.band[i][j].upper()
      print ''
      imax = plt.axes([xstart,ystart,xsize,ysize])
      f.display(cmap='gaia',sighigh=lensinfo.sighi[i][j],subimdef=imdef,
                subimcent=[lensinfo.ra[i][j],lensinfo.dec[i][j]],
                subimsize=[imsizen,imsizen])
   
   
      """ Label the plots """
      bandlab = '%s/%s' % (lensinfo.inst[i][j].upper(),band)
      plt.text(0.05,0.9,lensinfo.root[i],horizontalalignment='left',
               color='w',transform=imax.transAxes,fontsize=10)
      plt.text(0.95,0.9,bandlab,horizontalalignment='right',
               color='w',transform=imax.transAxes,fontsize=10)
   
      plt.axis('off')
      del f

#plt.show()
plt.savefig('multiband_test.pdf')
