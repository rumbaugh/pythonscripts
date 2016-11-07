"""
ccdredux.py - A library of functions to do various basic CCD image processing
              operations

High-level Functions:
   make_bias         - combines input bias or dark frames into one master file
   make_flat         - combines input flat frames into a master flat
   make_flat_files   - exactly like make_flat, but input is a file
                        containing a list of flat-field files rather than
                        an array of frame numbers
   median_combine    - does a generic median combination of a list of images
   apply_calib       - applies calibration corrections to input files
   hdr_offsets       - uses header information in the input files to derive 
                        initial guesses for the offsets between dithered 
                        exposures on a field.
   xcorr_offsets     - does a fft-based cross-correlation to get the initial
                        estimate of the offsets between dithered exposures
                        on a field

Low-level Functions:
   define_trimsec    - sets the section of an image to trim
   sigma_clip        - sigma-clips a data array
"""

import pyfits as pf
import scipy as sp
import numpy as n
import imfuncs as imf
import wcs, coords
from math import cos,sin,pi,sqrt,atan2

#-----------------------------------------------------------------------

def sigma_clip(data,nsig=3.,verbose=False):
   # Only compute outputs for data values that are numbers
   if verbose:
      print " sigma_clip: Full size of data       = %d" % data.size
      print " sigma_clip: Number of finite values = %d" % \
          data[n.isfinite(data)].size
   d = data[n.isfinite(data)].flatten()
   avg = d.mean()
   std = d.std()
   
   delta = 1
   while delta:
      size = d.size
      d = d[abs(d-avg)<nsig*std]
      avg = d.mean()
      std = d.std()
      delta = size-d.size
   return avg,std

#-----------------------------------------------------------------------

def set_param_array(hdulen,inval):
   """
   Converts an input parameter value, which may have been passed as an
   integer or float, into an array so that it can be used in general 
   image-processing functions.
   
   Inputs:
      hdulen - the number of HDUs in the associated fits file
      inval  - the value of the variable

   Output:
      valarr - the input value as an array, if necessary
   """

   if (isinstance(inval,int)) or (isinstance(inval,float)):
      if hdulen == 1:
         valarr = n.array([inval]).astype(int)
      else:
         valarr = n.zeros(hdulen-1).astype(int)
         for i in range(hdulen-1):
            valarr[i] = inval
   else:
      valarr = inval

   return valarr

#-----------------------------------------------------------------------

def define_trimsec(hdu,x1,x2,y1,y2):
   xmax = hdu.header['NAXIS1']
   ymax = hdu.header['NAXIS2']
   #xmax = pf.getval(fullfits,'NAXIS1')
   #ymax = pf.getval(fullfits,'NAXIS2')
   if x1 != 0:
      x1 = int(x1)
   if x2 != 0:
      x2 = int(x2)
   else:
      x2 = xmax
   #xtrim = slice(x1,x2,1)
   if y1 != 0:
      y1 = int(y1)
   if y2 != 0:
      y2 = int(y2)
   else:
      y2 = ymax
   #ytrim = slice(y1,y2,1)
   return x1,x2,y1,y2

#-----------------------------------------------------------------------

def divide_images(a,b,output):
   print "Dividing images: '%s' / '%s' = '%s'" % (a,b,output)
   try:
      hdua = pf.open(a)[0]
   except:
      hdua = pf.open(a,ignore_missing_end=True)[0]
   try:
      hdub = pf.open(b)[0]
   except:
      hdub = pf.open(b,ignore_missing_end=True)[0]
   hdu = pf.PrimaryHDU(hdua.data/hdub.data)
   hdu.writeto(output,output_verify='ignore')

def subtract_images(a,b,output,hexta=0,hextb=0):
   print("Subtracting images: '%s' - '%s' = '%s'" % (a,b,output))
   try:
      hdua = pf.open(a)[hexta]
   except:
      hdua = pf.open(a,ignore_missing_end=True)[hexta]
   try:
      hdub = pf.open(b)[hextb]
   except:
      hdub = pf.open(b,ignore_missing_end=True)[hextb]
   hdu = pf.PrimaryHDU(hdua.data - hdub.data)
   hdu.writeto(output,output_verify='ignore',clobber=True)

def add_images(a,b,output,preserve_header=0):
   print("Adding images: '%s' + '%s' = '%s'" % (a,b,output))
   try:
      hdua = pf.open(a)[0]
   except:
      hdua = pf.open(a,ignore_missing_end=True)[0]
   try:
      hdub = pf.open(b)[0]
   except:
      hdub = pf.open(b,ignore_missing_end=True)[0]

   if preserve_header == 1:
      hdu = pf.PrimaryHDU(hdua.data + hdub.data,header=hdua.header)
   elif preserve_header == 2:
      hdu = pf.PrimaryHDU(hdua.data + hdub.data,header=hdub.header)
   else:
      hdu = pf.PrimaryHDU(hdua.data + hdub.data)

   hdu.writeto(output,output_verify='ignore',clobber=True)

#-----------------------------------------------------------------------

def read_calfile(filename, file_description):
   print 'Reading in %s file: %s' % (file_description, filename)
   try:
      calhdulist = pf.open(filename)
   except:
      try:
         calhdulist = pf.open(filename,ignore_missing_end=True)
      except:
         print " ERROR: Requested %s file %s does not exist" % \
             (file_description, filename)
         print ""
         return -1
   return calhdulist

#-----------------------------------------------------------------------

def median_combine(input_files,output_file,x1=0,x2=0,y1=0,y2=0,
                   biasfile=None,gain=-1.0,normalize=False,zeromedian=False,
                   NaNmask=False,hdu0only=False):
   """ Given a list of input file names, this function will:
      1. Subtract a bias frame (if the optional biasfile parameter is set)
      2. Multiply by the gain, required to be in e-/ADU (if the optional
         gain parameter is set)
      3. Normalize the frame (if the optional normalize parameter is set)
      4. Subtract the median (if the optional zeromedian parameter is set)
      5. Median combine the resulting data
      6. Write the output to a file.
   """

   print "median_combine: Inputs:"
   print "-----------------------"
   if (biasfile != None):
      # NB: Should not need to trim the input bias frame, since it should
      #  have been created out of trimmed files.
      bias = pf.open(biasfile)
      bias.info()
   else:
      print "  bias frame: [No bias file]"

   print ""
   print "median_combine: Loading files"
   print "-----------------------------"
   files = []
   for filename in input_files:
      print " %s" % filename
      try:
         f = pf.open(filename)
      except:
         f = pf.open(filename,ignore_missing_end=True)
      files.append(f)

   # Use first file to check for multiple HDUs, and set some variables
   #  depending on the result
   print ""
   print "median_combine: Getting info on first file"
   print "------------------------------------------"
   files[0].info()
   hdulen = len(files[0])
   if (hdulen == 1) or (hdu0only == True):
      imhdu = n.arange(1)
   else:
      imhdu = n.arange(1,hdulen)
      phdu = pf.PrimaryHDU()
      phdu.header.add_comment(
         'This file contains %d image extensions' % (hdulen-1),after='extend')
      hdulist = pf.HDUList([phdu])

   # Set up the trim variables as numpy arrays, if they are not
   #  already in that format

   x1 = set_param_array(hdulen,x1)
   x2 = set_param_array(hdulen,x2)
   y1 = set_param_array(hdulen,y1)
   y2 = set_param_array(hdulen,y2)
   gain = set_param_array(hdulen,gain)

   # Start outer loop on the extension number.  For old-school FITS files
   #  there is only one extension, which is the main image.

   for j in imhdu:

      # Set trim section for this HDU in the input files, and use that
      #  to define the container for the stack of images

      k = j-1
      xt1,xt2,yt1,yt2 = define_trimsec((files[0])[j], \
                                          x1[k],x2[k],y1[k],y2[k])
      print ""
      print "median_combine: setting up stack for images (HDU %d)" % j
      print "----------------------------------------------------"
      print "Stack will have dimensions (%d,%d,%d)" \
          %(len(input_files),yt2-yt1,xt2-xt1)
      stack = n.zeros((len(input_files),yt2-yt1,xt2-xt1))

      # Inner loop is on the input files

      count = 0
      for i in range(len(input_files)):
         
         print " %s" % files[i].filename()

         # Process the data (bias and gain only), if desired
         if (biasfile == None):
            process_data(files[i],j,gain=gain[k],x1=x1[k],x2=x2[k], \
                            y1=y1[k],y2=y2[k])
         else:
            process_data(files[i],j,bias,gain=gain[k],x1=x1[k], \
                            x2=x2[k],y1=y1[k],y2=y2[k])

         tmpf = (files[i])[j].data

         # Normalize or set to zero median, if desired
         if(normalize == True):
            frame_med = n.median(tmpf,axis=None)
            print "    Normalizing %s by %f" % (files[i].filename(),frame_med)
            tmpf /= frame_med
         if(zeromedian == True):
            print "    Subtracting the median from %s" % files[i].filename()
            tmpf -= n.median(tmpf,axis=None)
         stack[count] = tmpf.copy()
         count += 1
      
      print ""

      # Actually form the median
      # Note that numpy median converts the data type to float64
      if(NaNmask == True):
         print "median_combine: Computing median frame using NaN masking"
         print "   Can take a while..."
         median = sp.stats.stats.nanmedian(stack,axis=0)
      else:
         print "median_combine: Computing median frame (can take a while)..."
         median = n.median(stack,axis=0)
      del stack

      # Save the median HDU
      # For multi-extension files, use pf.ImageHDU(median)
      if(hdulen == 1) or (hdu0only):
         phdu = pf.PrimaryHDU(median)
         hdulist = pf.HDUList([phdu])
      else:
         hdu = pf.ImageHDU(median)
         hdulist.append(hdu)

   # Write the output median file

   hdulist.writeto(output_file,output_verify='ignore',clobber=True)
   print "   ... Writing output to %s." % output_file

   # Clean up

   for i in range(len(input_files)):
      files[i].close()

#-----------------------------------------------------------------------

def apply_rough_wcs(hdu, pixscale, phdu=None):
   """
   Takes the RA and Dec pointing info from the fits header, along with
   a pixel scale, and converts that information into the standard WCS
   headers (CD1_1, etc.).

   For multi-extension FITS files may contain the RA and Dec information
   only in the PHDU and not in the image extension HDUs.  In this case,
   the function can be passed the PHDU via the optional phdu parameter.
   If phdu=None then the function will look for the RA and Dec header
   cards in the hdu that has been passed via the hdu parameter.
   """

   if phdu is None:
      phdu = hdu

   if pixscale>0.0:

      """ 
      Get the size of the data array.  Note that if this fails, the
      probable reason is that the selected hdu is the PHDU of a
      multiextension fits file.  In that case, don't set the CRPIXes
      """

      containsdata = True
      try:
         shape = hdu.data.shape
      except:
         containsdata = False

      """
      Read the header cards
      """
      wcsread = True
      try:
         ra  = wcs.ra2deg(hdu.header['RA'].strip())
      except:
         try:
            ra  = wcs.ra2deg(phdu.header['RA'].strip())
         except:
            print 'ERROR. Attempts to read RA header card failed.'
            print 'No wcs information'
            wcsread = False
      try:
         dec = wcs.dec2deg(hdu.header['DEC'].strip())
      except:
         try:
            dec = wcs.dec2deg(phdu.header['DEC'].strip())
         except:
            print 'ERROR. Attempts to read RA header card failed.'
            print 'No wcs information'
            wcsread = False

      """ Clean out old WCS and mosaic info """
      foo = hdu.header.ascardlist()[4:]
      for k in range(0,len(foo)):
         keyname = foo[k].key
         if keyname[0:4] == 'CD1_' or keyname[0:4] == 'CD2_' or \
                keyname[0:4] == 'PANE' or keyname[0:5] == 'CRVAL' \
                or keyname[0:5] == 'CRPIX' or keyname[0:5] == 'CDELT' \
                or keyname[0:5] == 'CTYPE' or keyname[0:5] == 'CUNIT' \
                or keyname[0:5] == 'CRDER' or keyname[0:5] == 'CSYER' \
                or keyname[0:7] == 'WCSNAME' or keyname[0:3] == 'ADC':
            del hdu.header[keyname]

      """ Apply the rough WCS info """
      if wcsread:
         hdu.header.update('CRVAL1',ra)
         hdu.header.update('CRVAL2',dec)
         if containsdata:
            hdu.header.update('CRPIX1',shape[1]/2.)
            hdu.header.update('CRPIX2',shape[0]/2.)
         hdu.header.update('CD1_1',-pixscale/3600.)
         hdu.header.update('CD1_2',0.)
         hdu.header.update('CD2_1',0.)
         hdu.header.update('CD2_2',pixscale/3600.)
         hdu.header.update('CTYPE1','RA---TAN')
         hdu.header.update('CTYPE2','DEC--TAN')
         hdu.header.update('EQUINOX',2000.0)
         hdu.header.update('RADECSYS','FK5')

#-----------------------------------------------------------------------

def make_bias(bias_frames, raw_prefix, rawdir="../Raw", rawext='fits',
      outfile="Bias.fits", x1=0,x2=0,y1=0,y2=0,hdu0only=False):
   """ 

       This function takes as input the frame numbers of the
       dark frames (either true darks, or bias frames) and
       median-combines them to create a master dark

       Required inputs:
        dark_frames - an array of frame numbers (e.g., [101,102,103,105])
        raw_prefix  - prefix for frame numbers (e.g., "lred0")

       Optional inputs:
        rawdir      - directory with raw files (default="../Raw")
        rawext      - extension for raw filenames (default = 'fits')
        outfile     - output filename (default="masterdark.fits")
        x1          - to set a trim section that is smaller than the full frame
        x2          - to set a trim section that is smaller than the full frame
        y1          - to set a trim section that is smaller than the full frame
        y2          - to set a trim section that is smaller than the full frame

   """

   # Make file list
   filenames = []
   for i in bias_frames:
      filenames.append('%s/%s%d.%s'%(rawdir,raw_prefix,i,rawext))

   # Call median_combine
   median_combine(filenames,outfile,x1,x2,y1,y2,hdu0only=hdu0only)

#-----------------------------------------------------------------------

def make_flat(flat_frames, raw_prefix, rawdir="../Raw", rawext='fits',
      outfile="Flat.fits", biasfile=None, gain=-1.0, normalize=True,
      x1=0, x2=0, y1=0, y2=0, framesig=0):
   """ 

       This function takes as input the frame numbers of the flat-field
       frames and median-combines them to create a master flat

       Required inputs:
        flat_frames - array of input frame numbers (e.g., [101,102,103,105])
        raw_prefix  - prefix for frame numbers (e.g., "lred0")

       Optional inputs:
        rawdir      - directory with raw files (default="../Raw")
        rawext      - extension for raw filenames (default = 'fits')
        outfile     - output filename (default="Flat.fits")
        biasfile    - input bias file to subtract before combining. 
                      (default=None)
        gain        - gain to divide by before combining (default=1.0)
        normalize   - normalize flat-field frames before combining?
                      (default = True)
        x1          - to set a trim section that is smaller than the full frame
        x2          - to set a trim section that is smaller than the full frame
        y1          - to set a trim section that is smaller than the full frame
        y2          - to set a trim section that is smaller than the full frame
        framesig    - number of significant digits in the frame numbers.  Keep
                      at the default (framesig=0) if all of the frame numbers
                      have the same number of significant digits (e.g., if all of
                      them are between 100 and 999).  However, if there is a
                      range of significant figures (e.g., frame numbers between
                      75 and 125), then set this to the higher number of figures.

   """

   # Make file list
   filenames = []
   for i in flat_frames:
      if framesig == 2:
         filenames.append('%s/%s%02d.%s'%(rawdir,raw_prefix,i,rawext))
      elif framesig == 3:
         filenames.append('%s/%s%03d.%s'%(rawdir,raw_prefix,i,rawext))
      elif framesig == 4:
         filenames.append('%s/%s%04d.%s'%(rawdir,raw_prefix,i,rawext))
      else:
         filenames.append('%s/%s%d.%s'%(rawdir,raw_prefix,i,rawext))

   # Call median_combine
   median_combine(filenames,outfile,biasfile=biasfile,gain=gain,
                  normalize=normalize,x1=x1,x2=x2,y1=y1,y2=y2)

#-----------------------------------------------------------------------

def make_flat_files(file_with_filelist, indir=None, 
                    outfile="Flat.fits", biasfile=None, gain=1.0, normalize=True,
                    x1=0, x2=0, y1=0, y2=0):
   """ 

       This function takes as input the frame numbers of the flat-field
       frames and median-combines them to create a master flat

       Required inputs:
        file_with_filelist - file containing list of input files 
          (e.g., "filelist.txt")

       Optional inputs:
        indir       - directory containing input files (default=None)
        outfile     - output filename (default="Flat.fits")
        biasfile    - input bias file to subtract before combining. 
                      (default=None)
        gain        - gain to divide by before combining (default=1.0)
        normalize   - normalize flat-field frames before combining?
                      (default = True)
        x1          - to set a trim section that is smaller than the full frame
        x2          - to set a trim section that is smaller than the full frame
        y1          - to set a trim section that is smaller than the full frame
        y2          - to set a trim section that is smaller than the full frame

   """

   # Extract base file names from input file
   basenames = open(file_with_filelist).read().split()

   # Make file list
   if indir is None:
      filenames = basenames
   else:
      filenames = []
      for i in basenames:
         filenames.append('%s/%s'%(indir,i))

   # Call median_combine
   median_combine(filenames,outfile,biasfile=biasfile,gain=gain,
    normalize=normalize,x1=x1,x2=x2,y1=y1,y2=y2)

#-----------------------------------------------------------------------

def make_fringe(fringe_frames, in_prefix, indir=None, inext='fits',
      outfile="Fringe.fits", biasfile=None, gain=1.0, normalize=False,
      zeromedian=True, x1=0, x2=0, y1=0, y2=0):
   """ 

       This function takes an array containing the flat-fielded frame 
       numbers and subtracts from each flat-fielded science frame
       its individual median.  The resulting data are median-combined
        to create a master fringe file

       Required inputs:
        fringe_frames - array of input frame numbers (e.g., [101,102,103,105])
        in_prefix     - prefix for frame numbers (e.g., "ff")

       Optional inputs:
        indir       - directory containing input files (default=None)
        inext       - extension for input filenames (default = 'fits')
        outfile     - output filename (default="Fringe.fits")
        biasfile    - input bias file to subtract before combining. 
                      (default=None)
        gain        - gain to divide by before combining (default=1.0)
        zeromedian  - subtract median from each science frame before combining?
                      (default = True)
        normalize   - normalize science frames before combining?
                      (default = False)
        x1          - to set a trim section that is smaller than the full frame
        x2          - to set a trim section that is smaller than the full frame
        y1          - to set a trim section that is smaller than the full frame
        y2          - to set a trim section that is smaller than the full frame

   """

   # Make file list
   if indir is None:
      indir = '.'
   filenames = []
   for i in fringe_frames:
      filenames.append('%s/%s%d.%s'%(indir,in_prefix,i,inext))

   # Call median_combine
   median_combine(filenames,outfile,biasfile=biasfile,gain=gain,
    normalize=normalize,zeromedian=zeromedian,x1=x1,x2=x2,y1=y1,y2=y2)

#-----------------------------------------------------------------------

def make_fringe_files(fringe_frames, in_prefix, indir=None, 
      outfile="Fringe.fits", biasfile=None, gain=1.0, normalize=False,
      zeromedian=True, x1=0, x2=0, y1=0, y2=0):
   """ 

       This function takes a list of the flat-fielded input files
       and subtracts their individual medians before median-combining
       them to create a master fringe file

       Required inputs:
        file_with_filelist - file containing list of input files 
          (e.g., "filelist.txt").  The input files should be
          flat-fielded science files.

       Optional inputs:
        indir       - directory containing input files (default=None)
        outfile     - output filename (default="Flat.fits")
        biasfile    - input bias file to subtract before combining. 
                      (default=None)
        gain        - gain to divide by before combining (default=1.0)
        zeromedian  - subtract median from each science frame before combining?
                      (default = True)
        normalize   - normalize science frames before combining?
                      (default = False)
        x1          - to set a trim section that is smaller than the full frame
        x2          - to set a trim section that is smaller than the full frame
        y1          - to set a trim section that is smaller than the full frame
        y2          - to set a trim section that is smaller than the full frame

   """

   # Extract base file names from input file
   basenames = open(file_with_filelist).read().split()

   # Make file list
   if indir is None:
      filenames = basenames
   else:
      filenames = []
      for i in basenames:
         filenames.append('%s/%s'%(indir,i))

   # Call median_combine
   median_combine(filenames,outfile,biasfile=biasfile,gain=gain,
    normalize=normalize,zeromedian=zeromedian,x1=x1,x2=x2,y1=y1,y2=y2)

#-----------------------------------------------------------------------

def process_data(hdulist, hdunum, bias=None, flat=None, fringe=None, 
                 darksky=None, gain=-1.0, texp_key=None, skysub=False,
                 flip=0, pixscale=0.0, x1=0, x2=0, y1=0, y2=0):

   """ This function applies calibration corrections to the passed HDU.  All
        of the calbration steps are by default turned off (keywords set to None).
        To apply a particular calibration step, set the appropriate keyword.
        The possible steps, along with their keywords are:

          Keyword     Calibration step
          ----------  ----------------------------------
          bias        Bias subtraction
          gain        Convert from ADU to electrons if set to value > 0
                      NB: Gain must be in e-/ADU
          flat        Flat-field correction
          fringe      Fringe subtraction
          darksky     Dark-sky flat correction
          skysub      Subtract mean sky level if keyword set to True
          texp_key    Divide by exposure time (set keyword to fits header name)
          flip        0 => no flip
                      1 => PFCam style (flip x then rotate -90), 
                      2 => P60 CCD13 style (not yet implemented)
                      3 => flip x-axis
          pixscale    If >0, apply a rough WCS using this pixel scale (RA and
                       Dec come from telescope pointing info in fits header)
      
       Required inputs:
        in_frames   - list of input file names
        in_prefix   - input prefix that will be replaced by the output prefix
        out_prefix  - prefix for the flattened output frames

       Optional inputs (in addition to those listed above in the keyword list):
        rawdir      - directory with raw files (default="../Raw")
        x1          - to set a trim section that is smaller than the full frame
        x2          - to set a trim section that is smaller than the full frame
        y1          - to set a trim section that is smaller than the full frame
        y2          - to set a trim section that is smaller than the full frame
   """

   # Make the name of the passed HDU more tractable
   tmp = hdulist[hdunum]
   phdu = hdulist[0]

   # Trim the data if requested
   xt1,xt2,yt1,yt2 = define_trimsec(tmp,x1,x2,y1,y2)
   tmp.data = tmp.data[yt1:yt2,xt1:xt2].astype(n.float64)
   if xt2-xt1 != tmp.header['naxis1'] or yt2-yt1 != tmp.header['naxis2']:
      tmp.header.update('trim',
         'Trim data section is [%d:%d,%d:%d] ([xrange,yrange])' %
         (xt1,xt2,yt1,yt2))
      print "   Trimmed data to section [xrange,yrange] [%d:%d,%d:%d]" \
          % (xt1,xt2,yt1,yt2)
   
   # Set up a string for use in header keywords
   if hdunum==0:
      hdustr = 'these data'
   else:
      hdustr = 'HDU %d' % hdunum
   
   # Bias-subtract if requested
   if bias is not None:
      biasdata = bias[hdunum].data
      tmp.data -=  biasdata
      biasmean = biasdata.mean()
      if (hdunum == 0):
         keystr = 'biassub'
      else:
         keystr = 'biassub'+str(hdunum)
      tmp.header.update(keystr,
         'Bias frame for %s is %s with mean %f' % 
                        (hdustr,bias.filename(),biasmean))
      print "   Subtracted bias frame %s" % bias.filename()
   
   # Convert to electrons if requested
   if gain>0:
      tmp.data *= gain
      tmp.header.update('gain', 1.0, 'Units are now electrons')
      if (hdunum == 0):
         keystr = 'gainorig'
      else:
         keystr = 'gainori'+str(hdunum)
      tmp.header.update(keystr, gain, 'Original gain for %s in e-/ADU' % hdustr,
         after='gain')
      tmp.header.update('bunit', 'Electrons', 
         'Converted from ADU in raw image')
      if (hdunum == 0):
         keystrb1 = 'binfo_1'
      else:
         keystrb1 = 'binfo'+str(hdunum)+'_1'
      keystrb1 = keystrb1.upper()
      tmp.header.update(keystrb1,
         'Units for %s changed from ADU to e- using gain=%6.3f e-/ADU' %
                        (hdustr,gain))
      print "   Converted units to e- using gain = %f" % gain
   
   # Divide by the exposure time if requested
   if texp_key is not None:
      texp_good = True
      try:
         texp = tmp.header[texp_key]
      except:
         try:
            texp = phdu.header[texp_key]
         except:
            print("")
            print("ERROR: No exposure time keyword called %s found in header"
                  % texp_key)
            print(" Setting texp = 1.0")
            texp_good = False
            texp = 1.0
      tmp.data /= texp
      if texp_good:
         if (hdunum == 0):
            keystr = 'binfo_2'
         else:
            keystr = 'binfo'+str(hdunum)+'_2'
         keystr = keystr.upper()
         tmp.header.update('gain', texp, 
                           'If units are e-/s then gain=t_exp')
         tmp.header.update('bunit','Electrons/sec','See %s header' % keystr,
                           keystr)
         tmp.header.update(keystr,
                           'Units for %s changed from e- to e-/s using texp=%7.2f'
                           % (hdustr,texp), after=keystrb1)
         print "   Converted units from e- to e-/sec using exposure time %7.2f" \
             % texp
      else:
         tmp.header.update('bunit','Electrons',
                           'There was an ERROR in converting to e-/s')
         tmp.header.update(keystr,'ERROR: Exposure time query failed.')
   
   # Apply the flat-field correction if requested
   if flat is not None:
      flatdata = flat[hdunum].data
      tmp.data /= flatdata
      flatmean = flatdata.mean()
      # Set up a bad pixel mask based on places where the flat frame = 0,
      #  since dividing by zero gives lots of problems
      zeromask = flatdata==0
      # Correct for any zero pixels in the flat-field frame
      tmp.data[zeromask] = 0
      if (hdunum == 0):
         keystr = 'flatcor'
      else:
         keystr = 'flatcor'+str(hdunum)
      tmp.header.update(keystr,
         'Flat field image for %s is %s with mean=%f' %
                        (hdustr,flat.filename(),flatmean))
      print "   Divided by flat-field image: %s" % flat.filename()
   
   # Apply the fringe correction if requested
   if fringe is not None:
      tmp.data -= fringe
      if (hdunum == 0):
         keystr = 'fringcor'
      else:
         keystr = 'frngcor'+str(hdunum)
      tmp.header.update(keystr,
         'Fringe image for %s is %s with mean=%f' % 
                        (hdustr,fringe.filename(),fringemean))
      print "   Subtracted fringe image: %s" % fringe.filename()
   
   # Apply the dark sky flat-field correction if requested
   if darksky is not None:
      tmp.data /= darksky
      # Correct for any zero pixels in the flat-field frame
      tmp.data[dszeromask] = 0
      if (hdunum == 0):
         keystr = 'darksky'
      else:
         keystr = 'darksky'+str(hdunum)
      tmp.header.update(keystr,
         'Dark-sky flat image for %s is %s with mean=%f' % 
         (hdustr,darksky.filename(),darkskymean))
      print "   Divided by dark-sky flat: %s" % darksky.filename()

   # Subtract the sky level if requested
   if skysub:
      m,s = sigma_clip(tmp.data)
      tmp.data -= m
      if (hdunum == 0):
         keystr = 'skysub'
      else:
         keystr = 'skysub'+str(hdunum)
      tmp.header.update(keystr,'For %s, subtracted mean sky level of %f' %
                        (hdustr,m))
      print '   Subtracted mean sky level of %f' % m
   
   # Flip if requested
   if flip == 1:
      tmp.data = tmp.data.T[::-1,::-1]
   if flip == 3:
      tmp.data = tmp.data[:,::-1]
   
   # Add a very rough WCS if requested
   if pixscale > 0.0:
      if hdunum>0:
         apply_rough_wcs(tmp,pixscale,hdulist[0])
      else:
         apply_rough_wcs(tmp,pixscale)

#   return hdu

#-----------------------------------------------------------------------

def apply_calib(in_frames, in_prefix, out_prefix, split=False,
        biasfile=None, flatfile=None, fringefile=None, darkskyfile=None,
        skysub=False, gain=-1.0, texp_key=None, flip=0, pixscale=0.0,
        rawdir="../Raw", rawext='fits', x1=0, x2=0, y1=0, y2=0):
   """ This function applies calibration corrections to the input files,
        which are designated by an array of frame numbers.  All of the
        calbration steps are by default turned off (keywords set to None).
        To apply a particular calibration step, set the appropriate keyword.
        The possible steps, along with their keywords are:

          Keyword     Calibration step
          ----------  ----------------------------------
          biasfile    Bias subtraction
          gain        Convert from ADU to electrons if set to value > 0
                      NB: gain must be in e-/ADU
          flatfile    Flat-field correction
          fringefile  Fringe subtraction
          darkskyfile Dark-sky flat correction
          skysub      Subtract mean sky level if keyword set to True
          texp_key    Divide by exposure time (set keyword to fits header name)
          flip        0 => no flip
                      1 => PFCam-style (flip x then rotate -90), 
                      2 => P60 CCD13 style (not yet implemented)
                      3 => flip x-axis
          pixscale    If >0, apply a rough WCS using this pixel scale (RA and
                       Dec come from telescope pointing info in fits header)
      
       Required inputs:
        in_frames   - list of input file names
        in_prefix   - input prefix that will be replaced by the output prefix
        out_prefix  - prefix for the flattened output frames

       Optional inputs (in addition to those listed above in the keyword list):
        rawdir      - directory with raw files (default="../Raw")
        rawext      - extension for raw filenames (default = 'fits')
        x1          - to set a trim section that is smaller than the full frame
        x2          - to set a trim section that is smaller than the full frame
        y1          - to set a trim section that is smaller than the full frame
        y2          - to set a trim section that is smaller than the full frame
   """

   # Read in calibration frames if they have been selected

   print ''
   if biasfile is not None:
      bias = read_calfile(biasfile,'bias')
      if (bias == -1):
         return
   else:
      bias = None

   # Put file error checking for flats, fringe, and darksky files

   if flatfile is not None:
      print 'Reading in flat-field file: %s' % flatfile
      flat = read_calfile(flatfile,'flat-field')
      if (flat == -1):
         return
   else:
      flat = None

   if fringefile is not None:
      print 'Reading in fringe file: %s' % fringefile
      fringe = pf.getdata(fringefile)
      fringemean = fringe.mean()
   else:
      fringe = None

   if darkskyfile is not None:
      print 'Reading in dark-sky flat file: %s' % darkskyfile
      darksky = pf.getdata(darkskyfile)
      darkskymean = darksky.mean()
      # Set up a bad pixel mask based on places where the flat frame = 0,
      #  since dividing by zero gives lots of problems
      dszeromask = darksky==0
   else:
      darksky = None

   print ""
   print "Processing files..."
   print "-------------------"

   write_one_output_file = True
   for i in in_frames:
      filename = '%s/%s%s.%s'%(rawdir,in_prefix,i,rawext)
      print "%s:" % filename

      # Open the input file.  Read in this way in order to be able to handle
      #  multi-extension fits files

      try:
         hdulist = pf.open(filename)
      except:
         hdulist = pf.open(filename,ignore_missing_end=True)

      # Set things up depending on whether there are extensions or not
      hdulen = len(hdulist)
      if hdulen == 1:
         imhdu = n.arange(1)
      else:
         imhdu = n.arange(1,hdulen)
         phdu = hdulist[0]
         phdu.header.add_comment(
            'This file contains %d image extensions' % (hdulen-1),after='extend')
         outhdulist = pf.HDUList([phdu])

      # Set up the trim and flip variables as numpy arrays, if they are not
      #  already in that format

      x1 = set_param_array(hdulen,x1)
      x2 = set_param_array(hdulen,x2)
      y1 = set_param_array(hdulen,y1)
      y2 = set_param_array(hdulen,y2)
      gain = set_param_array(hdulen,gain)
      flip = set_param_array(hdulen,flip)

      # Loop over fits extensions.  This loop assumes that if the input file
      #  is a multi-extension fits file, then HDU 0 does not contain any
      #  data, and only HDU's 1-N do.  Of course, if this is a straightforward
      #  fits file with no extensions, then HDU 0 will contain the data


      for j in imhdu:

         # Check that extension is valid (to be implemented)

         # Read in the data
         if hdulen>0:
            print " Processing image extension %d" % j

         # Process the data
         k = j-1
         process_data(hdulist,j,bias,flat,fringe,darksky,gain[k],texp_key,
                      skysub,flip[k],pixscale,x1[k],x2[k],y1[k],y2[k])

         # If multiple extensions, append the current extension to the HDU list
         if hdulen>1:
            if(split):
               write_one_output_file = False
               print " Splitting image extension %d to output file."
            else:
               outhdulist.append(hdulist[j])
         else:
            outhdulist = hdulist[j]

      # Write out final file
      # if hdulen>1:
      #    newhdu = pf.PrimaryHDU(tmp.data)
      #    foo = phdu.header.ascardlist()[4:]
      #    for k in range(0,len(foo)):
      #       if foo[k].key[0:5] != 'NAXIS' and foo[k].key[0:4] != 'PANE' and \
      #              foo[k].key != 'COMMENT':
      #          newhdu.header.update(foo[k].key,foo[k].value,foo[k].comment)
      # else:
      #    newhdu = tmp

      if(write_one_output_file):
         outname = '%s%s.fits'%(out_prefix,i)
         outname = outname.strip()
         outhdulist.writeto(outname,output_verify='ignore',clobber=True)
         print " Writing output file: %s" % outname

      hdulist.close()

# For LRIS B
#  x1 = [400, 51, 51, 400]
#  x2 = 1068
#  y1 = 775
#  y2 = 3200

#-----------------------------------------------------------------------

def hdr_offsets(files, pixscale=0, rakey=None, deckey=None, rot=None,
                oformat='pix', hext=0, verbose=True):
   """
   Uses the telescope-provided RA and Dec in the fits header to calculate
   an initial guess for the offsets between the input files.

   The default behavior for each file is to:
     1. Try to get the RA and Dec values from the CRPIXn header keywords
        and the pixel scale from the CDn_m keywords
     2. If step 1 fails for the RA and Dec, then try to get the RA and
        Dec information from the telescope pointing keywords defined by
        the rakey and deckey parameters.

   Inputs:
      files:     list of files for which to calculate offsets
      pixscale:  pixel scale of the input files, in arcsec/pix.  The default
                 value of 0 means that the pixel scale should be derived
                 from the CD matrix (e.g., CD1_1, etc.) in the fits header
      rakey:     name of the fits header keyword for RA. Default=None implies
                 that the default behavior described above is followed.
                 Set to a non-None value to override default
      deckey:    name of the fits header keyword for Dec.  Default=None implies
                 that the default behavior described above is followed.
                 Set to a non-None value to override default  
      rot:       Rotation of coordinate axes, North through East, in degrees.
                 Set this parameter to override default setting, which is the
                 rotation set by the CD matrix (if it exists in the fits header)
                 or 0.0 (if the CD matrix doesn't exist)
      oformat:   Output format for offsets.  Default is 'pix' for pixel offsets.
                 The only other option is 'arcsec' for offsets in arcseconds.
      hext:      HDU extension.  Default is hext=0, i.e., the primary HDU,
                 which is also the only HDU for many files
      verbose:   set to True for verbose output.  Default=True
   """

   # Initialize containers
   ra   = n.zeros(len(files))
   dec  = n.zeros(len(files))
   dx   = n.zeros(len(files))
   dy   = n.zeros(len(files))
   cdmatx = []

   # Print out some header information
   if verbose:
      print ""
      print \
          "   File                RA       Dec     xpxscl ypxscl rotation"
      print \
          "------------------- --------- --------- ------ ------ --------"
   

   # Read in the WCS information and process it
   count = 0
   for f in files:

      # Open the file and get the appropriate header info
      hdulist = imf.open_fits(f)
      if hdulist is None:
         return None,None
      hdr = hdulist[hext].header

      # Initialize the containers
      ra[count] = n.nan
      dec[count] = n.nan

      # Start with the default behavior of trying to read in the WCS info
      # If this is successful, set the parameter values based on the WCS info
      foundwcs = True
      badwcs = False
      try:
         wcsinfo = wcs.parse_header(hdr)
      except:
         foundwcs = False
      if foundwcs:
         cd = wcsinfo[2]
         if wcsinfo[4] == 'RA':
            ra[count] = wcsinfo[0][1]
         else:
            ra[count] = wcsinfo[0][0]
         if wcsinfo[4] == 'DEC':
            dec[count] = wcsinfo[0][1]
         else:
            dec[count] = wcsinfo[0][0]
      else:
         cd = n.array([[-1.0,0.0],[0.0,1.0]])

      # Calculate the RA and Dec based on whether the override parameters are set
      if rakey is not None:
         try:
            ra[count] = hdr[rakey]
         except:
            print ""
            print "ERROR: Could not read RA using keyword %s from file %s" % \
                (rakey,f)
            return None,None

      if deckey is not None:
         try:
            dec[count] = hdr[deckey]
         except:
            print ""
            print "ERROR: Could not read DEC using keyword %s from file %s" % \
                (deckey,f)
            return None,None

      # Set the CD matrix if pixscale is set to override. Default rotation = 0
      if pixscale > 0.:
         if rot is None:
            rot = 0.0
         cd = coords.rscale_to_cdmatrix(pixscale,rot,verbose=False)

      # Print out information about the input files
      if verbose:
         cd1,cd2,cr2 = coords.cdmatrix_to_rscale(cd)
         print "%-18s  %9.5f %+9.5f %+6.3f %+6.3f %+6.1f" \
             % (files[count],ra[count],dec[count],3600.*cd1,3600.*cd2,
                cr2*180./pi)

      # Finally, save the constructed CD matrix
      cdmatx.append(cd)

      count += 1
      del hdr
      del hdulist

   # Calculate the offsets
   dalpha = 3600. * ((ra - ra[0]) * cos(pi * dec[0]/180.))
   ddelta = 3600. * (dec - dec[0])
   if verbose:
      print ""
      print \
          "   File              da(asec) dd(asec)  dx(pix) dy(pix)"
      print \
          "-------------------  -------- --------  ------- -------"
   for i in range(dx.size):
      dx[i],dy[i] = coords.darcsec_to_dpix(dalpha[i],ddelta[i],cdmatx[i])
      if verbose:
         print "%-18s    %+7.3f  %+7.3f  %+7.2f %+7.2f" \
             % (files[i],dalpha[i],ddelta[i],dx[i],dy[i])

   # Return the offsets
   if oformat == 'arcsec':
      return dalpha,ddelta
   else:
      return dx,dy

#-----------------------------------------------------------------------

def plot_hdr_offsets(files, pixscale=0, rakey=None, deckey=None, rot=None, 
                     oformat='pix', hext=0, verbose=True):
   """
   Makes an initial calculation for the offsets between the input files based on
   information in the fits header cards (via hdr_offsets).  
   These offsets get plotted in the desired units, either pixels (the default)
   or arcseconds.

   Inputs:
      files:     list of files for which to calculate offsets
      pixscale:  pixel scale of the input files, in arcsec/pix.  The default
                 value of 0 means that the pixel scale should be derived
                 from the CD matrix (e.g., CD1_1, etc.) in the fits header
      rakey:     name of the fits header keyword for RA. Default=None implies
                 that the default behavior described in hdr_offsets is followed.
                 Set to a non-None value to override default
      deckey:    name of the fits header keyword for Dec.  Default=None implies
                 that the default behavior described in hdr_offsets is followed.
                 Set to a non-None value to override default           
      oformat:   Output format for offsets.  Default is 'pix' for pixel offsets.
                 The only other option is 'arcsec' for offsets in arcseconds.
      hext:      HDU extension.  Default is hext=0, i.e., the primary HDU,
                 which is also the only HDU for many files
      verbose:   set to True for verbose output.  Default=True
   """

   # Set up
   from matplotlib import pyplot as plt

   # Get the offsets
   dx,dy = hdr_offsets(files,pixscale,rakey,deckey,rot,oformat,hext,verbose)

   if dx is None or dy is None:
      print "ERROR. Failed to calculate offsets"
      return

   # Plot the offsets
   plt.scatter(dx,dy,marker='+')
   if oformat=='arcsec':
      plt.xlabel(r'd$\alpha$ (arcsec)')
      plt.ylabel(r'd$\delta$ (arcsec)')
   else:
      plt.xlabel('dx (pixels)')
      plt.ylabel('dy (pixels)')

#-----------------------------------------------------------------------------

def xcorr_offsets(infiles, hext=0):
   """
   Estimates the pixel shifts between dithered exposures on a field by doing
   a cross-correlation between the images.  The cross-correlation is done
   via fft.

   Inputs:
      infiles  - input files
      hext     - HDU of the data (default = 0)

   Output:
      xshift   - returned array of x offsets
      yshift   - returned array of y offsets
   """

   from scipy import fftpack

   """ Set the first file as the reference """
   ref = imf.Image(infiles[0])
   refdat = ref.hdu[hext].data.astype(n.float32)
   if refdat.shape[1]%2 == 0:
      xcent = refdat.shape[1] / 2
   else:
      xcent = (refdat.shape[1] - 1)/2
   if refdat.shape[0]%2 == 0:
      ycent = refdat.shape[0] / 2
   else:
      ycent = (refdat.shape[0] - 1)/2
   Fref = fftpack.fft2(refdat)

   """ Loop through the list of files """
   xshift = n.zeros(len(infiles))
   yshift = n.zeros(len(infiles))
   for i in range(1,len(infiles)):
      comp = imf.Image(infiles[i])
      compdat = comp.hdu[hext].data.astype(n.float32)
      Fprod = Fref * fftpack.fft2(compdat).conj()
      xcorr = fftpack.ifft2(Fprod)
      corrshift = fftpack.fftshift(xcorr)
      ymax,xmax = n.unravel_index(corrshift.argmax(),corrshift.shape)
      xshift[i] = xmax - xcent
      yshift[i] = ymax - ycent
      del compdat,Fprod,xcorr,corrshift
      #comp.close()
      del comp

   """ Clean up and return shifts """
   del refdat,Fref
   #ref.close()
   del ref
   return xshift,yshift
   
#------------------------------------------------------------------------------

def coadd_intshift(infiles, xshifts, yshifts, outfile, origsize=False, hext=0):

   """ For now assume that all of the inputs files have the same size """
   ref = imf.Image(infiles[0])
   ysize,xsize = ref.hdu[hext].data.shape
   tmpxshift = xshifts - xshifts.min()
   tmpyshift = yshifts - yshifts.min()
   if origsize:
      xmed = n.median(tmpxshift)
      ymed = n.median(tmpyshift)
   del ref

   """ Modify the size based on the shifts """
   dimx = xsize + int(n.ceil(xshifts.max() - xshifts.min()))
   dimy = ysize + int(n.ceil(yshifts.max() - yshifts.min()))

   """ Create containers for the output data and nexp info """
   odat = n.zeros((len(infiles),dimy,dimx))
   nexp = n.zeros((len(infiles),dimy,dimx))

   """ Load data into the containers """
   for i in range(len(infiles)):
      tmpim = imf.Image(infiles[i])
      tmpdat = tmpim.hdu[hext].data.copy()
      xstart = int(tmpxshift[i])
      xend = int(xstart + xsize)
      ystart = int(tmpyshift[i])
      yend = int(ystart + ysize)
      odat[i,ystart:yend,xstart:xend] = tmpdat.copy()
      nexp[i,ystart:yend,xstart:xend] = 1.
      del tmpdat,tmpim

   """ Take the weighted average """
   nexp[nexp==0] = 0.1
   outdat = odat.sum(axis=0) / nexp.sum(axis=0)
   if origsize:
      outdat = outdat[int(ymed):int(ymed+ysize),int(xmed):int(xmed+xsize)]
   pf.PrimaryHDU(outdat).writeto(outfile,clobber=True)
   del odat,nexp

#------------------------------------------------------------------------------

def fixpix_wht(datafile, whtfile, outfile=None, boxsize=11, datahdu=0,
               whthdu=0):
   """
   Replaces bad pixels in a fits image (datafile) based on the associated
   weight file.  Pixels for which the weight file is zero will be replaced
   by the median-filter value in the data file (i.e., a version of the
   data file is created where each pixel is replaced by the median pixel
   value in a boxsize x boxsize box centered on the pixel).
   If the optional outfile parameter is given, then the results are written
   out to a new file.  If not, then the data file is modified in place.
   """

   from scipy.ndimage import filters

   """ Open the input files """
   if outfile is None:
      dhdulist = imf.open_fits(datafile,'update')
      data = dhdulist[datahdu].data
   else:
      dhdulist = imf.open_fits(datafile)
      data = dhdulist[datahdu].data.copy()

   whdulist = imf.open_fits(whtfile)
   whtdat = whdulist[whthdu].data.copy()

   """ Median filter the data """
   meddat = filters.median_filter(data,boxsize)

   """ Replace the bad data (where whtdat==0) with the median-filtered values """
   mask = whtdat==0
   data[mask] = meddat[mask]
   dhdulist[datahdu].header.update('fixpix','ccdredux.fixpix_wht',
                                   'Ran code to fix bad pixels')

   """ Write output """
   if outfile is None:
      dhdulist.flush()
   else:
      pf.PrimaryHDU(data,dhdulist[datahdu].header).writeto(outfile)

#------------------------------------------------------------------------------

def fixpix_rms(datafile, rms_high=5., rms_low=5., rms_sigclip=3., outfile=None, 
               boxsize=11, datahdu=0, verbose=True):
   """
   Replaces bad pixels in a fits image (datafile), where the bad pixels are
   determined by:
    (1) running sigma_clip on the image data to get a clipped mean and sigma
    (2) finding pixels where the counts are (> mean + rms_high*sigma) or
        (< mean - rms_low*sigma)
   Bad pixels will be replaced by the median-filter value in the data file 
   (i.e., a version of the data file is created where each pixel is replaced by 
   the median pixel value in a boxsize x boxsize box centered on the pixel).
   If the optional outfile parameter is given, then the results are written
   out to a new file.  If not, then the data file is modified in place.
   """

   from scipy.ndimage import filters

   """ Open the input file """
   if outfile is None:
      dhdulist = imf.open_fits(datafile,'update')
      data = dhdulist[datahdu].data
   else:
      dhdulist = imf.open_fits(datafile)
      data = dhdulist[datahdu].data.copy()

   """ Sigma clip the input data and set the bad pixel mask """
   m,s = sigma_clip(data,verbose=verbose)
   mask = (data > m+rms_high*s) | (data < m-rms_low*s)

   """ Pass 1, replace bad pixels with the overall median value """
   data[mask] = n.median(data)

   """ Median filter the data """
   meddat = filters.median_filter(data,boxsize)

   """ Pass 2, replace the bad data with the median-filtered values """
   data[mask] = meddat[mask]

   """ Write output """
   if outfile is None:
      dhdulist.flush()
   else:
      pf.PrimaryHDU(data,dhdulist[datahdu].header).writeto(outfile)

