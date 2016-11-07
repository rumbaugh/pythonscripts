# make_color.py
# Takes 3 input images and makes a rgb output.

import sys,pyfits,scipy,pylab
from scipy import ndimage

def make_color(file1=None,file2=None,file3=None,bands=['F390W','F555W','F814W'],noweight=False,std_wfile_format=True,wfile1=None,wfile2=None,wfile3=None,outfile=None):

#file1 = sys.argv[1]
#file2 = sys.argv[2]
#file3 = sys.argv[3]
    data1,data2,data3 = None,None,None
    if file1 != None: data1 = pyfits.open(file1)[0].data.astype(scipy.float32)
    if file2 != None: data2 = pyfits.open(file2)[0].data.astype(scipy.float32)
    if file3 != None: data3 = pyfits.open(file3)[0].data.astype(scipy.float32)

# Load weight files
    if std_wfile_format:
        if file1 != None: wfile1 = file1.replace(".fits",".weight.fits")
        if file2 != None: wfile2 = file2.replace(".fits",".weight.fits")
        if file3 != None: wfile3 = file3.replace(".fits",".weight.fits")
    if noweight:
        if file1 != None: w1 = np.ones(np.shape(data1))
        if file2 != None: w2 = np.ones(np.shape(data2))
        if file3 != None: w3 = np.ones(np.shape(data3))
    else:
        if file1 != None: w1 = pyfits.open(wfile1)[0].data.astype(scipy.float32)
        if file2 != None: w2 = pyfits.open(wfile2)[0].data.astype(scipy.float32)
        if file3 != None: w3 = pyfits.open(wfile3)[0].data.astype(scipy.float32)
    if file1 == None: w1 = np.ones(np.shape(data2))
    if file2 == None: w2 = np.ones(np.shape(data3))
    if file3 == None: w3 = np.ones(np.shape(data1))

# Mask of good points in ALL THREE images
    w = w1*w2*w3

# The array bad masks non-image reasons. This is different than the weight
#   array w because w masks parts within the image interior that we would like
#   to keep in the output. The maximum_filter 'fills in' these bad pixels.
#   It looks like the wht map from multidrizzle or swarp is screwing up; I
#   don't know why these pixels are zero (a small number sure, but zero?!)....
    w1 = ndimage.maximum_filter(w1,11)
    w2 = ndimage.maximum_filter(w2,11)
    w3 = ndimage.maximum_filter(w3,11)
    bad = w1*w2*w3

    cond = w!=0
    cond2 = bad==0

# Perform statistics for scaling
    data_arr = np.array([data1,data2,data3])
    for data in data_arr:
        if data != None:
            tmp = data[cond]
            tmp.sort()
            min = tmp[tmp.size*0.1]
            max = tmp[tmp.size*0.995]
            data[data<min] = min
            data[data>max] = max
            data -= min
            data /= data[cond].max()
            data[cond2] = 0.

# Create the RGB array
    if data1 != None:
        data = scipy.empty((data1.shape[0],data1.shape[1],3))
    else:
        data = scipy.empty((data2.shape[0],data2.shape[1],3))
# Note that for some reason the Y-axis is flipped. The ::-1 syntax flips an
#   array about the axis.
    if data1 != None:
        data[:,:,0] = data1[::-1].copy()
    else:
        data[:,:,0] = np.zeros(np.shape(data1))
    if data2 != None:
        data[:,:,1] = data2[::-1].copy()
    else:
        data[:,:,1] = np.zeros(np.shape(data2))
    if data3 != None:
        data[:,:,2] = data3[::-1].copy()
    else:
        data[:,:,2] = np.zeros(np.shape(data3))
    pylab.imshow(data)
    pylab.show()
    if outfile != None: pylab.savefig(outfile)
