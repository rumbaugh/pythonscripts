import math as m
import numpy as np
import cmath
import sys

def nextpow2(num): 
    n = 2 
    while n < num: 
        n = n * 2 
    return n 

def bitrev(x): 
    N, x = len(x), x[:] 
    if N != nextpow2(N): 
        raise ValueError, 'N is not power of 2' 
    for i in range(N): 
        k, b, a = 0, N>>1, 1 
        while b >= a: 
            if b & i: k = k | a 
            if a & i: k = k | b 
            b, a = b>>1, a<<1 
        if i < k: # important not to swap back 
            x[i], x[k] = x[k], x[i] 
    return x 

def fft(fx,isign):
    nn = 0L
    nn = len(fx)
    if nn != nextpow2(nn):
        raise ValueError, 'Input array length should be a power of 2'
    n = nn << 1
    fx = bitrev(fx)
    x = np.zeros(2*nn)
    for i in range(0L,nn):
        x[2*i],x[2*i+1] = fx[i].real,fx[i].imag
    mmax = 2
    while(n > mmax):
        istep = mmax << 1
        theta = isign*(2*m.pi/mmax)
        wtemp = m.sin(0.5*theta)
        wpr = -2.0*wtemp**2
        wpi = m.sin(theta)
        wr = 1.0
        wi = 0.0
        for mi in range(1L,mmax,2):
            for i in range(mi,n+1,istep):
                i2 = i + mmax
                tempr = wr*x[i2-1]-wi*(x[i2])
                tempi = wr*x[i2]+wi*(x[i2-1])
                x[i2-1] = x[i-1] -tempr
                x[i2]=x[i]-tempi
                x[i-1] += tempr
                x[i] += tempi
            wtemp = wr
            wr = wr*wpr-wi*wpi+wr
            wi = wi*wpr+wtemp*wpi+wi
        mmax = istep
    xout = [0]*nn
    for i in range(0,nn): xout[i] = x[2*i]+x[2*i+1]*1j
    return xout

def fftS(fx):
#    if len(fx) != len(x): sys.exit("Inputs to FFT must be of same length")
    N = len(fx)
    if N == 1: return fx
    even = fftS(fx[0::2]) # slice notation, very convenient
    odd =  fftS(fx[1::2])
    M = N/2.0
    expBase = 2j*m.pi/N
    l = [even[k] + cmath.exp(expBase*k) * odd[k] for k in range(0L,int(M))]
    r = [even[k] - cmath.exp(expBase*k) * odd[k] for k in range(0L,int(M))]
 
    return l + r

def fft2D(datain,nn,isign):
    data = np.zeros(2*len(datain))
    for i in range(0L,len(datain)):
        data[2*i],data[2*i+1]=datain[i].real,datain[i].imag
    ntot = nn**2
    nprev = 1
    for idim in range(2,0,-1):
        n = nn
        nrem = ntot/(n*nprev)
        ip1 =nprev << 1
        ip2 = ip1 *n
        ip3 = ip2*nrem
        i2rev=1
        for i2 in range(1L,ip2+1,ip1):
            if i2 < i2rev:
                for i1 in range(i2,i2+ip1-1,2):
                    for i3 in range(i1,ip3+1,ip2):
                        data[i3-1],data[i3rev-1]=data[i3rev-1],data[i3-1]
            ibit=ip2>>1
            while ((ibit >= ip1) & (i2rev > ibit)):
                i2rev -= ibit
                ibit >>= 1
            i2rev += ibit
        ifp1=ip1
        while (ifp1 < ip2):
            ifp2=ifp1 << 1
            theta = isign*2*m.pi/(ifp2/(1.0*ip1))
            wtemp=m.sin(0.5*theta)
            wpr = -2.0*wtemp*wtemp
            wpi = m.sin(theta)
            wr,wi = 1.0,0.0
            for i3 in range(1L,ifp1+1,ip1):
                for i1 in range(i3,i3+ip1-1,2):
                    for i2 in range(i1,ip3+1,ifp2):
                        k1=i2
                        k2=k1+ifp1
                        tempr=wr*data[k2-1]-wi*data[k2]
                        tempi=wr*data[k2]+wi*data[k2-1]
                        data[k2-1]=data[k1-1]-tempr
                        data[k2]=data[k1]-tempi
                        data[k1-1]+= tempr
                        data[k1] += tempi
                wtemp=wr
                wr=wr*wpr-wi*wpi+wr
                wi *= wpr+wtemp+wpi
            ifp1=ifp2
        nprev*=n
        ic = 1
        if isign < 0: ic = 0
        for i in range(0L,len(datain)):
            datain[i]=data[2*i]+data[2*i+1]*1j*ic
    return datain

def Conv2d(xin,resp,nn):
    #Convolves xin with response function resp. However, resp must be 
    #the fourier transform of the response function. nn is the dimensions
    # of the arrays. They must be square. 
    xFT = fft2D(xin,nn,1)
    multArr = xFT*resp
    xOut1 = fft2D(multArr,nn,-1)
    xOut = np.zeros((nn,nn))
    for i in range(0,nn):
        for j in range(0,nn):
            xOut[i,j] = xOut1[i+nn*j]
    return xOut
