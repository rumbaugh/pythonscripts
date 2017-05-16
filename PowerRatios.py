import numpy as np
execfile('/home/rumbaugh/pythonscripts/calc_SF_Xray.py')

def a_m(data,m,R,xcen=None,ycen=None):
    xsize,ysize=np.shape(data)[-2],np.shape(data)[-1]
    if xcen==None:xcen=0.5*(xsize-1)
    if ycen==None:ycen=0.5*(ysize-1)
    xcoord,ycoord=np.arange(xsize)*np.ones(ysize).reshape((ysize,1))-xcen,np.ones(xsize)*np.arange(ysize).reshape((ysize,1))-ycen
    r=np.sqrt(xcoord**2+ycoord**2)
    phi=np.arctan(ycoord*1./xcoord)
    phi[(xcoord==0)&(ycoord>=0)]=np.pi/2
    phi[(xcoord==0)&(ycoord<0)]=-np.pi/2
    phi[xcoord<0]+=np.pi
    g=np.where((np.abs(xcoord)<R)&(np.abs(ycoord)<R))
    if len(np.shape(data))>2:
        out=np.zeros(np.shape(data)[0])
        for i in range(0,len(out)): out[i]=np.sum(data[i][g]*np.cos(m*phi[g])*r[g]**m)#*.492**2
        return out
    else:
        return np.sum(data[g]*np.cos(m*phi[g])*r[g]**m)#*.492**2

def b_m(data,m,R,xcen=None,ycen=None):
    xsize,ysize=np.shape(data)[-2],np.shape(data)[-1]
    if xcen==None:xcen=0.5*(xsize-1)
    if ycen==None:ycen=0.5*(ysize-1)
    xcoord,ycoord=np.arange(xsize)*np.ones(ysize).reshape((ysize,1))-xcen,np.ones(xsize)*np.arange(ysize).reshape((ysize,1))-ycen
    r=np.sqrt(xcoord**2+ycoord**2)
    phi=np.arctan(ycoord*1./xcoord)
    phi[(xcoord==0)&(ycoord>=0)]=np.pi/2
    phi[(xcoord==0)&(ycoord<0)]=-np.pi/2
    phi[xcoord<0]+=np.pi
    g=np.where((np.abs(xcoord)<R)&(np.abs(ycoord)<R))
    if len(np.shape(data))>2:
        out=np.zeros(np.shape(data)[0])
        for i in range(0,len(out)): out[i]=np.sum(data[i][g]*np.sin(m*phi[g])*r[g]**m)#*.492**2
        return out
    else:
        return np.sum(data[g]*np.sin(m*phi[g])*r[g]**m)#*.492**2

def PowerRatio(data,m,R,xcen=None,ycen=None):
    mask=np.zeros(np.shape(data)[-2:])
    xsize,ysize=np.shape(data)[-2],np.shape(data)[-1]
    xcoord,ycoord=np.arange(xsize)*np.ones(ysize).reshape((ysize,1)),np.ones(xsize)*np.arange(ysize).reshape((ysize,1))
    mask[np.sqrt((xcoord-0.5*(xsize-1))**2+(ycoord-0.5*(ysize-1))**2)<=R]=1
    data*=mask
    return (a_m(data,m,R,xcen,ycen)**2+b_m(data,m,R,xcen,ycen)**2)/(2.*m**2*R**(2*m))/(a_m(data,0,R,xcen,ycen)*np.log(R))**2

def PowerRatioError(data,m,R,r_smooth,xcen=None,ycen=None,ntrials=10):
    randcounts=np.random.poisson(data,((ntrials,)+np.shape(data)))
    randSF=np.zeros(np.shape(randcounts))
    for i in range(0,ntrials):
        randSF[i]=calc_SF_Xray(randcounts[i],'tophat',[r_smooth])
    randnoise=np.random.poisson(np.ones(np.shape(data))*np.sum(data),((ntrials,)+np.shape(data)))
    noiseSF=np.zeros(np.shape(randnoise))
    for i in range(0,ntrials):
        noiseSF[i]=calc_SF_Xray(randnoise[i],'tophat',[r_smooth])
    if np.shape(m)==():
        m=np.array([m])
    PLBs,PUBs=np.zeros(len(m)),np.zeros(len(m))
    for mcur,im in zip(m,np.arange(len(m))):
        Pcur=PowerRatio(randSF,mcur,R)
        PLBs[im],PUBs[im]=np.sort(Pcur)[int((0.5-0.5*0.682689492)*ntrials)],np.sort(Pcur)[int((0.5+0.5*0.682689492)*ntrials)]
        Pnoise=PowerRatio(noiseSF,mcur,R)
        PLBs[im]-=np.sort(Pnoise)[int((0.5-0.5*0.682689492)*ntrials)]
        PUBs[im]-=np.sort(Pnoise)[int((0.5+0.5*0.682689492)*ntrials)]
        print np.sort(Pnoise)[int((0.5-0.5*0.682689492)*ntrials)],np.sort(Pnoise)[int((0.5+0.5*0.682689492)*ntrials)]
    return PLBs,PUBs
        
