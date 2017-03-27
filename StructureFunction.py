import numpy as np
import matplotlib.pyplot as plt
import sys

def ConstructStructureFunctionArray(S,ltime,Serr=None,fullerror=True):
    Varr = np.zeros(0)
    tauarr = np.zeros(0)
    if len(S) != len(ltime): sys.exit("S and ltime don't have same dimensions.")
    if Serr != None: 
        Verr = np.zeros(0)
        if len(S) != len(Serr): sys.exit("S and Serr don't have same dimensions.")
    S1,S2=np.repeat(S,len(S)**2),np.tile(S,len(S)**2)
    t1,t2=np.repeat(ltime,len(ltime)**2),np.tile(ltime,len(ltime)**2)
    g12=np.where(t1>t2)[0]
    S1,S2,t1,t2=S1[g12],S2[g12],t1[g12],t2[g12]
    Varr,tauarr=(S1-S2)**2,t1-t2
    if Serr != None: 
        err1,err2=np.repeat(Serr,len(S)**2),np.tile(Serr,len(S)**2)
        err1,err2=err1[g12],err2[g12]
        if fullerror:
            Verr=2*np.abs(S1-S2)*np.sqrt(err1**2+err2**2)
        else:
            Verr=np.sqrt(err1**2+err2**2)
    as_tauarr = np.argsort(tauarr)
    SFarr = np.zeros((len(tauarr),2))
    if Serr != None: SFarr = np.zeros((len(tauarr),3))
    SFarr[:,0],SFarr[:,1] = tauarr[as_tauarr],Varr[as_tauarr]
    if Serr != None: SFarr[:,2] = Verr[as_tauarr]
    return SFarr

def CalcStructureFunction_IQR(S,ltime,nbins=10):
    Varr,tauarr = np.zeros(0),np.zeros(0)
    if len(S) != len(ltime): sys.exit("S and ltime don't have same dimensions.")
    S1,S2=np.repeat(S,len(S)**2),np.tile(S,len(S)**2)
    t1,t2=np.repeat(ltime,len(ltime)**2),np.tile(ltime,len(ltime)**2)
    g12=np.where(t1>t2)[0]
    S1,S2,t1,t2=S1[g12],S2[g12],t1[g12],t2[g12]
    Varr,tauarr=S1-S2,t1-t2
    as_tauarr = np.argsort(tauarr)
    SF_arr = np.zeros((len(tauarr),2))
    SF_arr[:,0],SF_arr[:,1] = tauarr[as_tauarr],Varr[as_tauarr]
    V0,V0days=np.nan,2
    while np.isnan(V0):
        SF0=np.sort(SF_arr[:,1][SF_arr[:,0]<V0days])
        V0=SF0[(3*len(SF0))/4]-SF0[len(SF0)/4]
        V0days*=2
    npairs = np.shape(SF_arr)[0]
    binsize = npairs/nbins
    binsize0 = binsize + np.mod(npairs,nbins)
    V_arr = np.zeros(nbins)
    tau_arr = np.zeros(nbins)
    V_temp=np.sort(SF_arr[:binsize0][:,1])
    V_arr[0] = V_temp[(3*len(V_temp))/4]-V_temp[len(V_temp)/4]
    tau_arr[0] = np.average(SF_arr[:binsize0][:,0])
    Varr_tmp,tauarr_tmp=np.sort(np.reshape(SF_arr[binsize0:][:,1],((nbins-1,binsize))),axis=1),np.reshape(SF_arr[binsize0:][:,0],((nbins-1,binsize)))
    V_arr[1:],tau_arr[1:]=Varr_tmp[:,(3*np.shape(Varr_tmp)[-1])/4]-Varr_tmp[:,np.shape(Varr_tmp)[-1]/4],np.average(tauarr_tmp,axis=1)
    return tau_arr,np.sqrt(0.549*(V_arr**2-V0**2))


def CalcStructureFunction_eq14(S,ltime,nbins=10):
    SF_arr=ConstructStructureFunctionArray(S,ltime)
    npairs = np.shape(SF_arr)[0]
    binsize = npairs/nbins
    binsize0 = binsize + np.mod(npairs,nbins)
    V_arr = np.zeros(nbins)
    tau_arr = np.zeros(nbins)
    V_arr[0] = np.average(SF_arr[:binsize0][:,1])
    tau_arr[0] = np.average(SF_arr[:binsize0][:,0])
    Varr_tmp,tauarr_tmp=np.reshape(SF_arr[binsize0:][:,1],((nbins-1,binsize))),np.reshape(SF_arr[binsize0:][:,0],((nbins-1,binsize)))
    V_arr[1:],tau_arr[1:]=np.average(Varr_tmp,axis=1),np.average(tauarr_tmp,axis=1)
    return tau_arr,V_arr

def CalcStructureFunction_eq15(S,ltime,Serr,nbins=10):
    SF_arr=ConstructStructureFunctionArray(S,ltime,Serr)
    SF_arr[:,2]*=SF_arr[:,2]
    npairs = np.shape(SF_arr)[0]
    binsize = npairs/nbins
    binsize0 = binsize + np.mod(npairs,nbins)
    V_arr = np.zeros(nbins)
    tau_arr = np.zeros(nbins)
    V_arr[0] = np.sqrt(np.average(SF_arr[:binsize0][:,1])+np.average(SF_arr[:binsize0][:,2]))
    tau_arr[0] = np.average(SF_arr[:binsize0][:,0])
    Varr_tmp,tauarr_tmp,errarr_tmp=np.reshape(SF_arr[binsize0:][:,1],((nbins-1,binsize))),np.reshape(SF_arr[binsize0:][:,0],((nbins-1,binsize))),np.reshape(SF_arr[binsize0:][:,2],((nbins-1,binsize)))
    V_arr[1:],tau_arr[1:]=np.sqrt(np.average(Varr_tmp,axis=1)+np.average(errarr_tmp,axis=1)),np.average(tauarr_tmp,axis=1)
    return tau_arr,V_arr

def CalcStructureFunction_eq16(S,ltime,Serr,nbins=10):
    SF_arr=ConstructStructureFunctionArray(S,ltime,Serr)
    SF_arr[:,2]*=SF_arr[:,2]
    npairs = np.shape(SF_arr)[0]
    binsize = npairs/nbins
    binsize0 = binsize + np.mod(npairs,nbins)
    V_arr = np.zeros(nbins)
    tau_arr = np.zeros(nbins)
    V_arr[0] = np.sqrt(0.5*np.pi*np.average(SF_arr[:binsize0][:,1])+np.average(SF_arr[:binsize0][:,2]))
    tau_arr[0] = np.average(SF_arr[:binsize0][:,0])
    Varr_tmp,tauarr_tmp,errarr_tmp=np.reshape(SF_arr[binsize0:][:,1],((nbins-1,binsize))),np.reshape(SF_arr[binsize0:][:,0],((nbins-1,binsize))),np.reshape(SF_arr[binsize0:][:,2],((nbins-1,binsize)))
    V_arr[1:],tau_arr[1:]=np.sqrt(0.5*np.pi*np.average(Varr_tmp,axis=1)+np.average(errarr_tmp,axis=1)),np.average(tauarr_tmp,axis=1)
    return tau_arr,V_arr

def CalcStructureFunction_eq17(S,ltime,nbins=10):
    SF_arr=ConstructStructureFunctionArray(S,ltime)
    npairs = np.shape(SF_arr)[0]
    binsize = npairs/nbins
    binsize0 = binsize + np.mod(npairs,nbins)
    V_arr = np.zeros(nbins)
    tau_arr = np.zeros(nbins)
    V_arr[0] = np.median(SF_arr[:binsize0][:,1])
    tau_arr[0] = np.average(SF_arr[:binsize0][:,0])
    Varr_tmp,tauarr_tmp=np.reshape(SF_arr[binsize0:][:,1],((nbins-1,binsize))),np.reshape(SF_arr[binsize0:][:,0],((nbins-1,binsize)))
    V_arr[1:],tau_arr[1:]=np.median(Varr_tmp,axis=1),np.average(tauarr_tmp,axis=1)
    return tau_arr,V_arr

def CalcStructureFunction_eq18(S,ltime,Serr,nbins=10):
    SF_arr=ConstructStructureFunctionArray(S,ltime,Serr)
    SF_arr[:,1],SF_arr[:,2]=np.sqrt(SF_arr[:,1]),np.sqrt(SF_arr[:,2])
    npairs = np.shape(SF_arr)[0]
    binsize = npairs/nbins
    binsize0 = binsize + np.mod(npairs,nbins)
    V_arr = np.zeros(nbins)
    tau_arr = np.zeros(nbins)
    V_arr[0] = np.average(np.sqrt(0.5*np.pi)*SF_arr[:binsize0][:,1]+SF_arr[:binsize0][:,2])
    tau_arr[0] = np.average(SF_arr[:binsize0][:,0])
    Varr_tmp,tauarr_tmp,errarr_tmp=np.reshape(SF_arr[binsize0:][:,1],((nbins-1,binsize))),np.reshape(SF_arr[binsize0:][:,0],((nbins-1,binsize))),np.reshape(SF_arr[binsize0:][:,2],((nbins-1,binsize)))
    V_arr[1:],tau_arr[1:]=np.average(np.sqrt(0.5*np.pi)*Varr_tmp+errarr_tmp,axis=1),np.average(tauarr_tmp,axis=1)
    return tau_arr,V_arr


    

def CalcStructureFunction(S,ltime,Serr=None,nbins=10,doplot=True,plotfile=None,ylimits=None,plotlog=True,output=False,output_errors=False,outfile=None,noiseSF=False,ntrials=None):
    SF_arr = ConstructStructureFunctionArray(S,ltime,Serr)
    npairs = np.shape(SF_arr)[0]
    if ((noiseSF) & (ntrials == None)):
        if Serr == None: sys.exit('Serr must be input for noiseSF_S85')
        noisevar = np.var(Serr)
        #next two lines must be in that order
        SF_arr[:,1] += 2*noisevar
        SF_arr[:,2] = np.sqrt(8*noisevar*SF_arr[:,1]/npairs)
    elif ntrials != None:
        nSFarrTot = np.zeros((ntrials,npairs))
        noisevar = np.var(Serr)
        for iteration in range(0,ntrials):
            shuff = np.copy(Serr)
            np.random.shuffle(shuff)
            Serrtmp = shuff+np.random.normal(0.,noisevar,len(Serr))
            SFtmp = ConstructStructureFunctionArray(Serrtmp,ltime)
            nSFarrTot[iteration] = SFtmp[:,1]
        noiseSF_arr = np.average(nSFarrTot,axis=0)
        SF_arr[:,1] -= noiseSF_arr
        SF_arr[:,2] = np.sqrt(np.average(noiseSF_arr**2))*np.ones(len(SF_arr[:,2]))
        print SF_arr[:,1],noiseSF_arr,SF_arr[:,2]
    binsize = npairs/nbins
    binsize0 = binsize + np.mod(npairs,nbins)
    V_arr = np.zeros(nbins)
    tau_arr = np.zeros(nbins)
    V_arr[0] = 0.5*np.average(SF_arr[:binsize0][:,1])
    tau_arr[0] = np.average(SF_arr[:binsize0][:,0])
    if Serr != None: 
        Verr_arr = np.zeros(nbins)
        Verr_arr[0] = 0.5*np.sqrt(np.sum(SF_arr[:binsize0][:,2]**2))/binsize0
        terr_arr = np.zeros(nbins)
        terr_arr[0] = np.std(SF_arr[:binsize0][:,0])
    for iv in range(1,nbins):
        V_arr[iv] = 0.5*np.average(SF_arr[binsize0+(iv-1)*binsize:binsize0+iv*binsize][:,1])
        tau_arr[iv] = np.average(SF_arr[binsize0+(iv-1)*binsize:binsize0+iv*binsize][:,0])
        if Serr != None: 
            Verr_arr[iv] = 0.5*np.sqrt(np.sum(SF_arr[binsize0+(iv-1)*binsize:binsize0+iv*binsize][:,2]**2))/binsize
            terr_arr[iv] = np.std(SF_arr[binsize0+(iv-1)*binsize:binsize0+iv*binsize][:,0])
    #V_arr[len(V_arr)-1] = 0.5*np.average(SF_arr[npairs-binsize0:][:,1])
    #tau_arr[len(V_arr)-1] = np.average(SF_arr[npairs-binsize0:][:,0])
    #for iv in range(0,nbins-1):
    #    V_arr[iv] = 0.5*np.average(SF_arr[(iv)*binsize:(iv+1)*binsize])
    #    tau_arr[iv] = np.average(SF_arr[(iv)*binsize:(iv+1)*binsize])
    if doplot:
        plt.figure(1)
        plt.clf()
        if plotlog:
            if Serr != None:
                plt.errorbar(np.log10(tau_arr),np.log10(V_arr),xerr=np.log10(np.e)*terr_arr/tau_arr,yerr=np.log10(np.e)*Verr_arr/V_arr,fmt='ro',lw=2,capsize=3,mew=1)
            plt.scatter(np.log10(tau_arr),np.log10(V_arr))
            plt.xlabel('log(' + r'$\tau$)')
            plt.ylabel('log(V)')
        else: 
            if Serr != None: plt.errorbar(tau_arr,V_arr,xerr=terr_arr,yerr=Verr_arr,fmt='ro',lw=2,capsize=3,mew=1)
            plt.scatter(tau_arr,V_arr)
        if ylimits != None: plt.ylim(0,ylimits)
        if plotfile != None: plt.savefig(plotfile)
    if outfile != None:
        SFFILE = open(outfile,'w')
        SFFILE.write('# V(tau) tau\n')
        for iv in range(0,len(V_arr)): SFFILE.write('%f %f\n'%(V_arr[iv],tau_arr[iv]))
        SFFILE.close()
    if output: 
        if output_errors:
            return tau_arr,V_arr,terr_arr,Verr_arr
        return tau_arr,V_arr 
