import numpy as np
import matplotlib.text as txt
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.backends.backend_pdf as bpdf
import itertools as it
import pyfits as py

def plot_flux(ax,fluxes,fluxerrs=None,label=None,curcol='k',bands=np.array(['g','r','i','z'])):
    bcens={'u': 3876.63943790537, 'g': 4841.83358196563, 'r': 6438/534828217, 'i': 7820.99282740933, 'z': 9172.34266385718, 'Y': 9877.80238651117}
    if len(fluxes)!=len(bands):
        print 'Lengths of fluxes and bands must be equal'
        return
    cens=np.zeros(len(bands))
    for ib,b in zip(np.arange(len(bands)),bands): cens[ib]=bcens[b] 
    if 'i' in bands:
        gi=np.where(bands=='i')[0][0]
        refpnt=fluxes[gi]
        if fluxerrs!=None:fluxerrs/=fluxes[gi]
        fluxes/=fluxes[gi]
    elif 'r' in bands:
        gi=np.where(bands=='r')[0][0]
        refpnt=fluxes[gi]
        if fluxerrs!=None:fluxerrs/=fluxes[gi]
        fluxes/=fluxes[gi]
    else:
        return
    if label==None:
        if fluxerrs!=None:ax.errorbar(cens,fluxes,yerr=fluxerrs,color=curcol,fmt=' ',lw=2,capsize=3,mew=1,zorder=3,label=None)
        ax.scatter(cens,fluxes,color=curcol,s=36,zorder=4)
    else:
        ax.scatter(cens,fluxes,color=curcol,s=36,label=label,zorder=4)
        if fluxerrs!=None:ax.errorbar(cens,fluxes,yerr=fluxerrs,color=curcol,fmt=' ',lw=2,capsize=3,mew=1,zorder=3,label=None)

def calc_flux(mjd,mag,magerr,cbands,band,mjdcen):
    gband=np.where(cbands==band)[0]
    magplot=mag[gband]
    magploterr=magerr[gband]
    g100=np.where((magplot<100)&(np.isnan(magplot)==False))[0]
    if len(g100)>1:
        gcen=np.argsort(np.abs(mjd[gband][g100]-mjdcen))[0]
        medmag=magplot[g100][gcen]
        mederr=magploterr[g100][gcen]
    else:
        medmag=magplot[g100]
        mederr=magploterr[g100]
    if np.shape(medmag)!=(): 
        if len(medmag)==0:
            return 0,0
        medmag=medmag[0]
        mederr=mederr[0]
    return 10**(medmag/-2.5),np.abs(np.log(10)/-2.5*10**(medmag/-2.5)*mederr)

def plot_band(ax,mjd,mag,magerr,cbands,band,connectpoints=True,nolabels=False,outlierarr=None,overridecolor=None):
    gband=np.where(cbands==band)[0]
    magplot=mag[gband]
    magploterr=magerr[gband]
    g100=np.where(magplot<100)[0]
    coldict={'g': 'green','r': 'red', 'i': 'magenta', 'z': 'blue', 'Y': 'cyan'}
    try:
        curcol=coldict[band]
    except KeyError:
        print '%s is not a valid band'%band
        return
    if overridecolor!=None: curcol=overridecolor
    gsort=np.argsort(mjd[gband][g100])
    if connectpoints:
        ax.plot(mjd[gband][g100][gsort],magplot[g100][gsort],color=curcol,lw=2,zorder=15)
    ax.errorbar(mjd[gband][g100][gsort],magplot[g100][gsort],yerr=magploterr[g100],color=curcol,fmt=' ',lw=2,capsize=3,mew=1,zorder=-30)
    if nolabels:
        ax.scatter(mjd[gband][g100],magplot[g100],color=curcol,zorder=20,edgecolor='k')
    else:
        ax.scatter(mjd[gband][g100],magplot[g100],color=curcol,label=band,zorder=20,edgecolor='k')
    try:
        for outlier_set in outlierarr.keys():
            gout=np.intersect1d(np.arange(len(mag))[outlierarr[outlier_set]['g']],gband[g100])
            ax.scatter(mjd[gout],mag[gout],color=outlierarr[outlier_set]['c'],marker=outlierarr[outlier_set]['marker'],lw=2,s=100,zorder=100,edgecolor='k')
    except:
        if np.shape(outlierarr)!=():
            gout=np.intersect1d(np.arange(len(mag))[outlierarr],gband[g100])
            ax.scatter(mjd[gout],mag[gout],color='r',marker='x',lw=3,s=50,zorder=50)
    #return


def plot_lightcurve(dbid,mjd,mag,magerr,bands,survey,trueredshift,DBdir,psfpage=None,fname=None,DESfname=None,connectpoints=True,specfile=None,WavLL=3000,WavUL=10500,outlierflag=0,zoominband=None,outlierarr=None):
    try:
        test1=mjd.index
        plotmacleod=True
        mjd,mag,magerr,bands,mjdmc,magmc,magerrmc,bandsmc=mjd[0],mag[0],magerr[0],bands[0],mjd[1],mag[1],magerr[1],bands[1]
    except:
        plotmacleod=False
    bcens={'u': 3876.63943790537, 'g': 4841.83358196563, 'r': 6438/534828217, 'i': 7820.99282740933, 'z': 9172.34266385718, 'Y': 9877.80238651117}
    VBfile='%s/VanderBerk_datafile1.txt'%DBdir
    crv=np.loadtxt(VBfile,skiprows=23)
    redshift=np.copy(trueredshift)
    if redshift<0:redshift=0
    gdes,gsdss,gposs=np.where(survey=='DES')[0],np.where(survey=='SDSS')[0],np.where(survey=='POSS')[0]
    POSSbands=np.array(['g','r','i'])
    if len(gposs)>0:
        POSSmagdict,POSSmagerrdict,POSSmjddict={b: mag[gposs][bands[gposs]==b] for b in POSSbands},{b: np.zeros(0) for b in POSSbands},{b: np.zeros(0) for b in POSSbands}
        for band in POSSbands:
            if len(POSSmagdict[band])>0:
                POSSmagdict[band]=np.array([np.median(POSSmagdict[band])])
                POSSmagerrdict[band]=np.array([np.mean(POSSmagdict[band])])
        if (len(POSSmagdict['g'])>0)&(len(POSSmagdict['r'])>0):
            mag[gposs][bands[gposs]=='g'],mag[gposs][bands[gposs]=='r']=mag[gposs][bands[gposs]=='g']+0.392*(POSSmagdict['g']-POSSmagdict['r'])-0.28, mag[gposs][bands[gposs]=='r'] +0.127*(POSSmagdict['g']-POSSmagdict['r'])+0.1
            magerr[gposs][bands[gposs]=='g'],magerr[gposs][bands[gposs]=='r']=np.sqrt(1.392**2*magerr[gposs][bands[gposs]=='g']**2+0.392**2*magerr[gposs][bands[gposs]=='r']**2),  np.sqrt(magerr[gposs][bands[gposs]=='r']**2+0.127**2*(magerr[gposs][bands[gposs]=='g']**2+magerr[gposs][bands[gposs]=='r']**2))
        else: 
            if len(POSSmagdict['g'])>0: mag[gposs][bands[gposs]=='g']=0
            if len(POSSmagdict['r'])>0: mag[gposs][bands[gposs]=='r']=0
        if (len(POSSmagdict['i'])>0)&(len(POSSmagdict['r'])>0):   
            mag[gposs][bands[gposs]=='i']=mag[gposs][bands[gposs]=='i']+0.27*(POSSmagdict['r']-POSSmagdict['i'])+0.32
            magerr[gposs][bands[gposs]=='i']=np.sqrt(magerr[gposs][bands[gposs]=='i']**2+0.27**2*(magerr[gposs][bands[gposs]=='r']**2+magerr[gposs][bands[gposs]=='i']**2))
        else:
            if len(POSSmagdict['i'])>0: mag[gposs][bands[gposs]=='i']=0
    bestdiff={b: {'diff': 0, 'ihi': 0, 'ilo': 0} for b in ['g','r','i','z']}
    #for b in ['g','r','i','z']:
    for b in ['g']:
        gb=np.where((bands==b)&(survey!='POSS'))[0]
        #if ((len(gsdss)>0)&(len(gdes)>0)):
        #    gsdssb,gdesb=np.where(bands[gsdss]==b)[0],np.where(bands[gdes]==b)[0]
        #    if ((len(gsdssb)>0)&(len(gdesb)>0)):
        if len(gb)>0:
            magpairs=np.zeros([len(gb)*len(gb),2])
            magpairs[:,1],magpairs[:,0]=np.repeat(mag[gb],len(gb)),np.tile(mag[gb],len(gb))
            magerrpairs=np.zeros([len(gb)*len(gb),2])
            magerrpairs[:,1],magerrpairs[:,0]=np.repeat(magerr[gb],len(gb)),np.tile(magerr[gb],len(gb))
            magdiffs,differrs=magpairs[:,0]-magpairs[:,1],np.max(magerrpairs,axis=1)#,np.sqrt(np.sum(magerrpairs**2,axis=1))
            ipairs=np.zeros([len(gb)*len(gb),2])
            ipairs[:,1],ipairs[:,0]=np.repeat(gb,len(gb)),np.tile(gb,len(gb))
            diffsigs=magdiffs/differrs
            gsig=np.where(differrs<0.15)[0]
            if len(gsig)>0:
                gmax=np.argsort(magdiffs[gsig])[-1]
                imax,imin=ipairs[:,1][gsig[gmax]],ipairs[:,0][gsig[gmax]]
                bestdiff[b]['diff']=np.max(magdiffs[gsig])
                bestdiff[b]['ihi'],bestdiff[b]['ilo']=imax,imin
    fig=plt.figure(1)
    fig.clf()
    ax3=plt.subplot2grid((2,10),(1,0),colspan=6)
    plt.rc('axes',linewidth=2)
    plt.fontsize = 14
    plt.tick_params(which='major',length=8,width=2,labelsize=14)
    plt.tick_params(which='minor',length=4,width=1.5,labelsize=14)
    plt.locator_params(nbins=4)
    if zoominband==None:
        totdiffs=np.zeros(4)
        for ib,b in zip(np.arange(4),['g','r','i','z']):
            totdiffs[ib]=bestdiff[b]['diff']
        ibest=np.argsort(totdiffs)[-1]
        bbest=['g','r','i','z'][ibest]
        gbbest=np.where(bands==bbest)[0]
        imax,imin=bestdiff[bbest]['ihi'],bestdiff[bbest]['ilo']
        maxfluxes,minfluxes=np.zeros(4),np.zeros(4)
        maxfluxerrs,minfluxerrs=np.zeros(4),np.zeros(4)
        for ib,b in zip(np.arange(4),['g','r','i','z']):
            gbt=np.where(bands==b)[0]
            maxfluxes[ib],maxfluxerrs[ib]=calc_flux(mjd,mag,magerr,bands,b,mjd[imax])
            minfluxes[ib],minfluxerrs[ib]=calc_flux(mjd,mag,magerr,bands,b,mjd[imin])
        if specfile!=None:
            try:
                shdu=py.open(specfile)
                specdata=shdu[1].data
                sflux,swav=specdata['flux'],10**(specdata['loglam'])
                s_closei=np.where(np.abs(swav-bcens['i'])<20)[0]
                smwid=10
                swav=swav[smwid:-smwid]
                normflux=sflux/np.mean(sflux[s_closei])
                smoothflux=[np.mean(normflux[x-smwid:x+smwid+1]) for x in np.arange(smwid,len(normflux)-smwid)]
                ax3.plot(swav,smoothflux,lw=1,color='magenta',zorder=1)
            except IOError:
                specfile=None
        v_closei=np.where(np.abs(crv[:,0]*(1.+redshift)-bcens['i'])<20)[0]
        gvrange=np.where((crv[:,0]*(1.+redshift)>WavLL)&(crv[:,0]*(1.+redshift)<WavUL))[0]
        if len(gvrange)>0:
            vmax=np.max(crv[:,1][gvrange]/np.mean(crv[:,1][v_closei]))
        else:
            vmax=np.max(crv[:,1])
        if trueredshift>0:ax3.plot(crv[:,0]*(1.+redshift),crv[:,1]/np.mean(crv[:,1][v_closei]),color='k',lw=1,zorder=2)
        plot_flux(ax3,maxfluxes,maxfluxerrs,label='Max',curcol='r')
        maxplot=np.max(maxfluxes)
        plot_flux(ax3,minfluxes,minfluxerrs,label='Min',curcol='b')
        if np.max(minfluxes)>maxplot:maxplot=np.max(minfluxes)
        ax3.set_xlabel('Wavelength (A)')
        ax3.set_ylabel('Flux (Arb. Units)')
        #ax3.legend(loc='lower right')
        plt.xlim(WavLL,WavUL)
        plt.ylim(-0.05,vmax*1.05)
        if trueredshift<=0: 
            plt.ylim(-0.05,1.15*maxplot)
    else:
        plot_band(ax3,mjd,mag,magerr,bands,zoominband,connectpoints=connectpoints,nolabels=False,outlierarr=outlierarr)
        if plotmacleod:
            plot_band(ax3,mjdmc,magmc,magerrmc,bandsmc,zoominband,connectpoints=connectpoints,nolabels=False,overridecolor='magenta')
        plt.axvline(mjd[imax],ls='dashed',lw=1,color='r')
        plt.axvline(mjd[imin],ls='dashed',lw=1,color='b')
        ylim=plt.ylim()
        if ylim[1]>30:
            ylim=(ylim[0],np.max(mag)+0.1)
        if ylim[1]>30: ylim=(ylim[0],30)
        if ylim[0]<15:
            ylim=(np.min(mag)-0.1,ylim[1])
        if ylim[0]<15: ylim=(15,ylim[1])
        plt.ylim(ylim[1],ylim[0])
        ax3.set_xlabel('MJD')
        ax3.set_ylabel('%s_PSF'%zoominband)
    ax1=plt.subplot2grid((2,10),(0,0),colspan=6)
    plt.rc('axes',linewidth=2)
    plt.fontsize = 14
    plt.tick_params(which='major',length=8,width=2,labelsize=14)
    plt.tick_params(which='minor',length=4,width=1.5,labelsize=14)
    plt.locator_params(nbins=4)
    for b in ['g','r','i','z','Y']:
        plot_band(ax1,mjd,mag,magerr,bands,b,connectpoints=connectpoints,nolabels=False)
    xlim=plt.xlim()
    plt.xlim(xlim[0],xlim[1]+0.33*(xlim[1]-xlim[0]))
    ylim=plt.ylim()
    plt.axvline(mjd[imax],ls='dashed',lw=1,color='r')
    plt.axvline(mjd[imin],ls='dashed',lw=1,color='b')
    if ylim[1]>30:
        ylim=(ylim[0],np.max(mag)+0.1)
    if ylim[1]>30: ylim=(ylim[0],30)
    if ylim[0]<15:
        ylim=(np.min(mag)-0.1,ylim[1])
    if ylim[0]<15: ylim=(15,ylim[1])
    plt.ylim(ylim[1],ylim[0])
    ax1.legend()
    xlim=plt.xlim()
    if outlierflag==1:
        ax1.text(0.5*(xlim[0]+xlim[1]),15./16*ylim[0]+ylim[1]/16.,'OUTLIER',color='r',horizontalalignment='center')
    elif outlierflag==2:
        ax1.text(0.5*(xlim[0]+xlim[1]),15./16*ylim[0]+ylim[1]/16.,'BAD PHOTOMETRY',color='r',horizontalalignment='center')
    ax1.set_xlabel('MJD')
    ax1.set_ylabel('Mag_PSF')
    if redshift>0:
        ax1.set_title('%s, z=%.4f'%(dbid,trueredshift))
    else:
        ax1.set_title(dbid)
    if not(DESfname in [None,'False']):
        ax4=plt.subplot2grid((2,10),(1,6),colspan=4,xticks=[],yticks=[])
        img4=mpimg.imread('%s/imagestamps/%s'%(DBdir,DESfname))
        ax4.imshow(img4)
    if len(gsdss)>0:
        SDSSfname='%s/imagestamps/%s_SDSScutout.jpeg'%(DBdir,dbid)
        try:
            img3=mpimg.imread(SDSSfname)
            ax3=plt.subplot2grid((2,10),(0,6),colspan=4,xticks=[],yticks=[])
            ax3.imshow(img3)
        except:
            pass
    if psfpage!=None:plt.savefig(psfpage,format='pdf')
    return

def DES2SDSS_gr(g,r):
    return (133625*g-9375*r-218)/124250.,(69.*g+925*r)/994.+516./62125.

def DES2SDSS_iz(i,z):
    return (-89*np.sqrt(-96000*i+96000*z+181561)+8000*z+37827)/8000.,(-17*np.sqrt(-96000*i+96000*z+181561)+24000*z+6731)/24000.

def plot_DB_lightcurves(DBIDs,outputfile,DBdir='/data2/rumbaugh/var_database/Y3A1',WavLL=3000,WavUL=10500,convertDESmags=False,outlierflags=None,zoominband=None,calc_outliers=False,outlier_window=300,outlier_thresh=0.5,plotmacleod=False,connectpoints=True,load_macleod=False,load_outliers=False):
    if np.shape(outlierflags)==(): outlierflags=np.zeros(len(DBIDs))
    hdu=py.open('%s/masterfile.fits'%DBdir)
    data=hdu[1].data
    psfpage=bpdf.PdfPages(outputfile)
    crdescutout=np.loadtxt('%s/imagestamps/DESimagestamp_index.dat'%DBdir,dtype={'names':('DBID','ra','dec','tile','fname'),'formats':('|S32','f8','f8','|S12','|S25')})

    crdb=np.loadtxt('%s/database_index.dat'%DBdir,dtype={'names':('DBID','CID','thingid','sdr7id','MQrownum','SPrownum','SDSSNAME'),'formats':('|S64','i8','i8','|S24','i8','i8','|S64')})
    crpmf1=np.loadtxt('%s/MQ_SDSS_y3a1_tidplatemjdfiber.csv'%DBdir,delimiter=',',skiprows=1,dtype={'names':('tid','plate','mjd','fiber'),'formats':('i8','i8','i8','i8')})
    crpmf2=np.loadtxt('%s/DR7_BH_SDSS_tidplatemjdfiber.csv'%DBdir,delimiter=',',skiprows=1,dtype={'names':('tid','plate','mjd','fiber'),'formats':('i8','i8','i8','i8')})
    crpmf=np.concatenate((crpmf1,crpmf2),axis=0)

    coldict={'g': 'green','r': 'red', 'i': 'magenta', 'z': 'blue', 'Y': 'cyan'}
    SDSSbands=np.array(['u','g','r','i','z'])
    SDSS_colnames={b:'%s_SDSS'%b for b in SDSSbands}
    POSSbands=np.array(['g','r','i'])

    for DBID,idb in zip(DBIDs,np.arange(len(DBIDs))):
        print DBID
        gdc=np.where(crdescutout['DBID']==DBID)[0]
        DESfname=None
        if len(gdc)>0:
            if crdescutout['fname'][gdc[0]]!='False':
                DESfname='%s.tif'%(crdescutout['fname'][gdc[0]])
        gdb=np.where(crdb['DBID']==DBID)[0][0]
        gmf=np.where(data['DatabaseID']==DBID)[0]
        if ((len(gmf)==0)&(DBID[:2]!='MQ')&(crdb['MQrownum'][gdb]>-1)):
            print len(gmf),DBID[:2],crdb['MQrownum'][gdb]
            gmf=np.where(data['DatabaseID']=='MQ%i'%crdb['MQrownum'][gdb])[0]
        if ((len(gmf)==0)&(DBID[:2]=='DR')&(crdb['SPrownum'][gdb]>-1)):
            gmf=np.where(data['DatabaseID']=='SDSSPOSS%i'%crdb['SPrownum'][gdb])[0]
        try:
            gmf=gmf[0]
            trueredshift=data['Redshift'][gmf]
        except:
            trueredshift=0
        redshift=np.copy(trueredshift)
        if redshift<0: redshift=0
        tid=crdb['thingid'][gdb]
        plate,pmf_mjd,fiber=-1,-1,-1
        if tid!=0:
            gpmf=np.where(crpmf['tid']==tid)[0]
            if len(gpmf)>0:
                gpmf=gpmf[0]
                plate,pmf_mjd,fiber=crpmf['plate'][gpmf],crpmf['mjd'][gpmf],crpmf['fiber'][gpmf]
        if (plate>-1)&(pmf_mjd>-1)&(fiber>-1):
            specfile='%s/spec/spec-%04i-%05i-%04i.fits'%(DBdir,plate,pmf_mjd,fiber)
        else:
            specfile=None
        cr=np.loadtxt('%s/%s/LC.tab'%(DBdir,DBID),dtype={'names':('DatabaseID','Survey','SurveyCoaddID','SurveyObjectID','RA','DEC','MJD','TAG','BAND','MAGTYPE','MAG','MAGERR','FLAG'),'formats':('|S64','|S20','|S20','|S20','f8','f8','f8','|S20','|S12','|S12','f8','f8','i8')},skiprows=1)
        if ((plotmacleod)|(load_macleod)):
            try:
                crmac=np.loadtxt('%s/%s/Macleod_LC.tab'%(DBdir,DBID),dtype={'names':('DatabaseID','RA','DEC','MJD','BAND','MAG','MAGERR','FLAG'),'formats':('|S24','f8','f8','f8','|S4','f8','f8','i8')})
                if not(load_macleod):crmac=crmac[(crmac['MAG']>0)&(crmac['MAG']<30)&(crmac['MAGERR']<5)]
                if load_outliers:
                    try:
                        crout=np.loadtxt('%s/%s/outliers.tab'%(DBdir,DBID),dtype='i8')
                    except:
                        crout=np.zeros(len(cr))
                    if load_macleod:
                        try:
                            croutmac=np.loadtxt('%s/%s/outliers_Macleod.tab'%(DBdir,DBID),dtype='i8')
                        except:
                            try:
                                croutmac=-np.ones(len(crmac))
                                crmac=crmac[croutmac>-1]
                                croutmac=croutmac[croutmac>-1]
                                outlier_arr=np.array(np.append(crout,croutmac),dtype='bool')
                            except IndexError:
                                croutmac=-np.ones(1)
                                crmac=crmac[croutmac>-1]
                                croutmac=croutmac[croutmac>-1]
                                outlier_arr=np.array(np.append(crout,croutmac),dtype='bool')
                            except TypeError:
                                croutmac=np.ones(0)
                                outlier_arr=np.array(crout,dtype='bool')
                            
                    else:
                        outlier_arr=np.array(crout,dtype='bool')
            except:
                crmac=None
                if load_outliers:
                    try:
                        crout=np.loadtxt('%s/%s/outliers.tab'%(DBdir,DBID),dtype='i8')
                    except:
                        crout=np.zeros(len(cr))
                    outlier_arr=np.array(crout,dtype='bool')
        else:
            crmac=None
            if load_outliers:
                try:
                    crout=np.loadtxt('%s/%s/outliers.tab'%(DBdir,DBID),dtype='i8')
                except:
                    crout=np.zeros(len(cr))
                outlier_arr=np.array(crout,dtype='bool')
        if ((calc_outliers)&(not(load_outliers))): 
            outlier_arr=np.zeros(len(cr),dtype='bool')
        if load_macleod:
            if crmac!=None:
                newcr=np.zeros((len(cr),),dtype={'names':('DatabaseID','RA','DEC','MJD','BAND','MAG','MAGERR','FLAG'),'formats':('|S24','f8','f8','f8','|S4','f8','f8','i8')})
                newcr['DatabaseID'],newcr['RA'],newcr['DEC'],newcr['MJD'],newcr['BAND'],newcr['MAG'],newcr['MAGERR'],newcr['FLAG']=cr['DatabaseID'],cr['RA'],cr['DEC'],cr['MJD'],cr['BAND'],cr['MAG'],cr['MAGERR'],cr['FLAG']
                cr=np.append(newcr,crmac)
        gorig=np.arange(len(cr))[(cr['MAG']>0)&(cr['MAG']<30)&(cr['MAGERR']<5)]
        cr=cr[(cr['MAG']>0)&(cr['MAG']<30)&(cr['MAGERR']<5)]
        gdes=np.where(cr['Survey']=='DES')[0]
        if ((len(gdes)>1)&(convertDESmags)):
            gg,gr=np.where(cr['BAND'][gdes]=='g')[0],np.where(cr['BAND'][gdes]=='r')[0]
            if (len(gg)>0)&(len(gr)>0):
                if len(gg)==1:
                    medg=cr['MAG'][gdes][gg]
                else:
                    medg=np.median(cr['MAG'][gdes][gg])
                if len(gr)==1:
                    medr=cr['MAG'][gdes][gr]
                else:
                    medr=np.median(cr['MAG'][gdes][gr])
                newg,dum1=DES2SDSS_gr(cr['MAG'][gdes][gg],medr)
                dum2,newr=DES2SDSS_gr(medg,cr['MAG'][gdes][gr])
                cr['MAG'][gdes[gg]],cr['MAG'][gdes[gr]]=newg,newr
            gi,gz=np.where(cr['BAND'][gdes]=='i')[0],np.where(cr['BAND'][gdes]=='z')[0]
            if (len(gi)>0)&(len(gz)>0):
                if len(gi)==1:
                    medi=cr['MAG'][gdes][gi]
                else:
                    medi=np.median(cr['MAG'][gdes][gi])
                if len(gz)==1:
                    medz=cr['MAG'][gdes][gz]
                else:
                    medz=np.median(cr['MAG'][gdes][gz])
                newi,dum1=DES2SDSS_iz(cr['MAG'][gdes][gi],medz)
                dum2,newz=DES2SDSS_iz(medi,cr['MAG'][gdes][gz])
                cr['MAG'][gdes[gi]],cr['MAG'][gdes[gz]]=newi,newz
        mjd,mag,magerr,bands,survey=cr['MJD'],cr['MAG'],cr['MAGERR'],cr['BAND'],cr['Survey']
        if ((calc_outliers)&(not(load_macleod))):
            gb=np.where(bands=='g')[0]
            if np.shape(outlier_window)!=():
                outliercolorarr=['red','cyan','green','blue','orange','magenta']
                outliermarkerarr=['x','+','d','o','*']
                outlier_arr={}
                for w,iw in zip(outlier_window,np.arange(len(outlier_window))):
                    outlier_arr[w]= {'g': np.zeros(len(mag),dtype='bool'), 'c': outliercolorarr[iw], 'marker': outliermarkerarr[iw]}
                    for ipt in np.arange(len(gb)):
                        gthresh=np.where(np.abs(mjd[gb]-mjd[gb[ipt]])<w)[0]
                        if len(gthresh)>2:
                            outlier_arr[w]['g'][gb[ipt]]= np.abs(np.median(mag[gb[gthresh]])-mag[gb[ipt]]) > outlier_thresh
                        elif len(gthresh)==2:
                            outlier_arr[w]['g'][gb[ipt]]= np.abs(np.median(mag[gb[gthresh]])-mag[gb[ipt]]) > outlier_thresh
            else:
                for ipt in np.arange(len(gb)):
                    gthresh=np.where(np.abs(mjd[gb]-mjd[gb[ipt]])<outlier_window)[0]
                    if len(gthresh)>2:
                        outlier_arr[gorig[gb[ipt]]]= np.abs(np.median(mag[gb[gthresh]])-mag[gb[ipt]]) > outlier_thresh
                    elif len(gthresh)==2:
                        outlier_arr[gorig[gb[ipt]]]= np.abs(np.median(mag[gb[gthresh]])-mag[gb[ipt]]) > outlier_thresh
                np.savetxt('%s/%s/outliers.tab'%(DBdir,DBID),outlier_arr)
                outlier_arr=outlier_arr[gorig]
        elif not(load_macleod):
            outlier_arr=None
        if ((np.shape(crmac)!=())&(not(load_macleod))):
            mjd,mag,magerr,bands=(mjd,crmac['MJD']),(mag,crmac['MAG']),(magerr,crmac['MAGERR']),(bands,crmac['BAND'])
        plot_lightcurve(DBID,mjd,mag,magerr,bands,survey,trueredshift,DBdir,psfpage,specfile=specfile,DESfname=DESfname,WavLL=WavLL,WavUL=WavUL,outlierflag=outlierflags[idb],zoominband=zoominband,outlierarr=outlier_arr,connectpoints=connectpoints)
    psfpage.close()

def plot_DBID(DBID,DBdir='/data2/rumbaugh/var_database/Y3A1',WavLL=3000,WavUL=10500,convertDESmags=False):
    hdu=py.open('%s/masterfile.fits'%DBdir)
    data=hdu[1].data
    crdescutout=np.loadtxt('%s/imagestamps/DESimagestamp_index.dat'%DBdir,dtype={'names':('DBID','ra','dec','tile','fname'),'formats':('|S32','f8','f8','|S12','|S25')})
    crdb=np.loadtxt('%s/database_index.dat'%DBdir,dtype={'names':('DBID','CID','thingid','sdr7id','MQrownum','SPrownum','SDSSNAME'),'formats':('|S64','i8','i8','|S24','i8','i8','|S64')})
    crpmf1=np.loadtxt('%s/MQ_SDSS_y3a1_tidplatemjdfiber.csv'%DBdir,delimiter=',',skiprows=1,dtype={'names':('tid','plate','mjd','fiber'),'formats':('i8','i8','i8','i8')})
    crpmf2=np.loadtxt('%s/DR7_BH_SDSS_tidplatemjdfiber.csv'%DBdir,delimiter=',',skiprows=1,dtype={'names':('tid','plate','mjd','fiber'),'formats':('i8','i8','i8','i8')})
    crpmf=np.concatenate((crpmf1,crpmf2),axis=0)

    coldict={'g': 'green','r': 'red', 'i': 'magenta', 'z': 'blue', 'Y': 'cyan'}
    SDSSbands=np.array(['u','g','r','i','z'])
    SDSS_colnames={b:'%s_SDSS'%b for b in SDSSbands}
    POSSbands=np.array(['g','r','i'])
    gdc=np.where(crdescutout['DBID']==DBID)[0]
    DESfname=None
    if len(gdc)>0:
        if crdescutout['fname'][gdc[0]]!='False':
            DESfname='%s.tif'%(crdescutout['fname'][gdc[0]])
    try:
        gmf=np.where(data['DatabaseID']==DBID)[0][0]
        trueredshift=data['Redshift'][gmf]
    except:
        trueredshift=0
    redshift=np.copy(trueredshift)
    if redshift<0: redshift=0
    gdb=np.where(crdb['DBID']==DBID)[0][0]
    tid=crdb['thingid'][gdb]
    plate,pmf_mjd,fiber=-1,-1,-1
    if tid!=0:
        gpmf=np.where(crpmf['tid']==tid)[0]
        if len(gpmf)>0:
            gpmf=gpmf[0]
            plate,pmf_mjd,fiber=crpmf['plate'][gpmf],crpmf['mjd'][gpmf],crpmf['fiber'][gpmf]
    if (plate>-1)&(pmf_mjd>-1)&(fiber>-1):
        specfile='%s/spec/spec-%04i-%05i-%04i.fits'%(DBdir,plate,pmf_mjd,fiber)
    else:
        specfile=None
    cr=np.loadtxt('%s/%s/LC.tab'%(DBdir,DBID),dtype={'names':('DatabaseID','Survey','SurveyCoaddID','SurveyObjectID','RA','DEC','MJD','TAG','BAND','MAGTYPE','MAG','MAGERR','FLAG'),'formats':('|S64','|S20','|S20','|S20','f8','f8','f8','|S20','|S12','|S12','f8','f8','i8')},skiprows=1)
    cr=cr[(cr['MAG']>0)&(cr['MAG']<30)&(cr['MAGERR']<5)]
    gdes=np.where(cr['Survey']=='DES')[0]
    if ((len(gdes)>1)&(convertDESmags)):
        gg,gr=np.where(cr['BAND'][gdes]=='g')[0],np.where(cr['BAND'][gdes]=='r')[0]
        if (len(gg)>0)&(len(gr)>0):
            if len(gg)==1:
                medg=cr['MAG'][gdes][gg]
            else:
                medg=np.median(cr['MAG'][gdes][gg])
            if len(gr)==1:
                medr=cr['MAG'][gdes][gr]
            else:
                medr=np.median(cr['MAG'][gdes][gr])
            newg,dum1=DES2SDSS_gr(cr['MAG'][gdes][gg],medr)
            dum2,newr=DES2SDSS_gr(medg,cr['MAG'][gdes][gr])
            cr['MAG'][gdes[gg]],cr['MAG'][gdes[gr]]=newg,newr
        gi,gz=np.where(cr['BAND'][gdes]=='i')[0],np.where(cr['BAND'][gdes]=='z')[0]
        if (len(gi)>0)&(len(gz)>0):
            if len(gi)==1:
                medi=cr['MAG'][gdes][gi]
            else:
                medi=np.median(cr['MAG'][gdes][gi])
            if len(gz)==1:
                medz=cr['MAG'][gdes][gz]
            else:
                medz=np.median(cr['MAG'][gdes][gz])
            newi,dum1=DES2SDSS_iz(cr['MAG'][gdes][gi],medz)
            dum2,newz=DES2SDSS_iz(medi,cr['MAG'][gdes][gz])
            cr['MAG'][gdes[gi]],cr['MAG'][gdes[gz]]=newi,newz
    mjd,mag,magerr,bands,survey=cr['MJD'],cr['MAG'],cr['MAGERR'],cr['BAND'],cr['Survey']
    plot_lightcurve(DBID,mjd,mag,magerr,bands,survey,trueredshift,DBdir,specfile=specfile,DESfname=DESfname,WavLL=WavLL,WavUL=WavUL)
    plt.show()

def plot_CLQ_candidates(magdrop,outputfile,DBdir='/data2/rumbaugh/var_database/Y3A1',WavLL=3000,WavUL=10500,convertDESmags=False,outputDBIDs=None):
    crdrop=np.loadtxt('%s/max_mag_drop.dat'%DBdir,dtype={'names':('DBID','maxdiff'),'formats':('|S128','f8')},skiprows=1)
    good_dbids=crdrop['DBID'][crdrop['maxdiff']>=magdrop]
    plot_DB_lightcurves(good_dbids,outputfile,DBdir=DBdir,WavLL=WavLL,WavUL=WavUL,convertDESmags=convertDESmags)
    if outputDBIDs!=None:np.savetxt(outputDBIDs,good_dbids,fmt='%s')
