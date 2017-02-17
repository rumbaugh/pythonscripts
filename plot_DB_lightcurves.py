import numpy as np
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
        if fluxerrs!=None:ax.errorbar(cens,fluxes,yerr=fluxerrs,color=curcol,fmt='ro',lw=2,capsize=3,mew=1,zorder=3,label=None)
        ax.scatter(cens,fluxes,color=curcol,s=36,zorder=4)
    else:
        ax.scatter(cens,fluxes,color=curcol,s=36,label=label,zorder=4)
        if fluxerrs!=None:ax.errorbar(cens,fluxes,yerr=fluxerrs,color=curcol,fmt='ro',lw=2,capsize=3,mew=1,zorder=3,label=None)

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

def plot_band(ax,mjd,mag,magerr,cbands,band,connectpoints=True,nolabels=False):
    gband=np.where(cbands==band)[0]
    magplot=mag[gband]
    magploterr=magerr[gband]
    g100=np.where(magplot<100)[0]
    try:
        curcol=coldict[band]
    except KeyError:
        print '%s is not a valid band'%band
        return
    if connectpoints:
        gsort=np.argsort(mjd[gband][g100])
        ax.plot(mjd[gband][g100][gsort],magplot[g100][gsort],color=curcol,lw=2)
    ax.errorbar(mjd[gband][g100],magplot[g100],yerr=magploterr[g100],color=curcol,fmt='ro',lw=2,capsize=3,mew=1)
    if nolabels:
        ax.scatter(mjd[gband][g100],magplot[g100],color=curcol)
    else:
        ax.scatter(mjd[gband][g100],magplot[g100],color=curcol,label=band)
    #return


def plot_lightcurve(dbid,mjd,mag,magerr,bands,survey,trueredshift,DBdir,fname=None,DESfname=None,connectpoints=True,specfile=None,WavLL=3000,WavUL=10500):
    bcens={'u': 3876.63943790537, 'g': 4841.83358196563, 'r': 6438/534828217, 'i': 7820.99282740933, 'z': 9172.34266385718, 'Y': 9877.80238651117}
    VBfile='%s/VanderBerk_datafile1.txt'%DBdir
    crv=np.loadtxt(VBfile,skiprows=23)
    redshift=np.copy(trueredshift)
    if redshift<0:redshift=0
    gdes,gsdss,gposs=np.where(survey=='DES')[0],np.where(survey=='SDSS')[0],np.where(survey=='POSS')[0]
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
            mag[gposs][bands[gposs]=='g'],mag[gposs][bands[gposs]=='r']=np.zeros(0),np.zeros(0)
        if (len(POSSmagdict['i'])>0)&(len(POSSmagdict['r'])>0):   
            mag[gposs][bands[gposs]=='i']=mag[gposs][bands[gposs]=='i']+0.27*(POSSmagdict['r']-POSSmagdict['i'])+0.32
            magerr[gposs][bands[gposs]=='i']=np.sqrt(magerr[gposs][bands[gposs]=='i']**2+0.27**2*(magerr[gposs][bands[gposs]=='r']**2+magerr[gposs][bands[gposs]=='i']**2))
        else:
            mag[gposs][bands[gposs]=='i']=np.zeros(0)
    bestdiff={b: {'diff': 0, 'ihi': 0, 'ilo': 0} for b in ['g','r','i','z']}
    for b in ['g','r','i','z']:
        gb=np.where(bands==b)[0]
        if len(gb)>1:
            combis=np.array(list(it.combinations(np.arange(len(mag[gb])),2)))
            i1,i2=combis[:,0],combis[:,1]
            sigma=np.abs((mag[gb][i1]-mag[gb][i2])/np.sqrt(magerr[gb][i1]**2+magerr[gb][i2]**2))
            totdiffstmp=np.abs(mag[gb][i1]-mag[gb][i2])
            ggooddiff=np.where((sigma>=3)&(mag[gb][i1]>1)&(mag[gb][i1]<30)&(mag[gb][i2]>1)&(mag[gb][i2]<30))[0]
            if len(ggooddiff)>0:
                bestdiff[b]['diff']=np.max(totdiffstmp[ggooddiff])
                gsort=np.argsort(totdiffstmp[ggooddiff])
                gsortis=np.argsort([mag[gb][i1][ggooddiff][gsort[-1]],mag[gb][i2][ggooddiff][gsort[-1]]])
                imax,imin=[i1,i2][gsortis[0]][ggooddiff][gsort[-1]],[i1,i2][gsortis[1]][ggooddiff][gsort[-1]]
                bestdiff[b]['ihi'],bestdiff[b]['ilo']=imax,imin
    fig=plt.figure(1)
    fig.clf()
    ax3=plt.subplot2grid((2,10),(1,0),colspan=6)
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
        maxfluxes[ib],maxfluxerrs[ib]=calc_flux(mjd,mag,magerr,bands,b,mjd[gbbest][imax])
        minfluxes[ib],minfluxerrs[ib]=calc_flux(mjd,mag,magerr,bands,b,mjd[gbbest][imin])
    plt.rc('axes',linewidth=2)
    plt.fontsize = 14
    plt.tick_params(which='major',length=8,width=2,labelsize=14)
    plt.tick_params(which='minor',length=4,width=1.5,labelsize=14)
    plt.locator_params(nbins=4)
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
    survmax,survmin=survey[gbbest][imax],survey[gbbest][imin]
    if ((not('DES' in [survmax,survmin]))&('DES' in survey[bands==bbest])):
        sortmjdcens=np.argsort([mjd[gbbest][imax],mjd[gbbest][imin]])
        iDESex=np.argsort(mag[(bands==bbest)&(survey=='DES')])[-sortmjdcens[0]]
        DESexfluxes,DESexfluxerrs=np.zeros(4),np.zeros(4)
        for ib,b in zip(np.arange(4),['g','r','i','z']):
            gbt=np.where(bands==b)[0]
            DESexfluxes[ib],DESexfluxerrs[ib]=calc_flux(mjd,mag,magerr,bands,b,mjd[(bands==bbest)&(survey=='DES')][iDESex])
        plot_flux(ax3,DESexfluxes,DESexfluxerrs,label='DES %s'%['Max','Min'][-sortmjdcens[0]],curcol='cyan')
        if np.max(DESexfluxes)>maxplot:maxplot=np.max(DESexfluxes)
    #ax3.legend(loc='lower right')
    plt.xlim(WavLL,WavUL)
    plt.ylim(-0.05,vmax*1.05)
    if trueredshift<=0: 
        plt.ylim(-0.05,1.15*maxplot)
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
    plt.axvline(mjd[gbbest][imax],ls='dashed',lw=1,color='r')
    plt.axvline(mjd[gbbest][imin],ls='dashed',lw=1,color='b')
    if ((not('DES' in [survmax,survmin]))&('DES' in survey[bands==bbest])):plt.axvline(mjd[gbbest][iDESex],ls='dashed',lw=1,color='cyan')
    if ylim[1]>30:
        ylim=(ylim[0],np.max(mag)+0.1)
    if ylim[1]>30: ylim=(ylim[0],30)
    if ylim[0]<15:
        ylim=(np.min(mag)-0.1,ylim[1])
    if ylim[0]<15: ylim=(15,ylim[1])
    plt.ylim(ylim[1],ylim[0])
    ax1.legend()
    ax1.set_xlabel('MJD')
    ax1.set_ylabel('Mag_PSF')
    if redshift>0:
        ax1.set_title('%i, z=%.4f'%(dbid,trueredshift))
    else:
        ax1.set_title(dbid)
    if not(DESfname in [None,'False']):
        ax4=plt.subplot2grid((2,10),(1,6),colspan=4,xticks=[],yticks=[])
        img4=mpimg.imread('%s/imagestamps/%s'%(DBdir/DESfname))
        ax4.imshow(img4)
    if len(gsdss)>0:
        SDSSfname='%s/imagestamps/%s_SDSScutout.jpeg'%(DBdir,dbid)
        try:
            img3=mpimg.imread(SDSSfname)
            ax3=plt.subplot2grid((2,10),(0,6),colspan=4,xticks=[],yticks=[])
            ax3.imshow(img3)
        except:
            pass
    plt.savefig(psfpdf,format='pdf')
    return

def DES2SDSS_gr(g,r):
    return (133625*g-9375*r-218)/124250.,(69.*g+925*r)/994.+516./62125.

def DES2SDSS_iz(i,z):
    return (-89*np.sqrt(-96000*i+96000*z+181561)+8000*z+37827)/8000.,(-17*np.sqrt(-96000*i+96000*z+181561)+24000*z+6731)/24000.

def plot_DB_lightcurves(DBIDs,outputfile,DBdir='/data2/rumbaugh/var_database/Y3A1',WavLL=3000,WavUL=10500,convertDESmags=False):
    hdu=py.open('%s/masterfile.fits'%DBdir)
    data=hdu[1].data
    psfpdf=bpdf.PdfPages(outputfile)

    crdescutout=np.loadtxt('%s/imagestamps/DESimagestamp_index.dat'%DBdir,dtype={'names':('DBID','ra','dec','tile','fname'),'formats':('|S32','f8','f8','|S12','|S25')})

    crdb=np.loadtxt('%s/database_index.dat'%DBdir,dtype={'names':('DBID','CID','thingid','sdr7id','MQrownum','SP_rownum','SDSSNAME'),'formats':('|S64','i8','i8','|S24','i8','i8','|S64')})
    crpmf1=np.loadtxt('%s/MQ_SDSS_y3a1_tidplatemjdfiber.csv'%DBdir,delimiter=',',skiprows=1,dtype={'names':('tid','plate','mjd','fiber'),'formats':('i8','i8','i8','i8')})
    crpmf2=np.loadtxt('%s/DR7_BH_SDSS_tidplatemjdfiber.csv'%DBdir,delimiter=',',skiprows=1,dtype={'names':('tid','plate','mjd','fiber'),'formats':('i8','i8','i8','i8')})
    crpmf=np.concatenate((crpmf1,crpmf2),axis=0)

    coldict={'g': 'green','r': 'red', 'i': 'magenta', 'z': 'blue', 'Y': 'cyan'}
    SDSSbands=np.array(['u','g','r','i','z'])
    SDSS_colnames={b:'%s_SDSS'%b for b in SDSSbands}
    POSSbands=np.array(['g','r','i'])

    for DBID in DBIDs:
        gdc=np.where(crdescutout['DBID']==DBID)[0]
        if len(gdc)>0:
            if crdescutout['fname'][gdc[0]]!='False':
                DESfname='%s.tif'%(crdescutout['fname'][gdc[0]])
        gmf=np.where(data['DatabaseID']==DBID)[0][0]
        trueredshift=data['Redshift'][gmf]
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
    psfpdf.close()

def plot_CLQ_candidates(magdrop,outputfile,DBdir='/data2/rumbaugh/var_database/Y3A1',WavLL=3000,WavUL=10500,convertDESmags=False):
    crdrop=np.loadtxt('%s/max_mag_drop.dat'%DBdir,dtype={'names':('DBID','maxdiff'),'formats':('|S128','f8')},skiprows=1)
    good_dbids=crdrop['DBID'][crdrop['maxdiff']>=magdrop]
    plot_DB_lightcurves(good_dbids,outputfile,DBdir=DBdir,WavLL=WavLL,WavUL=WavUL,convertDESmags=convertDESmags)
