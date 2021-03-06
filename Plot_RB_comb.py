import numpy as np
import spec_simple as ss
import matplotlib.pyplot as plt

def plot_RB_comb(mask,slit,side,line=1,redux_dir='/mnt/data2/rumbaugh/LRIS/2011_01/reduced',smooth=7.,output=False):
    if line == None:
        pass
    else:
        tdict,tdict2 = {x: ' Line %i'%(x) for x in range(1,10)},{x: '_line%i'%(x) for x in range(1,10)}
        tdict[1],tdict2[1] = '',''
        bfile,rfile,title = '%s/%s/spec_output/outspec.%s_%s_blue_%s_coadd_bgsub.dat'%(redux_dir,mask,mask,slit,side),'%s/%s/spec_output/outspec.%s_%s_red_%s_coadd_bgsub.dat'%(redux_dir,mask,mask,slit,side),'%s %s %s%s'%(mask,slit,side,tdict[line])
        wb,fb,vb = ss.read_spectrum(bfile,line=line)
        wr,fr,vr = ss.read_spectrum(rfile,line=line)
        m_overlap = 0.5*(np.max(wb)+np.min(wr))
        gb,gr = np.where(wb <= m_overlap)[0],np.where(wr > m_overlap)[0]
        w,f,v = np.append(wb[gb],wr[gr]),np.append(fb[gb],fr[gr]),np.append(vb[gb],vr[gr])
        if smooth == None:
            ss.plot_spectrum_array(w,f,title=title,clear=True)
        else:
            if output:
                w,f,v = ss.smooth_boxcar(None,smooth,varwt=True,title=title,line=line,output=output,clear=True,w_in=w,f_in=f,v_in=v)
            else:
                ss.smooth_boxcar(None,smooth,varwt=True,title=title,line=line,output=False,clear=True,w_in=w,f_in=f,v_in=v)
        plt.xlim(np.min(wb)-0.025*(np.max(wr)-np.min(wb)),np.max(wr)+0.025*(np.max(wr)-np.min(wb)))
    if output:
        return w,f,v
    else:
        return
