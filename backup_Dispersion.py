#these functions calculate the dispersions as defined in Pelt et al. (1996,1998)
#and used in Fassnacht et al. 1999
import time
import numpy as np
import sys

def D_2_old(A,B,A_t,B_t,A_err,B_err,mu,tau):
    G_ref = np.append(np.zeros(len(A),dtype='int'),np.ones(len(B),dtype='int'))
    comb_flux = np.append(A,B*mu)
    comb_err = np.append(A_err,mu*B_err)
    comb_t = np.append(A_t,B_t+tau)
    comb_t_argsort = np.argsort(comb_t)
    sum_W_k_G_k = 0.0
    D_2_top = 0.0
    for iwg in range(0,len(comb_t)-1):
        if G_ref[comb_t_argsort[iwg]] + G_ref[comb_t_argsort[iwg+1]] == 1:
            errt = 1.0/(comb_err[comb_t_argsort[iwg]]**2+comb_err[comb_t_argsort[iwg+1]]**2)
            sum_W_k_G_k += errt
            D_2_top += errt*(comb_flux[comb_t_argsort[iwg+1]]-comb_flux[comb_t_argsort[iwg]])**2
    return D_2_top/sum_W_k_G_k

def D_2(A,B,A_t,B_t,A_err,B_err,mu,tau):
    G_ref = np.append(np.zeros(len(A),dtype='int'),np.ones(len(B),dtype='int'))
    comb_flux = np.append(A,B*mu)
    comb_err = np.append(A_err,B_err)
    comb_t = np.append(A_t,B_t+tau)
    comb_t_argsort = np.argsort(comb_t)
    sum_W_k_G_k = 0.0
    D_2_top = 0.0
    G_ref_check = G_ref[comb_t_argsort] - np.roll(G_ref[comb_t_argsort],-1)
    gref = np.where(G_ref_check[:len(G_ref_check)-1] != 0)
    gref = gref[0]
    for iwg in range(0,len(gref)):
        errt = 1.0/(comb_err[comb_t_argsort[gref[iwg]]]**2+comb_err[comb_t_argsort[gref[iwg]+1]]**2)
        sum_W_k_G_k += errt
        D_2_top += errt*(comb_flux[comb_t_argsort[gref[iwg]]]-comb_flux[comb_t_argsort[gref[iwg]+1]])**2
    return D_2_top/sum_W_k_G_k

#The difference between D_2 and D_2b is that in D_2 the B flux is multiplied
#by mu while in D_2b the A flux is divided by mu. There is a similar
#difference between D_4_2 and D_4_2b
def D_2b(A,B,A_t,B_t,A_err,B_err,mu,tau):
    G_ref = np.append(np.zeros(len(A),dtype='int'),np.ones(len(B),dtype='int'))
    comb_flux = np.append(A/mu,B)
    comb_err = np.append(A_err,B_err)
    comb_t = np.append(A_t,B_t+tau)
    comb_t_argsort = np.argsort(comb_t)
    sum_W_k_G_k = 0.0
    D_2_top = 0.0
    G_ref_check = G_ref[comb_t_argsort] - np.roll(G_ref[comb_t_argsort],-1)
    gref = np.where(G_ref_check[:len(G_ref_check)-1] != 0)
    gref = gref[0]
    for iwg in range(0,len(gref)):
        errt = 1.0/(comb_err[comb_t_argsort[gref[iwg]]]**2+comb_err[comb_t_argsort[gref[iwg]+1]]**2)
        sum_W_k_G_k += errt
        D_2_top += errt*(comb_flux[comb_t_argsort[gref[iwg]]]-comb_flux[comb_t_argsort[gref[iwg]+1]])**2
    return D_2_top/sum_W_k_G_k

def D_4_2_old(A,B,A_t,B_t,A_err,B_err,mu,tau,delta):
    G_ref = np.append(np.zeros(len(A),dtype='int'),np.ones(len(B),dtype='int'))
    comb_flux = np.append(A,B*mu)
    comb_err = np.append(A_err,B_err)
    comb_t = np.append(A_t,B_t+tau)
    comb_t_argsort = np.argsort(comb_t)
    sum_bottom,sum_top = 0.0,0.0
    for iwg in range(0,len(comb_t)-1):
        for iwg2 in range(iwg+1,len(comb_t)):
            if ((G_ref[comb_t_argsort[iwg]] + G_ref[comb_t_argsort[iwg2]] == 1) & (np.fabs(comb_t[comb_t_argsort[iwg]] - comb_t[comb_t_argsort[iwg2]]) <= delta)):
                errt = 1.0/(comb_err[comb_t_argsort[iwg]]**2+comb_err[comb_t_argsort[iwg2]]**2)
                sum_bottom += (1-np.fabs(comb_t[comb_t_argsort[iwg]] - comb_t[comb_t_argsort[iwg+1]])/delta)*errt
                sum_top += (1-np.fabs(comb_t[comb_t_argsort[iwg]] - comb_t[comb_t_argsort[iwg2]])/delta)*(comb_flux[comb_t_argsort[iwg2]]-comb_flux[comb_t_argsort[iwg]])**2*errt
    if ((sum_bottom == 0) | (sum_top == 0)):
        return 999999999.
    else:
        if sum_top*1.0/sum_bottom == 0.0:
            return 999999.
        else:
            return sum_top*1.0/sum_bottom

def D_4_2(A,B,A_t,B_t,A_err,B_err,mu,tau,delta):
    G_ref = np.append(np.zeros(len(A),dtype='int'),np.ones(len(B),dtype='int'))
    comb_flux = np.append(A,B*mu)
    comb_err = np.append(A_err,B_err)
    comb_t = np.append(A_t,B_t+tau)
    comb_t_argsort = np.argsort(comb_t)
    sum_bottom,sum_top = 0.0,0.0
    for iwg in range(0,len(comb_t)-1):
        giwg = np.where((comb_t[comb_t_argsort[iwg]] - comb_t[comb_t_argsort[iwg+1:]] <= delta) & (comb_t[comb_t_argsort[iwg]] - comb_t[comb_t_argsort[iwg+1:]] >= -1*delta))
        giwg = giwg[0]
        for iwg2 in range(0,len(giwg)):
            if ((G_ref[comb_t_argsort[iwg]] + G_ref[comb_t_argsort[iwg+1+giwg[iwg2]]] == 1)):
                errt = 1.0/(comb_err[comb_t_argsort[iwg]]**2+comb_err[comb_t_argsort[iwg+1+giwg[iwg2]]]**2)
                fabst = (1-np.fabs(comb_t[comb_t_argsort[iwg]] - comb_t[comb_t_argsort[iwg+1+giwg[iwg2]]])/delta)
                sum_bottom += fabst*errt
                sum_top += fabst*(comb_flux[comb_t_argsort[iwg+1+giwg[iwg2]]]-comb_flux[comb_t_argsort[iwg]])**2*errt
    if ((sum_bottom == 0) | (sum_top == 0)):
        return 999999999.
    else:
        if sum_top*1.0/sum_bottom == 0.0:
            return 999999.
        else:
            return sum_top*1.0/sum_bottom

def D_4_2b(A,B,A_t,B_t,A_err,B_err,mu,tau,delta):
    G_ref = np.append(np.zeros(len(A),dtype='int'),np.ones(len(B),dtype='int'))
    comb_flux = np.append(A/mu,B)
    comb_err = np.append(A_err,B_err)
    comb_t = np.append(A_t,B_t+tau)
    comb_t_argsort = np.argsort(comb_t)
    sum_bottom,sum_top = 0.0,0.0
    for iwg in range(0,len(comb_t)-1):
        giwg = np.where((G_ref[comb_t_argsort[iwg]] + G_ref[comb_t_argsort[iwg+1:]] == 1) & (np.fabs(comb_t[comb_t_argsort[iwg]] - comb_t[comb_t_argsort[iwg+1:]]) <= delta))[0]
        #for iwg2 in range(0,len(giwg)):
        #    errt = 1.0/(comb_err[comb_t_argsort[iwg]]**2+comb_err[comb_t_argsort[iwg+1+giwg[iwg2]]]**2)
        #    fabst = (1-np.fabs(comb_t[comb_t_argsort[iwg]] - comb_t[comb_t_argsort[iwg+1+giwg[iwg2]]])/delta)
        #    sum_bottom += fabst*errt
        #    sum_top += fabst*(comb_flux[comb_t_argsort[iwg+1+giwg[iwg2]]]-comb_flux[comb_t_argsort[iwg]])**2*errt
        errt = 1.0/(comb_err[comb_t_argsort[iwg]]**2+comb_err[comb_t_argsort[iwg+1+giwg]]**2)
        fabst = (1-np.fabs(comb_t[comb_t_argsort[iwg]] - comb_t[comb_t_argsort[iwg+1+giwg]])/delta)
        sum_bottom += (fabst*errt).sum()
        sum_top += (fabst*(comb_flux[comb_t_argsort[iwg+1+giwg]]-comb_flux[comb_t_argsort[iwg]])**2*errt).sum()
    if ((sum_bottom == 0) | (sum_top == 0)):
        return 999999999.
    else:
        if sum_top*1.0/sum_bottom == 0.0:
            return 999999.
        else:
            return sum_top*1.0/sum_bottom

def calc_disp_delay(A,B,A_t,B_t,A_err,B_err,maxtime,timestep,minmu,maxmu,mustep,disp_type,delta=3.5,output=1,print_times=False,disparray=False,dispmatrix=False,mintime=None,inner50=True,simplemuerr=False,ALAG=31.5,use_overlap_mean=False,outfile=None,verbose=True):
    #calculates the dispersion on a grid of tau (time delay) and mu values
    #output = 1 returns just the minimum delay and mu
    #output = 2 also returns the minimum dispersion
    #output = 3 also returns dispersion matrix
    #output = 4 also returns mu0
    #if disparray=True, an array
    #if outfile is set, dispersion values for each mu,tau pair are written
    #to a text file
    maxtimestep=maxtime
    if outfile != None:
        FILE = open(outfile,'w')
        FILE.write('#tau mu disp\n')
    B_err_t = B_err.copy()
    disparrout,disparrmu,disparrtime = np.zeros(0),np.zeros(0),np.zeros(0)
    start_t = time.time()
    mindisp,mintau,minmu_out = -99,0.,0.
    if use_overlap_mean:
        galag,gblag = np.where((A_t > np.min(B_t)+ALAG) & (A_t < np.max(B_t) + ALAG)),np.where((B_t + ALAG > np.min(A_t)) & (B_t + ALAG < np.max(A_t)))
        galag,gblag = galag[0],gblag[0]
        mu0 = np.mean(A[galag])/np.mean(B[gblag])
        mu0err = mu0*np.sqrt(np.sum(A_err[galag]*A_err[galag])/np.mean(A[galag])/np.mean(A[galag])/((len(galag))**2)+np.sum(B_err[gblag]*B_err[gblag])/np.mean(B[gblag])/np.mean(B[gblag])/((len(gblag))**2))
    elif inner50:
        mu0 = np.mean(A[np.floor(len(A)/4)-1:np.ceil(3*len(A)/4)])/np.mean(B[np.floor(len(A)/4)-1:np.ceil(3*len(A)/4)])
        mu0err = mu0*np.sqrt(np.sum((A_err[np.floor(len(A)/4)-1:np.ceil(3*len(A)/4)])**2)/((np.mean(A[np.floor(len(A)/4)-1:np.ceil(3*len(A)/4)]))**2)+np.sum((B_err[np.floor(len(A)/4)-1:np.ceil(3*len(A)/4)])**2)/((np.mean(B[np.floor(len(A)/4)-1:np.ceil(3*len(A)/4)]))**2))/(np.ceil(3*len(A)/4)-np.floor(len(A)/4))
    else:
        mu0 = np.mean(A)*1.0/np.mean(B)
        mu0err = mu0*np.sqrt(np.sum(A_err*A_err)/np.mean(A)/np.mean(A)+np.sum(B_err*B_err)/np.mean(B)/np.mean(B))/len(A)
    if ((disp_type != 'D_2b') & (disp_type != 'D_2') & (disp_type != 'D_4_2b') & (disp_type != 'D_4_2')): sys.exit("disp_type must be either 'D_2' pr 'D_4_2' for calc_disp_delay")
    tau = -1.0*maxtimestep
    #if ((mintime != None) & (mintime > tau) & (mintime < -1.0*tau)): tau = mintime
    if ((mintime != None)): 
        tau = mintime
    basetime = tau
    st2 = time.time()
    while tau <= maxtimestep:
        #Figuring out where the overlap between A and shifted B curves are
        #Then, finds the flux ratio between them
        #if tau > 0:
        #    gmu0A,gmu0B = np.where(A_t >= tau+A_t.min()),np.where(B_t <= B_t.max()-tau)
        #    gmu0A,gmu0B = gmu0A[0],gmu0B[0]
        #    Acmu0,Bcmu0 = A[gmu0A],B[gmu0B]
        #    Acmu0err,Bcmu0err = A_err[gmu0A],B_err[gmu0B]
        #elif tau < 0:
        #    gmu0B,gmu0A = np.where(B_t >= tau+B_t.min()),np.where(A_t <= A_t.max()-tau)
        #    gmu0A,gmu0B = gmu0A[0],gmu0B[0]
        #    Acmu0,Bcmu0 = A[gmu0A],B[gmu0B]
        #    Acmu0err,Bcmu0err = A_err[gmu0A],B_err[gmu0B]
        #else:
        #    Acmu0,Bcmu0 = A,B
        #    Acmu0err,Bcmu0err = A_err,B_err
        #mu0 = np.mean(Acmu0)*1.0/np.mean(Bcmu0)
        #mu0err = mu0*np.sqrt(np.sum(Acmu0err*Acmu0err)/np.mean(Acmu0)/np.mean(Acmu0)+np.sum(Bcmu0err*Bcmu0err)/np.mean(Bcmu0)/np.mean(Bcmu0))/len(Acmu0)
        mu = minmu*mu0
        muerr = minmu*mu0err
        mindisp2,minmu2 = -99,0.
        disparrtmp = np.zeros(0)
        while mu <= maxmu*mu0:
            if simplemuerr:
                if ((disp_type == 'D_2b') | (disp_type == 'D_4_2b')):
                    A_err_t = A_err/mu
                    B_err_t = B_err
                else:
                    A_err_t = A_err
                    B_err_t = B_err*mu
            else:
                for ibet in range(0,len(B_err)): B_err_t[ibet] = np.sqrt(muerr**2/mu/mu+(B_err[ibet])**2/B[ibet]/B[ibet])*B[ibet]*mu
            if disp_type == 'D_2b': 
                if not simplemuerr: sys.exit('D_2b option only works with simplemuerr')
                D_tmp = D_2b(A,B,A_t,B_t,A_err_t,B_err,mu,tau)
            if disp_type == 'D_4_2b': 
                if not simplemuerr: sys.exit('D_4_2b option only works with simplemuerr')
                D_tmp = D_4_2b(A,B,A_t,B_t,A_err_t,B_err,mu,tau,delta)
            if disp_type == 'D_2': D_tmp = D_2(A,B,A_t,B_t,A_err,B_err_t,mu,tau)
            if disp_type == 'D_4_2': D_tmp = D_4_2(A,B,A_t,B_t,A_err,B_err_t,mu,tau,delta)
            if ((D_tmp < mindisp) | (mindisp == -99)): mindisp,mintau,minmu_out = D_tmp,tau,mu
            if ((D_tmp < mindisp2) | (mindisp2 == -99)): mindisp2,minmu2 = D_tmp,mu
            disparrtmp = np.append(disparrtmp,D_tmp)
            if outfile != None:
                #print D_tmp
                FILE.write('%f %f %E\n'%(tau,mu,D_tmp))
            mu += mustep*mu0
            muerr += mustep*mu0err
        if (((tau == -1*maxtimestep) & (mintime == None)) | (tau == mintime)):
            dispmatrixout = np.array([disparrtmp])
        else:
            dispmatrixout = np.append(dispmatrixout,np.array([disparrtmp]),axis=0)
        disparrout,disparrmu,disparrtime = np.append(disparrout,mindisp2),np.append(disparrmu,minmu2),np.append(disparrtime,tau)
        if verbose:
            if tau < (basetime + 0.5*timestep): 
                t_mu_1 = time.time()
                print "Initial ETA: %i seconds"%(int((maxtimestep-basetime)/timestep*(t_mu_1-st2)-(t_mu_1-st2)))
        tau += timestep
        if verbose:
            if ((tau >= (basetime + 0.25*(maxtimestep-basetime))) & (tau < (basetime + 0.25*(maxtimestep-basetime) + timestep))): 
                t_25 = time.time()
                print "25%% Done - ETA: %i seconds"%(int(3*(t_25-start_t)))
            if ((tau >= (basetime + 0.50*(maxtimestep-basetime))) & (tau < (basetime + 0.50*(maxtimestep-basetime) + timestep))): 
                t_50 = time.time()
                print "50%% Done - ETA: %i seconds"%(int((t_50-start_t)))
            if ((tau >= (basetime + 0.75*(maxtimestep-basetime))) & (tau < (basetime + 0.75*(maxtimestep-basetime) + timestep))): 
                t_75 = time.time()
                print "75%% Done - ETA: %i seconds"%(int((t_75-start_t)/3.0))
    if verbose: 
        if time.time()-start_t<2: print 'Total time elapsed: %f seconds'%(time.time()-start_t)
        else: print 'Total time elapsed: %i seconds'%(time.time()-start_t)
    if outfile != None: FILE.close()
    if dispmatrix:
        if output == 3:
            return dispmatrixout,mintau,minmu_out,mindisp
        elif output == 4:
            return dispmatrixout,mintau,minmu_out,mindisp,mu0
        else:
            return dispmatrixout
    elif disparray:
        return disparrout,disparrmu,disparrtime
    elif output == 1:
        return mintau,minmu_out
    else:
        return mintau,minmu_out,mindisp

    

def calc_disp_delay_test(A,B,A_t,B_t,A_err,B_err,maxtime,timestep,minmu,maxmu,mustep,disp_type,delta=3.5,output=1,print_times=False,disparray=False,dispmatrix=False,mintime=None,inner50=True,simplemuerr=False):
    #calculates the dispersion on a grid of tau (time delay) and mu values
    #output = 1 returns just the minimum delay and mu
    #output = 2 also returns the minimum dispersion
    #if disparray=True, an array
    B_err_t = B_err.copy()
    n_timepts,n_mupts = 2*maxtime/timestep+1,(maxmu-minmu)/mustep+1
    disparrout,disparrmu,disparrtime = np.zeros(n_timepts),np.zeros(n_timepts),np.zeros(n_timepts)
    dispmatrixout = np.zeros((n_timepts,n_mupts))
    start_t = time.time()
    mindisp,mintau,minmu_out = -99,0.,0.
    if inner50:
        mu0 = np.mean(A[np.floor(len(A)/4)-1:np.ceil(3*len(A)/4)])/np.mean(B[np.floor(len(A)/4)-1:np.ceil(3*len(A)/4)])
        mu0err = mu0*np.sqrt(np.sum((A_err[np.floor(len(A)/4)-1:np.ceil(3*len(A)/4)])**2)/((np.mean(A[np.floor(len(A)/4)-1:np.ceil(3*len(A)/4)]))**2)+np.sum((B_err[np.floor(len(A)/4)-1:np.ceil(3*len(A)/4)])**2)/((np.mean(B[np.floor(len(A)/4)-1:np.ceil(3*len(A)/4)]))**2))/(np.ceil(3*len(A)/4)-np.floor(len(A)/4))
    else:
        mu0 = np.mean(A)*1.0/np.mean(B)
        mu0err = mu0*np.sqrt(np.sum(A_err*A_err)/np.mean(A)/np.mean(A)+np.sum(B_err*B_err)/np.mean(B)/np.mean(B))/len(A)
    if ((disp_type != 'D_2b') & (disp_type != 'D_2') & (disp_type != 'D_4_2')): sys.exit("disp_type must be either 'D_2' pr 'D_4_2' for calc_disp_delay")
    for itau in range(0,n_timepts):
        tau = -1*maxtime+itau*timestep
        muerr = minmu*mu0err
        mindisp2,minmu2 = -99,0.
        for imu in range(0,n_mupts):
            mu = minmu*mu0+imu*mustep*mu0
            if simplemuerr:
                if (disp_type == 'D_2b'):
                    A_err_t = A_err/mu
                    B_err_t = B_err
                else:
                    A_err_t = A_err
                    B_err_t = B_err*mu
            else:
                for ibet in range(0,len(B_err)): B_err_t[ibet] = np.sqrt(muerr**2/mu/mu+(B_err[ibet])**2/B[ibet]/B[ibet])*B[ibet]*mu
            if disp_type == 'D_2b': 
                if not simplemuerr: sys.exit('D_2b option only works with simplemuerr')
                D_tmp = D_2b(A,B,A_t,B_t,A_err_t,B_err,mu,tau)
            if disp_type == 'D_2': D_tmp = D_2(A,B,A_t,B_t,A_err,B_err_t,mu,tau)
            if disp_type == 'D_4_2': D_tmp = D_4_2(A,B,A_t,B_t,A_err,B_err_t,mu,tau,delta)
            if ((D_tmp < mindisp) | (mindisp == -99)): mindisp,mintau,minmu_out = D_tmp,tau,mu
            if ((D_tmp < mindisp2) | (mindisp2 == -99)): mindisp2,minmu2 = D_tmp,mu
            if dispmatrix: dispmatrixout[int((tau+maxtime)/timestep)][int((mu-minmu*mu0)/mustep/mu0)] = D_tmp
            muerr += mustep*mu0err
        if disparray: 
            itau = int((tau+maxtime)/timestep)
            disparrout[itau],disparrmu[itau],disparrtime[itau] = mindisp2,minmu2,tau
        if tau < (-1.0*maxtimestep + 0.5*timestep): 
            t_mu_1 = time.time()
            print "Initial ETA: %i seconds"%(int((maxtimestep-basetime)/timestep/(t_mu_1-start_t)-(t_mu_1-start_t)))
        if ((tau >= (-1.0*maxtimestep + 0.25*(maxtimestep-basetime))) & (tau < (-1.0*maxtimestep + 0.25*(maxtimestep-basetime) + timestep))): 
            t_25 = time.time()
            print "25%% Done - ETA: %i seconds"%(int(3*(t_25-start_t)))
        if ((tau >= (-1.0*maxtimestep + 0.50*(maxtimestep-basetime))) & (tau < (-1.0*maxtimestep + 0.50*(maxtimestep-basetime) + timestep))): 
            t_50 = time.time()
            print "50%% Done - ETA: %i seconds"%(int((t_50-start_t)))
        if ((tau >= (-1.0*maxtimestep + 0.75*(maxtimestep-basetime))) & (tau < (-1.0*maxtimestep + 0.75*(maxtimestep-basetime) + timestep))): 
            t_75 = time.time()
            print "75%% Done - ETA: %i seconds"%(int((t_75-start_t)/3.0))
    if dispmatrix:
        return dispmatrixout
    elif disparray:
        return disparrout,disparrmu,disparrtime
    elif output == 1:
        return mintau,minmu_out
    else:
        return mintau,minmu_out,mindisp

    
