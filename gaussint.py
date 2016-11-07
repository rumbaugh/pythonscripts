import math as m
import numpy as np
import sys

#guassian cdf table - prob of being within N sigma
# 1  -  0.682689492137
# 2  -  0.954499736104
# 3  -  0.997300203937
# 4  -  0.999936657516 
# 5  -  0.999999426697
# 6  -  0.999999998027 	
gausscdftable = np.zeros(7)
#this table is probability of being within N sigma of mean (0 is 0)
gausscdftable[1] = 0.682689492137
gausscdftable[2] = 0.954499736104
gausscdftable[3] = 0.997300203937
gausscdftable[4] = 0.999936657516 
gausscdftable[5] = 0.999999426697
gausscdftable[6] = 0.999999998027 	


#gaussian amplitude for pdf
#amppdf = 1.0/m.sqrt(2*m.pi*variance)
def setgaussamppdf(var):
    if var <= 0: sys.exit("setgaussamppdf variance must be positive")
    amppdf = 1.0/m.sqrt(2*m.pi*var)
    return amppdf

def gausseval(point,center,var,amp):
    if var <= 0: sys.exit("gausseval variance must be positive")
    geval = amp*m.exp(-0.5*(mu-point)**2/var)
    return geval

def gaussint(loB,upB,step,center,var,amp):
    if var <= 0: sys.exit("gaussint variance must be positive")
    if upB <= loB: sys.exit("gaussint upper bound must be higher than lower bound (%f,%f)"%(loB,upB))
    if step > (upB-loB): step = (upB-loB)/2.0
    if step <= 0: sys.exit("gaussint step size must be positive")
    gintsum = 0.0
    giX = loB + 0.5*step
    gintsum += step*gausseval(giX,center,var,amp)
    while (giX + 1.5*step) <= upB:
        giX += step
        gintsum += step*gausseval(giX,center,var,amp)
    gintsum += (upB - (giX + 0.5*step))*gausseval(0.5*(upB + (giX + 0.5*step)),center,var,amp)
    return gintsum

def gausspdfbin(loB=0,upB=1,step=0.001,center=0.5,var=0.5,loBeqNI=False,upBeqI=False):
    gtemptab = gausscdftable
    #gtemptab[0] = 0.5
    if var <= 0: sys.exit("gausspdfbin variance must be positive")
    std = m.sqrt(var)
    if ((loBeqNI) & (upBeqI)): return 1.0
    amp = setgaussamppdf(var)
    gpdfsum = 0.0
    if loBeqNI: 
        if upB < center:
            sigs = m.ceil((center-upB)/m.sqrt(var)+0.01)
            if sigs > 6: sys.exit("gausspdfbin bounds must be within 6 sigma of mean")
            gpdfsum = 0.5*(1-gtemptab[int(sigs)])
            loB = center-sigs*std
        if upB >= center:
            sigs = m.floor((upB-center)/m.sqrt(var)-0.01)
            if sigs > 6: sys.exit("gausspdfbin bounds must be within 6 sigma of mean")
            if sigs >= 0: 
                gpdfsum = 0.5*gtemptab[int(sigs)] + 0.5
            else:
                gpdfsum = 0.5*(1-gtemptab[int(-1*sigs)])
            loB = center+sigs*std
    if upBeqI:
        if loB > center:
            sigs = m.ceil((loB-center)/m.sqrt(var)+0.01)
            if sigs > 6: sys.exit("gausspdfbin bounds must be within 6 sigma of mean")
            gpdfsum = 0.5*(1-gtemptab[int(sigs)])
            upB = center+sigs*std
        if loB <= center:
            sigs = m.floor((center-loB)/m.sqrt(var)-0.01)
            if sigs > 6: sys.exit("gausspdfbin bounds must be within 6 sigma of mean")
            if sigs >= 0: 
                gpdfsum = 0.5*gtemptab[int(sigs)] + 0.5
            else:
                gpdfsum = 0.5*(1-gtemptab[int(-1*sigs)])
            upB = center-sigs*std
    pdfbin = gaussint(loB,upB,step,center,var,amp) + gpdfsum
    return pdfbin



