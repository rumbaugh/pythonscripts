import math as m
import numpy as np

def Interpolate(func,x,intpos):
    g = np.where(x == intpos)
    if len(g) > 0: 
        g = g[0]
        if len(g) > 0:
            return func[g]
    insertpos = np.searchsorted(x,intpos)
    if insertpos == 0:
#        intval = 0.0
        intval = func[0]+(func[0]-func[1])*(x[0]-intpos)/(1.0*x[1]-x[0])
    elif insertpos == len(x):
#        intval = 0.0
        intval = func[len(x)-1]+(func[len(x)-1]-func[len(x)-2])*(intpos-x[len(x)-1])/(x[len(x)-1]-x[len(x)-2])
    else:
        intval = func[insertpos-1]+(func[insertpos]-func[insertpos-1])*(intpos-x[insertpos-1])/(x[insertpos]-x[insertpos-1])
    return intval

def KCorCalc(FE,FEwav,flux,galwav,z):
    USum = 0.0
    LSum = 0.0
    galwav_prop = galwav*(1.0+z)
    flux_prop = flux*(1.0+z)
#    LBind = np.searchsorted(galwav,4770)
#    UBind = np.searchsorted(galwav,6480)
#    LBindL = np.searchsorted(galwav,4770/(1+z))
#    UBindL = np.searchsorted(galwav,6480/(1+z))
#    dlamda = 0.5*(galwav[LBind+1]-galwav[LBind])
#    for i in range(LBind,UBind+1):
    dlamda = 0.5*(FEwav[1]-FEwav[0])
    for i in range(0,len(FEwav)):
        Uflux_temp = Interpolate(flux,galwav,FEwav[i])
        Uintegrand = FEwav[i]*FE[i]*Uflux_temp*dlamda
#        FEtemp = Interpolate(FE,FEwav,galwav[i])
#        Uintegrand = FEtemp*flux[i]*galwav[i]*dlamda
        USum += Uintegrand
#        if i == UBind-1:
#            dlamda = 0.5*(galwav[UBind]-galwav[UBind-1])
#        elif i < UBind-1:
#            dlamda = 0.5*(galwav[i+2]-galwav[i])
#        if i == len(FEwav)-2:
#            dlamda = 0.5*(FEwav[len(FEwav)-1]-FEwav[len(FEwav)-2])
#        elif i < len(FEwav)-1:
#            dlamda = 0.5*(FEwav[i+2]-FEwav[i])
#    dlamda = 0.5*(galwav[LBindL+1]-galwav[LBindL])
#    for i in range(LBindL,UBindL+1):
        Lflux_temp = Interpolate(flux_prop,galwav_prop,FEwav[i])
        Lintegrand = FEwav[i]*FE[i]*Lflux_temp*dlamda
#        FEtemp = Interpolate(FE,FEwav,galwav[i]*(1.0+z))
#        Lintegrand = FEtemp*flux[i]*galwav[i]*dlamda
        LSum += Lintegrand
        if i == len(FEwav)-2:
            dlamda = 0.5*(FEwav[len(FEwav)-1]-FEwav[len(FEwav)-2])
        elif i < len(FEwav)-1:
            dlamda = 0.5*(FEwav[i+2]-FEwav[i])
#        if i == UBindL-1:
#            dlamda = 0.5*(galwav[UBindL]-galwav[UBindL-1])
#        elif i < UBindL-1:
#            dlamda = 0.5*(galwav[i+2]-galwav[i])
#    print USum/((1.0+z)*LSum)
#    print -2.5*m.log10(USum/((1.0+z)*LSum))
#    return 2.5*m.log10(USum/((1.0+z)*LSum))
#    return -2.5*m.log10(USum*(1.0+z)/LSum)
    return 2.5*m.log10((1+z)*USum/LSum)

crV = read_file("/home/rumbaugh/Download/V.dat")
crE = read_file("/home/rumbaugh/Download/E.flux.dat")
crS = read_file("/home/rumbaugh/Download/Sbc.flux.dat")

Vwav = get_colvals(crV,'col1')
VFE = get_colvals(crV,'col2')
Ewav = get_colvals(crE,'col1')
Eflux = get_colvals(crE,'col2')
Swav = get_colvals(crS,'col1')
Sflux = get_colvals(crS,'col2')

crCos = read_file("/home/rumbaugh/cosmocalcoutput.nh.txt")
zs = get_colvals(crCos,'col1')
EKcors = np.zeros(len(zs))
SbcKcors = np.zeros(len(zs))
FILE = open("ast267.kcors.dat","w")
j = 0
while zs[j] <= 1.0:
    EKcors[j] = KCorCalc(VFE,Vwav,Eflux,Ewav,zs[j])
    SbcKcors[j] = KCorCalc(VFE,Vwav,Sflux,Swav,zs[j])
    FILE.write(str(zs[j]) + " " + str(EKcors[j]) + " " + str(SbcKcors[j])+ "\n")
    j += 1
FILE.close()
