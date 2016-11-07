import math as m
import numpy as np
import sys
def LinReg(x,y,z=np.zeros(0),wx=np.zeros(0),wy=np.zeros(0),wz=np.zeros(0),numvar=2,A=-99,B=-99,C=-99,err=False):
    #Performs linear regression to fit two or three variable data to the equation
    # y = A + Bx (2 vars) or z = A + Bx + Cy (3 vars). If numvar = 2, only two 
    # variables are fit. If A, B, or C are set (to something other than -99), they
    # are assumed fixed. Currently, only C can be fixed. Weights can be set for y
    # for 2 variables, and for either z or both z and y (only if C is fixed) for 3 
    # variables. Errors are output for A, B, and C if err is set to True. 
    #len(x),len(y),len(z),len(wx),len(wy),len(wz)
    #print 'a'
    wxset = 0
    wyset = 0
    wzset = 0
    numcons = 0
    if A == -99: numcons += 1
    if B == -99: numcons += 1
    if ((C == -99) & (numvar > 2)): numcons += 1
    if ((numcons != 1) & (numcons != 2) & (numcons != 3)): sys.exit("Can only solve for 1, 2, or 3 coefficients")
    if ((numvar != 2) & (numvar != 3)): sys.exit("Can only have two or three variables")
    if numcons > numvar: sys.exit("Can't have more coefficients than variables")
    if ((A != -99) | (B != -99)): sys.exit("Currently can't freeze A or B")
    if ((len(wx) > 0) & (len(wx) != len(x))): "x weights must be same length as x"
    if ((len(wy) > 0) & (len(wy) != len(y))): "y weights must be same length as y"
    if ((len(wz) > 0) & (len(wz) != len(z))): "z weights must be same length as z"
    if len(wx) > 0:
        gwxz = np.where(wx == 0)
        wxset = 1
        if len(gwxz[0]) > 0.1: sys.exit("All x weights must be non-zero")
    if len(wy) > 0:
        gwyz = np.where(wy == 0)
        wyset = 1
        if len(gwyz[0]) > 0.1: sys.exit("All y weights must be non-zero")
    if len(wz) > 0:
        gwzz = np.where(wz == 0)
        wzset = 1
        if len(gwzz[0]) > 0.1: sys.exit("All z weights must be non-zero")
    if len(wx) == 0: wx = np.ones(len(x))
    if len(wy) == 0: wy = np.ones(len(y))
    if len(wz) == 0: wz = np.ones(len(z))
    if wxset != 0:
        if numvar == 2: 
            sys.exit("Only weights on y can be set when using 2 variables")
        elif numvar == 3: 
            #print wxset,wyset,wzset,A,B,C
            if ((C == -99) | (A != -99) | (B != -99)):
                sys.exit("For 3 variables, weights can be set for z only, or for y and z with C frozen")
        else:
            sys.exit("Shouldn't be able to get here")
    if wyset != 0:
        if ((numvar == 3) & (C == -99)): sys.exit("Weights for y can only be set if C is frozen")
        # Calculating these weights is major problem. I am not entirely sure
        # what is the best way to do it
        #if ((numvar == 3) & (C != -99)): wz = 1.0/(C*C*1.0/(wy)+1.0/wz)
        if ((numvar == 3) & (C != -99)): wz = wz*wy
    errx,erry,errz = wx**-2,wy**-2,wz**-2
    varx,vary,varz = 1.0/wx,1.0/wy,1.0/wz
    if wxset == 0: varx = np.zeros(len(wx))
    if wyset == 0: vary = np.zeros(len(wy))
    if wzset == 0: varz = np.zeros(len(wz))
    if numvar == 2:
        Sx = np.sum(wy*x)
        Sxx = np.sum(wy*x*x)
        Sy = np.sum(wy*y)
        Sxy = np.sum(wy*x*y)
        Sw = np.sum(wy)
        Del = Sw*Sxx-Sx*Sx
        A = (Sxx*Sy-Sx*Sxy)/Del
        B = (Sw*Sxy-Sx*Sy)/Del
        #print wy,Sx,Sxx,Sy,Sxy,Sw,Del,A,B
        if not err:
            return A,B
        else:
            varSx = np.sum(wy*wy*varx)
            varSxx = np.sum(4*wy*wy*x*varx)
            varSy = np.sum(wy*wy*vary)
            varSxy = np.sum(wy*wy*(y*y*varx+x*x*vary))
            varDel = Sw*Sw*varSxx+4*Sx*Sx*varSx
            sigy = np.sqrt(1./(len(y)-2)*np.sum((y-A-B*x)**2))
            errA = np.sqrt(Sxx/Del)*sigy
            errB = np.sqrt(Sw/Del)*sigy
            return A,B,errA,errB
    elif numvar == 3:
        Sx = np.sum(wz*x)
        Sy = np.sum(wz*y)
        Sz = np.sum(wz*z)
        Sxx = np.sum(wz*x*x)
        Syy = np.sum(wz*y*y)
        Sxy = np.sum(wz*x*y)
        Sxz = np.sum(wz*x*z)
        Syz = np.sum(wz*y*z)
        Sw = np.sum(wz)
        if C == -99:
            Del = Sw*(Sxy*Sxy-Sxx*Syy)-2*Sx*Sxy*Sy+Sxx*Sy*Sy+Sx*Sx*Syy
            A = (Sxz*(Sx*Syy-Sxy*Sy)+Sxx*(Sy*Syz-Syy*Sz)+Sxy*(Sxy*Sz-Sx*Syz))/Del
            B = (Sw*(Sxy*Syz-Sxz*Syy)+Sx*(Syy*Sz-Sy*Syz)+Sy*(Sxz*Sy-Sxy*Sz))/Del
            C = (Sw*(Sxy*Sxz-Sxx*Syz)+Sx*(Sx*Syz-Sxy*Sz)+Sy*(Sxx*Sz-Sx*Sxz))/Del
            if not err:
                return A,B,C
            else:
                varSx = np.sum(wz*wz*varx)
                varSy = np.sum(wz*wz*vary)
                varSz = np.sum(wz*wz*varz)
                varSxx = np.sum(4*wz*wz*x*x*varx)
                varSyy = np.sum(4*wz*wz*y*y*vary)
                varSxy = np.sum(wz*wz*(x*x*vary+y*y*varx))
                varSxz = np.sum(wz*wz*(x*x*varz+z*z*varx))
                varSyz = np.sum(wz*wz*(z*z*vary+y*y*varz))
                varDel = np.sum(Sw*Sw*(4*Sxy*Sxy*varSxy+Sxx*Sxx*varSyy+Syy*Syy*varSxx)+4*(Sx*Sx*Sxy*Sxy*varSy+Sxy*Sxy*Sy*Sy*varSx+Sx*Sx*Sy*Sy*varSxx)+Sy*Sy*(4*Sxx*Sxx*varSy+Sy*Sy*varSxx)+Sx*Sx*(4*Syy*Syy*varSx+Sx*Sx*VarSyy))
                errA = (1.0/m.fabs(Del))*m.sqrt((Sx*Syy-Sxy*Sy)**2*varSxz+Sxz*Sxz*(Sx*Sx*varSyy+varSx*Syy*Syy+Sxy*Sxy*varSy+Sy*Sy*varSxy)+varSxx*(Sy*Syz-Syy*Sz)**2+Sxx*Sxx*(Sy*Sy*varSyz+varSy*Syz*Syz+Sz*Sz*varSyy+Syy*Syy*varSz)+varSxy*(Sx*Syz)**2+Sxy*Sxy*(Sxy*Sxy*varSz+Sx*Sx*varSyz+varSx*Syz*Syz)+4*Sxy*Sxy*Sz*Sz*varSxy+A*A*varDel)
                errB = (1.0/m.fabs(Del))*m.sqrt(Sw*Sw*(Sxy*Sxy*varSyz+varSxy*Syz*Syz+Sxz*Sxz*varSyy+varSxz*Syy*Syy)+varSx*(Syy*Sz-Sy*Syz)**2+Sx*Sx*(Syy*Syy*varSz+Sz*Sz*varSyy+Sy*Sy*varSyz+varSy*Syz*Syz)+Sy*Sy*(Sxy*Sxy*varSz+Sz*Sz*varSxy)+varSy*Sxy*Sxy*Sz*Sz+varSxz*Sy*Sy*Sy*Sy+4*Sy*Sy*Sxz*Sxz*varSy+A*A*varDel)
                errC = (1.0/m.fabs(Del))*m.sqrt(Sw*Sw*(Sxy*Sxy*varSxz+Sxz*Sxz*varSxy+Sxx*Sxx*varSyz+Syz*Syz+varSxx)+varSy*(Sxx*Sz-Sz*Sxz)**2+Sy*Sy*(Sxx*Sxx*varSz+varSxx*Sz*Sz+Sx*Sx*varSxz+varSx*Sxz*Sxz)+varSz*(Sxy*Sxy*Sz*Sz)+Sx*Sx*(Sx*Sx*varSyzSxy*Sxy*varSz+varSxy*Sz*Sz)+4*Sx*Sx*Syz*Syz*varSx+A*A*varDel)
                return A,B,C,errA,errB,errC
        elif ((C != -99) & (A == -99) & (B == -99)):
            #print varx,vary,varz
            A = (Sz-Sx*Sxz*1.0/Sxx+C*(Sxy*Sx*1.0/Sxx-Sy))*1.0/(Sw-Sx*Sx*1.0/Sxx)
            B = (C*(Sw*Sxy-Sx*Sy)+Sz*Sx-Sw*Sxz)/(Sx*Sx-Sw*Sxx)
            if not err:
                return A,B
            else:
                varSx = np.sum(wz*wz*varx)
                varSy = np.sum(wz*wz*vary)
                varSz = np.sum(wz*wz*varz)
                varSxx = np.sum(4*wz*wz*x*x*varx)
                varSyy = np.sum(4*wz*wz*y*y*vary)
                varSxy = np.sum(wz*wz*(x*x*vary+y*y*varx))
                varSxz = np.sum(wz*wz*(x*x*varz+z*z*varx))
                errA = (1.0/m.fabs(Sw-Sx*Sx*1.0/Sxx))*m.sqrt(varSz+varSx*Sxz*Sxz*1.0/Sxx/Sxx+Sx*Sx*varSxz*1.0/Sxx/Sxx+Sx*Sx*Sxz*Sxz*varSxx*(Sxx**-4)+C*C*(varSxy*Sx*Sx*1.0/Sxx/Sxx+Sxy*Sxy*varSx*1.0/Sxx/Sxx+Sxy*Sxy*Sx*Sx*varSxx*Sxx**-4+varSy)+A*A*(varSxx*(Sx*1.0/Sxx)**4+4*varSx*(Sx*1.0/Sxx)**2))
                errB = (1.0/Sx)*m.sqrt(varSx*B*B+varSz+C*C*varSy+Sw*Sw*errA*errA)
                return A,B,errA,errB
        else:
            sys.exit("Shouldn't be able to get here - 2")
            return 0
    else:
        sys.exit("Shouldn't be able to get here - 3")
        return 0

