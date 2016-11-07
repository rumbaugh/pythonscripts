import numpy as np
import math as m
import time
import sys

def ML_find_a(a,isderiv):
    if isderiv == 0:
        return N/a-Sum-N*m.log(b)/(b**a)
    else:
        return -N/(a*a)+N*m.log(b)*m.log(b)*b**(-a)

try:
    func
except NameError:
    func = ML_find_a

def NewtonRaphson(init,func,tol):
    xtemp = init
    err = 999.9
    while err > tol:
        xtemplast = xtemp
        xtemp = xtemplast - func(xtemplast,0)/func(xtemplast,0)
        err = m.fabs((xtemp-xtemplast)/xtemp)
    return xtemp
