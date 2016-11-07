import numpy as np
import os
import math as m
import sys
import time

def LQFit(X,Y):
    #Fits y = a + bx
    N = len(X)
    if len(Y) != N: sys.exit("LQ input arrays must have equal lengths")
    Del = N*sum(X**2)-(sum(X))**2
    a = (sum(Y)*sum(X**2)-sum(X)*sum(X*Y))/Del
    b = (N*sum(X*Y)-sum(X)*sum(y))/Del
    return a,b

def LQFitYErr(X,Y,Yerr):
    #Fits y = a + bx and returns error in a and b, assuming the error
    #in X is approximately zero
    N = len(X)
    if len(Y) != N: sys.exit("LQ input arrays must have equal lengths")
    if len(Yerr) != N: sys.exit("LQ input arrays must have equal lengths")
    Del = N*sum(X**2)-(sum(X))**2
    a = (sum(Y)*sum(X**2)-sum(X)*sum(X*Y))/Del
    b = (N*sum(X*Y)-sum(X)*sum(Y))/Del
    da = m.sqrt((sum(Yerr**2)*(sum(X**2))**2+(sum(X))**2*sum((X*Yerr)**2))/(Del**2))
    db = m.sqrt((N**2*sum((X*Yerr)**2)+(sum(X))**2*sum(Yerr**2))/(Del**2))
    return a,b,da,db

