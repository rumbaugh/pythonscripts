import numpy as np
import sys

EPS1 = 0.001
EPS2 = 10.0**-8

def probks(alam):
    fac = 2.0
    psum = 0.0
    termbf = 0.0
    a2 = -2.0*alam*alam
    for j in range(1,101):
        term = fac*np.exp(a2*j*j)
        psum += term
        if ((j>3) & ((np.fabs(term) <= EPS1*termbf) | (np.fabs(term) <= EPS2*psum))): 
            return psum
        fac *= -1
        termbf = np.fabs(term)
    return 1.0

def KStest(arr1,arr2):
    #Inputs should not be CDFs
    maxdif = 0.0
    sort1,sort2 = np.sort(arr1),np.sort(arr2)
    for iks in range(0,len(sort1)):
        g1 = np.where(sort2 < sort1[iks])
        frac1 = np.fabs(iks/(1.0*len(sort1))-len(g1[0])/(1.0*len(sort2)))
        if frac1 > maxdif: maxdif = frac1
    for iks in range(0,len(sort2)):
        g2 = np.where(sort1 < sort2[iks])
        frac2 = np.fabs(iks/(1.0*len(sort2))-len(g2[0])/(1.0*len(sort1)))
        if frac2 > maxdif: maxdif = frac2
    Ne = len(arr1)*len(arr2)/(1.0*len(arr1)+len(arr2))
    prob = probks(maxdif*(np.sqrt(Ne) + 0.12 + 0.11/np.sqrt(Ne)))
    return maxdif,prob

def KStest1var(arr1,cdf1):
    #values of cdf1 are cumulative probabilities for corresponding values of arr1 
    #based on the distribution being tested as source of arr1
    if len(cdf1) != len(arr1): sys.exit("length of input array and cdf1 must be the same for 1 variable KS test")
    maxdif = 0.0
    sort1args = np.argsort(arr1)
    sort1 = np.sort(arr1)
    for iks in range(0,len(sort1)):
        if iks < len(sort1)-1: 
            frac = np.fabs((iks+1)/(1.0*len(sort1))-cdf1[sort1args[iks]])
            if frac > maxdif: maxdif = frac
        frac = np.fabs(iks/(1.0*len(sort1))-cdf1[sort1args[iks]])
        if frac > maxdif: maxdif = frac
    Ne = len(arr1)
    prob = probks(maxdif*(np.sqrt(Ne) + 0.12 + 0.11/np.sqrt(Ne)))
    return maxdif,prob
