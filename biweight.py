import numpy as np
import math as m
import sys

def MAD(data):
    #calculates the median absolute deviation
    tmed = np.median(data)
    ADtemp = data-tmed
    for iAD in range(0,len(data)): ADtemp[iAD] = m.fabs(ADtemp[iAD])
    return np.median(ADtemp)

def biweight_loc(data,biweight_loc_tuning_con=6.0):
    #calculates the biweight location estimator of data
    tmed = np.median(data)
    tMAD = MAD(data)
    u_i = (data-tmed)/(biweight_loc_tuning_con*tMAD)
    abs_u_i = u_i
    for i_u_i in range(0,len(u_i)): abs_u_i[i_u_i] = m.fabs(u_i[i_u_i])
    g_abs = np.where(abs_u_i < 1)
    g_abs = g_abs[0]
    C_BI = tmed + (np.sum((data[g_abs]-tmed)*(1-u_i[g_abs]**2)**2)/(np.sum((1-u_i[g_abs]**2)**2)))
    return C_BI

def biweight_scale(data,biweight_scale_tuning_con=9.0):
    #calculates the biweight scale estimator of data
    #note that this stuff is from Beer et al. (1990)
    tmed = np.median(data)
    tMAD = MAD(data)
    u_i = (data-tmed)/(biweight_scale_tuning_con*tMAD)
    abs_u_i = u_i
    for i_u_i in range(0,len(u_i)): abs_u_i[i_u_i] = m.fabs(u_i[i_u_i])
    g_abs = np.where(abs_u_i < 1)
    g_abs = g_abs[0]
    S_BI = m.sqrt(len(data))*m.sqrt(np.sum((data[g_abs]-tmed)**2*(1-u_i[g_abs]**2)**4))/m.fabs(np.sum((1-u_i[g_abs]**2)*(1-5*u_i[g_abs]**2)))
    return S_BI
    
    
