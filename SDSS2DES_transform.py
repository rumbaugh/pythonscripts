import numpy as np

def SDSS2DES_mag(u_sdss,g_sdss,r_sdss,i_sdss,z_sdss):
    u_des = u_sdss - 0.479 + 0.466*(g_sdss-r_sdss) - 0.350*(g_sdss-r_sdss)**2 
    g_des = g_sdss + 0.001 - 0.075*(g_sdss-r_sdss) 
    r_des = r_sdss - 0.009 - 0.069*(g_sdss-r_sdss) 
    i_des = i_sdss + 0.014 - 0.214*(i_sdss-z_sdss) - 0.096*(i_sdss-z_sdss)**2
    z_des = z_sdss + 0.022 - 0.068*(i_sdss-z_sdss) 
    Y_des = z_sdss + 0.045 - 0.306*(i_sdss-z_sdss)
    return u_des,g_des,r_des,i_des,z_des,Y_des
