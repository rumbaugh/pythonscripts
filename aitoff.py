import numpy as np

def aitoff(ra,dec):
    psi,lam=dec*np.pi/180.,(ra-180)*np.pi/180.
    alpha=np.arccos(np.cos(psi)*np.cos(0.5*lam))
    x,y=(2*np.cos(psi)*np.sin(0.5*lam))/np.sinc(alpha/np.pi),np.sin(psi)/np.sinc(alpha/np.pi)
    return x,y
