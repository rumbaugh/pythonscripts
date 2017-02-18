import numpy as np

def SphDist(RAi,Deci,Raf,Decf):
	#Distance in arcminutes. Uses Haversine Formula
	dist = 2.0*np.arcsin(np.sqrt((np.sin(0.5*np.pi*(Decf-Deci)/180))**2+np.cos(np.pi*Decf/180.0)*np.cos(np.pi*Deci/180.0)*(np.sin(0.5*np.pi*(Raf-RAi)/180))**2))
	dist *= 180.0*60/np.pi
	return dist
