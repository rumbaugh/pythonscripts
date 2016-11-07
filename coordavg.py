import numpy as np
import math as m
import sys

infile = "posse_coords.txt"
cr = read_file(infile)
RA = get_colvals(cr,'col2')
DEC = get_colvals(cr,'col1')

length = len(RA)
X = np.zeros(length)
Y = np.zeros(length)
Z = np.zeros(length)

theta = np.zeros(length)
phi = np.zeros(length)

for i in range(0,length):
    #setting R = 1
    theta[i] = RA[i]*m.pi/180.0
    phi[i] = m.pi/2.0 - DEC[i]*m.pi/180.0
    X[i] = m.cos(theta[i])*m.sin(phi[i])
    Y[i] = m.sin(theta[i])*m.sin(phi[i])
    Z[i] = m.cos(phi[i])

aX = sum(X)/length
aY = sum(Y)/length
aZ = sum(Z)/length

aR = m.sqrt(aX*aX + aY*aY + aZ*aZ)
aTheta = m.atan(aY/aX)
aPhi = m.acos(aZ/aR)

# ArcTan only gives out numbers in a range from -90 to 90 degrees
# Right here we check to see if this is giving us problems, and if it
# is, we fix it
Xchk = aR*m.cos(aTheta)*m.sin(aPhi)
if m.fabs(Xchk - aX) > 0.1:
    aTheta -= m.pi
    Xchk2 = aR*m.cos(aTheta)*m.sin(aPhi)
    if m.fabs(Xchk2 - aX) > 0.1:
        sys.error("Degeneracy problems")

aRA = 180*aTheta/m.pi
aDEC = 180*(m.pi/2.0-aPhi)/m.pi
#print aRA,aDEC

if aRA > 180:
    EorW = "W"
    aRA = 360 - aRA
elif aRA < 0 and aRA >= -180: 
    EorW = "W"
    aRA *= -1
elif aRA < -180: 
    EorW = "W"
    aRA += 360
else:
    EorW = "E"

if aDEC > 90:
    NorS = "S"
    aDEC = 180 - aDEC
elif aDEC < 0:
    NorS = "S"
    aDEC *= -1
else:
    NorS = "N"


print "Average Position is " + str(aDEC) + " " + NorS + ", "+ str(aRA) + " " + EorW
