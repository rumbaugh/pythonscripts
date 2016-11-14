import numpy as np
import easyaccess as ea
import shapely.geometry as geo

con_op=ea.connect(section='desoper')

q2='SELECT f.RA_CENT,f.DEC_CENT,f.RACMAX,f.RACMIN,f.DECCMAX,f.DECCMIN,f.RAC1,f.RAC2,f.RAC3,f.RAC4,f.DECC1,f.DECC2,f.DECC3,f.DECC4,f.ID,d.TILENAME FROM RUMBAUGH.Y1A1_TILENAMES_LIST d,FELIPE.COADDTILE_GEOM_NEW f WHERE d.tilename=f.tilename'

TDF=con_op.query_to_pandas(q2)

Tpolys=np.zeros(len(TDF),dtype='object')
altpolys=np.zeros(0,dtype='object')
altpolynames=np.zeros(0,dtype='|S30')
altpolynums=np.zeros(0,dtype='i8')
altpolyracens,altpolydeccens=np.zeros(0,dtype='f8'),np.zeros(0,dtype='f8')
for i in range(0,len(TDF)):
    Tpolys[i]=geo.Polygon([(TDF['RAC1'][i],TDF['DECC1'][i]),(TDF['RAC2'][i],TDF['DECC2'][i]),(TDF['RAC3'][i],TDF['DECC3'][i]),(TDF['RAC4'][i],TDF['DECC4'][i])])
    if ((TDF['RACMAX'][i]<TDF['RACMIN'][i]) | (TDF['RACMIN'][i]<0) | (TDF['RACMAX'][i]>360)):
        ractmps=np.array([TDF['RAC1'][i],TDF['RAC2'][i],TDF['RAC3'][i],TDF['RAC4'][i]])
        ractmps[ractmps>300]-=360
        Tpolys[i]=geo.Polygon([(ractmps[0],TDF['DECC1'][i]),(ractmps[1],TDF['DECC2'][i]),(ractmps[2],TDF['DECC3'][i]),(ractmps[3],TDF['DECC4'][i])])
        ractmps=np.array([TDF['RAC1'][i],TDF['RAC2'][i],TDF['RAC3'][i],TDF['RAC4'][i]])
        ractmps[ractmps<300]+=360
        altpolys=np.append(altpolys,geo.Polygon([(ractmps[0],TDF['DECC1'][i]),(ractmps[1],TDF['DECC2'][i]),(ractmps[2],TDF['DECC3'][i]),(ractmps[3],TDF['DECC4'][i])]))
        altpolynames=np.append(altpolynames,TDF['TILENAME'][i])
        altpolynums=np.append(altpolynames,i)
        altpolydeccens=np.append(altpolydeccens,TDF['DEC_CENT'][i])
        if TDF['RA_CENT'][i]>300: TDF['RA_CENT'][i]-=360
        altpolyracens=np.append(altpolyracens,TDF['RA_CENT'][i]+360)
            

def find_tile(rap,decp):
    tilenameout='None'
    tilenamematches=''
    curp=geo.Point(rap,decp)
    gct=np.where((np.abs(decp-TDF['DEC_CENT'])<2.)&(np.abs(rap-TDF['RA_CENT'])*np.cos(decp*np.pi/180)<2.))[0]
    gcta=np.where((np.abs(decp-altpolydeccens)<2.)&(np.abs(rap-altpolyracens)*np.cos(decp*np.pi/180)<2.))[0]
    for i in range(0,len(gct)):
        if Tpolys[gct[i]].contains(curp): tilenamematches='%s;%s'%(tilenamematches,TDF['TILENAME'][gct[i]])
    for i in range(0,len(gcta)):
        if altpolys[gcta[i]].contains(curp): tilenamematches='%s;%s'%(tilenamematches,altpolynames[i])
    if tilenamematches!='': 
        if tilenamematches[0]==';': tilenamematches=tilenamematches[1:]
        tilenameout=tilenamematches
    return tilenameout
            
