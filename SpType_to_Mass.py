# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 09:15:45 2020

@author: Sophia Stasevic

Outputs text file with the following columns: RA, Dec, Effective Temp, Mass, Spectral Type

Use SpType_logT.txt for clusters closer to 3.5 Myrs old, and SpType_logT.txt for clusters
closer to 13 Myrs old
"""

import numpy as np
from astropy.table import Table 
import math


def OBASort():
    
    stars=Table.read('NGC2362_McSwain_Unique_2MASS.txt', format='ascii') #table containing star data
    stype=np.loadtxt('SpType_logT.txt',dtype=str) #table of conversions between spectral type and teff
    iso=np.loadtxt('MIST_iso_5_Myrs.cmd') #MIST iscochroen for cluster age containing relevant spectroscopic band information
    
    #if Teff is known but not SpT:
    """
    #teff=np.log10(np.array(stars['Teff'])) #or
    teff=np.array(stars['logT'])
    """
    
    #for stars with spectral type but not Teff:
    spt=np.array(stars['SpT'])
    #spt=StringSplit(spt) #if the spectral type numbers are decimals
    teff=np.empty(len(spt))
    data = ["" for x in range(len(spt))]
    for i in range(0,len(data)):
        if len(spt[i])<3:
            data[i]= '{star}V'.format(star=spt[i]) #assumes all stars with no specified luminosity class are class V
        else:
            data[i]=spt[i]
        same_type=np.where(stype[:,0]==data[i])[0]
        if same_type.size<1:
            data[i]=ManualSType(spt,stype,i)
            same_type=np.where(stype[:,0]==data[i])[0]
        same_type=int(same_type)
        teff[i]=stype[same_type][1]
    
    AbsMag=MagConversion(stars)
    
    mass=MagTeffInterpol(AbsMag,teff,iso) #teff known
    final=np.column_stack((np.array(stars['_RAJ2000']),np.array(stars['_DEJ2000']),mass,teff,stars['SpT'])) #when SpT is known
    
    
    """
    mass=MagInterpol(AbsMag,iso) #teff unknown
    teff=mass[1]
    mass=mass[0]
    spt=TefftoSpT(teff)
    final=np.column_stack((np.array(stars['_RA']),np.array(stars['_DE']),mass,teff,spt)) #when SpT isn't known
    
    """
    
    np.savetxt('NGC2362_McSwain_2MASS_Masses.txt',final, delimiter=" ", fmt="%s")          

def ManualSType(spt,stype,i):
    
    print('The spectral type is {type}'.format(type=spt[i]))
    s=input('Please input the new spectral type: ')
    
    return s

def MagConversion(stars):
    
    import SpinRatesAreHere as spin 
    cluster = spin.NGC2362()
    cluster.ClusterInfo()
    
    #conversion from Av depends on what spectroscopic bands data is available for
    med_Av=0.31
    mj=stars['Jmag']
    #Av=3.2*stars['E_B-V_']
    #Av=stars['Av']
    Av=med_Av
    #for i in range(len(Av)):
        #if Av[i]==0:
            #Av[i]=med_Av
    
    #BP_RP=stars['BP-RP'] #if BP-RP is known but only median A is known --> only for calculating G band extinction
    #AG=makeAvGaia(BP_RP,med_Av)   
    Aj=0.807*Av
    
    M=np.array(mj-5*np.log10(cluster.dist/10)-Aj)
    #M=np.array(mj-5*np.log10(100/stars['Plx'])-Aj) #if parallax for individual stars is known
  
    return M

def MagTeffInterpol(AbsMag,iso,teff):
    
    #column in isochrone corresponding to the band the absolute magnitude was calculated for
    upperVmag=np.where(iso[:,14]>-4.4)[0] #want the point where the isochrone curves back on itself
    lowerVmag=np.where(iso[:,14]<=-4.4)[0]
    
    mass=np.empty(len(AbsMag))
    for i in range(0,len(AbsMag)):
        if AbsMag[i]>-4.4:
            mass[i]=np.interp(teff[i],iso[:,4][upperVmag],iso[:,2][upperVmag])
        else:
            mass[i]=np.interp(teff[i],iso[:,4][lowerVmag],iso[:,2][lowerVmag])

    return mass

def MagInterpol(AbsMag,iso):
    
    sort=np.argsort(iso[:,14])
    magsort=iso[:,14][sort]
    masssort=iso[:,2][sort]
    teffsort=iso[:,4][sort]
    
    upperVmag=np.where(magsort>-4.4)[0] #want the point where the isochrone curves back on itself
    lowerVmag=np.where(magsort<=-4.4)[0]
    
    mass=np.empty(len(AbsMag))
    teff=np.empty(len(AbsMag))
    for i in range(0,len(AbsMag)):
        if AbsMag[i]>-4.4:
            teff[i]=np.interp(AbsMag[i],magsort[upperVmag],teffsort[upperVmag])
            mass[i]=np.interp(AbsMag[i],magsort[upperVmag],masssort[upperVmag])
        else:
            teff[i]=np.interp(AbsMag[i],magsort[lowerVmag],teffsort[lowerVmag])
            mass[i]=np.interp(AbsMag[i],magsort[lowerVmag],masssort[lowerVmag])
            
    return mass, teff

def makeAvGaia(BPRP,A0): #function written by Julia Roquette
        """
        Define the parameters needed for building Gaia reddening vectors
        https://www.aanda.org/articles/aa/abs/2018/08/aa32843-18/aa32843-18.html
        """
        cG=[0.9761,-0.1704,0.0086,0.0011,-0.0438, 0.0013, 0.0099]
        cBP=[1.1517, -0.0871, -0.0333, 0.0173, -0.0230, 0.0006, 0.0043] 
        cRP=[0.6104,-0.0170,-0.0026, -0.0017, -0.0078, 0.00005, 0.0006]                 
        def AvGaia(cX=cG,BPRP=BPRP,A0=A0):
            """
            cX is the array containing the c values for a certain X band
            """
            return A0*(cX[0]+cX[1]*BPRP+cX[2]*(BPRP**2)+cX[3]*(BPRP**3)+cX[4]*A0+cX[5]*A0*A0+cX[6]*BPRP*A0)              
        star_AG=AvGaia(cX=cG,BPRP=BPRP,A0=A0)
        
        return star_AG

#rounds the number of the spectral type to the nearest integer/.5    
def StringSplit(data):
    
    for i in range(len(data)):
        a=data[i][0]
        b=[data[i][1],data[i][2],data[i][3]] #extend to number of decimal places
        x=float(''.join(b))
        if x>9.2:
            x=9
        if a=='O' or (a=='B' and x<2.8): #for Currie conversion, change to just O types
            y=str(round(x*2)/2) #rounds to nearest 0.5
            if int(y[2])!=5:
                y=str(round(x))
        else:
            y=str(round(x))
        z=[a,y]
        data[i]=(''.join(z))
        
    return data

def TefftoSpT(teff):
    
    stype=np.loadtxt('MatchSpT.txt',dtype='str')
    
    spt=["" for x in range(len(teff))]
    table_spt=np.empty(len(stype[:,1]))
    
    for i in range(len(stype[:,1])):
        table_spt[i]=float(stype[i][1])
    
    for i in range(len(teff)):
        nearest=(np.abs(table_spt-teff[i])).argmin()
        spt[i]=stype[nearest][0]
    
    return spt
   
OBASort()
#makeAvGaia(BPRP=0.52925,A0=2.3)
