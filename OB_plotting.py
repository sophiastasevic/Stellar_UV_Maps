# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 12:39:22 2020

@author: Sophia Stasevic

Relations from Parravano et al. (2003)

Data needed for calculating FUV/EUV/MS lifetime:
    -'ra' 
    -'dec'
    -'mass' (M_Sun)
    -'spt' (spectral type)
Text file output from SpType_to_Mass.py will need to be edited so that the 
first row is: '#ra   dec   teff   mass   spt'

If using SaveOutput(), need to edit output so that the first row is:
'#ra   dec   FUV   EUV   MS_lifetime'
    -FUV and EUV in L_Sun
    -MS lifetime in Myr
    
"""
"""
==========================    PACKAGES + IMPORTS    ================================
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table 
import astropy.coordinates as ac
import math
import seaborn as sns

import matplotlib as mpl
import matplotlib.ticker as ticker
import matplotlib.cm as cm

import SpinRatesAreHere as spin 
import matplotlib.patheffects as fx

"""
============================	GLOBAL VARIABLES	=================================
"""
name_input='NGC2264' #for ease of plot titles and saving
cluster = spin.NGC2264()

disked=np.where(cluster.Disk>0)[0]
diskless=np.where(cluster.Disk<0)[0]
unknown=np.where(cluster.Disk==0)[0]

#file containing at least ra, dec, spectral type, and mass of each massive star
data=Table.read('{name}_OBA_full_data_UV.txt'.format(name=name_input), format='ascii') 

letter=np.array(data['spt'])
number=np.empty(len(data['spt']))

for i in range(len(letter)):
    number[i]=int(letter[i][1])
    letter[i]=letter[i][0]
    
O_type=np.where(letter=='O')[0]
B_type=np.where(letter=='B')[0]
A_type=np.where(letter=='A')[0]
early=np.where(number<=5)[0]
late=np.where(number>5)[0]

O_stars=np.array(list(set(O_type)))
B_stars_early=np.array(list(set(B_type)&set(early)))
B_stars_late=np.array(list(set(B_type)&set(late)))
OB_stars=[*O_stars,*B_stars_early]
A_stars=np.array(list(set(A_type)))

#chi=np.where(data['ra']>35.25) #ra split between hPer and chiPer

L_Sun = 3.486e26 #(W/m^2)
parsec = 3.086e+16 #(m)

"""
============================	FUNCTIONS	=================================
"""

def CoordMap():

    plt.scatter(cluster.RA[disked], cluster.Dec[disked], color='k', s=4, alpha=0.6, label='Disk')
    plt.scatter(cluster.RA[diskless], cluster.Dec[diskless], color='r', s=4, alpha=0.6, label='Diskless')

    plt.scatter(data['ra'][OB_stars],data['dec'][OB_stars], color='k', s=20, marker='*', label='OB stars',zorder=10)

    #plt.legend(loc='lower right')    
    #plt.xlabel('Right Ascension [deg]')
    #plt.ylabel('Declination [deg]')
    #plt.title('RA-Dec Map of {name}'.format(name=name_input))
    

def MS_lifetime(mass_list):
    
    length=len(mass_list)
    
    lifetime=np.empty(length)
    for i in range(0,length):
        if mass_list[i]>=1.2 and mass_list[i]<3:
            lifetime[i]=(7.65*pow(10,3))*pow(mass_list[i],(-2.8))
        elif mass_list[i]>=3 and mass_list[i]<6:
            lifetime[i]=(4.73*pow(10,3))*pow(mass_list[i],(-2.36))
        elif mass_list[i]>=6 and mass_list[i]<9:
            lifetime[i]=(2.76*pow(10,3))*pow(mass_list[i],(-2.06))
        elif mass_list[i]>=9 and mass_list[i]<12:
            lifetime[i]=(1.59*pow(10,3))*pow(mass_list[i],(-1.81))
        elif mass_list[i]>=12 and mass_list[i]<120:
            lifetime[i]=(7.6*pow(10,3))*pow(mass_list[i],(-1.57))+2.3
        else:
            lifetime[i]=0
    
    return lifetime

def EUV_luminosity(mass_list):
    
    length=len(mass_list)
    
    EUV_lum=np.empty(length) 
    for i in range(0,length):
        if mass_list[i]>=5 and mass_list[i]<7:
            EUV_lum[i]=(2.32*pow(10,34))*pow(mass_list[i],(11.5))
        elif mass_list[i]>=7 and mass_list[i]<12:
            EUV_lum[i]=(3.69*pow(10,36))*pow(mass_list[i],(8.87))
        elif mass_list[i]>=12 and mass_list[i]<20:
            EUV_lum[i]=(4.8*pow(10,38))*pow(mass_list[i],(7.85))
        elif mass_list[i]>=20 and mass_list[i]<30:
            EUV_lum[i]=(3.12*pow(10,41))*pow(mass_list[i],(4.91))
        elif mass_list[i]>=30 and mass_list[i]<40:
            EUV_lum[i]=(2.8*pow(10,44))*pow(mass_list[i],(2.91))
        elif mass_list[i]>=40 and mass_list[i]<60:
            EUV_lum[i]=(3.49*pow(10,45))*pow(mass_list[i],(2.23))
        elif mass_list[i]>=60 and mass_list[i]<120:
            EUV_lum[i]=(2.39*pow(10,46))*pow(mass_list[i],(1.76))
        else:
            EUV_lum[i]=0
            
    return EUV_lum

def FUV_luminosity(mass_list):
    
    length=len(mass_list)
    
    H2_lum=np.empty(length)
    for i in range(0,length):
        mass_list[i]=float(mass_list[i])
        if mass_list[i]>=1.8 and mass_list[i]<3:
            H2_lum[i]=(1.98*pow(10,-14))*pow(mass_list[i],(26.6))
        elif mass_list[i]>=3 and mass_list[i]<4:
            H2_lum[i]=(2.86*pow(10,-8))*pow(mass_list[i],(13.7))
        elif mass_list[i]>=4 and mass_list[i]<6:
            H2_lum[i]=(1.35*pow(10,-4))*pow(mass_list[i],(7.61))
        elif mass_list[i]>=6 and mass_list[i]<9:
            H2_lum[i]=(1.1*pow(10,-2))*pow(mass_list[i],(5.13))
        elif mass_list[i]>=9 and mass_list[i]<12:
            H2_lum[i]=(1.07*pow(10,-1))*pow(mass_list[i],(4.09))
        elif mass_list[i]>=12 and mass_list[i]<15:
            H2_lum[i]=(5.47*pow(10,-1))*pow(mass_list[i],(3.43))
        elif mass_list[i]>=15 and mass_list[i]<30:
            H2_lum[i]=(9.07*pow(10,0))*pow(mass_list[i],(2.39))
        elif mass_list[i]>=30 and mass_list[i]<120:
            H2_lum[i]=(9.91*pow(10,1))*pow(mass_list[i],(1.69))
        else:
            H2_lum[i]=0
    
    FUV_H2_lum=np.empty(length) #luminosity FUV-H2
    for i in range(0,length):
        mass_list[i]=float(mass_list[i])
        if mass_list[i]>=1.8 and mass_list[i]<2:
            FUV_H2_lum[i]=(2.77*pow(10,-4))*pow(mass_list[i],(11.8))
        elif mass_list[i]>=2 and mass_list[i]<2.5:
            FUV_H2_lum[i]=(1.88*pow(10,-3))*pow(mass_list[i],(9.03))
        elif mass_list[i]>=2.5 and mass_list[i]<3:
            FUV_H2_lum[i]=(1.19*pow(10,-2))*pow(mass_list[i],(7.03))
        elif mass_list[i]>=3 and mass_list[i]<6:
            FUV_H2_lum[i]=(1.47*pow(10,-1))*pow(mass_list[i],(4.76))
        elif mass_list[i]>=6 and mass_list[i]<9:
            FUV_H2_lum[i]=(8.22*pow(10,-1))*pow(mass_list[i],(3.78))
        elif mass_list[i]>=9 and mass_list[i]<12:
            FUV_H2_lum[i]=(2.29*pow(10,0))*pow(mass_list[i],(3.31))
        elif mass_list[i]>=12 and mass_list[i]<30:
            FUV_H2_lum[i]=(2.7*pow(10,1))*pow(mass_list[i],(2.32))
        elif mass_list[i]>=30 and mass_list[i]<120:
            FUV_H2_lum[i]=(3.99*pow(10,2))*pow(mass_list[i],(1.54))
        else:
            FUV_H2_lum[i]=0

    return FUV_H2_lum+H2_lum


def IncidentFlux():

    #coords of massvie stars
    cluster.ClusterInfo()
    RA=data['ra']#[B_stars_late]
    DE=data['dec']#[B_stars_late]
    dist=float(cluster.dist)

    #cat. used for NGC2362 to convert to degrees
    #cat = ac.SkyCoord(cluster.RA, cluster.Dec, unit="deg")
    
    #coords of low mass stars 
    RA_ref=cluster.RA #cat.ra.deg*15 #
    DE_ref=cluster.Dec #cat.dec.deg #
    radius= Distance(RA_ref, DE_ref, RA, DE, dist) #for flux at the point of low mass stars; only needed to calculate flux for the first time
    
    step_RA=(max(RA_ref)-min(RA_ref))/50
    step_DE=(max(DE_ref)-min(DE_ref))/50
    RA_bin=np.arange(min(RA_ref),max(RA_ref)+step_RA,step_RA).tolist()
    DE_bin=np.arange(min(DE_ref),max(DE_ref)+step_DE,step_DE).tolist()
    coord=CoordGrid(RA_bin,DE_bin,RA,DE,dist) #for the UV map
    
    #luminosity=data['FUV']#[B_stars_late] #if FUV is already known
    
    #for calculating flux from mass:
    mass_list= MassEstimate(data)
    luminosity= FUV_luminosity(mass_list)
    
    #for calculating flux of low mass stars 
    flux=np.empty((len(RA_ref),len(RA)))
    for i in range(0,len(RA_ref)):
        for n in range(0,len(RA)):
            flux[i][n]= luminosity[n]/(4*math.pi*pow(radius[i][n],2))
    
    #coordinate grid for massive star plots
    flux_coord=np.empty((len(RA),len(DE_bin),len(RA_bin)))
    for i in range(0,len(RA)):
        for n in range(0,len(DE_bin)):
            for m in range(0,len(RA_bin)):
                flux_coord[i][n][m]=luminosity[i]/(4*math.pi*pow(coord[i][n][m]+step_RA/10,2)) #buffer to avoid going to infinity when r=0

    #LowMassUV(RA_ref,DE_ref) #if FUV of low mass stars is known
    #LowMassUV(RA_ref,DE_ref,flux,radius) #if FUV is not known
    UVMap(RA_bin,DE_bin,flux_coord,coord,step_RA,step_DE)

#calculates seperation for every low mass star from each massive star
def Distance(RA_ref, DE_ref, RA, DE, dist):
    
    ang_sep=np.empty((len(RA_ref),len(RA)))
    radius=np.empty((len(RA_ref),len(RA)))
    RA_dif=np.empty((len(RA_ref),len(RA)))
    DE_dif=np.empty((len(RA_ref),len(RA)))  
    
    for i in range(0,len(RA_ref)):
        for n in range(0,len(RA)):
            RA_dif[i][n]=RA_ref[i]-RA[n]
            DE_dif[i][n]=DE_ref[i]-DE[n]
        for n in range(0,len(RA)):
            ang_sep[i][n]=math.sqrt(pow(RA_dif[i][n],2)+pow(DE_dif[i][n],2))
            radius[i][n]=math.tan(math.radians(ang_sep[i][n]/2))*dist*2 
    
    return radius

   
def MassEstimate(data):
    
    mass=data['mass']#[B_stars_late] 
    return mass

#creates a 2D array of periodic RA and DE seperations from each massive star
def CoordGrid(RA_bin,DE_bin,RA,DE,dist):

    coord=np.empty((len(RA),len(DE_bin),len(RA_bin)))
    RA_dif=np.empty((len(RA),len(RA_bin)))
    DE_dif=np.empty((len(DE),len(DE_bin)))
    for i in range(0,len(RA)):
        for n in range(0,len(RA_bin)):
            RA_dif[i][n]=RA[i]-RA_bin[n]
        for m in range(0,len(DE_bin)):
            DE_dif[i][m]=DE[i]-DE_bin[m]
        for n in range(0,len(DE_bin)):
            for m in range(0,len(RA_bin)):
                coord[i][n][m]=math.sqrt(pow(RA_dif[i][m],2)+pow(DE_dif[i][n],2))
                coord[i][n][m]=math.tan(math.radians(coord[i][n][m]/2))*dist*2
           
    return coord
    
def UVMap(RA_bin,DE_bin,flux_coord,coord,step_RA,step_DE):
    
    #for overplotting, need to comment out lines 297, 303, and 304
    
    #ClusterMaps() #for overplotting contours on density/ disk fraciton maps

    for i in range(len(RA_bin)): #so that ticks lie in the centre of each data pixel
        RA_bin[i]=RA_bin[i]-(step_RA/2)
    for i in range(len(DE_bin)):
        DE_bin[i]=DE_bin[i]-(step_DE/2)
    x=RA_bin
    y=DE_bin
    
    flux_total=flux_coord[0]
    for i in range(1,len(flux_coord)):
        flux_total=flux_total+flux_coord[i]
    z=np.log10(flux_total/3.98) #converson from Lo/pc^2 to G0
    
    pc= plt.pcolor(x,y,z,cmap = plt.cm.Spectral_r,norm=mpl.colors.Normalize(vmin=np.log10(1.7),vmax=5)) #1.7 is G0 in solar neighbourhood
    
    cmax=np.amax(z)
    print(cmax)
    levels=np.log10([50,100,500,1000,3000,30000])
    contour = plt.contour(z,levels,origin='lower', linewidths=1.5,colors=['b','lime','yellow','orange','r','purple'],extent=(min(x),max(x)+step_RA,min(y),max(y)+step_DE))
    cb=plt.colorbar(pc, orientation="vertical", pad=0.01,aspect=15)
    cb.ax.set_ylabel('FUV flux [log G$_0$]', rotation=270,linespacing=5,fontsize=10,labelpad=20)
    
    #fmtd=DistanceLabels(RA_bin,z,coord,levels) #plots the radius from the star at which the flux level is --> only possible with single stars
    G_0=[50,100,500,1000,3000,30000]
    fmt={}
    string="G$_{0}$"
    strs=["{}{}".format(i,string) for i in G_0]
    for l, s in zip(levels, strs):
        fmt[l] = s
        
    #plt.clabel(contour,levels,inline=True, manual=False, colors = 'k', fmt=fmtd, fontsize=10) #distance labels 
    plt.clabel(contour,levels,inline=True, manual=False, colors='k' , fmt=fmt, fontsize=10)
    plt.xlabel('Right Ascension [deg]')
    plt.ylabel('Declination [deg]')
    #CoordMap() #plots stars on top of UV map
    plt.axis([max(x),min(x),min(y),max(y)])
    plt.title('Incident FUV flux map of {name}'.format(name=name_input))
    plt.savefig('{name}_low_mass_FUV_Map_.png'.format(name=name_input))

def LowMassUV(RA_ref,DE_ref,flux,radius): #(RA_ref,DE_ref) if FUV is known
    
    flux_total=flux[:,0]
    for i in range(len(flux_total)):
        for n in range(len(flux[0])):
            flux_total[i]=flux_total[i]+flux[i][n]
    
    z=np.log10(flux_total/3.98) #conversion from Lo/pc^2 to Go
    
    #z=cluster.OB_FUV-np.log10(3.98) #if FUV is known
    x=RA_ref
    y=DE_ref
    
    final=np.column_stack((x,y,z))
    np.savetxt('{name}_low_mass_fuv.txt'.format(name=name_input),final, delimiter=" ", fmt="%s")
    
    def CoordPlot(x,y,z):
        a=plt.scatter(x,y,c=z,alpha=0.6,cmap = plt.cm.Spectral_r,norm=mpl.colors.Normalize(vmin=np.log10(1.7),vmax=5))
        cb=plt.colorbar(a, orientation="vertical", pad=0.01,aspect=15)
        cb.ax.set_ylabel('FUV flux [log G$_0$]', rotation=270,linespacing=5,fontsize=10,labelpad=20)
        
        plt.xlabel('Right Ascension [deg]')
        plt.ylabel('Declination [deg]')
        plt.axis([max(x),min(x),min(y),max(y)])
        plt.title('Incident FUV flux of {name} on Low Mass Stars'.format(name=name_input))
        
    CoordPlot(x,y,z)
    
def DistanceLabels(RA_bin,z,coord,levels):
   
    length=int(len(RA_bin)/2)
    sort=np.argsort(z[length])
    radius=np.interp(levels,z[length][sort],coord[0][length][sort])
    
    from decimal import Decimal
    for i in range(len(radius)):    
        radius[i]=Decimal(radius[i])
        radius[i]=round(radius[i],1)
    fmt = {}
    string="pc"
    strs=["{}{}".format(i,string) for i in radius]
    for l, s in zip(levels, strs):
        fmt[l] = s
    
    return fmt
   
#cat. used for NGC2362 as coords in table aren't in degrees    
def ClusterMaps():
    
    fig,ax = plt.subplots()
    #cat = ac.SkyCoord(cluster.RA, cluster.Dec, unit="deg")
    #cat_disked = ac.SkyCoord(cluster.RA[disked], cluster.Dec[disked], unit="deg")

    x = cluster.RA #cat.ra.deg*15
    y = cluster.Dec #cat.dec.deg
    x_disk = cluster.RA[disked] #cat_disked.ra.deg*15
    y_disk = cluster.Dec[disked] #cat_disked.dec.deg
    
    y_range=max(y)-min(y)
    x_range=max(x)-min(x)
    y_step=y_range/10
    x_step=x_range/10
    y_bin=np.arange(min(y),max(y),y_step).tolist()
    x_bin=np.arange(min(x),max(x),x_step).tolist()
    
    def DensityMap(x,y,x_bin,y_bin):
        
        plt.hist2d(x,y,(x_bin,y_bin),cmap=plt.cm.PuBuGn)
        plt.colorbar(pad=0.01, aspect=15, label='Number of Stars')
        plt.xlabel('Right Ascension [deg]')
        plt.ylabel('Declination [deg]')
        plt.title('Density Map of ')
        ax.invert_xaxis()
    
    def DiskFractMap(x,y,x_bin,y_bin):
        
        a=plt.hist2d(x,y,(x_bin,y_bin),cmap=plt.cm.PuBuGn)
        b=plt.hist2d(x_disk,y_disk,(x_bin,y_bin),cmap=plt.cm.PuBuGn)
        c=100*(b[0]/a[0])
        ra=a[1]+x_step/2
        dec=a[2]+y_step/2
        
        ra.round(decimals=2)
        
        dfract= plt.pcolor(ra,dec,c,cmap = plt.cm.PuBuGn,norm=mpl.colors.Normalize(vmin=0,vmax=50))
        
        plt.xlabel('Right Ascension [deg]')
        plt.ylabel('Declination [deg]')
        plt.title('Disk Fraction Map of (x < 5 stars)')
        ax.invert_xaxis()
        plt.axis([max(ra),min(ra),min(dec),max(dec)])
        
        cb2=plt.colorbar(dfract, orientation="vertical", pad=0.01,aspect=15)
        cb2.ax.set_ylabel('Disk Fraction [%]', rotation=90,linespacing=5,fontsize=10,labelpad=10)
        
        #label bins with x if number of stars is below a certain value
        for n in range(0,len(ra)-1):
            for m in range(0,len(dec)-1):
                if a[0][m,n]<5:
                   ax.text(ra[n]+x_step/1.85,dec[m]+y_step/2.3,'x')       
                ax.text(ra[n]+x_step/1.5,dec[m]+y_step/2.5,'%.0f' %a[0][m,n])
        
        #label bins with number of stars in them
        #for y in range(0,c.shape[0]):
         #   for x in range(c.shape[1]):
          #      ax.text(x + 0.5, y + 0.5, '%.4f' % c[y, x],
           #          horizontalalignment='center',
            #         verticalalignment='center',
             #        )
    
    DensityMap(x,y,x_bin,y_bin)
    #DiskFractMap(x,y,x_bin,y_bin)
    

#for saving FUV/EUV/MS lifetime to a file
def SaveOutput():
     
    mass_list= MassEstimate(data)
    FUV = FUV_luminosity(mass_list)
    EUV = EUV_luminosity(mass_list)
    lt = MS_lifetime(mass_list)
    
    final=np.column_stack((np.array(data['ra']),np.array(data['dec']),FUV,EUV,lt))
    np.savetxt('{name}_uvs.txt'.format(name=name_input),final, delimiter=" ", fmt="%s")
    

"""
================================================================================
"""
IncidentFlux()   
#ClusterMaps()
#CoordMap()
SaveOutput()      