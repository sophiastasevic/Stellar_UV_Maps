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
name_input='test' #for ease of plot titles and saving

#file containing position and uv data 
data_table=Table.read('tau_d_fuv_Cluster_Kroupa_288Msun_r_eff_2_1_plummer_8OB.fits', format='fits') 

#table contains both high (OB) mass stars and very low mass stars --> want to separate the massive stars for the full map
massive=np.where(np.isfinite(data_table['L_FUV']))[0] 
data=data_table[massive]
cluster=data_table

L_Sun = 3.486e26 #(W/m^2)
parsec = 3.086e+16 #(m)
erg = 1e-7 #(J)

Calc_emitted_FUV = False #set to False if the FUV of each massive star is already calculated
Calc_incident_FUV = False #set to False if FUV flux at position of each low mass star is already calculated

"""
============================	FUNCTIONS	=================================
"""

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
    x=data['x_pc']
    y=data['y_pc']    

    #coords of low mass stars 
    x_ref=cluster['x_pc']
    y_ref=cluster['y_pc']
        
    step_x=(max(x_ref)-min(x_ref))/100
    step_y=(max(y_ref)-min(y_ref))/100
    x_bin=np.arange(min(x_ref),max(x_ref)+step_x/2,step_x).tolist()
    y_bin=np.arange(min(y_ref),max(y_ref)+step_y,step_y).tolist()
        
    coord=CoordGrid(x_bin,y_bin,x,y) #creates grid of x and Dec in the range of the low mass stars
        
    if Calc_emitted_FUV == True:
        mass_list= MassEstimate(data)
        luminosity= FUV_luminosity(mass_list)
    else:
        luminosity=data['L_FUV']/L_Sun*erg #conversion from cgs to Lo
        
    #calculates incident FUV flux from the massive stars at each point in the coordinate grid
    flux_coord=np.empty((len(x),len(y_bin),len(x_bin)))
    for i in range(0,len(x)):
        for n in range(0,len(y_bin)):
            for m in range(0,len(x_bin)):
                flux_coord[i][n][m]=luminosity[i]/(4*math.pi*pow(coord[i][n][m]+step_x/10,2)) #buffer to avoid going to infinity when r=0

    contour_z=UVMap(x_bin,y_bin,flux_coord,coord,step_x,step_y)
    
    if Calc_incident_FUV == True:
        radius= Distance(x_ref, y_ref, x, y)
        flux=np.empty((len(x_ref),len(x)))
        for i in range(0,len(x_ref)):
            for n in range(0,len(x)):
                flux[i][n]= luminosity[n]/(4*math.pi*pow(radius[i][n],2))
        LowMassUV2(x_ref,y_ref,flux,radius,contour_z)
    else:
        LowMassUV(x_ref,y_ref,contour_z)

#calculates seperation for every low mass star from each massive star
def Distance(x_ref, y_ref, x, y):
    
    radius=np.empty((len(x_ref),len(x)))
    x_dif=np.empty((len(x_ref),len(x)))
    y_dif=np.empty((len(x_ref),len(x)))  
    
    for i in range(0,len(x_ref)):
        for n in range(0,len(x)):
            x_dif[i][n]=x_ref[i]-x[n]
            y_dif[i][n]=y_ref[i]-y[n]
        for n in range(0,len(x)):
            radius[i][n]=math.sqrt(pow(x_dif[i][n],2)+pow(y_dif[i][n],2))
    
    return radius
   
def MassEstimate(data):
    
    mass=data['Mass_Msun']
    return mass

#creates a 2D array of coordinates for binning of the full FUV flux map of massive stars
def CoordGrid(x_bin,y_bin,x,y):

    coord=np.empty((len(x),len(y_bin),len(x_bin)))
    x_dif=np.empty((len(x),len(x_bin)))
    y_dif=np.empty((len(y),len(y_bin)))
    for i in range(0,len(x)):
        for n in range(0,len(x_bin)):
            x_dif[i][n]=x[i]-x_bin[n]
        for m in range(0,len(y_bin)):
            y_dif[i][m]=y[i]-y_bin[m]
        for n in range(0,len(y_bin)):
            for m in range(0,len(x_bin)):
                coord[i][n][m]=math.sqrt(pow(x_dif[i][m],2)+pow(y_dif[i][n],2))
                
    return coord
    
def UVMap(x_bin,y_bin,flux_coord,coord,step_x,step_y):
    
    x=x_bin
    y=y_bin
    
    flux_total=flux_coord[0]
    for i in range(1,len(flux_coord)):
        flux_total=flux_total+flux_coord[i]
    z=np.log10(flux_total/3.98) #converson from Lo/pc^2 to G0
    
    pc= plt.pcolor(x,y,z,cmap = plt.cm.Spectral_r,norm=mpl.colors.Normalize(vmin=np.log10(1.7),vmax=5)) #1.7 is G0 in solar neighbourhood
    
    levels=np.log10([50,100,500,1000,3000,30000])
    contour = plt.contour(z,levels,origin='lower', linewidths=1.5,colors=['b','lime','yellow','orange','r','purple'],extent=(min(x),max(x),min(y),max(y)))
    cb=plt.colorbar(pc, orientation="vertical", pad=0.01,aspect=15)
    cb.ax.set_ylabel('FUV flux [log G$_0$]', rotation=270,linespacing=5,fontsize=10,labelpad=20)
    
    #fmtd=DistanceLabels(x_bin,z,coord,levels) #plots the radius from the star at which the flux level is --> only possible with single stars
    G_0=[50,100,500,1000,3000,30000]
    fmt={}
    string="G$_{0}$"
    strs=["{}{}".format(i,string) for i in G_0]
    for l, s in zip(levels, strs):
        fmt[l] = s
        
    #plt.clabel(contour,levels,inline=True, manual=False, colors = 'k', fmt=fmtd, fontsize=10) #distance labels 
    plt.clabel(contour,levels,inline=True, manual=False, colors='k' , fmt=fmt, fontsize=10)
    plt.xlabel('x [pc]')
    plt.ylabel('y [pc]')
    #CoordMap() #plots massive stars on top of UV map
    plt.axis([max(x),min(x),min(y),max(y)])
    plt.title('2D Incident FUV flux map of {name}'.format(name=name_input))
    plt.show()
    #plt.savefig('{name}_Full_FUV_Map.png'.format(name=name_input))
    
    return z

def LowMassUV(x_ref,y_ref,contour_z):
    
    z=cluster['F_FUV_2d']
    x=x_ref
    y=y_ref
    
    LowMassPlot(x,y,z,contour_z)
    
def LowMassUV2(x_ref,y_ref,flux,radius,contour_z):
    
    flux_total=flux[:,0]
    for i in range(len(flux_total)):
        for n in range(len(flux[0])):
            flux_total[i]=flux_total[i]+flux[i][n]
    
    z=np.log10(flux_total/3.98) #conversion from Lo/pc^2 to Go
    
    x=x_ref
    y=y_ref
    
    final=np.column_stack((x,y,np.log10(flux_total)))
    np.savetxt('{name}_low_mass_fuv.txt'.format(name=name_input),final, delimiter=" ", fmt="%s")
    
    LowMassPlot(x,y,z,contour_z)
    
def LowMassPlot(x,y,z,contour_z):
    
    a=plt.scatter(x,y,c=z,alpha=0.6,cmap = plt.cm.Spectral_r,norm=mpl.colors.Normalize(vmin=np.log10(1.7),vmax=5))
    cb=plt.colorbar(a, orientation="vertical", pad=0.01,aspect=15)
    cb.ax.set_ylabel('FUV flux [log G$_0$]', rotation=270,linespacing=5,fontsize=10,labelpad=20)
     
    levels=np.log10([50,100,500,1000,3000,30000])
    contour = plt.contour(contour_z,levels,origin='lower', linewidths=1.5,colors=['b','lime','yellow','orange','r','purple'],extent=(min(x),max(x),min(y),max(y)))
    
    #fmtd=DistanceLabels(x_bin,z,coord,levels) #plots the radius from the star at which the flux level is --> only possible with single stars
    G_0=[50,100,500,1000,3000,30000]
    fmt={}
    string="G$_{0}$"
    strs=["{}{}".format(i,string) for i in G_0]
    for l, s in zip(levels, strs):
        fmt[l] = s
        
    #plt.clabel(contour,levels,inline=True, manual=False, colors = 'k', fmt=fmtd, fontsize=10) #distance labels 
    plt.clabel(contour,levels,inline=True, manual=False, colors='k' , fmt=fmt, fontsize=10)
    
    plt.xlabel('x [pc]')
    plt.ylabel('y [pc]')
    plt.axis([max(x),min(x),min(y),max(y)])
    plt.title('Incident FUV flux of {name} at Low Mass Stars'.format(name=name_input))
    plt.show()
    #plt.savefig('{name}_Low_Mass_Star_FUV_Map.png'.format(name=name_input))
      
def DistanceLabels(x_bin,z,coord,levels):
   
    length=int(len(x_bin)/2)
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
#CoordMap()
#SaveOutput()      