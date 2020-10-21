"""
==========================    PACKAGES + IMPORTS    ================================
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table 
import astropy.coordinates as ac
import math

import matplotlib as mpl
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import matplotlib.patheffects as fx

"""
============================	GLOBAL VARIABLES	=================================
"""

#file containing at least ra, dec, spectral type, and mass of each massive star
data=Table.read('single_star_data.txt', format='ascii') 
spt_list=['O3','O5','O7','B1','B3','B5','O9']

L_Sun = 3.486e26 #(W/m^2)
parsec = 3.086e+16 #(m)

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

    #position of massive stars
    RA=data['ra']
    DE=data['dec']
    dist=100
    
    #range for the coordinate grid
    RA_ref=[5,15]
    DE_ref=[5,15]
    
    step_RA=(max(RA_ref)-min(RA_ref))/50
    step_DE=(max(DE_ref)-min(DE_ref))/50
    RA_bin=np.arange(min(RA_ref),max(RA_ref)+step_RA,step_RA).tolist()
    DE_bin=np.arange(min(DE_ref),max(DE_ref)+step_DE,step_DE).tolist()
    coord=CoordGrid(RA_bin,DE_bin,RA,DE,dist) #for the UV map
    
    mass_list= data['mass']
    luminosity= FUV_luminosity(mass_list)
    
    #coordinate grid for massive star plots
    flux_coord=np.empty((len(RA),len(DE_bin),len(RA_bin)))
    for i in range(0,len(RA)):
        for n in range(0,len(DE_bin)):
            for m in range(0,len(RA_bin)):
                flux_coord[i][n][m]=luminosity[i]/(4*math.pi*pow(coord[i][n][m]+step_RA/10,2)) #buffer to avoid going to infinity when r=0

    UVMap(RA_bin,DE_bin,flux_coord,coord,step_RA,step_DE)

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
    
    for i in range(len(RA_bin)): #so that ticks lie in the centre of each data pixel
        RA_bin[i]=RA_bin[i]-(step_RA/2)
    for i in range(len(DE_bin)):
        DE_bin[i]=DE_bin[i]-(step_DE/2)
    x=RA_bin
    y=DE_bin
    
    fig, axes = plt.subplots(3,2, figsize=(10,11))
    axes_list=[axes[0,0],axes[0,1],axes[1,0],axes[1,1],axes[2,0],axes[2,1]]
    for i in range(6):
        ax=axes_list[i]
        plt.subplot(ax)
        z=np.log10(flux_coord[i]/3.98) #converson from Lo/pc^2 to G0

        pc= plt.pcolor(x,y,z,cmap = plt.cm.Spectral_r,norm=mpl.colors.Normalize(vmin=np.log10(1.7),vmax=5)) #1.7 is G0 in solar neighbourhood
        
        levels=np.log10([50,100,500,1000,3000,30000])
        contour = plt.contour(z,levels,origin='lower', linewidths=1.5,colors=['b','lime','yellow','orange','r','purple'],extent=(min(x),max(x)+step_RA,min(y),max(y)+step_DE))
        cb=plt.colorbar(pc, orientation="vertical", pad=0.01,aspect=15)
        cb.ax.set_ylabel('FUV flux [log G$_0$]', rotation=270,linespacing=5,fontsize=14,labelpad=20)
        cb.ax.tick_params(labelsize=12)
        G_0=[50,100,500,1000,3000,30000]
        fmt={}
        string="G$_{0}$"
        strs=["{}{}".format(i,string) for i in G_0]
        for l, s in zip(levels, strs):
            fmt[l] = s
        
        plt.clabel(contour,levels,inline=True, manual=False, colors='k' , fmt=fmt, fontsize=12)
        plt.xlabel('Right Ascension [deg]',fontsize=14)
        plt.ylabel('Declination [deg]',fontsize=14)
        #CoordMap() #plots stars on top of UV map
        plt.axis([max(x),min(x),min(y),max(y)])
        plt.tick_params(axis='both', labelsize=12)
        plt.title('Incident FUV flux map of {name} star at 100pc'.format(name=spt_list[i]),fontsize=14)
        fig.tight_layout(pad=0.9)
    plt.savefig('sing_star_FUV_100pc.png'.format(name=spt_list[i]))

        
IncidentFlux()