# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 14:38:46 2020

Function for calculating rolling disk fraction
"""

import matplotlib.pyplot as plt
import numpy as np

import warnings
warnings.simplefilter("ignore")

import SpinRatesAreHere as spin 

name_input = 'USco'

if name_input =='NGC2264':
    cluster = spin.NGC2264()

    
elif name_input == 'USco':
    cluster = spin.USco()
    
elif name_input == 'NGC2362':
    cluster = spin.NGC2362()
    
elif name_input == 'NGC6530':
    cluster = spin.NGC6530()
    
else:
    print('Invalid name')
    exit()
    
mass=np.where(np.isfinite(cluster.Mass))[0]
prot=np.where(np.isfinite(cluster.Prot))[0]
disk=np.where(np.isfinite(cluster.Disk))[0]
has_disk_mass_prot=np.array(list(set(mass)&set(disk)&set(prot))) #only want stars with mass, period, and disk information

def RollingDiskFract(cluster):
    
    prot_sort=cluster.Prot[has_disk_mass_prot]
    disk_sort=cluster.Disk[has_disk_mass_prot]
    sort=np.argsort(prot_sort)
    prot_sort=prot_sort[sort] #sorts array in order of increasing period
    disk_sort=disk_sort[sort] #sorts disk information to be in the same order as the period sort
    data_length=int(len(prot_sort))
    window=int(data_length*0.25) #number of stars disk fraction will be calculated from at each point
    length=data_length-(window+1)
    
    y=np.empty(length)
    x=np.empty(length)
    for i in range(0,(length)):
        n_disk=0
        n_diskless=0
        n_unknown=0
        for n in range(i,(i+window)):
            if disk_sort[n]>0:
                n_disk=n_disk+1
            elif disk_sort[n]<0:
                n_diskless=n_diskless+1
            else:
                n_unknown=n_unknown+1
        x[i]=prot_sort[i]
        y[i]=100*n_disk/(n_disk+n_diskless+n_unknown)
   
    #a lot of the slower rotators are cut off due to the window size, added an extra section to try and keep more data in the plot
    window_slow=int(data_length*0.15) #must be lower than window
    s=np.empty(window-(window_slow+1))
    t=np.empty(window-(window_slow+1))
    

    for p in range(0,window-(window_slow+1)): 
        q_disk=0
        q_diskless=0
        q_unknown=0
        for q in range(p,(p+window_slow)):
            if disk_sort[length+q]>0:
                q_disk=q_disk+1
            elif disk_sort[length+q]<0:
                q_diskless=q_diskless+1
            else:
                q_unknown=q_unknown+1
        s[p]=prot_sort[length+p]
        t[p]=100*q_disk/(q_disk+q_diskless+q_unknown)
        
    w=[*x,*s] #joins the two sets of data together so that the curves can be overplotted to fix discontinuity between curves
    z=[*y,*t]

    plt.ylabel('Disk Fraction [%]',fontsize=14)
    plt.xlabel('Period [days]',fontsize=14)
    plt.title('Rolling Disk fraction-Period plot for {name}'.format(name=name_input),fontsize=14)
    plt.plot(w,z,color='b',label='window size = 25%')
    plt.plot(s,t,color='r', label= 'window size = 15%')
    plt.legend(loc='lower right',fontsize=14)
    plt.tick_params( axis='both',labelsize=12)
    
    plt.savefig('{name}_rolling_disk_fract_prot.png'.format(name=name_input))
    plt.show()
    
RollingDiskFract(cluster)