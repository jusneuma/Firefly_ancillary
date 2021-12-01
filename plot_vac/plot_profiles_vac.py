#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  5 09:58:00 2020

@author: Justus Neumann jusneuma.astro@gmail.com
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib.gridspec as gridspec
import sys, os

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
mpl.rcParams['xtick.labelsize']=12
mpl.rcParams['ytick.labelsize']=12
mpl.rcParams['font.size']=14.0
mpl.rcParams['xtick.major.size']=8
mpl.rcParams['xtick.minor.size']=4
mpl.rcParams['ytick.major.size']=8
mpl.rcParams['ytick.minor.size']=4
mpl.rcParams['xtick.major.width']=1.2
mpl.rcParams['xtick.minor.width']=1
mpl.rcParams['ytick.major.width']=1.2
mpl.rcParams['ytick.minor.width']=1
mpl.rcParams['xtick.minor.visible']=True
mpl.rcParams['ytick.minor.visible']=True
mpl.rcParams['xtick.direction']='in'
mpl.rcParams['ytick.direction']='in'
mpl.rcParams['ytick.right']=True
mpl.rcParams['xtick.top']=True
mpl.rcParams['axes.labelsize']='medium'


loc_vac = os.getcwd()

# read the vac file
vac = fits.open(os.path.join(loc_vac,'manga-firefly-v3_1_1-miles.fits'))  

props=['LW_AGE_VORONOI','MW_AGE_VORONOI','LW_Z_VORONOI','MW_Z_VORONOI','SURFACE_MASS_DENSITY_VORONOI','E(B_V)_VORONOI']

dict = {'LW_AGE_VORONOI':r'$\boldsymbol{\mathrm{log (Age_{\ LW}\ [Gyr])}}$',\
        'MW_AGE_VORONOI':r'$\boldsymbol{\mathrm{log (Age_{\ MW}\ [Gyr])}}$',\
        'LW_Z_VORONOI':r'$\boldsymbol{\mathrm{[Z/H]_{\ LW}}}$',\
        'MW_Z_VORONOI':r'$\boldsymbol{\mathrm{[Z/H]_{\ MW}}}$',\
        'SURFACE_MASS_DENSITY_VORONOI':r'$\boldsymbol{\mathrm{log\ \Sigma_\star\ [M_\odot\,kpc^{-1}]}}$',\
        'E(B_V)_VORONOI':r'$\boldsymbol{\mathrm{E_\text{B-V}}}$'}

def plot_prop(prop,ax,plate,ifu):
    basic = vac[1].data
    galid = basic['PLATEIFU']==str(plate)+'-'+str(ifu)
    if prop=='E(B_V)_VORONOI':
        prop_data = vac[prop].data[galid,:][0]
    elif (prop=='LW_AGE_VORONOI') or (prop=='MW_AGE_VORONOI'):
        prop_data = vac[prop].data[galid,:,0][0]
    else:
        prop_data = vac[prop].data[galid,:,0][0]
    radius = vac[4].data[galid,:,3][0]
    mdata = np.compress(~np.logical_or(prop_data<-10,radius<0),prop_data)
    mradius = np.compress(~np.logical_or(prop_data<-10,radius<0),radius)
    
    ax.scatter(mradius,mdata,c='0.5',s=1) 

    bins = np.linspace(0,1.5,16)
    idx  = np.digitize(mradius,bins)
    running_median = [np.median(mdata[idx==k]) for k in range(16)]

    ax.plot(bins-0.05,running_median,'go-')
    ax.set_ylabel(dict[prop],fontsize=10)

def main(plate,ifu):
    fig1 = plt.figure(figsize=(4,5))
    gs = gridspec.GridSpec(6,1, height_ratios=[1,1,1,1,1,1], width_ratios=[1])
    gs.update(left=0.1, right=0.9, bottom=0.07, top=0.93,hspace=0,wspace=0)
    
    plt.suptitle(str(plate)+'-'+str(ifu))
    
    for i in range(6):
        prop=props[i]
        ax = fig1.add_subplot(gs[i,0])
        plot_prop(prop,ax,plate,ifu)
        
        if i!=5:
            ax.set_xticklabels([])
        if i==5:
            ax.set_xlabel(r'$\boldsymbol{\mathrm{R/R_e}}$')
    
    plt.show()

    
print(sys.argv[1],sys.argv[2])
if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2])



