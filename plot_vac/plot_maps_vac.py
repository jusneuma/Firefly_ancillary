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
from matplotlib.colorbar import Colorbar
from matplotlib.ticker import FuncFormatter
import sys
import os

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


def minus_formatter(x, pos):
    return '{:.2f}'.format(x).replace(u'\u22122','-') 

loc_vac = os.getcwd()

# read the vac file
vac = fits.open(os.path.join(loc_vac,'manga-firefly-v3_1_1-miles.fits'))  


props=['LW_AGE_VORONOI','MW_AGE_VORONOI','LW_Z_VORONOI','MW_Z_VORONOI','SURFACE_MASS_DENSITY_VORONOI','E(B_V)_VORONOI']

dict = {'LW_AGE_VORONOI':r'$\boldsymbol{\mathrm{Age_{\ LW}\ [Gyr]}}$',\
        'MW_AGE_VORONOI':r'$\boldsymbol{\mathrm{Age_{\ MW}\ [Gyr]}}$',\
        'LW_Z_VORONOI':r'$\boldsymbol{\mathrm{[Z/H]_{\ LW}}}$',\
        'MW_Z_VORONOI':r'$\boldsymbol{\mathrm{[Z/H]_{\ MW}}}$',\
        'SURFACE_MASS_DENSITY_VORONOI':r'$\boldsymbol{\mathrm{log\ \Sigma_\star\ [M_\odot\,kpc^{-1}]}}$',\
        'E(B_V)_VORONOI':r'$\boldsymbol{\mathrm{E_\text{B-V}}}$'}

def plot_prop(prop,ax,cbax,plate,ifu):
    binid = vac[5].data
    basic = vac[1].data
    galid = basic['PLATEIFU']==str(plate)+'-'+str(ifu)
    if prop=='E(B_V)_VORONOI':
        prop_data = vac[prop].data[galid,:][0]
    elif (prop=='LW_AGE_VORONOI') or (prop=='MW_AGE_VORONOI'):
        prop_data = 10**vac[prop].data[galid,:,0][0]
    else:
        prop_data = vac[prop].data[galid,:,0][0]
    bin1d = vac[4].data[galid,:,0][0]
    image_sz = 80
    maps = np.zeros((image_sz,image_sz))-99
    for i in range(image_sz):
        for j in range(image_sz):
            idbin = (bin1d==binid[galid,i,j])
            if len(bin1d[idbin])==1:
                maps[i,j] = prop_data[idbin]

    R=vac[4].data[galid,:,3].flatten()
    c_idx = np.where(R==np.min(R[R>0]))[0]
    central = bin1d[c_idx]
    cx, cy = np.where(binid[galid,:,:]==central[0])[1][0], np.where(binid[galid,:,:]==central[0])[2][0]
    cx_off, cy_off = vac[4].data[galid,c_idx,1][0], vac[4].data[galid,c_idx,2][0]
	
    extent=[cx_off-cx*0.5,cx_off+(image_sz-cx)*0.5,cy_off-cy*0.5,cy_off+(image_sz-cy)*0.5]

    # only show the spaxels with non-empty values
    masked_array = np.ma.array(maps,mask=(maps<-10))
    
    # plot the map 
    mapl = ax.imshow(masked_array,interpolation='nearest',cmap='RdYlBu_r',origin='lower',extent=extent)
    
    cb = Colorbar(ax = cbax, mappable = mapl, orientation = 'vertical', ticklocation = 'right')
    cb.set_label(dict[prop], labelpad=10)
    cb.ax.tick_params(width=1)
    
    cb.ax.get_yaxis().set_major_formatter(FuncFormatter(minus_formatter))


def main(plate,ifu):
    fig1 = plt.figure(figsize=(4,5))
    gs = gridspec.GridSpec(3,5, height_ratios=[1,1,1], width_ratios=[1,0.05,0.3,1,0.05])
    gs.update(left=0.1, right=0.9, bottom=0.07, top=0.93,hspace=0.05,wspace=0)
    
    plt.suptitle(str(plate)+'-'+str(ifu))
    
    for i in range(6):
        prop=props[i]
        ax = fig1.add_subplot(gs[i//2,3*(i%2)])
        cbax = fig1.add_subplot(gs[i//2,3*(i%2)+1])
        plot_prop(prop,ax,cbax,plate,ifu)
        
        if i//2<2:
            ax.set_xticklabels([])
        if i%2>0:
            ax.set_yticklabels([])
        if i//2==2:
            ax.set_xlabel(r'$\boldsymbol{\mathrm{\Delta\alpha\ [arcsec]}}$')
        if i%2==0:
            ax.set_ylabel(r'$\boldsymbol{\mathrm{\Delta\delta\ [arcsec]}}$')
    
    plt.show()


    
print(sys.argv[1],sys.argv[2])
if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2])



