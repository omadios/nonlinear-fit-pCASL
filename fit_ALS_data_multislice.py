#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Nonlinear fit on multi PLD pCASL data. Will perform the fitting on a ROI (no voxel wise implementation, better to use FSL basil in that case)

@author: Marcello Venzi
"""

import numpy as np
import nibabel as nb
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import os
import sys

#%%

def parse_cmdln():
    parser=argparse.ArgumentParser()
    
    parser.add_argument("-i","--input", help="filename of input 4D mPLD ASL file with n PLD (already M0 calibrated) ", required=True)
    parser.add_argument("-m","--mask", help="filename of mask file containing ROI of interest", required=True)
    
    parser.add_argument("-out","--out_path", help="Output pathname for saving the results", default=os.getcwd())
    parser.add_argument("-dt","--slice_dt", help="slice delta",type=float,default=0.04)
    parser.add_argument("-l","--labeling_time", help="pCASL labeling duration (s)",type=float,default=1.8)
    parser.add_argument("--PLDs", help="list of PLDs in seconds",type=float,default=np.array((2.05,2.30,2.55,2.8,3.05,3.3,3.55,3.8,4.05,4.3,4.55,4.8)))
    parser.add_argument("-s","--save_pic", help="Save picture with the fit (1 or 0)",type=int,default=1)
    
    
    parser.parse_args()
    args=parser.parse_args()

    # check that input filenames have been specified
    
    if not os.path.isfile(args.input):
        raise Exception("Input ASL file does not exist")
    if not os.path.isfile(args.mask):
        raise Exception("Input mask does not exist")
    
    return args


#%%

def pCASL(t,BAT,CBF):
    
    pCASLcurve=np.zeros((t.shape[0],1))
    
    #Fixed parameters
    M0b=1
    alpha=0.85
    global tau
    #tau=1.8
    T1blood=1.65
    T1tissue=1.3
    lambd=0.9
    
    T1app=1/(1/T1tissue + CBF/lambd)
    
    for i in np.arange(0,len(t)):
        if t[i] < BAT:
            pCASLcurve[i,0]=0
        elif t[i] >= BAT and t[i] <= BAT + tau:
            pCASLcurve[i,0]=2*M0b*CBF*T1app*alpha*np.exp(-BAT/T1blood)*(1-np.exp((BAT-t[i])/T1app))
        elif t[i] > BAT + tau:
            pCASLcurve[i,0]=2*M0b*CBF*T1app*alpha*np.exp(-BAT/T1blood)*np.exp((-t[i]+tau+BAT)/T1app)*(1-np.exp(-tau/T1app))    
            
    return pCASLcurve.squeeze()



#%%
    
args=parse_cmdln()

tau = args.labeling_time
TIs=args.PLDs
slicedt=args.slice_dt



#%% 
    
#load with nb the average ASL file and the 
    
img_asl=nb.load(args.input)
img_mask=nb.load(args.mask)

#%%

img_asl_data=img_asl.get_data()
img_mask_data=img_mask.get_data()

#%%


zsum=np.sum(img_mask_data[:,:,:],axis=0).sum(axis=0) #sum over xy for each z -> zslices

zlist=[z for z in range(len(zsum)) if zsum[z] > 0]

meants=np.zeros((len(zlist),TIs.shape[0]))
deltaTI=np.zeros((len(zlist),TIs.shape[0]))

weighTS=np.zeros((len(zlist),TIs.shape[0]))

print(img_asl_data.shape)

for z in range(0,len(zlist)):
    c=np.nonzero(img_mask_data[:,:,zlist[z]])
    meants[z,:]=np.mean(img_asl_data[c[0],c[1],zlist[z],:],axis=0)
    deltaTI[z,:]=slicedt*zlist[z]+TIs #correct for slice timing acquisition
    weighTS[z,:]=weighTS[z,:]+zsum[z] #add weight proportional to the number of voxels in the slice
    #popt, pcov = curve_fit( pCASL, deltaTI[z,:].flatten() , meants[z,:].flatten())
    #print("slice " + str(z) + ' AAT = ' + str(popt[0]) + ' perf = ' + str(popt[1]*6000))
       
#%%
    
tswave=np.average(meants,weights = zsum[zlist], axis=0)
tiwave=np.average(deltaTI,weights = zsum[zlist], axis=0)

popt, pcov = curve_fit( pCASL, tiwave, tswave, bounds = ([0.2,0.0001],[4,0.03]))

  
if args.save_pic==1:
    
    plt.figure()
    plt.plot(np.arange(0,8,0.25), pCASL(np.arange(0,8,0.25), *popt), 'g--')
    plt.plot(tiwave, pCASL(tiwave, *popt), 'g-', label='nlfit: AAT=%5.3f, CBF=%5.3f' % tuple(np.array((popt[0],6000*popt[1]))))
    plt.plot(tiwave, tswave, 'bo', label='data')
    plt.xlabel('t(s)')
    plt.ylabel('deltaM/M0')
    plt.legend()
    name=str(args.out_path) + '/' + str.split(Path(args.input).stem,'.')[0]+ '.' +  str.split(Path(args.mask).stem,'.')[1] + '.png'
    plt.savefig(name)

print(str(np.around(popt[0],decimals=2)) + ',' + str(np.around(6000*popt[1],decimals=2)))
sys.exit(0)