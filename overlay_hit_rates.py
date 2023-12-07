#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 14:27:54 2023

@author: brpp
"""
from brian2 import *

import os

from scipy import signal
from scipy.stats import norm
import math
Z = norm.ppf

import pylab as plt
import numpy as np

liste_target_time=[300*msecond,400*msecond,500*msecond,600*msecond,700*msecond,800*msecond,900*msecond,1000*msecond,1100*msecond,1200*msecond,1300*msecond,1400*msecond,1500*msecond,1600*msecond,1700*msecond]
liste_target_time+=[350*msecond,450*msecond,550*msecond,650*msecond,750*msecond,850*msecond,950*msecond,1050*msecond,1150*msecond,1250*msecond,1350*msecond,1450*msecond,1550*msecond,1650*msecond]

liste_target_time+=[325*msecond,425*msecond,525*msecond,625*msecond,725*msecond,825*msecond,925*msecond,1025*msecond,1125*msecond,1225*msecond,1325*msecond,1425*msecond,1525*msecond,1625*msecond,1725*msecond]
liste_target_time+=[375*msecond,475*msecond,575*msecond,675*msecond,775*msecond,875*msecond,975*msecond,1075*msecond,1175*msecond,1275*msecond,1375*msecond,1475*msecond,1575*msecond,1675*msecond]

order=[0,29,15,44]
for i in range(1,14):
    order+=[0+i,29+i,15+i,44+i]
order+=[14,43]

liste_target_time=[liste_target_time[i] for i in order]

prefix = 'jitter' #'gPul' # 
suffix = '' #'mScm-2' # 
j_list = [0.75, 1.0] #[1, 2, 3, 4, 5] #
description = 'Jitter During Poor Theta Phase' #'mdPul Input Conductance'
short_description = 'jitter' # 'mdPul Input'

# prefix = 'gPul_offset' # 
# suffix = 'mScm-2' # 
# j_list = [0, 2.5, 5, 7.5, 10] #
# description = 'mdPul Input Conductance (Offset from 5/2.5 mScm-2 good/bad)'
# short_description = 'mdPul Input Offset'

# prefix = 'j_offset'
# suffix = 'uAcm-2'
# j_list = [0, -5, -10, -15]
# description = 'FEF & LIP Tonic Input (Offset, uA/cm^2)' #'mdPul Input Conductance (Offset from 5/2.5 mScm-2 good/bad)'
# short_description = 'FEF & LIP Tonic Offset'

# prefix = 'jRSFEFvm'
# suffix = 'uAcm-2'
# j_list = [54, 50, 45, 40, 36, 31, 27]
# description = 'FEFvm RS Tonic Input (uA/cm^2)' #'mdPul Input Conductance (Offset from 5/2.5 mScm-2 good/bad)'
# short_description = 'FEFvm RS Tonic'

r = linspace(0, 1, len(j_list)).tolist()
g = linspace(1, 0, len(j_list)).tolist()
b = ones(shape(r)).tolist()
colors = list(zip(r, g, b))

hit_rates = [[] for j in j_list]
spectra = [[] for j in j_list]

for i,j in enumerate(j_list):
    path="simulation_results/"+prefix+"_"+str(j)+suffix+"/figures/" #" gPul_offset_"+str(j)+"mScm-2/figures/"
    with open(path+'hit_rates.txt', 'r') as file:
        these_hit_rates = [float(line.strip()) for line in file if line] # file.readlines()
    hit_rates[i] = these_hit_rates
    f,Spectrum=signal.periodogram(these_hit_rates, 40,'flattop', scaling='spectrum')
    spectra[i] = Spectrum
    
mean_hit_rate = [mean(hr) for hr in hit_rates]

figure()
plot(j_list, mean_hit_rate,'ko--', markersize=12)
xticks(j_list)
xlabel(description) #)
ylabel('mean hit rate')

figure(figsize=(10,6))
# imshow(hit_rates) #liste_target_time, j_list, hit_rates)
for i,j in enumerate(j_list):
    plot(liste_target_time, hit_rates[i], color=colors[i], label=str(j)+suffix, alpha=0.75, linewidth=2)
xlabel('target time (s)')
ylabel('hit rate')
legend(loc='best', frameon=False, title=short_description)

figure()
for i,j in enumerate(j_list):
    plot(f, spectra[i], color=colors[i], label=str(j)+suffix, alpha = 0.75, linewidth=2.5)
xlabel('frequency (Hz)')
ylabel('power')
legend(loc='upper right', frameon=False, title=short_description)

t_suffixes = ['mean_hr','hr','hr_spectrum']
path="simulation_results/"
t_stem = prefix+'_'+str(min(j_list))+'to'+str(max(j_list))+'_'
if not os.path.exists(path):
    os.mkdir(path)
for i in get_fignums():
    current_fig=figure(i)
    current_fig.savefig(path+t_stem+t_suffixes[i-1]+'.eps', format='eps')

