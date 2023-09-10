#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 14:27:54 2023

@author: brpp
"""
from brian2 import *

import os

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

colors = [(0, 1, 1), (.25, .75, 1), (.5, .5, 1), (.75, .25, 1), (1, 0, 1)]#['tab:purple','blue','tab:green','yellow','tab:orange','red', 'k']
j_list = [45, 40, 36, 31, 27]

hit_rates = [[] for j in j_list]

for i,j in enumerate(j_list):
    path="simulation_results/jRSFEFvm_"+str(j)+"uAcm-2/figures/"
    with open(path+'hit_rates.txt', 'r') as file:
        these_hit_rates = [float(line.strip()) for line in file if line] # file.readlines()
    hit_rates[i] = these_hit_rates

figure()
for i,j in enumerate(j_list):
    plot(liste_target_time, hit_rates[i], color=colors[i], label=str(j)+"uAcm-2")
xlabel('target time (s)')
ylabel('hit rate')
legend(loc='upper right')