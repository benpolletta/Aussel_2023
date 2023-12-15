# -*- coding: utf-8 -*-
"""
Created on Wed May 27 08:13:18 2020

@author: amelie
"""

from brian2 import *

import os

from scipy.stats import norm
import math
Z = norm.ppf

import pylab as plt
import numpy as np
from scipy.stats import chi2
from scipy.stats import spearmanr


#name of the folder containing all the simulation results
# res_folder='uncued'
# res_folder='uncued_2'
# res_folder='uncued_012'
# res_folder='uncued_015_with_LFP'
# res_folder='uncued_015_1601'
# res_folder='un_010_0302'
# res_folder='un_015_0302'
# res_folder='un_020_0602'
# res_folder='uncued_lowFEF_with_LFP'
# res_folder='uncued_lowFEF_jRS25'
# res_folder='cued_015'
# res_folder='cued_010_new'
# res_folder='cued_015_more_input'
# res_folder='uncued_010_2402'
# res_folder='uncued_015_2402'
# res_folder='cued_015_2402'
# res_folder='cued_010_2402'
# res_folder='uncued_withFEFtoLIP_0103'
# res_folder='uncued_withmdPul_0203'
# res_folder='uncued_latinh_1003'
# res_folder='uncued_latinh_more_1303'
# res_folder='lessgranLIP_1303'
# res_folder='uncued_real_1503'
# res_folder='lessgranLIP_1603'
# res_folder='lessgranLIP_2003'
# res_folder='sameobject_2203'
# res_folder='sameobject_3003'
# res_folder='with_noise_2004'
# res_folder='with_noise10Hz_2404'
# res_folder='with_noise_0_3_2804'
# res_folder='with_noise_0_6_2804'
# res_folder='with_noise_0_10_2804'
# res_folder='with_noise_0_3_triple_0405'
# res_folder='poisson2ms'
res_folder = sys.argv[1] #"Jbegin_50_Jend_50"

res_folder = "sim"
for i,arg in enumerate(sys.argv):
    if "=" in arg:
        key, value = arg.split("=")
                 
        if key == "location":
            res_folder +=f"_{value}"
        else:
            res_folder += f"_{key}{value}"

theta_freq = 4
#"jRSFEFvm_27uAcm-2" #'theta_'+str(theta_freq)+'Hz'
theta_period = 1000/theta_freq

num_theta_bins = 10
theta_bin_ms = floor(theta_period/num_theta_bins)
theta_bin_degrees = floor(360/num_theta_bins)
theta_bin_limits = arange(300, 1825, theta_bin_ms)
theta_bin_centers = theta_bin_limits[:-1] + theta_bin_degrees/2

hit_bin_length = 25
hit_bin_limits = arange(300, 1825, hit_bin_length)
hit_bin_centers = hit_bin_limits[:-1]+hit_bin_length/2
total_hit_bins = len(hit_bin_centers)


#Number of simulations for each target time
N=50
# N=25

def read_raster_times(filename):
    all_times=[]
    raster_t=open(filename,'r')
    time_str=[]
    for line in raster_t:
        time_str=line.split(',')[:-1]
    for line in time_str:
        if line[-2]=='m':
            time=float(line[:-2])*msecond
        elif line[-2]=='u':
            time=float(line[:-2])*usecond
        else :
            time=float(line[:-1])*second
        all_times.append(time)
#    print(time_str[5],time_str[-3])
#    print(all_times[5],all_times[-3])
    raster_t.close()
    return all_times

def read_raster_times_and_indexes(file_t,file_i):
    all_times=[]
    raster_t=open(file_t,'r')
    time_str=[]
    for line in raster_t:
        time_str=line.split(',')[:-1]

    for line in time_str:
        if line[-2]=='m':
            time=float(line[:-2])*msecond
        elif line[-2]=='u':
            time=float(line[:-2])*usecond
        else :
            time=float(line[:-1])*second
        all_times.append(time)
    raster_t.close()
        
    all_i=[]
    raster_i=open(file_i,'r')
    for line in raster_i:
        all_i=line.split(',')[:-1]
    all_i=[int(i) for i in all_i]
    raster_i.close()
    return array(all_times),array(all_i)    

def plot_raster(file_t,file_i):
    all_times,all_i=read_raster_times_and_indexes(file_t,file_i)
    
    figure()
    plot(all_times,all_i,'r.')

    return     

def plot_sim_number(nsim):
    base=res_folder+"/results_"+str(nsim)
    RSLIP_t_name=base+"/raster_LIP RS_t.txt"
    RSLIP_i_name=base+"/raster_LIP RS_i.txt"
    R1_t,R1_i=read_raster_times_and_indexes(RSLIP_t_name,RSLIP_i_name)
    FSLIP_t_name=base+"/raster_LIP FS_t.txt"
    FSLIP_i_name=base+"/raster_LIP FS_i.txt"
    R2_t,R2_i=read_raster_times_and_indexes(FSLIP_t_name,FSLIP_i_name)
    SILIP_t_name=base+"/raster_LIP SI_t.txt"
    SILIP_i_name=base+"/raster_LIP SI_i.txt"
    R3_t,R3_i=read_raster_times_and_indexes(SILIP_t_name,SILIP_i_name)
    IBLIP_t_name=base+"/raster_LIP IB_t.txt"
    IBLIP_i_name=base+"/raster_LIP IB_i.txt"
    R4_t,R4_i=read_raster_times_and_indexes(IBLIP_t_name,IBLIP_i_name)
    RSgLIP_t_name=base+"/raster_LIP RS gran_t.txt"
    RSgLIP_i_name=base+"/raster_LIP RS gran_i.txt"
    R5_t,R5_i=read_raster_times_and_indexes(RSgLIP_t_name,RSgLIP_i_name)
    FSgLIP_t_name=base+"/raster_LIP FS gran_t.txt"
    FSgLIP_i_name=base+"/raster_LIP FS gran_i.txt"
    R6_t,R6_i=read_raster_times_and_indexes(FSgLIP_t_name,FSgLIP_i_name)
    SIdLIP_t_name=base+"/raster_LIP SI deep_t.txt"
    SIdLIP_i_name=base+"/raster_LIP SI deep_i.txt"
    R7_t,R7_i=read_raster_times_and_indexes(SIdLIP_t_name,SIdLIP_i_name)
 
    RSvm_t_name=base+"/raster_FEF RS vm_t.txt"
    RSvm_i_name=base+"/raster_FEF RS vm_i.txt"
    R8_t,R8_i=read_raster_times_and_indexes(RSvm_t_name,RSvm_i_name)
    # SI2vm_t_name=base+"/raster_FEF SI2 vm_t.txt"
    # SI2vm_i_name=base+"/raster_FEF SI2 vm_i.txt"
    # R9_t,R9_i=read_raster_times_and_indexes(SI2vm_t_name,SI2vm_i_name)
    SI1vm_t_name=base+"/raster_FEF SOM vm_t.txt"
    SI1vm_i_name=base+"/raster_FEF SOM vm_i.txt"
    R10_t,R10_i=read_raster_times_and_indexes(SI1vm_t_name,SI1vm_i_name)
                        
    RSv_t_name=base+"/raster_FEF RS v_t.txt"
    RSv_i_name=base+"/raster_FEF RS v_i.txt"
    R11_t,R11_i=read_raster_times_and_indexes(RSv_t_name,RSv_i_name)
    FSv_t_name=base+"/raster_FEF FS v_t.txt"
    FSv_i_name=base+"/raster_FEF FS v_i.txt"
    R12_t,R12_i=read_raster_times_and_indexes(FSv_t_name,FSv_i_name)
    VIPv_t_name=base+"/raster_FEF VIP v_t.txt"
    VIPv_i_name=base+"/raster_FEF VIP v_i.txt"
    R13_t,R13_i=read_raster_times_and_indexes(VIPv_t_name,VIPv_i_name)
    SIv_t_name=base+"/raster_FEF SI v_t.txt"
    SIv_i_name=base+"/raster_FEF SI v_i.txt"
    R14_t,R14_i=read_raster_times_and_indexes(SIv_t_name,SIv_i_name)
                              
    RSm_t_name=base+"/raster_FEF RS m_t.txt"
    RSm_i_name=base+"/raster_FEF RS m_i.txt"
    mon_RS_t,mon_RS_i=read_raster_times_and_indexes(RSm_t_name,RSm_i_name)
    
    
    runtime=2*second
    figure()
    plot(R1_t,R1_i+140,'r.',label='RS cells')
    plot(R2_t,R2_i+120,'b.',label='FS cells')
    plot(R3_t,R3_i+100,'g.',label='SI cells')
    plot(R5_t,R5_i+70,'.',label='Granular RS',color='C1')
    plot(R6_t,R6_i+50,'c.',label='Granular FS')
    plot(R4_t,R4_i+20,'m.',label='IB cells')
    plot(R7_t,R7_i,'.',label='Deep SI',color='lime')
    xlim(0,runtime/second)
    legend(loc='upper left')
    
    #FEF Plots    
    figure(figsize=(10,4))
    subplot(131)
    title('Visual Neurons')
    plot(R11_t,R11_i+0,'r.',label='RS')
    plot(R12_t,R12_i+20,'b.',label='FS')
    plot(R13_t,R13_i+40,'k.',label='VIP')
    plot(R14_t,R14_i+60,'g.',label='SI')
    xlim(0,runtime/second)
    legend(loc='upper left')   
    
    subplot(132)
    title('Visual-Motor Neurons')
    plot(R10_t,R10_i+0,'g.',label='SOM')
    plot(R8_t,R8_i+20,'r.',label='RS')
    # plot(R9_t,R9_i+40,'.',label='SI 2',color='lime')
    xlim(0,runtime/second)
    legend(loc='upper left') 
    
    subplot(133)
    title('Motor Neurons')
    plot(mon_RS_t,mon_RS_i+0,'r.',label='RS')
    xlim(0,runtime/second)
    legend(loc='upper left')     
    return

def get_N_spikes(times,target_time):
    if len(times)==0:
        return 0
    if len(where(times>target_time+100*msecond)[0])>0: #100ms or 200ms?
        return where(times>target_time+100*msecond)[0][0]-where(times>target_time)[0][0]
    else :
        if len(where(times<target_time)[0])>0:
#        return len(where(times>target_time)[0])   
            return len(times)-where(times<target_time)[0][-1]
        else : 
            return len(times)

def false_positive(times,target_time):
    if len(times)==0:
        return 0,0,0
    if len(where(times>target_time+100*msecond)[0])>0:
        times_modif=hstack((times[:where(times>target_time)[0][0]],times[where(times>target_time+100*msecond)[0][0]:]))
        times_modif=times_modif*second
    else :
        times_modif=times[where(times<target_time)]  
    
    n_spikes_per_bin=[]
#    for t_deb in [100*ms*i for i in range(20)]:
#        if len(where(times>t_deb+100*msecond)[0])>0:
#            n_spikes_per_bin.append(where(times>t_deb+100*msecond)[0][0]-where(times>t_deb)[0][0])
#        else :
#            n_spikes_per_bin.append(len(where(times>t_deb)[0]))
    for t_deb in [25*ms*(i+20) for i in range(60)]:
        if len(where(times_modif>t_deb+25*msecond)[0])>0:
            n_spikes_per_bin.append(where(times_modif>t_deb+25*msecond)[0][0]-where(times_modif>t_deb)[0][0])
        else :
            n_spikes_per_bin.append(len(where(times_modif>t_deb)[0]))
    n_spikes_per_bin=array(n_spikes_per_bin)
    good_theta=array([[0+10*i,1+10*i,2+10*i,3+10*i,4+10*i] for i in range(6)]).flatten()
    bad_theta=array([[5+10*i,6+10*i,7+10*i,8+10*i,9+10*i] for i in range(6)]).flatten()
    n_spikes_good_theta=n_spikes_per_bin[good_theta]
    n_spikes_bad_theta=n_spikes_per_bin[bad_theta]
    return min(1,sum(n_spikes_per_bin>10)),min(1,sum(n_spikes_good_theta>10)),min(1,sum(n_spikes_bad_theta>10))

def get_hits_per_bin(times,target_time):
    if len(times)==0:
        return [0,0,0,0]
    if len(where(times>target_time+100*msecond)[0])>0:
        times_modif=times[where(times>target_time)[0][0]:where(times>target_time+100*msecond)[0][0]]
        # times_modif=times_modif*second
    elif len(where(times>target_time)[0])>0:
        times_modif=times[where(times>target_time)[0][0]:where(times>target_time)[0][-1]]
    else :
        return [0,0,0,0]
    
    n_spikes_per_bin=[]
    for t_deb in [target_time+i*hit_bin_length*ms for i in range(4)]:
        if len(where(times_modif>t_deb+hit_bin_length*msecond)[0])>0:
            n_spikes_per_bin.append(where(times_modif>t_deb+hit_bin_length*msecond)[0][0]-where(times_modif>t_deb)[0][0])
        else :
            n_spikes_per_bin.append(len(where(times_modif>t_deb)[0]))
    n_spikes_per_bin=array(n_spikes_per_bin)
    hits_per_bin=[(n_spikes>10) for n_spikes in n_spikes_per_bin]
    return hits_per_bin


liste_target_time=[300*msecond,400*msecond,500*msecond,600*msecond,700*msecond,800*msecond,900*msecond,1000*msecond,1100*msecond,1200*msecond,1300*msecond,1400*msecond,1500*msecond,1600*msecond,1700*msecond]
liste_target_time+=[350*msecond,450*msecond,550*msecond,650*msecond,750*msecond,850*msecond,950*msecond,1050*msecond,1150*msecond,1250*msecond,1350*msecond,1450*msecond,1550*msecond,1650*msecond]

liste_target_time+=[325*msecond,425*msecond,525*msecond,625*msecond,725*msecond,825*msecond,925*msecond,1025*msecond,1125*msecond,1225*msecond,1325*msecond,1425*msecond,1525*msecond,1625*msecond,1725*msecond]
liste_target_time+=[375*msecond,475*msecond,575*msecond,675*msecond,775*msecond,875*msecond,975*msecond,1075*msecond,1175*msecond,1275*msecond,1375*msecond,1475*msecond,1575*msecond,1675*msecond]

liste_simus=[]
for t in liste_target_time:
    liste_simus+=[t]*N

liste_simus=[[liste_simus[i],i] for i in range(len(liste_simus))]
simus_pas_faites=[]

number_of_spikes=[[] for i in range(len(liste_target_time))]
phase_rm=[[] for i in range(len(liste_target_time))]
time_last_rm=[[] for i in range(len(liste_target_time))]
number_spikes_vm=[[] for i in range(len(liste_target_time))]
number_spikes_v=[[] for i in range(len(liste_target_time))]
phase_lip_gamma=[[] for i in range(len(liste_target_time))]
#phase_lip_beta=[[] for i in range(len(liste_target_time))]
number_spikes_no_target=[[] for i in range(len(liste_target_time))]
number_false_positives=[[] for i in range(len(liste_target_time))]
number_fp_good=[[] for i in range(len(liste_target_time))]
number_fp_bad=[[] for i in range(len(liste_target_time))]
#number_fp_per_bin=[[] for i in range(80)]
freq_lip=[[] for i in range(len(liste_target_time))]

# Nspikes_LIP_RSsup=[[] for i in range(len(liste_target_time))]
# Nspikes_FEFvm_RS=[[] for i in range(len(liste_target_time))]
Nspikes_FEFvm_SOM=[[] for i in range(len(liste_target_time))]
Nspikes_FEFvis_RS=[[] for i in range(len(liste_target_time))]
Nspikes_FEFvis_SOM=[[] for i in range(len(liste_target_time))]
Nspikes_FEFvis_VIP=[[] for i in range(len(liste_target_time))]

N_hits_per_bin=[0]*total_hit_bins #number of 25ms windows from the beginning of first target window time (300ms) to end of last target window time (1825ms)
SOM_per_bin=[0]*total_hit_bins

    
for simu in liste_simus:
    target_time,n_sim=simu
    i_target_time=liste_target_time.index(target_time)
    try :
        raster_path = "simulation_results/"+res_folder+"/results_"+str(n_sim+1)
        raster_t_name = raster_path + "/raster_FEF RS m_t.txt"
        times=array(read_raster_times(raster_t_name))*second
        number_of_spikes[i_target_time].append((get_N_spikes(times,target_time)))
        number_spikes_no_target[i_target_time].append(len(times)-number_of_spikes[i_target_time][-1])
        fp,fp_good,fp_bad=false_positive(times,target_time)
        number_false_positives[i_target_time].append(fp)
        number_fp_good[i_target_time].append(fp_good)
        number_fp_bad[i_target_time].append(fp_bad)
        
        hits_per_bin=get_hits_per_bin(times,target_time)
        for k in range(4):
            N_hits_per_bin[int((target_time-300*ms)/hit_bin_length/ms)+k-1]+=hits_per_bin[k]
        
        
        raster_t_name2=raster_path+"/raster_FEF RS vm_t.txt"
        times2=array(read_raster_times(raster_t_name2))*second
    
        before_rm=times2[where(times2<(target_time+50*msecond))[0][-1]]
        try:
            after_rm=times2[where(times2>(target_time+50*msecond))[0][0]]
        except:
            after_rm=(target_time+50*msecond)
        phase=(target_time+50*msecond-before_rm)/(after_rm-before_rm)*360
    
        phase_rm[i_target_time].append(phase)
        time_last_rm[i_target_time].append(target_time-before_rm)
        number_spikes_vm[i_target_time].append(get_N_spikes(times2,target_time))
    
        raster_t_name3=raster_path+"/raster_LIP RS_t.txt"
        times3=array(read_raster_times(raster_t_name3))*second
        number_spikes_v[i_target_time].append(get_N_spikes(times3,target_time))
    
        before_rs=times3[where(times3<(target_time+50*msecond))[0][-1]]
        try :
            after_rs=times3[where(times3>(target_time+50*msecond))[0][0]]
        except :
            after_rs=target_time+50*msecond
            
        phase=((target_time+50*msecond)-before_rs)/(after_rs-before_rs)*360
        phase_lip_gamma[i_target_time].append(phase)
        freq_lip[i_target_time].append(1//(after_rs-before_rs))
        
        raster_t_name4=raster_path+"/raster_FEF RS v_t.txt"
        times4=array(read_raster_times(raster_t_name4))*second
        Nspikes_FEFvis_RS[i_target_time].append(get_N_spikes(times4,target_time))
        
        raster_t_name4=raster_path+"/raster_FEF SI v_t.txt"
        times4=array(read_raster_times(raster_t_name4))*second
        Nspikes_FEFvis_SOM[i_target_time].append(get_N_spikes(times4,target_time))
        
        raster_t_name4=raster_path+"/raster_FEF VIP v_t.txt"
        times4=array(read_raster_times(raster_t_name4))*second
        Nspikes_FEFvis_VIP[i_target_time].append(get_N_spikes(times4,target_time))
        
        raster_t_nameSOM=raster_path+"/raster_FEF SOM vm_t.txt"
        timesSOM=array(read_raster_times(raster_t_nameSOM))*second
        Nspikes_FEFvm_SOM[i_target_time].append(get_N_spikes(times4,target_time))
        
        SOM_per_target=get_hits_per_bin(times4,target_time)
        for k in range(4):
            SOM_per_bin[int((target_time-300*ms)/hit_bin_length/ms)+k-1]+=SOM_per_target[k]
        
    except Exception as e:
        print(e)
        simus_pas_faites.append(n_sim)

if len(simus_pas_faites)>0:    
    print(str(len(simus_pas_faites))+' simulations were not processed.')    
    print('These simulations were not processed: '+str(simus_pas_faites))
    

close('all')

mean_spikes=[mean(Ns) for Ns in number_of_spikes]
std_spikes=[std(Ns) for Ns in number_of_spikes]
mean_spikes_v=[mean(Ns) for Ns in number_spikes_v]
std_spikes_v=[std(Ns) for Ns in number_spikes_v]
mean_spikes_vm=[mean(Ns) for Ns in number_spikes_vm]
std_spikes_vm=[std(Ns) for Ns in number_spikes_vm]
mean_spikes_no_target=[mean(Ns) for Ns in number_spikes_no_target]
std_spikes_no_target=[std(Ns) for Ns in number_spikes_no_target]
mean_spikes_FEFvRS=[mean(Ns) for Ns in Nspikes_FEFvis_RS]
std_spikes_FEFvRS=[std(Ns) for Ns in Nspikes_FEFvis_RS]
mean_spikes_FEFvSOM=[mean(Ns) for Ns in Nspikes_FEFvis_SOM]
std_spikes_FEFvSOM=[std(Ns) for Ns in Nspikes_FEFvis_SOM]
mean_spikes_FEFvVIP=[mean(Ns) for Ns in Nspikes_FEFvis_VIP]
std_spikes_FEFvVIP=[std(Ns) for Ns in Nspikes_FEFvis_VIP]

order=[0,29,15,44]
for i in range(1,14):
    order+=[0+i,29+i,15+i,44+i]
order+=[14,43]

number_of_spikes=[number_of_spikes[i] for i in order]
mean_spikes=[mean_spikes[i] for i in order]
std_spikes=[std_spikes[i] for i in order]
liste_target_time=[liste_target_time[i] for i in order]
phase_rm=[phase_rm[i] for i in order]
time_last_rm=[time_last_rm[i] for i in order]
phase_lip_gamma=[phase_lip_gamma[i] for i in order]
freq_lip=[freq_lip[i] for i in order]
#phase_lip_beta=[phase_lip_beta[i] for i in order]
mean_spikes_v=[mean_spikes_v[i] for i in order]
std_spikes_v=[std_spikes_v[i] for i in order]
mean_spikes_vm=[mean_spikes_vm[i] for i in order]
std_spikes_vm=[std_spikes_vm[i] for i in order]
number_spikes_no_target=[number_spikes_no_target[i] for i in order]
mean_spikes_no_target=[mean_spikes_no_target[i] for i in order]
std_spikes_no_target=[std_spikes_no_target[i] for i in order]

mean_spikes_FEFvRS=[mean_spikes_FEFvRS[i] for i in order]
std_spikes_FEFvRS=[std_spikes_FEFvRS[i] for i in order]
mean_spikes_FEFvSOM=[mean_spikes_FEFvSOM[i] for i in order]
std_spikes_FEFvSOM=[std_spikes_FEFvSOM[i] for i in order]
mean_spikes_FEFvVIP=[mean_spikes_FEFvVIP[i] for i in order]
std_spikes_FEFvVIP=[std_spikes_FEFvVIP[i] for i in order]

number_false_positives=[number_false_positives[i] for i in order]
number_fp_good=[number_fp_good[i] for i in order]
number_fp_bad=[number_fp_bad[i] for i in order]

threshold = 10
#threshold = 12
hit_rates=[sum(array(Ns)>threshold)/len(Ns) for Ns in number_of_spikes]

figure()
errorbar(liste_target_time, mean_spikes, yerr=std_spikes, fmt='o')
plot(liste_target_time, mean_spikes,'b')
xlabel('Target time')
ylabel('Number of decision cells spikes')
title(res_folder)

figure()
plot(liste_target_time,hit_rates)
xlabel('Target time')
ylabel('Hit rates')
title(res_folder)

figure()
errorbar(liste_target_time, mean_spikes_vm, yerr=std_spikes_vm, fmt='o')
plot(liste_target_time, mean_spikes_vm,'b')
xlabel('Target time')
ylabel('Number of FEFvm RS cells spikes')
title(res_folder)

figure()
errorbar(liste_target_time, mean_spikes_v, yerr=std_spikes_v, fmt='o')
plot(liste_target_time, mean_spikes_v,'b')
xlabel('Target time')
ylabel('Number of LIPsup RS cells spikes')
title(res_folder)

figure()
errorbar(liste_target_time, mean_spikes_FEFvRS, yerr=std_spikes_FEFvRS, fmt='o')
plot(liste_target_time, mean_spikes_FEFvRS,'b')
xlabel('Target time')
ylabel('Number of FEFv RS cells spikes')

figure()
errorbar(liste_target_time, mean_spikes_FEFvSOM, yerr=std_spikes_FEFvSOM, fmt='o')
plot(liste_target_time, mean_spikes_FEFvSOM,'b')
xlabel('Target time')
ylabel('Number of FEFv SOM cells spikes')

figure()
errorbar(liste_target_time, mean_spikes_FEFvVIP, yerr=std_spikes_FEFvVIP, fmt='o')
plot(liste_target_time, mean_spikes_FEFvVIP,'b')
xlabel('Target time')
ylabel('Number of FEFv VIP cells spikes')
title(res_folder)

figure()
plot(liste_target_time, mean_spikes_FEFvRS-mean(mean_spikes_FEFvRS),'r',label='FEFv RS')
plot(liste_target_time, mean_spikes_FEFvSOM-mean(mean_spikes_FEFvSOM),'g',label='FEFv SOM')
plot(liste_target_time, mean_spikes_FEFvVIP-mean(mean_spikes_FEFvVIP),'k',label='FEFv VIP')
xlabel('Target time')
ylabel('Number of spikes during target presentation - mean ')
legend()
title(res_folder)

from scipy import signal

f,Spectrum=signal.periodogram(hit_rates, 40,'flattop', scaling='spectrum')

figure()
plot(f,Spectrum)
xlabel('Frequency (Hz)')
ylabel('Hit rates Power spectrum')
title(res_folder)

f,Spectrum=signal.periodogram(mean_spikes, 40,'flattop', scaling='spectrum')

figure()
plot(f,Spectrum)
xlabel('Frequency (Hz)')
ylabel('Power spectrum of the number of decision cells spikes')
title(res_folder)

#    theta_phase=[((i+2)*36)%360 for i in range(len(liste_target_time))]
# theta_phase=[((i+4)*(theta_period/25))%360 for i in range(len(liste_target_time))]
theta_phase=[360*((i+50*ms)%(theta_period*ms)/(theta_period*ms)) for i in liste_target_time]
phase_indices = floor(theta_phase/theta_bin_degrees)
phase_indices = phase_indices.astype(int)
hr_per_theta=[[] for i in range(num_theta_bins)]
for k in range(len(hit_rates)):
    hr_per_theta[phase_indices[k]].append(hit_rates[k])
hr_per_theta=[mean(h) for h in hr_per_theta]

hr_per_increasing_theta=hr_per_theta#[6:10]+hr_per_theta[:6]
increasing_theta_phase=[(i+.5)*theta_bin_degrees for i in arange(num_theta_bins)]#theta_phase#[6:10]+theta_phase[:6]


figure()
plot(increasing_theta_phase, hr_per_increasing_theta,'b')
xlabel('Theta phase')
ylabel('Hit rates')
title(res_folder)

# #If some simulations are missing and there are "nan"s at the end :
# f,Spectrum=signal.periodogram(hit_rates[:-3], 40,'flattop', scaling='spectrum')

# figure()
# plot(f,Spectrum)
# xlabel('Frequency (Hz)')
# ylabel('Hit rates Power spectrum')


    
#theta_phase=[((i+2)*36)%360 for i in range(len(liste_target_time))]
# theta_phase=[((i+4)*36)%360 for i in range(len(liste_target_time))]
# good_theta=where(array(theta_phase)<=180)[0]
# bad_theta=where(array(theta_phase)>180)[0]

#print(good_theta,bad_theta)


# from scipy.stats import spearmanr

# bin_size=30 #size of bin for phase in degrees   
# flat_gamma_phase=array(phase_lip_gamma)[good_theta].flatten()
# flat_n_spikes=array(number_of_spikes)[good_theta].flatten()
# phase_binned=[int(i//bin_size) for i in flat_gamma_phase]
# spike_phase_binned_good=[[] for i in range(int(360//bin_size))]
# for i in range(len(flat_gamma_phase)):
#     spike_phase_binned_good[phase_binned[i]].append(flat_n_spikes[i])
    
# mean_phase_spikes=[mean(i) for i in spike_phase_binned_good]
# std_phase_spikes=[std(i) for i in spike_phase_binned_good]



# #theta_phase=[((i+4)*36)%360 for i in range(len(liste_target_time))]
# #theta_phase=[((i+2)*36)%360 for i in range(len(liste_target_time))]
# theta_phase=[((i-2)*36)%360 for i in range(len(liste_target_time))]
# good_theta=where(array(theta_phase)<=180)[0]
# bad_theta=where(array(theta_phase)>180)[0]

# bin_size=30 #size of bin for phase in degrees  
# flat_gamma_phase=array(phase_lip_gamma)[bad_theta].flatten()
# flat_n_spikes=array(number_of_spikes)[bad_theta].flatten()
# #flat_n_spikes=array(hit_rates)[bad_theta].flatten()
# phase_binned=[int(i//bin_size) for i in flat_gamma_phase]
# spike_phase_binned_bad=[[] for i in range(int(360//bin_size))]
# for i in range(len(flat_gamma_phase)):
#     spike_phase_binned_bad[phase_binned[i]].append(flat_n_spikes[i])
# #    spike_phase_binned_bad[phase_binned[i]].append(flat_n_spikes[i//50])
    
# mean_phase_spikes_bad=[mean(i) for i in spike_phase_binned_bad]
# std_phase_spikes_bad=[std(i) for i in spike_phase_binned_bad]

# figure()
# errorbar(list(range(0,360,bin_size)), mean_phase_spikes_bad, yerr=std_phase_spikes_bad, fmt='o')
# plot(list(range(0,360,bin_size)), mean_phase_spikes_bad,'b')
# xlabel('LIP Beta phase')
# ylabel('Number of decision cells spikes')
# ylim(-5,20)
# tight_layout()


# print('Correlation of n_spikes to LIP gamma in good theta phase :')
# print(spearmanr(array(number_of_spikes)[good_theta].flatten(), array(phase_lip_gamma)[good_theta].flatten()))
# print('Correlation of n_spikes  to LIP beta in bad theta phase :')
# print(spearmanr(array(number_of_spikes)[bad_theta].flatten(), array(phase_lip_gamma)[bad_theta].flatten()))

# bin_size=30 #size of bin for phase in degrees  
# flat_gamma_phase=array(phase_lip_gamma)[bad_theta].flatten()
# flat_n_spikes=array(number_of_spikes)[bad_theta].flatten()
# #flat_n_spikes=array(hit_rates)[bad_theta].flatten()
# phase_binned=[int(i//bin_size) for i in flat_gamma_phase]
# spike_phase_binned_bad=[[] for i in range(int(360//bin_size))]
# for i in range(len(flat_gamma_phase)):
#     spike_phase_binned_bad[phase_binned[i]].append(flat_n_spikes[i])
# #    spike_phase_binned_bad[phase_binned[i]].append(flat_n_spikes[i//50])
    

# hit_rates_binned=[len(where(array(spikes)>10)[0])/len(spikes) for spikes in spike_phase_binned_bad]

# figure()
# plot(list(range(0,360,bin_size)), hit_rates_binned, 'o-')
# xlabel('LIP Beta phase')
# ylabel('Hit rate')
# ylim(0,1)
# tight_layout()
    
#Plotting total number of hits per 25ms time bin
#This is not a hit rate
bin_time = [((i+.5)*hit_bin_length+300)*ms for i in range(total_hit_bins)]
figure()
plot(bin_time, N_hits_per_bin)
xlabel('Time (ms)')
ylabel('Hits per 25ms time bin')
title(res_folder)

f,Spectrum=signal.periodogram(N_hits_per_bin, 40,'flattop', scaling='spectrum')

figure()
plot(f,Spectrum)
xlabel('Frequency (Hz)')
ylabel('Hit rates per 25ms time bin, power spectrum')

#Plotting total hit rates per 1/20th period time bin as function of theta phase
phase_indices = floor((360*(hit_bin_centers%theta_period)/theta_period)/theta_bin_degrees)

#To make it simpler we will only consider the time bins that constitute a whole number of theta cycles
cycle_start_bin = np.where(np.diff(phase_indices)<0)[0][0]+1#floor((theta_period-300)/hit_bin_length)+1
cycle_end_bin = np.where(np.diff(phase_indices)<0)[0][-1]+1
N_bins = [sum(phase_indices[cycle_start_bin:cycle_end_bin] == k) for k in arange(num_theta_bins)]

N_hits_per_bin_per_theta=[0]*num_theta_bins
k=cycle_start_bin;
for hits in N_hits_per_bin[cycle_start_bin:cycle_end_bin]:
    N_hits_per_bin_per_theta[k]+=hits
    k=(k+1)%num_theta_bins
N_hits_per_bin_per_theta = divide(N_hits_per_bin_per_theta, N_bins)

figure()
plot([theta_bin_degrees*k for k in range(num_theta_bins)]-theta_bin_degrees/2,N_hits_per_bin_per_theta)
xlabel('Theta phase (degree)')
ylabel('Hits per 25ms time bin')

#Plotting total number of hits per 25ms time bin
#This is not a hit rate
figure()
plot(bin_time,SOM_per_bin, 'g', label='SOM spikes')
plot(bin_time,N_hits_per_bin, 'k', label='hits')
xlabel('Time (ms)')
ylabel('SOM spikes per 25ms time bin')


f,Spectrum=signal.periodogram(SOM_per_bin, 40,'flattop', scaling='spectrum')

figure()
plot(f,Spectrum)
xlabel('Frequency (Hz)')
ylabel('SOM spikes per 25ms time bin, power spectrum')

#Plotting total hit rates per 25ms time bin as function of theta phase
#to make it simpler we will only consider the time bins from 375ms to 1625ms so we will have 5 full theta cycles
SOM_per_bin_per_theta=[0]*num_theta_bins
k=cycle_start_bin;
for hits in SOM_per_bin[cycle_start_bin:cycle_end_bin]:
    SOM_per_bin_per_theta[k]+=hits
    k=(k+1)%num_theta_bins
SOM_per_bin_per_theta = divide(SOM_per_bin_per_theta, N_bins)

figure()
plot([theta_bin_degrees*k for k in range(num_theta_bins)]+theta_bin_degrees/2,SOM_per_bin_per_theta, 'g', label='SOM spikes')
plot([theta_bin_degrees*k for k in range(num_theta_bins)]+theta_bin_degrees/2,N_hits_per_bin_per_theta, 'k', label='hits')
xlabel('Theta phase (degree)')
ylabel('Mean over 25ms time bins')

titles = ['decision_cells', 'hit_rates', 'FEFvm_RS', 'LIPsup_RS', 'FEFv_RS', 'FEFv_SOM', 'FEFv_VIP', 'FEFv_RS_SOM_VIP', 'hit_rate_spectrum', 'decision_cell_spectrum', 'hit_rate_by_theta_phase', 'hit_rate_by_25ms_bin', 'hit_rate_by_25ms_bin_spectrum', 'hit_rate_25ms_bin_by_theta_phase', 'SOM_25ms_bin', 'SOM_25ms_bin_spectrum', 'SOM_25ms_bin_by_theta_phase']

path="simulation_results/"+res_folder+'/figures/'
if not os.path.exists(path):
    os.mkdir(path)
for i in get_fignums():
    current_fig=figure(i)
    title(res_folder)
    current_fig.savefig(path+titles[i-1]+'.png')
    
with open(path+'hit_rates.txt', 'w') as file:
    file.writelines("%f\n" % hr for hr in hit_rates)
    
    
# with open(path+'hit_rates.txt', 'r') as file:
#     other_hit_rates = file.readlines()

    
