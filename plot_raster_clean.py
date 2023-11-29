# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 09:43:39 2022

@author: ameli
"""

from brian2 import *

#Indicate which simulation to plot the raster from :
# folder_name='results_noVIP_noSOM_highgLIPFEF_015_2' #name of the folder containing all the simulation results
folder_name='results_uncued_latinh_1003' #name of the folder containing all the simulation results
folder_name='results_uncued_latinh_more_1303'
folder_name='results_lessgranLIP_1303'
folder_name='results_lessgranLIP_1603'
folder_name='results_uncued_real_1503'
folder_name='results_sameobject_2203'
folder_name='results_noise0_0905'
folder_name='results_noise0_1005'
folder_name='results_noise0_1205'
folder_name='results_noise0_1705'
folder_name='results_same_object_5tonicinh'
folder_name='results_same_object_45tonicinh'
foldecdr_name='results_same_object_35tonicinh_LIPRSonly'
folder_name='results_cued_with_VIP'
folder_name='results_cued_with_VIP_v2'
folder_name='results_JRSdec_30'
folder_name='results_JRSdec_10'
folder_name='results_JRSdec_0_bis'
folder_name='results_JRSdec_m10'
# folder_name='results_1sim'
# folder_name='results_VIP_everywhere_JRSdec_0'
# folder_name='results_1sim_J20'
folder_name='results_VIP_everywhere_JRSdec_10'
folder_name='results_VIP_everywhere_JRSdec_ramp_2050'
folder_name='results_1sim_Jramp'
folder_name='results_VIP_everywhere_JRSdec_ramp4020_noise10'
folder_name='results_VIP_everywhere_JRSdec_ramp4020_noise40'

folder_name='results_VIP_everywhere_uncued_JRSdec_ramp6050'

folder_name='results_VIP_everywhere_g2_JRSdec_50_noise10'

folder_name='results_VIP_everywhere_g2_nomdPul'
folder_name='results_VIP_everywhere_g0_nomdPul'

folder_name='results_VIP_everywhere_g0_noise80'

folder_name='results_VIP_everywhere_g0_uncued_noise50'
folder_name='results_VIP_everywhere_g0_uncued_noise10'

folder_name='results_VIP_everywhere_g0_noise100'
folder_name='simulation_results/jitter_0.5'



# folder_name='results_VIP_everywhere_g2_JRSdec_ramp6040_noise0'
# folder_name='results_VIP_everywhere_uncued_JRSdec_ramp8060'
# folder_name='results_VIP_everywhere_uncued_JRSdec_ramp8060_bis' #was actually 6040
# folder_name='results_VIP_everywhere_uncued_JRSdec_ramp8060_ter'
# folder_name='results_VIP_everywhere_uncued_JRSdec_ramp8060_quad'
# folder_name='results_VIP_everywhere_uncued_JRSdec_ramp8060_quint'
# target_time=650*msecond
target_time=500*msecond
tSI=20*msecond #use only if the set of simulations changed the time constant of inhibition
tFS=5*msecond #use only if the set of simulations changed the time constant of inhibition
n=0 #50 simulations are run with each set of parameter, indicate here which one to plot (between 0 and 49)


def read_raster_times_and_indexes(file_t,file_i):
    all_times=[]
    raster_t=open(file_t,'r')
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

def plot_sim_number(nsim,folder=folder_name):
    base=folder+"/results_"+str(nsim)
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
    try:
        SI2vm_t_name=base+"/raster_FEF SI2 vm_t.txt"
        SI2vm_i_name=base+"/raster_FEF SI2 vm_i.txt"
        R9_t,R9_i=read_raster_times_and_indexes(SI2vm_t_name,SI2vm_i_name)
    except:
        SI2vm_t_name=base+"/raster_FEF SOM vm_t.txt"
        SI2vm_i_name=base+"/raster_FEF SOM vm_i.txt"
        R9_t,R9_i=read_raster_times_and_indexes(SI2vm_t_name,SI2vm_i_name)
    try:
        SI1vm_t_name=base+"/raster_FEF SI1 vm_t.txt"
        SI1vm_i_name=base+"/raster_FEF SI1 vm_i.txt"
        R10_t,R10_i=read_raster_times_and_indexes(SI1vm_t_name,SI1vm_i_name)
    except:
        pass
                        
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
    
    try:
        VIPLIP_t_name=base+"/raster_LIP VIP_t.txt"
        VIPLIP_i_name=base+"/raster_LIP VIP_i.txt"
        VIP_LIP_t,VIP_LIP_i=read_raster_times_and_indexes(VIPLIP_t_name,VIPLIP_i_name)   
    except:
        pass
    
    runtime=2*second
    
    tgood_b=[0,]
    tgood_e=[0.125,]
    tbad_b=[0.125,]
    tbad_e=[0.25,]
    while tbad_e[-1]*second<runtime:
        tgood_b.append(tgood_b[-1]+0.25)
        tbad_b.append(tbad_b[-1]+0.25)
        tgood_e.append(tgood_e[-1]+0.25)
        tbad_e.append(tbad_e[-1]+0.25)
        
#    runtime=1*second
#    figure(figsize=(4,9))
    figure(figsize=(8,9))
    hlines([425]*len(tgood_b),tgood_b,tgood_e,'C0')
    hlines([425]*len(tbad_b),tbad_b,tbad_e,'C1')
    up=200
    plot(R1_t,R1_i+140+up,'r.',label='RS cells')
    plot(R2_t,R2_i+120+up,'b.',label='FS cells')
    plot(R3_t,R3_i+100+up,'g.',label='SI cells')
    try:
        plot(VIP_LIP_t,VIP_LIP_i+100+up,'k.',label='VIP cells')
    except:
        pass
    plot([0.2,runtime/second],[95+up,95+up],'k--')
    plot(R5_t,R5_i+70+up,'r.',label='Granular RS')
    plot(R6_t,R6_i+50+up,'b.',label='Granular FS')
    plot([0.2,runtime/second],[45+up,45+up],'k--')
    plot(R4_t,R4_i+20+up,'m.',label='IB cells')
    plot(R7_t,R7_i+up,'g.',label='Deep SI')
    xlim(0.2,runtime/second)
    plot([0.2,runtime/second],[up-10,up-10],'k')
    
    up=40
#    title('FEF Visual Neurons')
    plot(R11_t,R11_i+0+up,'r.',label='RS')
    plot(R12_t,R12_i+20+up,'b.',label='FS')
    plot(R13_t,R13_i+40+up,'k.',label='VIP')
    plot(R14_t,R14_i+60+up,'g.',label='SOM')
    xlim(0.2,runtime/second)
    plot([0.2,runtime/second],[up-10,up-10],'k')

    up=140
#    title('FEF Visual-Motor Neurons')
    hlines([185]*len(tgood_b),tgood_b,tgood_e,'C0')
    hlines([185]*len(tbad_b),tbad_b,tbad_e,'C1')
    plot(R8_t,R8_i+20+up,'r.',label='RS')
    plot(R9_t,R9_i+0+up,'g.',label='SOM')
    try:
        plot(R10_t,R10_i+0+up,'k.',label='VIP')
    except:
        pass
    xlim(0.2,runtime/second)
    plot([0.2,runtime/second],[up-10,up-10],'k')

#    title('FEF Decision cells')
    plot(mon_RS_t,mon_RS_i+0,'r.',label='RS')
    xlim(0.2,runtime/second)
    xlabel('Time (s)',fontsize=12)
    ylabel('Neuron index',fontsize=12)
    xticks(fontsize=12)
    yticks(fontsize=12)
    # ylim(-5,425)
    ylim(-5,430)
    
    tight_layout()    
    subplots_adjust(bottom=0.1)
    
    return
    

# close('all')

#indicate the list of time constants for inhibition here if they were changed in your set of simulation
liste_t_SI=[20*ms]
liste_t_FS=[5*ms]
    
i_tFS=liste_t_FS.index(tFS)
i_tSI=liste_t_SI.index(tSI)

liste_target_time=[300*msecond,400*msecond,500*msecond,600*msecond,700*msecond,800*msecond,900*msecond,1000*msecond,1100*msecond,1200*msecond,1300*msecond,1400*msecond,1500*msecond,1600*msecond,1700*msecond]
liste_target_time+=[350*msecond,450*msecond,550*msecond,650*msecond,750*msecond,850*msecond,950*msecond,1050*msecond,1150*msecond,1250*msecond,1350*msecond,1450*msecond,1550*msecond,1650*msecond]
liste_target_time+=[325*msecond,425*msecond,525*msecond,625*msecond,725*msecond,825*msecond,925*msecond,1025*msecond,1125*msecond,1225*msecond,1325*msecond,1425*msecond,1525*msecond,1625*msecond,1725*msecond]
liste_target_time+=[375*msecond,475*msecond,575*msecond,675*msecond,775*msecond,875*msecond,975*msecond,1075*msecond,1175*msecond,1275*msecond,1375*msecond,1475*msecond,1575*msecond,1675*msecond]

i_target_time=liste_target_time.index(target_time)

N=50
    
liste_simus=[]
for t in liste_target_time:
    liste_simus+=[t]*N
liste_simus=[[i,j,k] for i in liste_simus for j in liste_t_SI for k in liste_t_FS]

nsim=i_target_time*N*len(liste_t_SI)*len(liste_t_FS)+n*len(liste_t_SI)*len(liste_t_FS)+i_tSI*len(liste_t_FS)+i_tFS

nsim=1052

print(nsim)
print(liste_simus[nsim])

plot_sim_number(nsim)
