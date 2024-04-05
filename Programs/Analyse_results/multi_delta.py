import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import os

###=============================================================
###=============================================================

path_only_sig='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/angles/only_sig'
path_sig_bg='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/angles/bg_sig'
path_elec='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/angles/electronic'

path_sph3_Fe='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/angles/sph3/Fe'
path_sph3_N='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/angles/sph3/N'
path_sph3_P='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/angles/sph3/P'



real_theta=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_theta',header=None,sep='\s+')[0]) #number of PMT
real_phi=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_phi',header=None,sep='\s+')[0]) #number of PMT

real_phi_sph3_Fe=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/angels_phi_Fe',header=None,sep='\s+')[0])
real_phi_sph3_N=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/angels_phi_N',header=None,sep='\s+')[0])
real_phi_sph3_P=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/angels_phi_P',header=None,sep='\s+')[0])

##=== Open directories

name_files_only_sig = glob.glob(os.path.join(path_only_sig,'*'))
name_files_sig_bg = glob.glob(os.path.join(path_sig_bg,'*'))
name_files_elec = glob.glob(os.path.join(path_elec,'*'))

name_files_sph3_Fe = glob.glob(os.path.join(path_sph3_Fe,'*'))
name_files_sph3_N = glob.glob(os.path.join(path_sph3_N,'*'))
name_files_sph3_P = glob.glob(os.path.join(path_sph3_P,'*'))

###=============================================================
###=============================================================


def f_delta(theta_fit,phi_fit,theta_real,phi_real):
    cx = np.sin(theta_real) * np.cos(phi_real) #theta - истинные углы
    cy = np.sin(theta_real) * np.sin(phi_real)
    cz = np.cos(theta_real)
    cx_fit = np.sin(theta_fit) * np.cos(phi_fit) #theta' -полученные углы моделированием
    cy_fit = np.sin(theta_fit) * np.sin(phi_fit)
    cz_fit = np.cos(theta_fit)
    return (np.arccos(cx *cx_fit + cy * cy_fit + cz * cz_fit))*180/np.pi

deltas_bar=np.arange(0,20+0.5,0.5)

def graph(delta):
    delta_perc=[]
    for i in range(0,len(deltas_bar)-1):
        # delta_df=pd.DataFrame({'delta':delta_max})
        delta_perc.append(len(delta[ (delta>deltas_bar[i]) & (delta<deltas_bar[i+1])])*100/len(delta))
    return delta_perc


def results(name_files):
    delta_multiple=[]
    delta_total=[]
    d0=[]
    ratio=[]
    nbr_pmt=[]
    for files in range(0,len(name_files)):
        # print(name_files[files])
        result_angle=pd.read_csv(name_files[files], header=None, skiprows=[0],sep='\s+')  
        real_phi_analyse=real_phi[np.array(result_angle[7].astype('int')-1)]
        real_theta_analyse=real_theta[np.array(result_angle[7].astype('int')-1)]
        delta_max=f_delta(list(result_angle[0]),list(result_angle[1]),real_theta_analyse,real_phi_analyse)
        percent_delta=graph(delta_max)
    
        delta_total.append(delta_max)
        delta_multiple.append(percent_delta)
        nbr_pmt.append(result_angle[4])
        d0.append(result_angle[5])
        ratio.append(result_angle[6])
                
    return delta_total,delta_multiple,d0,ratio,nbr_pmt
 
    
def ratio_all(delta,ratio):
    delta_ratio=pd.DataFrame({'delta':delta,'ratio':ratio})
    delta_ratio['ratio']=delta_ratio['ratio']. round(0)

    moy_delta=delta_ratio.groupby(delta_ratio['ratio']).mean()
    std_delta=delta_ratio.groupby(delta_ratio['ratio']).std()
    
    x=list(moy_delta.index)
    mean=list(moy_delta['delta'])
    std=list(std_delta['delta'])
    
    return x,mean,std

#=====================================================================================
#=====================================================================================

def results_sph3_Fe(name_files):
    delta_multiple=[]
    delta_total=[]
    d0=[]
    ratio=[]
    for files in range(0,len(name_files)):
        result_angle=pd.read_csv(name_files[files], header=None, skiprows=[0],sep='\s+')  
        real_phi_analyse=real_phi_sph3_Fe[np.array(result_angle[7].astype('int'))]
        real_theta_analyse=[15/180*np.pi]*len(real_phi_analyse)
        
        delta_max=f_delta(list(result_angle[0]),list(result_angle[1]),real_theta_analyse,real_phi_analyse)
        
        percent_delta=graph(delta_max)
    
        delta_total.append(delta_max)
        delta_multiple.append(percent_delta)
        d0.append(result_angle[5])
        ratio.append(result_angle[6])
        
    return delta_total,delta_multiple,d0,ratio

def results_sph3_N(name_files):
    delta_multiple=[]
    delta_total=[]
    d0=[]
    ratio=[]
    for files in range(0,len(name_files)):
        result_angle=pd.read_csv(name_files[files], header=None, skiprows=[0],sep='\s+')  
        real_phi_analyse=real_phi_sph3_N[np.array(result_angle[7].astype('int'))]
        real_theta_analyse=[15/180*np.pi]*len(real_phi_analyse)
        
        delta_max=f_delta(list(result_angle[0]),list(result_angle[1]),real_theta_analyse,real_phi_analyse)
        
        percent_delta=graph(delta_max)
    
        delta_total.append(delta_max)
        delta_multiple.append(percent_delta)
        d0.append(result_angle[5])
        ratio.append(result_angle[6])
        
    return delta_total,delta_multiple,d0,ratio

def results_sph3_P(name_files):
    delta_multiple=[]
    delta_total=[]
    d0=[]
    ratio=[]
    for files in range(0,len(name_files)):
        result_angle=pd.read_csv(name_files[files], header=None, skiprows=[0],sep='\s+')  
        real_phi_analyse=real_phi_sph3_P[np.array(result_angle[7].astype('int'))]
        real_theta_analyse=[15/180*np.pi]*len(real_phi_analyse)
        
        delta_max=f_delta(list(result_angle[0]),list(result_angle[1]),real_theta_analyse,real_phi_analyse)
        
        percent_delta=graph(delta_max)
    
        delta_total.append(delta_max)
        delta_multiple.append(percent_delta)
        d0.append(result_angle[5])
        ratio.append(result_angle[6])
        
    return delta_total,delta_multiple,d0,ratio

###=============================================================
###=============================================================

#ONLY SIG

all_results_only_sig=results(name_files_only_sig)
delta_total,delta_multiple,d0=all_results_only_sig[0],all_results_only_sig[1],all_results_only_sig[2]

#SIG AND BG

all_results_sig_bg=results(name_files_sig_bg)
delta_total_sig_bg,delta_multiple_sig_bg,d0_bg,ratio_bg=all_results_sig_bg[0],all_results_sig_bg[1],all_results_sig_bg[2],all_results_sig_bg[3]

#ELECTRONIC

all_results_elec=results(name_files_elec)
delta_total_elec,delta_multiple_elec,d0_elec,ratio_elec,nbr_pmt_elec=all_results_elec[0],all_results_elec[1],all_results_elec[2],all_results_elec[3],all_results_elec[4]

#SPH3

all_results_sph3_N=results_sph3_N(name_files_sph3_N)
delta_total_sph3_N,delta_multiple_sph3_N,d0_sph3,ratio_sph3=all_results_sph3_N[0],all_results_sph3_N[1],all_results_sph3_N[2],all_results_sph3_N[3]

all_results_sph3_Fe=results_sph3_Fe(name_files_sph3_Fe)
delta_total_sph3_Fe,delta_multiple_sph3_Fe,d0_sph3,ratio_sph3=all_results_sph3_Fe[0],all_results_sph3_Fe[1],all_results_sph3_N[2],all_results_sph3_Fe[3]

all_results_sph3_P=results_sph3_P(name_files_sph3_P)
delta_total_sph3_P,delta_multiple_sph3_P,d0_sph3,ratio_sph3=all_results_sph3_P[0],all_results_sph3_P[1],all_results_sph3_P[2],all_results_sph3_P[3]




###=============================================================
###=============================================================

##==== Graphics

params = {"xtick.direction": "in", "ytick.direction": "in"}
plt.rcParams.update(params)

#===================

#ONLY SIG

#name_files_only_sig

#500 m 30 PeV

# N_500m_30PeV,p_500m_30PeV,Fe_500m_30PeV=-4,7,5
# percentage_save_500m_30=np.mean([len(delta_total[N_500m_30PeV])*100/6000,len(delta_total[p_500m_30PeV])*100/6000,len(delta_total[Fe_500m_30PeV])*100/6000])

# xx=deltas_bar[:10]+0.25
# plt.bar(xx,delta_multiple[Fe_500m_30PeV][:10],width=0.5,edgecolor="black",alpha=0.5,label='Fe,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total[Fe_500m_30PeV]),np.std(delta_total[Fe_500m_30PeV])),color='blue')
# plt.bar(xx,delta_multiple[N_500m_30PeV][:10],width=0.5,edgecolor="black",alpha=0.5,label='N,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total[N_500m_30PeV]),np.std(delta_total[N_500m_30PeV])),color='purple')
# plt.bar(xx,delta_multiple[p_500m_30PeV][:10],width=0.5,edgecolor="black",alpha=0.5,label='p,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total[p_500m_30PeV]),np.std(delta_total[p_500m_30PeV])),color='purple')
# plt.legend(title='500 м 30 ПэВ',title_fontsize=12)
# plt.xlabel('$\delta$, $^\circ$',fontsize=14)
# plt.ylabel('% от всей выборки',fontsize=14)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.ylim(0,35)
# plt.show()


# #500 m 10 PeV
# N_500m_10PeV,p_500m_10PeV,Fe_500m_10PeV=0,3,-2
# percentage_save_500m_10=np.mean([len(delta_total[N_500m_10PeV])*100/6000,len(delta_total[p_500m_10PeV])*100/6000,len(delta_total[Fe_500m_10PeV])*100/6000])

# xx=deltas_bar[:10]+0.25
# plt.bar(xx,delta_multiple[Fe_500m_10PeV][:10],width=0.5,edgecolor="black",alpha=0.5,label='Fe,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total[Fe_500m_10PeV]),np.std(delta_total[Fe_500m_10PeV])),color='blue')
# plt.bar(xx,delta_multiple[N_500m_10PeV][:10],width=0.5,edgecolor="black",alpha=0.5,label='N,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total[N_500m_10PeV]),np.std(delta_total[N_500m_10PeV])),color='purple')
# plt.bar(xx,delta_multiple[p_500m_10PeV][:10],width=0.5,edgecolor="black",alpha=0.5,label='p,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total[p_500m_10PeV]),np.std(delta_total[p_500m_10PeV])),color='red')
# plt.legend(title='500 м 10 ПэВ',title_fontsize=12)
# plt.xlabel('$\delta$, $^\circ$',fontsize=14)
# plt.ylabel('% от всей выборки',fontsize=14)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.ylim(0,40)
# plt.show()  
    

# #900 m 30 PeV

# N_900m_30PeV,p_900m_30PeV,Fe_900m_30PeV=1,2,-1
# percentage_save_900m_30=np.mean([len(delta_total[N_900m_30PeV])*100/6000,len(delta_total[p_900m_30PeV])*100/6000,len(delta_total[Fe_900m_30PeV])*100/6000])
# xx=deltas_bar[:10]+0.25
# plt.bar(xx,delta_multiple[Fe_900m_30PeV][:10],width=0.5,edgecolor="black",alpha=0.5,label='Fe,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total[Fe_900m_30PeV]),np.std(delta_total[Fe_900m_30PeV])),color='blue')
# plt.bar(xx,delta_multiple[N_900m_30PeV][:10],width=0.5,edgecolor="black",alpha=0.5,label='N,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total[N_900m_30PeV]),np.std(delta_total[N_900m_30PeV])),color='purple')
# plt.bar(xx,delta_multiple[p_900m_30PeV][:10],width=0.5,edgecolor="black",alpha=0.5,label='p,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total[p_900m_30PeV]),np.std(delta_total[p_900m_30PeV])),color='red')
# plt.legend(title='900 м 30 ПэВ',title_fontsize=12)
# plt.xlabel('$\delta$, $^\circ$',fontsize=14)
# plt.ylabel('% от всей выборки',fontsize=14)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.ylim(0,40)
# plt.show() 
    

# #900 m 10 PeV

# N_900m_10PeV,p_900m_10PeV,Fe_900m_10PeV=-3,-6,-8
# percentage_save_900m_10=np.mean([len(delta_total[N_900m_10PeV])*100/6000,len(delta_total[p_900m_10PeV])*100/6000,len(delta_total[Fe_900m_10PeV])*100/6000])
# cutting=20
# xx=deltas_bar[:cutting]+0.25
# plt.bar(xx,delta_multiple[Fe_900m_10PeV][:cutting],width=0.5,edgecolor="black",alpha=0.5,label='Fe,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total[Fe_900m_10PeV]),np.std(delta_total[Fe_900m_10PeV])),color='blue')
# plt.bar(xx,delta_multiple[N_900m_10PeV][:cutting],width=0.5,edgecolor="black",alpha=0.5,label='N,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total[N_900m_10PeV]),np.std(delta_total[N_900m_10PeV])),color='purple')
# plt.bar(xx,delta_multiple[p_900m_10PeV][:cutting],width=0.5,edgecolor="black",alpha=0.5,label='p,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total[p_900m_10PeV]),np.std(delta_total[p_900m_10PeV])),color='red')
# plt.legend(title='900 м 10 ПэВ',title_fontsize=12)
# plt.xlabel('$\delta$, $^\circ$',fontsize=14)
# plt.ylabel('% от всей выборки',fontsize=14)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.ylim(0,35)
# plt.show()  
    
      
#===================

    
#SIG & BG

#name_files_sig_bg

#500 m 30 PeV

N_500m_30PeV_bg,p_500m_30PeV_bg,Fe_500m_30PeV_bg=2,1,4
xx=deltas_bar[:10]+0.25
plt.bar(xx,delta_multiple_sig_bg[Fe_500m_30PeV_bg][:10],width=0.5,edgecolor="black",alpha=0.5,label='Fe,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_sig_bg[Fe_500m_30PeV_bg]),np.std(delta_total_sig_bg[Fe_500m_30PeV_bg])),color='blue')
plt.bar(xx,delta_multiple_sig_bg[N_500m_30PeV_bg][:10],width=0.5,edgecolor="black",alpha=0.5,label='N,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_sig_bg[N_500m_30PeV_bg]),np.std(delta_total_sig_bg[N_500m_30PeV_bg])),color='turquoise')
plt.bar(xx,delta_multiple_sig_bg[p_500m_30PeV_bg][:10],width=0.5,edgecolor="black",alpha=0.5,label='p,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_sig_bg[p_500m_30PeV_bg]),np.std(delta_total_sig_bg[p_500m_30PeV_bg])),color='green')
plt.legend(title='500 м 30 ПэВ',title_fontsize=12)
plt.xlabel('$\delta$, $^\circ$',fontsize=14)
plt.ylabel('% от всей выборки',fontsize=14)
plt.ylim(0,35)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()


# #900 m 10 PeV

N_900m_10PeV_bg,p_900m_10PeV_bg,Fe_900m_10PeV_bg=-3,0,-1
xx=deltas_bar[:20]+0.25
plt.bar(xx,delta_multiple_sig_bg[Fe_900m_10PeV_bg][:20],width=0.5,edgecolor="black",alpha=0.5,label='Fe,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_sig_bg[Fe_900m_10PeV_bg]),np.std(delta_total_sig_bg[Fe_900m_10PeV_bg])),color='blue')
plt.bar(xx,delta_multiple_sig_bg[N_900m_10PeV_bg][:20],width=0.5,edgecolor="black",alpha=0.5,label='N,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_sig_bg[N_900m_10PeV_bg]),np.std(delta_total_sig_bg[N_900m_10PeV_bg])),color='turquoise')
plt.bar(xx,delta_multiple_sig_bg[p_900m_10PeV_bg][:20],width=0.5,edgecolor="black",alpha=0.5,label='p,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_sig_bg[p_900m_10PeV_bg]),np.std(delta_total_sig_bg[p_900m_10PeV_bg])),color='green')
plt.legend(title='900 м 10 ПэВ',title_fontsize=12)
plt.xlabel('$\delta$, $^\circ$',fontsize=14)
plt.ylabel('% от всей выборки',fontsize=14)
plt.ylim(0,35)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show() 



#===================
    
#ELECTRONIC

#name_files_elec

#500 m 30 PeV

# N_500m_30PeV_elec,p_500m_30PeV_elec,Fe_500m_30PeV_elec=0,3,2 #before
# N_500m_30PeV_elec,p_500m_30PeV_elec,Fe_500m_30PeV_elec=2,1,-1 #lot
# percentage_save_500m_30_elec=np.mean([len(delta_total[N_500m_30PeV_elec])*100/6000,len(delta_total[p_500m_30PeV_elec])*100/6000,len(delta_total[Fe_500m_30PeV_elec])*100/6000])

# xx=deltas_bar[:10]+0.25
# plt.bar(xx,delta_multiple_elec[Fe_500m_30PeV_elec][:10],width=0.5,edgecolor="black",alpha=0.5,label='Fe,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_elec[Fe_500m_30PeV_elec]),np.std(delta_total_elec[Fe_500m_30PeV_elec])),color='orange')
# plt.bar(xx,delta_multiple_elec[N_500m_30PeV_elec][:10],width=0.5,edgecolor="black",alpha=0.5,label='N,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_elec[N_500m_30PeV_elec]),np.std(delta_total_elec[N_500m_30PeV_elec])),color='red')
# plt.bar(xx,delta_multiple_elec[p_500m_30PeV_elec][:10],width=0.5,edgecolor="black",alpha=0.5,label='p,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_elec[p_500m_30PeV_elec]),np.std(delta_total_elec[p_500m_30PeV_elec])),color='yellow')
# plt.legend(title='500 м 30 ПэВ',title_fontsize=12)
# plt.xlabel('$\delta$, $^\circ$',fontsize=14)
# plt.ylabel('% от всей выборки',fontsize=14)
# plt.ylim(0,35)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.show()


# # 900 m 10 PeV

# N_900m_10PeV_elec,p_900m_10PeV_elec,Fe_900m_10PeV_elec=-3,0,-2
# percentage_save_900m_10_elec=np.mean([len(delta_total[N_900m_10PeV_elec])*100/6000,len(delta_total[p_900m_10PeV_elec])*100/6000,len(delta_total[Fe_900m_10PeV_elec])*100/6000])
# xx=deltas_bar[:20]+0.25
# plt.bar(xx,delta_multiple_elec[Fe_900m_10PeV_elec][:20],width=0.5,edgecolor="black",alpha=0.5,label='Fe,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_elec[Fe_900m_10PeV_elec]),np.std(delta_total_elec[Fe_900m_10PeV_elec])),color='orange')
# plt.bar(xx,delta_multiple_elec[N_900m_10PeV_elec][:20],width=0.5,edgecolor="black",alpha=0.5,label='N,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_elec[N_900m_10PeV_elec]),np.std(delta_total_elec[N_900m_10PeV_elec])),color='red')
# plt.bar(xx,delta_multiple_elec[p_900m_10PeV_elec][:20],width=0.5,edgecolor="black",alpha=0.5,label='p,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_elec[p_900m_10PeV_elec]),np.std(delta_total_elec[p_900m_10PeV_elec])),color='yellow')
# plt.legend(title='900 м 10 ПэВ',title_fontsize=12)
# plt.xlabel('$\delta$, $^\circ$',fontsize=14)
# plt.ylabel('% от всей выборки',fontsize=14)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.ylim(0,35)
# plt.show() 


#УСИКИ

# plt.scatter([1+0.5,7+0.5,26+0.5],[np.mean(delta_total_sig_bg[p_900m_10PeV_bg]),np.mean(delta_total_sig_bg[N_900m_10PeV_bg]),np.mean(delta_total_sig_bg[Fe_900m_10PeV_bg])],label='Сигнал с фоном',color='darkorange')
# plt.errorbar([1+0.5,7+0.5,26+0.5],[np.mean(delta_total_sig_bg[p_900m_10PeV_bg]),np.mean(delta_total_sig_bg[N_900m_10PeV_bg]),np.mean(delta_total_sig_bg[Fe_900m_10PeV_bg])], yerr =[np.std(delta_total_sig_bg[p_900m_10PeV_bg]),np.std(delta_total_sig_bg[N_900m_10PeV_bg]),np.std(delta_total_sig_bg[Fe_900m_10PeV_bg])],fmt ='o',markersize=1, capsize=5,color='darkorange')

# plt.scatter([1,7,26],[np.mean(delta_total[p_900m_10PeV]),np.mean(delta_total[N_900m_10PeV]),np.mean(delta_total[Fe_900m_10PeV])],label='Чистый сигнал',color='blue')
# plt.errorbar([1,7,26],[np.mean(delta_total[p_900m_10PeV]),np.mean(delta_total[N_900m_10PeV]),np.mean(delta_total[Fe_900m_10PeV])], yerr =[np.std(delta_total[p_900m_10PeV]),np.std(delta_total[N_900m_10PeV]),np.std(delta_total[Fe_900m_10PeV])],fmt ='o',markersize=1, capsize=5,color='blue')

# plt.ylabel('$\delta$, $^\circ$')
# plt.xlabel('Z')
# plt.legend(title='900 м 10 ПэВ')
# plt.show() 


#===================
    
#SPHERE-3

##name_files_sph3_N
##name_files_sph3_P
##name_files_sph3_Fe


# #500 m 10 PeV

# N_500m_10PeV_sph3,p_500m_10PeV_sph3,Fe_500m_10PeV_sph3=0,1,1
# cutting=10
# xx=deltas_bar[:cutting]+0.25
# plt.bar(xx,delta_multiple_sph3_Fe[Fe_500m_10PeV_sph3][:cutting],width=0.5,edgecolor="black",alpha=0.5,label='Fe,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_sph3_Fe[Fe_500m_10PeV_sph3]),np.std(delta_total_sph3_Fe[Fe_500m_10PeV_sph3])),color='greenyellow')
# plt.bar(xx,delta_multiple_sph3_N[N_500m_10PeV_sph3][:cutting],width=0.5,edgecolor="black",alpha=0.5,label='N,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_sph3_N[N_500m_10PeV_sph3]),np.std(delta_total_sph3_N[N_500m_10PeV_sph3])),color='blue')
# plt.bar(xx,delta_multiple_sph3_P[p_500m_10PeV_sph3][:cutting],width=0.5,edgecolor="black",alpha=0.5,label='p,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_sph3_P[p_500m_10PeV_sph3]),np.std(delta_total_sph3_P[p_500m_10PeV_sph3])),color='darkturquoise')
# plt.legend(title='SPHERE-3 \n500 м 10 ПэВ')
# plt.xlabel('$\delta$, $^\circ$',fontsize=14)
# plt.ylabel('% от всей выборки',fontsize=14)
# plt.ylim(0,35)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.show()

# #1000 m 10 PeV

# N_500m_10PeV_sph3,p_500m_10PeV_sph3,Fe_500m_10PeV_sph3=1,0,0
# cutting=10
# xx=deltas_bar[:cutting]+0.25
# plt.bar(xx,delta_multiple_sph3_Fe[Fe_500m_10PeV_sph3][:cutting],width=0.5,edgecolor="black",alpha=0.5,label='Fe,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_sph3_Fe[Fe_500m_10PeV_sph3]),np.std(delta_total_sph3_Fe[Fe_500m_10PeV_sph3])),color='greenyellow')
# plt.bar(xx,delta_multiple_sph3_N[N_500m_10PeV_sph3][:cutting],width=0.5,edgecolor="black",alpha=0.5,label='N,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_sph3_N[N_500m_10PeV_sph3]),np.std(delta_total_sph3_N[N_500m_10PeV_sph3])),color='blue')
# plt.bar(xx,delta_multiple_sph3_P[p_500m_10PeV_sph3][:cutting],width=0.5,edgecolor="black",alpha=0.5,label='p,  <$\delta$>={:.2f}±{:.2f}$^\circ$'.format(np.mean(delta_total_sph3_P[p_500m_10PeV_sph3]),np.std(delta_total_sph3_P[p_500m_10PeV_sph3])),color='darkturquoise')
# plt.legend(title='SPHERE-3 \n1000 м 10 ПэВ')
# plt.xlabel('$\delta$, $^\circ$',fontsize=14)
# plt.ylabel('% от всей выборки',fontsize=14)
# plt.ylim(0,35)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.show()
