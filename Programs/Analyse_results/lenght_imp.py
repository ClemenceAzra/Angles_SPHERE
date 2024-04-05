# считает длину импульса всех файлов



import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob

###=============================================================

## Choose the nuclei

# nuclei='Fe'
# nuclei='N'
nuclei='P'

# Choose the altitude
H=500 #m, 900 

# Choose the energy
En=10 #PeV, 10

integral=0.5

# Choose the analyse
analyse='only_sig' #elec, bg
# analyse='electronic' #elec, bg
telescope='sph3'


###=============================================================
###=============================================================


#OPEN MULTI FILES

#SPHERE-3 only sig
if telescope=='sph3':
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/lenght_impulse/sph3'

elif telescope=='sph2':
    if analyse=='only_sig':
        path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/lenght_impulse/only_sig'
    elif analyse=='electronic':
        path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/lenght_impulse/electronic'


##=== Open files
name_files = glob.glob(os.path.join(path,'*'))

N_res=[]
t_res=[]
total_lenght=[]
for files in range(0,len(name_files)):
    print(name_files[files])
    result_lenght=pd.read_csv(name_files[files], header=None, skiprows=[0],sep='\s+')  
    all_res=list(abs(result_lenght[0]))
    bins_same_lenght=np.histogram(all_res,np.arange(0,1700,50))
    
    N_res.append(bins_same_lenght[0])
    t_res.append(bins_same_lenght[1:][0][1:])
    total_lenght.append(all_res)


#_______________________________________

#SPHERE-2 only sig

# plt.bar(t_res[0],N_res[0]*100/len(total_lenght[0]),alpha=0.6,label='500 м 10 ПэВ p <t>={:.0f} нс'.format(np.mean(total_lenght[0])),width=50)
# plt.bar(t_res[3],N_res[3]*100/len(total_lenght[3]),alpha=0.6,label='500 м 30 ПэВ p <t>={:.0f} нс'.format(np.mean(total_lenght[3])),width=50)
# plt.bar(t_res[2],N_res[2]*100/len(total_lenght[2]),alpha=0.6,label='500 м 30 ПэВ Fe <t>={:.0f} нс'.format(np.mean(total_lenght[2])),width=50)

# plt.bar(t_res[1],N_res[1]*100/len(total_lenght[1]),alpha=0.6,label='900 м 10 ПэВ p <t>={:.0f} нс'.format(np.mean(total_lenght[1])),width=50)
# plt.bar(t_res[4],N_res[4]*100/len(total_lenght[4]),alpha=0.6,label='900 м 30 ПэВ p <t>={:.0f} нс'.format(np.mean(total_lenght[4])),width=50)

#_______________________________________

#SPHERE-3 only sig
if telescope=='sph2' and analyse=='electronic': #elec, bg

    # 500 m
    N_500,P_500,Fe_500= 1,-2,-3
    plt.bar(t_res[P_500],N_res[P_500]*100/len(total_lenght[P_500]),alpha=0.6,label='500 м 30 ПэВ p  <t>={:.0f} нс'.format(np.mean(total_lenght[P_500])),width=50,color='aqua')
    plt.bar(t_res[N_500],N_res[N_500]*100/len(total_lenght[N_500]),alpha=0.6,label='500 м 30 ПэВ N  <t>={:.0f} нс'.format(np.mean(total_lenght[N_500])),width=50,color='green')
    plt.bar(t_res[Fe_500],N_res[Fe_500]*100/len(total_lenght[Fe_500]),alpha=0.6,label='500 м 30 ПэВ Fe <t>={:.0f} нс'.format(np.mean(total_lenght[Fe_500])),width=50,color='dodgerblue')
    
    #1000 m
    N_900,P_900,Fe_900=0,-1,2
    plt.bar(t_res[P_500],N_res[P_900]*100/len(total_lenght[P_900]),alpha=0.6,label='900 м 10 ПэВ p  <t>={:.0f} нс'.format(np.mean(total_lenght[P_900])),width=50,color='salmon')
    plt.bar(t_res[N_500],N_res[N_900]*100/len(total_lenght[N_900]),alpha=0.6,label='900 м 10 ПэВ N  <t>={:.0f} нс'.format(np.mean(total_lenght[N_900])),width=50,color='crimson')
    plt.bar(t_res[Fe_500],N_res[Fe_900]*100/len(total_lenght[Fe_900]),alpha=0.6,label='900 м 10 ПэВ Fe <t>={:.0f} нс'.format(np.mean(total_lenght[Fe_900])),width=50,color='purple')
    
    plt.xlabel('Длина импульса, нс',fontsize=14)
    plt.legend(title='СФЕРА-2, электронный сигнал',fontsize=8,title_fontsize=11)
    plt.ylabel('% от выборки',fontsize=14)
    plt.ylim(0,45)
    plt.xlim(0,1600)
    # plt.savefig('foo.pdf')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()


if telescope=='sph2' and analyse=='only_sig': #elec, bg

    # 500 m
    N_500,P_500,Fe_500=-1,-4,-7
    plt.bar(t_res[P_500]-50,N_res[P_500]*100/len(total_lenght[P_500]),alpha=0.5,label='500 м 30 ПэВ p  <t>={:.0f} нс'.format(np.mean(total_lenght[P_500])),width=50,color='aqua')
    # plt.bar(t_res[dop_500]-50,N_res[dop_500]*100/len(total_lenght[dop_500]),alpha=0.5,label='500 м 10 ПэВ p  <t>={:.0f} нс'.format(np.mean(total_lenght[dop_500])),width=50,color='orange')
    plt.bar(t_res[N_500]-50,N_res[N_500]*100/len(total_lenght[N_500]),alpha=0.5,label='500 м 30 ПэВ N  <t>={:.0f} нс'.format(np.mean(total_lenght[N_500])),width=50,color='green')
    plt.bar(t_res[Fe_500]-50,N_res[Fe_500]*100/len(total_lenght[Fe_500]),alpha=0.5,label='500 м 30 ПэВ Fe <t>={:.0f} нс'.format(np.mean(total_lenght[Fe_500])),width=50,color='dodgerblue')

    #1000 m
    N_900,P_900,Fe_900=-2,-3,-8
    plt.bar(t_res[P_900]-50,N_res[P_900]*100/len(total_lenght[P_900]),alpha=0.5,label='900 м 10 ПэВ p  <t>={:.0f} нс'.format(np.mean(total_lenght[P_900])),width=50,color='salmon')
    # plt.bar(t_res[dop_900]-50,N_res[dop_900]*100/len(total_lenght[dop_900]),alpha=0.5,label='900 м 30 ПэВ p  <t>={:.0f} нс'.format(np.mean(total_lenght[dop_900])),width=50,color='orange')
    plt.bar(t_res[N_900]-50,N_res[N_900]*100/len(total_lenght[N_900]),alpha=0.5,label='900 м 10 ПэВ N  <t>={:.0f} нс'.format(np.mean(total_lenght[N_900])),width=50,color='crimson')
    plt.bar(t_res[Fe_900]-50,N_res[Fe_900]*100/len(total_lenght[Fe_900]),alpha=0.5,label='900 м 10 ПэВ Fe <t>={:.0f} нс'.format(np.mean(total_lenght[Fe_900])),width=50,color='purple')
    
    plt.xlabel('Длина импульса, нс',fontsize=14)
    plt.ylabel('% от выборки',fontsize=14)
    plt.legend(title='СФРЕА-2, чистый сигнал',fontsize=8,title_fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylim(0,45)
    plt.xlim(0,1600)
    # plt.savefig('foo.pdf')
    plt.show()
    
    N_500,P_500,Fe_500=-1,-4,-7
    plt.bar(t_res[P_500]-50,N_res[P_500]*100/len(total_lenght[P_500]),alpha=0.5,label='500 м 30 ПэВ p  <t>={:.0f} нс'.format(np.mean(total_lenght[P_500])),width=50,color='aqua')
    plt.bar(t_res[N_500]-50,N_res[N_500]*100/len(total_lenght[N_500]),alpha=0.5,label='500 м 30 ПэВ N  <t>={:.0f} нс'.format(np.mean(total_lenght[N_500])),width=50,color='green')
    plt.bar(t_res[Fe_500]-50,N_res[Fe_500]*100/len(total_lenght[Fe_500]),alpha=0.5,label='500 м 30 ПэВ Fe <t>={:.0f} нс'.format(np.mean(total_lenght[Fe_500])),width=50,color='dodgerblue')

    #1000 m
    N_900,P_900,Fe_900=1,6,2
    plt.bar(t_res[P_900]-50,N_res[P_900]*100/len(total_lenght[P_900]),alpha=0.5,label='500 м 10 ПэВ p  <t>={:.0f} нс'.format(np.mean(total_lenght[P_900])),width=50,color='orange')
    # plt.bar(t_res[dop_900]-50,N_res[dop_900]*100/len(total_lenght[dop_900]),alpha=0.5,label='900 м 30 ПэВ p  <t>={:.0f} нс'.format(np.mean(total_lenght[dop_900])),width=50,color='orange')
    plt.bar(t_res[N_900]-50,N_res[N_900]*100/len(total_lenght[N_900]),alpha=0.5,label='500 м 10 ПэВ N  <t>={:.0f} нс'.format(np.mean(total_lenght[N_900])),width=50,color='red')
    plt.bar(t_res[Fe_900]-50,N_res[Fe_900]*100/len(total_lenght[Fe_900]),alpha=0.5,label='500 м 10 ПэВ Fe <t>={:.0f} нс'.format(np.mean(total_lenght[Fe_900])),width=50,color='yellow')
    
    plt.xlabel('Длина импульса, нс')
    plt.legend(title='СФРЕА-2, чистый сигнал',fontsize=7,title_fontsize=7)
    plt.ylabel('% от выборки')
    plt.ylim(0,45)
    plt.xlim(0,1600)
    # plt.savefig('foo.pdf')
    plt.show()

if telescope=='sph3': #elec, bg

    # 500 m
    N_500,P_500,Fe_500=-1,-2,0
    plt.bar(t_res[P_500]-50,N_res[P_500]*100/len(total_lenght[P_500]),alpha=0.5,label='500 м 10 ПэВ p  <t>={:.0f} нс'.format(np.mean(total_lenght[P_500])),width=50,color='aqua')
    plt.bar(t_res[N_500]-50,N_res[N_500]*100/len(total_lenght[N_500]),alpha=0.5,label='500 м 10 ПэВ N  <t>={:.0f} нс'.format(np.mean(total_lenght[N_500])),width=50,color='green')
    plt.bar(t_res[Fe_500]-50,N_res[Fe_500]*100/len(total_lenght[Fe_500]),alpha=0.5,label='500 м 10 ПэВ Fe <t>={:.0f} нс'.format(np.mean(total_lenght[Fe_500])),width=50,color='dodgerblue')

    #1000 m
    N_500,P_500,Fe_500=2,-3,1
    plt.bar(t_res[P_500]-50,N_res[P_500]*100/len(total_lenght[P_500]),alpha=0.5,label='1000 м 10 ПэВ p  <t>={:.0f} нс'.format(np.mean(total_lenght[P_500])),width=50,color='salmon')
    plt.bar(t_res[N_500]-50,N_res[N_500]*100/len(total_lenght[N_500]),alpha=0.5,label='1000 м 10 ПэВ N  <t>={:.0f} нс'.format(np.mean(total_lenght[N_500])),width=50,color='crimson')
    plt.bar(t_res[Fe_500]-50,N_res[Fe_500]*100/len(total_lenght[Fe_500]),alpha=0.5,label='1000 м 10 ПэВ Fe <t>={:.0f} нс'.format(np.mean(total_lenght[Fe_500])),width=50,color='purple')

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
   
    plt.xlabel('Длина импульса, нс',fontsize=14)
    plt.legend(title='СФЕРА-3, чистый сигнал',fontsize=8,title_fontsize=11)
    plt.ylabel('% от выборки',fontsize=14)
    # plt.savefig('foo.pdf')
  
    plt.show()
    
    
# kualaloumpour_elec=[np.mean(total_lenght[P_500]),np.mean(total_lenght[N_500]),np.mean(total_lenght[Fe_500])]
# kualaloumpour_elec_900=[np.mean(total_lenght[P_900]),np.mean(total_lenght[N_900]),np.mean(total_lenght[Fe_900])]

# plt.scatter([1,7,26],[np.mean(total_lenght[P_500]),np.mean(total_lenght[N_500]),np.mean(total_lenght[Fe_500])],label='500 м 30 ПэВ',color='aqua' ) 
# plt.plot([1,7,26],[np.mean(total_lenght[P_500]),np.mean(total_lenght[N_500]),np.mean(total_lenght[Fe_500])],label='Чистый сигнал',color='aqua' ,linestyle='--') 
# plt.errorbar([1,7,26],[np.mean(total_lenght[P_500]),np.mean(total_lenght[N_500]),np.mean(total_lenght[Fe_500])], yerr = [np.std(total_lenght[P_500]),np.std(total_lenght[N_500]),np.std(total_lenght[Fe_500])],fmt ='o',markersize=3, capsize=4,color='aqua')


# plt.scatter([1,7,26],[np.mean(total_lenght[P_900]),np.mean(total_lenght[N_900]),np.mean(total_lenght[Fe_900])],label='900 м 10 ПэВ',color='darkred') 
# plt.plot([1,7,26],[np.mean(total_lenght[P_900]),np.mean(total_lenght[N_900]),np.mean(total_lenght[Fe_900])],label='Чистый сигнал',color='darkred',linestyle='--') 
# # plt.errorbar(x_d0_max,delta_max_mean, yerr = delta_max_std,fmt ='o',markersize=3, capsize=4,color='red')

# plt.scatter([1,7,26],kualaloumpour_elec,label='500 м 30 ПэВ',color='aqua' ) 
# plt.plot([1,7,26],kualaloumpour_elec,label='Электронный сигнал',color='aqua' ,linestyle=':') 
# plt.scatter([1,7,26],kualaloumpour_elec_900[0],label='900 м 10 ПэВ',color='darkred') 
# plt.plot([1,7,26],kualaloumpour_elec_900[0],label='Электронный сигнал',color='darkred',linestyle=':') 

# plt.legend(fontsize=8,loc='upper center')
# plt.ylabel('<t>, нс')
# plt.xlabel('Z')
# plt.show()
    
    
    

