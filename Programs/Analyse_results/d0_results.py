import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

###=============================================================
###=============================================================

## !!! CHANGE ONLY HERE

# nuclei='Fe'
# nuclei='N'
nuclei='P'

# Choose the altitude
H=500 #m, 900 

# Choose the energy
En=30 #PeV, 10

integral=0.5

# Choose the analyse
analyse='only_sig' #elec, bg
# analyse='electronic' #elec, bg
telescope='sph2'

## !!! END CHANGE ONLY HERE

###=============================================================
###=============================================================

#==================================================================
if telescope!='sph3':
    all_val_x0_y0=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/d0/{}/{}_d0_{}m_{}PeV_{}_pmt.csv'.format(analyse,analyse,H,En,nuclei),header=None,sep='\s+')
    vision=180
    in_legend='СФЕРА-2'
    if H==500:
        vision=180
        max_x_axis=180
    elif H==900:
        vision=310
        max_x_axis=330

if telescope=='sph3':
    all_val_x0_y0=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/d0/{}/{}_d0_{}m_{}PeV_{}_sph3_pmt.csv'.format(telescope,analyse,H,En,nuclei),header=None,sep='\s+')
    in_legend='СФЕРА-3'
    if H==500:
        vision=165
        max_x_axis=180
    elif H==1000:
        vision=330
        max_x_axis=330



if analyse=='electronic':
    
    files=list(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/files/{}/{}_files_{}m_{}PeV_{}_pmt.csv'.format(analyse,analyse,H,En,nuclei),header=None)[0])
    order_elec=[ int(files[i][-8:-5] + files[i][-3:]) for i in range(0,len(files))]
    all_val_x0_y0['order_elec']=order_elec
    all_val_x0_y0=all_val_x0_y0.sort_values(by=['order_elec'])

    real_d0=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/d0/{}/d0_real_{}m_{}PeV_{}.csv'.format(analyse,H,En,nuclei),header=None,sep='\s+')
    real_d0=real_d0.sort_values(by=[2])

    all_val_x0_y0[0]=real_d0[0]
    all_val_x0_y0[1]=real_d0[1]
    
    in_legend='Электронный сигнал'
    
# elif analyse=='only_sig' and telescope!='sph3':
#     in_legend='Чистый сигнал'

    
#Extract values
  
x0_real=np.array(all_val_x0_y0[0])
y0_real=np.array(all_val_x0_y0[1])

x0_max=np.array(all_val_x0_y0[2])
y0_max=np.array(all_val_x0_y0[3])

x0_grav_center=np.array(all_val_x0_y0[4]) 
y0_grav_center=np.array(all_val_x0_y0[5]) 


#ERROR


def err_d0(results):
    d0_real=np.sqrt(x0_real**2+y0_real**2)
    d0_max=np.sqrt(x0_max**2+y0_max**2)
    d0_grav_center=np.sqrt(x0_grav_center**2+y0_grav_center**2)

    delta_d0_max=abs(d0_real-d0_max)
    delta_d0_grav_center=abs(d0_real-d0_grav_center)
    
    d0=pd.DataFrame({'d0_real':d0_real,'d0_max':d0_max,'d0_grav_center':d0_grav_center,'delta_max':delta_d0_max,'delta_grav':delta_d0_grav_center}).sort_values(by='d0_real')
    lenght=np.arange(0,vision+10,10)
    x_d0_grav=[]
    x_d0_max=[]
    delta_grav_mean=[]
    delta_grav_std=[]
    delta_max_mean=[]
    delta_max_std=[]
    for i in range(0,len(lenght)-1):
        delta_max=list(d0['delta_max'][ (d0['d0_real']>lenght[i]) & (d0['d0_real']<lenght[i+1]) ])
        delta_grav=list(d0['delta_grav'][ (d0['d0_real']>lenght[i]) & (d0['d0_real']<lenght[i+1]) ])
        if len(delta_grav)>2:
            x_d0_grav.append(lenght[i])
            delta_grav_mean.append(np.mean(delta_grav))
            delta_grav_std.append(np.std(delta_grav))
        if len(delta_max)>2:
            x_d0_max.append(lenght[i])    
            delta_max_mean.append(np.mean(delta_max))
            delta_max_std.append(np.std(delta_max))

    return delta_grav_mean, delta_max_mean, delta_grav_std,delta_max_std,x_d0_grav,x_d0_max

#============

results_d0=err_d0(all_val_x0_y0)
delta_grav_mean,delta_max_mean,delta_grav_std,delta_max_std,x_d0_grav,x_d0_max=results_d0[0],results_d0[1],results_d0[2],results_d0[3],results_d0[4],results_d0[5]

#Normalize capsize to min =0

verify_if_below_0_max=np.array(delta_max_mean)-np.array(delta_max_std)
verify_if_below_0_grav=np.array(delta_grav_mean)-np.array(delta_grav_std)

df_verify_below_0=pd.DataFrame({'delta_max_std':delta_max_std,'verify_if_below_0_max':verify_if_below_0_max,'delta_grav_std':delta_grav_std,'verify_if_below_0_grav':verify_if_below_0_grav})

#=====================================================================

plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in"})


plt.scatter(x_d0_max,delta_max_mean,color='red',label='По макс. ФЭУ  <{:.0f}>±{:.0f}м'.format(np.mean(delta_max_mean),np.mean(delta_max_std)), edgecolor ='black',s=15)
plt.scatter(x_d0_grav,delta_grav_mean,color='green',label='По центру тяжести <{:.0f}>±{:.0f}м'.format(np.mean(delta_grav_mean),np.mean(delta_grav_std)), edgecolor ='black',s=15)

plt.errorbar(x_d0_max,delta_max_mean, yerr = delta_max_std,fmt ='o',markersize=3, capsize=4,color='red')
plt.errorbar(x_d0_grav,delta_grav_mean, yerr = delta_grav_std,fmt ='o',markersize=3, capsize=4,color='green')
plt.ylabel(' < |∆d0| >, м',fontsize=14)
plt.xlabel('real d0, м',fontsize=14)
plt.legend(title='{}, {}ПэВ, {} м, {}'.format(nuclei,En,H,in_legend), loc='upper center',fontsize=12,title_fontsize=14)
plt.xlim(-5,max_x_axis)
plt.ylim(0,60)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()


