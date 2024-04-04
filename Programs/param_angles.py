import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

###=============================================================

## Choose the nuclei

# nuclei='Fe'
# nuclei='N'
nuclei='P'

# Choose the altitude
H=500 #m, 900 

# Choose the energy
En=30 #PeV, 10

# Choose the analyse
analyse='only_sig' #elec, bg

###=============================================================
###=============================================================

##=== Open files

results=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/param_angles/only_sig_param_t_500m_30PeV_P.csv'.format(En,analyse,H,En,nuclei),header=None,sep='\s+')

files=list(pd.read_csv('//Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/param_angles/only_sig_files_t_500m_30PeV_P.csv'.format(En,analyse,H,En,nuclei),header=None,sep='\s+')[0]) 
        
# theta_real=15/180*np.pi
real_theta=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_theta',header=None,sep='\s+')[0]) #number of PMT
real_phi=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_phi',header=None,sep='\s+')[0]) #number of PMT

theta_real=[real_theta[int(files[i][-8:-5])-1] for i in range(0,len(files))]#10-20º P, N
phi_real=[real_phi[int(files[i][-8:-5])-1] for i in range(0,len(files))]#10-20º P, N

###=============================================================
###=============================================================

results['theta_r'], results['phi_r']=theta_real, phi_real

results=results.sort_values(by='theta_r')
results['theta_r']=results['theta_r']*180/np.pi
moy_param=results.groupby(results['theta_r']).mean()

moy_param['max_a']=list(results.groupby(results['theta_r']).max()[0])
moy_param['min_a']=list(results.groupby(results['theta_r']).min()[0])

moy_param['max_b']=list(results.groupby(results['theta_r']).max()[1])
moy_param['min_b']=list(results.groupby(results['theta_r']).min()[1])

moy_param['max_c']=list(results.groupby(results['theta_r']).max()[2])
moy_param['min_c']=list(results.groupby(results['theta_r']).min()[2])

#===========================

x=np.arange(-100,100)
approx_10deg=moy_param.iloc[0][0]*x**2+moy_param.iloc[0][1]*x+moy_param.iloc[0][2]
approx_12_5deg=moy_param.iloc[19][0]*x**2+moy_param.iloc[19][1]*x+moy_param.iloc[19][2]
approx_15deg=moy_param.iloc[30][0]*x**2+moy_param.iloc[30][1]*x+moy_param.iloc[30][2]
approx_17_5deg=moy_param.iloc[41][0]*x**2+moy_param.iloc[41][1]*x+moy_param.iloc[41][2]
approx_19_5deg=moy_param.iloc[-1][0]*x**2+moy_param.iloc[-1][1]*x+moy_param.iloc[-1][2]


###=============================================================
###=============================================================


##==== Graphics

params = {"xtick.direction": "in", "ytick.direction": "in"}
plt.rcParams.update(params)


plt.plot(x,approx_10deg,label='~10º,    a={:.4f}, b={:.2f}, c={:.2f}'.format(moy_param.iloc[0][0],moy_param.iloc[0][1],moy_param.iloc[0][2]))
plt.plot(x,approx_12_5deg,label='~12.5º, a={:.4f}, b={:.2f}, c={:.2f}'.format(moy_param.iloc[19][0],moy_param.iloc[19][1],moy_param.iloc[19][2]))
plt.plot(x,approx_15deg,label='~15º,    a={:.4f}, b={:.2f}, c={:.2f}'.format(moy_param.iloc[30][0],moy_param.iloc[30][1],moy_param.iloc[30][2]))
plt.plot(x,approx_17_5deg,label='~17.5º, a={:.4f}, b={:.2f}, c={:.2f}'.format(moy_param.iloc[41][0],moy_param.iloc[41][1],moy_param.iloc[41][2]))
plt.plot(x,approx_19_5deg,label='~19.5º, a={:.4f}, b={:.2f}, c={:.2f}'.format(moy_param.iloc[-1][0],moy_param.iloc[-1][1],moy_param.iloc[-1][2]))
plt.xlabel('d, m')
plt.ylabel('t, ns')
plt.legend(title='err a = ±0.0012 / err b = ±0.15 / err c = ±10')
plt.show()

