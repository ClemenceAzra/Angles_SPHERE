import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

###=============================================================
###=============================================================

## !!! CHANGE ONLY HERE

# nuclei='Fe'
# nuclei='N'
nuclei='Fe'

# Choose the altitude
H=900 #m, 900 

# Choose the energy
En=10 #PeV, 10

integral=0.5

# Choose the analyse
analyse='only_sig' #elec, bg
analyse='electronic' #elec, bg
# analyse='sph3'
telescope='sph2'

## !!! END CHANGE ONLY HERE

###=============================================================
###=============================================================

##=== Open files

#Real angles
if analyse!='sph3':

    real_theta=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_theta',header=None,sep='\s+')[0]) #number of PMT
    real_phi=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_phi',header=None,sep='\s+')[0]) #number of PMT
    
    #Obtain results
    
    results=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/angles/{}/{}_total_{}m_{}PeV_{}_pmt.csv'.format(analyse,analyse,H,En,nuclei),header=None,sep='\s+')
    files=list(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/files/{}/{}_files_{}m_{}PeV_{}_pmt.csv'.format(analyse,analyse,H,En,nuclei),header=None)[0])
    
    theta=results[0]
    phi=results[1]
    
    nmbr_pmt=results[4]
    d0=results[5]
    ratio=results[6]
    lenght=results[6]
    
    max_impuse_summed=results[9]


#SPH 3
if analyse=='sph3':

    results=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/angles/sph3/{}/only_sig_total_{}m_{}PeV_{}_sph3_pmt.csv'.format(nuclei,H,En,nuclei),header=None,sep='\s+')
    files=list(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/ALL_RESULTS/files/sph3/{}/only_sig_files_{}m_{}PeV_{}_sph3_pmt.csv'.format(nuclei,H,En,nuclei),header=None)[0])

    
    theta=results[0]
    phi=results[1]
    
    nmbr_pmt=results[4]
    d0=results[5]
    max_impuse_summed=results[9]
    lenght=results[6]
    
    # theta_real=15/180*np.pi
    real_phi_sph3_Fe=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/angels_phi_{}'.format(nuclei),header=None,sep='\s+')[0])
    theta_real=15*np.pi/180


###=============================================================
###=============================================================

##=== Calculation of errors

#Real angles for each obtain files
if analyse!='sph3':  
    theta_real=[real_theta[int(files[i][-8:-5])-1] for i in range(0,len(files))]
    phi_real=[real_phi[int(files[i][-8:-5])-1] for i in range(0,len(files))]

elif analyse=='sph3':  
    phi_real=[real_phi_sph3_Fe[int(files[i][-8:-5])] for i in range(0,len(files))]


#Functions

#DELTA

def f_delta(theta_fit,phi_fit,theta_real,phi_real):
    cx = np.sin(theta_real) * np.cos(phi_real) #theta - истинные углы
    cy = np.sin(theta_real) * np.sin(phi_real)
    cz = np.cos(theta_real)
    cx_fit = np.sin(theta_fit) * np.cos(phi_fit) #theta' -полученные углы моделированием
    cy_fit = np.sin(theta_fit) * np.sin(phi_fit)
    cz_fit = np.cos(theta_fit)
    return (np.arccos(cx *cx_fit + cy * cy_fit + cz * cz_fit))*180/np.pi

#DEPENDANCES

def analyse_dep(results,colonne,step):
    table=results[[colonne,'delta_max']] #dataframe of nbr PMT (col 1) and delta (col 2)
    lenght=np.arange(results[colonne].min(),results[colonne].max()+5,step)
    x_step=[]
    y_mean=[]
    y_std=[]
    for i in range(0,len(lenght)-1):
        y=table[ (table[colonne]>lenght[i]) & (table[colonne]<lenght[i+1]) ]['delta_max']
        if len(y)>2:
            x_step.append(lenght[i])
            y_mean.append(y.mean())
            y_std.append(table[ (table[colonne]>lenght[i]) & (table[colonne]<lenght[i+1]) ]['delta_max'].std())
    return  x_step, y_mean, y_std

#==========

delta_max=f_delta(theta,phi,theta_real,phi_real)

print(np.mean(delta_max))
print(nuclei)

#=================================
#=================================


results['delta_max']=delta_max

results['real_th']=np.array(theta_real)*180/np.pi

results_limits=results[ (results[5]>100) & (results['delta_max']>15)]

bad_files=np.array(files)[np.where(delta_max>10)]

perc_below_3=100-len(bad_files)*100/len(files)

# np.savetxt('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/bad_files_elec/900m_{}_10PeV'.format(nuclei),bad_files,fmt='%s')

###=============================================================
###=============================================================

## RESULTS OF FUNCTION

#Dependance: delta(nbr PMT)

delta_pmt=analyse_dep(results,4,2)
x_delta_pmt,y_delta_pmt,y_delta_pmt_std=delta_pmt[0],delta_pmt[1],delta_pmt[2]

#Dependance: delta(lenght of the signal)
if telescope=='sph3':
    nbr=6
if H==900 and En==10 and analyse=='electronic':
    nbr=10
else:
    nbr=5

delta_lenght=analyse_dep(results,nbr,50)
x_delta_lenght,y_delta_lenght,y_delta_lenght_std=delta_lenght[0],delta_lenght[1],delta_lenght[2]

#Dependance: delta(ratio)
if analyse!='only_sig' and analyse!='sph3':
    delta_ratio=analyse_dep(results,6,2)
    x_delta_ratio,y_delta_ratio,y_delta_ratio_std=delta_ratio[0],delta_ratio[1],delta_ratio[2]

#Dependance: delta(d0)

delta_d0=analyse_dep(results,8,20)
x_delta_d0,y_delta_d0,y_delta_d0_std=delta_d0[0],delta_d0[1],delta_d0[2]

#Dependance:(max sig)

delta_max_sig=analyse_dep(results,9,100)
x_delta_max_sig,y_delta_max_sig,y_delta_max_sig_std=delta_max_sig[0],delta_max_sig[1],delta_max_sig[2]

#Dependance:(real theta)

delta_theta=analyse_dep(results,'real_th',1)
x_delta_theta,y_delta_theta,y_delta_theta_std=delta_theta[0],delta_theta[1],delta_theta[2]



###=============================================================
###=============================================================

##==== GRAPHICS

params = {"xtick.direction": "in", "ytick.direction": "in"} #ticks in box
plt.rcParams.update(params)

####===== 


# plt.hist(delta_max,label='{} \n< $\delta$ > ={:.2f} \n$\delta$ max={:.2f}º'.format(nuclei,np.mean(delta_max),max(delta_max)))
# plt.xlabel('$\delta$')
# plt.legend(title='H={} m, E={} PeV  \nby max'.format(H,En))
# # plt.yscale('log')
# plt.show()

plt.scatter(nmbr_pmt,delta_max,s=1)
plt.ylabel('$\delta$')
plt.xlabel('Число ФЭУ')
plt.show()


# plt.scatter(d0,delta_max,s=1)
# plt.ylabel('$\delta$, º')
# plt.xlabel('d0, м')
# plt.show()

# plt.scatter(ratio,delta_max,s=1)
# plt.ylabel('$\delta$, º')
# plt.xlabel('ratio')
# plt.show()

plt.scatter(lenght,delta_max,s=1)
plt.ylabel('$\delta$, º')
plt.xlabel('lenght, ns')
plt.show()


plt.scatter(max_impuse_summed,delta_max,s=1)
plt.ylabel('$\delta$, º')
plt.xlabel('max')
plt.show()

###=============================================================

## DEPENDANCES

#Dependance: delta(nbr PMT)

plt.scatter(x_delta_pmt,y_delta_pmt,color='grey')
plt.errorbar(x_delta_pmt, y_delta_pmt, yerr = y_delta_pmt_std,fmt ='o',markersize=0.1, capsize=3,color='grey')
plt.ylabel('<$\delta$>, º',fontsize=14)
plt.xlabel('Число ФЭУ',fontsize=14)
plt.ylim(0,5)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

# Dependance: delta(lenght of the signal)

plt.scatter(x_delta_lenght,y_delta_lenght,color='dodgerblue')
plt.errorbar(x_delta_lenght, y_delta_lenght, yerr = y_delta_lenght_std,fmt ='o',markersize=0.1, capsize=3,color='dodgerblue')
plt.ylabel('<$\delta$>, º',fontsize=14)
plt.xlabel('Длина импульса, нс ',fontsize=14)
plt.ylim(0,5)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

if analyse!='only_sig' and analyse!='sph3':

    # Dependance: delta(ratio)
    
    plt.scatter(x_delta_ratio,y_delta_ratio,color='green')
    plt.errorbar(x_delta_ratio, y_delta_ratio, yerr = y_delta_ratio_std,fmt ='o',markersize=0.1, capsize=3,color='green')
    plt.ylabel('<$\delta$>, º',fontsize=14)
    plt.xlabel('Отношение макс. сигнал / фон ',fontsize=14)
    plt.ylim(0,5)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()

#Dependance: delta(d0)

plt.scatter(x_delta_d0,y_delta_d0,color='orange')
plt.errorbar(x_delta_d0, y_delta_d0, yerr = y_delta_d0_std,fmt ='o',markersize=0.1, capsize=3,color='orange')
plt.ylabel('<$\delta$>, º',fontsize=14)
plt.xlabel('Ось ливня, м',fontsize=14)
plt.ylim(0,5)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

#Dependance: delta(max sig)

plt.scatter(x_delta_max_sig,y_delta_max_sig,color='green')
plt.errorbar(x_delta_max_sig, y_delta_max_sig, yerr = y_delta_max_sig_std,fmt ='o',markersize=0.1, capsize=3,color='green')
plt.ylabel('<$\delta$>, º',fontsize=14)
plt.xlabel('Максимум сигнала, фот',fontsize=14)
plt.ylim(0,5)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.show()

#Dependance: delta(theta)
if analyse!='sph3':
    plt.scatter(x_delta_theta,y_delta_theta,color='purple')
    plt.errorbar(x_delta_theta, y_delta_theta, yerr = y_delta_theta_std,fmt ='o',markersize=0.1, capsize=3,color='purple')
    plt.ylabel('<$\delta$>, º',fontsize=14)
    plt.xlabel('$\Theta$r, º',fontsize=14)
    plt.ylim(0,5)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()


