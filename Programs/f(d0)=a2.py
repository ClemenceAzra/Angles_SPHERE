import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

results=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/choose_a2_coulissante/choose_a0_500m_P_30PeV_res_10-20deg.csv',header=None,sep='\s+').sort_values(by=[2])    

#==== Angles ~ 10º, ~ 20º to analyze odd thing with a2
files=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/choose_a2_coulissante/choose_a0_500m_P_30PeV_files_10-20deg.csv',header=None,sep='\s+')[0]    


real_theta=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_theta',header=None,sep='\s+')[0]) #number of PMT


# num_ev= [ real_theta[int(files[i][-8:-5])] for i in range(0,len(files)) ]
theta_real=[real_theta[int(files[i][-8:-5])-1]*180/np.pi for i in range(0,len(files))]#10-20º P, N

new_pd=pd.DataFrame({'a2':list(results[1]),'d0':list(results[2]),'theta_real':theta_real})

# a2=list(new_pd['a2'][ (new_pd['theta_real']<20) & (new_pd['theta_real']>19) ])
# d0=list(new_pd['d0'][ (new_pd['theta_real']<20) & (new_pd['theta_real']>19) ])

# delta, a2, d0, mean monter, sum

a2=results[1]
d0=results[2]

a2_d0=pd.DataFrame({'a2':a2,'d0':d0}).sort_values(by=['d0'])
a2_d0['d0']=a2_d0['d0']. round(0)

moy_a2=a2_d0.groupby(a2_d0['d0']).mean()

#Approximation

a,b,c=np.polyfit(list(moy_a2.index)[1:],list(moy_a2['a2'])[1:],2)
# a,b,c=np.polyfit(d0,a2,2)

xx=np.arange(0,max(list(moy_a2.index)),1)
yy=a*xx**2+b*xx+c


# d0_axis=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/500m_P_30PeV_sig_total_d0.csv',header=None,sep='\s+').sort_values(by=[2])    
#d0 real, optimal number of pmt, min err d0, ∑ phot event
print(results[0].mean())

###==================================================================

#========= GRAPHICS

plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in"})

# results.plot.scatter(x=d0,y=a2,s=0.7)
# plt.plot(results[2][results[2]>40],yy,color='red', label='a=-6 *10^-6 \nb=1.6*10^-3')
plt.scatter(d0,a2,s=0.4)
# plt.plot(xx,yy,color='red')
plt.xlabel('d0, m')
plt.ylabel('a0')
# plt.yscale('log')
# plt.legend(title='19-20 deg')
plt.show()

#========= 
plt.scatter(list(moy_a2.index),list(moy_a2['a2']),s=1)
plt.plot(xx,yy,color='red',label='a={:.9f} \nb={:.7f} \nc={:.4f}'.format(a,b,c))
plt.ylabel('mean a0')
plt.grid(True)
# plt.xlim(0,None)
# plt.ylim(0,None)
plt.legend()
plt.xlabel('d0, m')
plt.show()

# plt.scatter(d0_axis[3],d0_axis[1],s=0.6)
# plt.xlabel('∑ photons')
# plt.ylabel('N PMT')
# plt.show()

# plt.scatter(d0_axis[0],d0_axis[1],s=0.6)
# plt.xlabel('d0 real')
# plt.ylabel('N PMT')
# plt.show()