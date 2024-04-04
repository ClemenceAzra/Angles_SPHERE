%%time
##=== IMPORT MODULES

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
import glob
import os
import random

###=========================================================================

##=== !!! CHANGE ONLY HERE !!!

#Type of telescope
name_telesc='SPHERE-2' #SPHERE-2, SPHERE-3

#Altitude of the telescope
H=500 #500, 900, 1000 m

#Type of analyse
# type_analyse='only_sig'
# type_analyse='sig_&_bg'
type_analyse='electronic'

#Energy
En=10

###=========================================================================

##=== FONCTIONS

#===============

#Path of the directory of the events

#Only signal

if type_analyse!='electronic':
    
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/2_data_signal/only_sig_{}PeV_P_{}m'.format(En,H) #Signal only

    def read_mosaic_hits_file(path, step = 12.5):
        
        a = pd.read_csv(path, header=None, skiprows=[0],sep='\s+')     
        a[4]=a[4]-min(a[4]) # номер столбца со временем - 4
       
        pd_event=pd.DataFrame({'pmt_ev':list(a[0]),'t_ev':list(a[4])}) #creation of DataFrame to facilite calculations
        q = np.zeros((109, 1020))  # empty array109 PMT * 1020 cells
        cells_here=np.arange(0,1021*step,step)
     
        shift = random.randint(400, 600) #first bin of the signal
        
        for i in num_pmt:
            impulse = list((pd_event.query(f'pmt_ev=={i}'))['t_ev']) #impulse in the i pmt
            n_photons_ev,t_photon_evv=np.histogram(impulse,bins=cells_here) #number and time of photons in each cell of the i pmt
            n_phot_shift=[0]*shift+n_photons_ev[:-shift].tolist() #same (n) with shift
            q[i]=n_phot_shift

#Signal with background

            if type_analyse=='sig_&_bg':
                data_photons_bg = pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/for_BG-distribution_photons/{}m'.format(H))

                q[i]=n_phot_shift+np.random.choice(data_photons_bg['{}'.format(i)].dropna(),size=1020)
                
        return pd.DataFrame(q.T)
#Electronics

if type_analyse=='electronic':
    
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/1_data_electronics/electronics_{}PeV_P_{}m'.format(En,H) #Electronic signal
    
    def read_electonic_event_file(path):
        event = pd.read_csv(path, header=None, sep='\s+', skiprows=[0],on_bad_lines='skip')
        for pmt in range(event.shape[1]):
            event.iloc[0::2, pmt] -= event.iloc[0:400:2, pmt].mean()
            event.iloc[1::2, pmt] -= event.iloc[1:400:2, pmt].mean()
        return event


name_files = glob.glob(os.path.join(path,'*'))

#== Start values for multistart

def theta_initt(theta_init):
    return theta_init*len(theta_init)

#== Fonction of the approximation of the front

def tau_fnc(theta,phi,a0,a1,a2):
    x_casc=((np.cos(theta)*np.cos(phi))*(x_front-x0_front)+(np.cos(theta)*np.sin(phi))*(y_front-y0_front)) #координаты x фотонов в системе ливня
    y_casc=(-np.sin(phi)*(x_front-x0_front)+np.cos(phi)*(y_front-y0_front)) #координаты y фотонов в системе ливня
    z_casc=((np.sin(theta)*np.cos(phi))*(x_front-x0_front)+(np.sin(theta)*np.sin(phi))*(y_front-y0_front)) #координаты z фотонов в системе ливня
    R=np.sqrt(x_casc**2+y_casc**2) # расстояние фотонов от центра оси в системе ливня
    return a0+a1*R+a2*R**2+z_casc/c_ns

###=========================================================================

##Coordinates of PMT in mosaic

#PMT on mosaic: x (mm),y (mm), num
if name_telesc=='SPHERE-2':
    x_y_mos=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_background/Data/Data_init/mosaic_pmt_coords_df.csv')

if name_telesc=='SPHERE-3':
    pixel_data = pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/SPHERE3_pixel_data_A.dat',header=None, names=['x','y','z','nx','ny','nz','theta','phi'], sep='\s+') #number of PMT
    pixel_data['segment'] = pixel_data.index // 7 #== add num segment
    pixel_data['pixel'] =  pixel_data.index % 7 #== add num pixel

###=========================================================================

##=== CONSTANT VALUES

cells=np.arange(0,1020*12.5,12.5) #1020 cells with 12.5 ns each
num_pmt=np.arange(0,109) #Num of PMT from 1 to 109 

if name_telesc=='SPHERE-2':
    if H==500:
        a=0.0005444285453203233 #a, b for different H --> find in another work
        b=0.9051146902370498
    elif H==900:
        a=0.0008982897531551649
        b=1.6455812746636922
       
if name_telesc=='SPHERE-3':
    b=(1.1314175328971585)/(1000/H) #If H=500 - a --> a/2


c_ns=0.3 #speed of light m/ns

part_of_int=0.5 #part of integral 

#=====FNC

#a0 fnc
# a_a0=-0.630484 #more in the article "Алгоритм восстановления направления прихода ШАЛ для телескопа СФЕРА"
# b_a0=53.541113661

#Multistart
theta_init=np.array([0.05,0,1,0.15, 0.2,0.25, 0.3,0.35]).tolist() #in rad
phi_init=np.arange(0.5,6.5,0.5) #in rad

theta_initt=theta_initt(theta_init)
phi_initt=phi_init.tolist()*len(theta_init)


#=====DELTA

# real_theta_deg=0.1745*180/np.pi
# real_phi_deg=3.0564*180/np.pi

theta_real=0.3277
phi_real=3.0564

# # theta_real=0.2175
# # phi_real=2.982

# # theta_real=0.1745
# # phi_real=3.0564

##=== PRELIMINARY CALCULATIONS

#==Translation x,y from mos to snow

if name_telesc=='SPHERE-2':
    x_snow=-b*x_y_mos['x']
    y_snow=-b*x_y_mos['y']

if name_telesc=='SPHERE-3':
    x_snow=-b*np.array(pixel_data.groupby('segment').mean()['x'])
    y_snow=-b*np.array(pixel_data.groupby('segment').mean()['y'])

#==Translation t from mos to snow (path from mos to snow)
  
t_path=((np.sqrt(H**2+(np.sqrt(np.array(x_snow)**2+np.array(y_snow)**2))**2))/c_ns) #path in ns

###=========================================================================

theta_res=[]
phi_res=[]
d0_res=[]
sum_phot_res=[]
files_res=[]

deltaa=[]
sum_pmt=[]

#===
for files in range(0,len(name_files)):

    #Open file
    if type_analyse!='electronic':
        event = read_mosaic_hits_file(name_files[files])
    if type_analyse=='electronic':
        event = read_electonic_event_file(name_files[files])
    
    if len(event)<1020:
        continue
    
    ###========================================================================
    
    
    #Search the impulse in the отклик
    
    pulse = event.sum(axis=1) # total impulse in 1020 cells
    
    if max(pulse[:-1])<max(pulse[0:400])*2:
        continue
    
    imax = pulse[:-1].idxmax() #index of the max impulse --> delete the last cell -> > max of impulse
    start, stop = imax, imax
    
    maxpart = 1
    bgd = maxpart * max(pulse[0:400])
    while pulse[start] > bgd:
        start -= 1
    while pulse[stop] > bgd:
        stop += 1
    
    margin_left,margin_right= 30, 30 
    
    stop  += margin_right
    start -= margin_left
    
    event_diapason=event[start:stop] #diapason of the impulse - n photons
    cells_diapason=cells[start:stop] #diapason of the impulse - time
    
    ###============================================
    
    #Save the PMT and bins in with no BG 
    
    saved_pmt=[]
    t_front_reconstr=[]
    x_front_reconstr=[]
    y_front_reconstr=[]
    number_photons_i=[]
    
    impulse_each_pmt=[]
    cells_each_pmt=[]
    
    for i in range(0,event_diapason.shape[1]):
       #  event_pmt=event_diapason[i]
        
       #  #Condition 1 : delete the PMT with max signal < max BG
       #  if max(event_pmt)<max(event[i][0:400]):
       #      continue
        
       #  imax_pmt = event_pmt.idxmax() #index of the max impulse
       #  start_pmt, stop_pmt = imax_pmt, imax_pmt

       #  bgd_pmt = max(event[i][0:400])
       #  while event_pmt[start_pmt] > bgd_pmt and start_pmt>event_pmt.index[0]:
       #      start_pmt -= 1  
            
       #  while event_pmt[stop_pmt] > bgd_pmt and stop_pmt<event_pmt.index[-1]:
       #      stop_pmt += 1
            
       #  # #Impulse in the finded diapason
       #  impulse=event_pmt[start_pmt-min(event_pmt.index):stop_pmt-min(event_pmt.index)]
       #  cells_impulse=cells_diapason[start_pmt-min(event_pmt.index):stop_pmt-min(event_pmt.index)]

       #  #Condition 2: delete bins < max BG
       #  impulse_filtrate=impulse[impulse>max(event[i][0:400])] 
       #  cells_impulse_filtrate=cells_impulse[impulse>max(event[i][0:400])] 

            
       # #Condition 3 : delete the PMT with 0 bins
       #  if len(cells_impulse_filtrate)<1: #if number of bins in PMT less than 4 - low signal
       #      continue
        
       #  #Aoment of the arrival of shower in the PMT
       #  t_arrival=cells_impulse_filtrate[len(impulse_filtrate.cumsum()[impulse_filtrate.cumsum() < impulse_filtrate.sum()*part_of_int])-1]

       #  #SAVE DATA
        
       #  t_front_reconstr.append(t_arrival-t_path[i]) #t - obtain front
       #  x_front_reconstr.append(x_snow[i]) #x - obtain front
       #  y_front_reconstr.append(y_snow[i]) #y - obtain front
        
       #  saved_pmt.append(num_pmt[i]) # Nº of the saved PMT with signal
       #  number_photons_i.append(impulse_filtrate.sum()) #number of signal photons in each PMT
       event_pmt=event_diapason[i]
       
       event_window=event_pmt.rolling(window=4).sum().dropna()
       noise_widnow=event[i][0:400].rolling(window=4).sum().dropna()
       # event_pmt.append(num_pmt[i])
       #___________________

       #Condition 1 for PMT: delete the PMT with max signal < max BG

       event_pmt=event_window
       if max(event_pmt)<max(noise_widnow):
           continue
       
       #___________________
       
       imax_pmt = event_pmt.idxmax() #index of the max impulse
       start_pmt, stop_pmt = imax_pmt, imax_pmt

       bgd_pmt = max(noise_widnow)
       while event_pmt[start_pmt] > bgd_pmt and start_pmt>event_pmt.index[0]:
           start_pmt -= 1  
           
       while event_pmt[stop_pmt] > bgd_pmt and stop_pmt<event_pmt.index[-1]:
           stop_pmt += 1
           
       # #Impulse in the finded diapason
       impulse=event_pmt[start_pmt-min(event_pmt.index):stop_pmt-min(event_pmt.index)]
       cells_impulse=cells_diapason[start_pmt-min(event_pmt.index):stop_pmt-min(event_pmt.index)]

       #___________________

       #Condition 2 for bins: delete bins < max BG
       
       impulse_filtrate=impulse[impulse>max(noise_widnow)] 
       cells_impulse_filtrate=cells_impulse[impulse>max(noise_widnow)] 

       #___________________
           
      #Condition 3 for bins: delete the PMT with 0 bins
      
       if len(cells_impulse_filtrate)<1: 
           continue
       
       #___________________

       #Aoment of the arrival of shower in the PMT
       t_arrival=cells_impulse_filtrate[len(impulse_filtrate.cumsum()[impulse_filtrate.cumsum() < impulse_filtrate.sum()*part_of_int])-1]

       #___________________

       #SAVE DATA
       
       t_front_reconstr.append(t_arrival-t_path[i]) #t - obtain front
       x_front_reconstr.append(x_snow[i]) #x - obtain front
       y_front_reconstr.append(y_snow[i]) #y - obtain front
       
       saved_pmt.append(num_pmt[i]) # Nº of the saved PMT with signal
       number_photons_i.append(impulse_filtrate.sum()) #number of signal photons in each PMT

       cells_impulse_filtrate_each_pmt.append(cells_impulse_filtrate) #to see image
       impulse_filtrate_each_pmt.append(impulse_filtrate) #to see image

       event_window_each_pmt.append(event_window)
       noise_window_each_pmt.append(noise_widnow)
       
       event_each_pmt.append(event_diapason[i])

    ###============================================
    
    ###============================================
    
    x0_front=x_front_reconstr[number_photons_i.index(max(number_photons_i))] 
    y0_front=y_front_reconstr[number_photons_i.index(max(number_photons_i))]
    t0_front=t_front_reconstr[number_photons_i.index(max(number_photons_i))]

    d0_from_center=np.sqrt(x0_front**2+y0_front**2)

    d_x_to_x_max=np.sqrt( (np.array(x_front_reconstr)-x0_front)**2  + (np.array(y_front_reconstr)-y0_front)**2  )

    if d0_from_center>140 and d0_from_center<180:
        max_d=50 

    if d0_from_center<140:
        max_d=100
    else:
        # print('no this')
        continue
        
    x_circle_around_max=np.array(x_front_reconstr)[np.where(d_x_to_x_max<max_d)[0]]
    y_circle_around_max=np.array(y_front_reconstr)[np.where(d_x_to_x_max<max_d)[0]]
    t_circle_around_max=np.array(t_front_reconstr)[np.where(d_x_to_x_max<max_d)[0]]
    number_photons_around_max=np.array(number_photons_i)[np.where(d_x_to_x_max<max_d)[0]]

    x0=np.sum(x_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max)
    y0=np.sum(y_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max)
    t0=np.sum(t_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max)

    d0=np.sqrt(x0**2+y0**2)
    
    #============== Analyse t in ascending order

    t_fnc_reconstr_d0=np.array(t_front_reconstr)-t0

    df_analyse_t=pd.DataFrame({'x':x_front_reconstr,'y':y_front_reconstr,'t':t_fnc_reconstr_d0,'N_phot':number_photons_i}).sort_values(by=['t'])
       
    a_t,b_t,c_t=np.polyfit(np.arange(0,len(t_fnc_reconstr_d0)),df_analyse_t['t'],2,w=df_analyse_t['N_phot'])
    t_fnc_reconstr_approx=a_t*np.arange(0,len(t_fnc_reconstr_d0))**2+b_t*np.arange(0,len(t_fnc_reconstr_d0))+c_t
       
       
    df_analyse_t['delta_t']=abs(t_fnc_reconstr_approx-df_analyse_t['t']) #add delta t
       
    x_front_after_t=df_analyse_t['x'][df_analyse_t['delta_t']<df_analyse_t['delta_t'].mean()*3]
    y_front_after_t=df_analyse_t['y'][df_analyse_t['delta_t']<df_analyse_t['delta_t'].mean()*3]
    t_front_after_t=df_analyse_t['t'][df_analyse_t['delta_t']<df_analyse_t['delta_t'].mean()*3]

    # #=======

    t_sinus=np.array(t_front_after_t.sort_index())
    x_sinus=np.array(x_front_after_t.sort_index())
    y_sinus=np.array(y_front_after_t.sort_index())

    diff_t_sinus=abs(t_sinus[:-1]-t_sinus[1:])

    diff_t_each=pd.DataFrame({'diff_t':diff_t_sinus}).rolling(window=2).mean().dropna()
    diff_t_each.loc[len(diff_t_each)+1], diff_t_each.loc[0] =abs(t_sinus[-1]-t_sinus[-2]), abs(t_sinus[0]-t_sinus[1])
    diff_t_each.sort_index(inplace=True)

    df_analyse_tt=pd.DataFrame({'x':x_sinus,'y':y_sinus,'t':t_sinus,'delta_t':list(diff_t_each['diff_t'])})

    x_front=df_analyse_tt['x'][ df_analyse_tt['delta_t']<df_analyse_tt['delta_t'].mean()*2]
    y_front=df_analyse_tt['y'][ df_analyse_tt['delta_t']<df_analyse_tt['delta_t'].mean()*2]
    t_fnc=df_analyse_tt['t'][ df_analyse_tt['delta_t']<df_analyse_tt['delta_t'].mean()*2]

    # ============== Minuit
    

    # a_a0,b_a0,c_a0= -0.001177952671304051,-0.3904769724822963, 44.96860327440242
    a_a0,b_a0,c_a0=  -0.0014842947511517158,-0.3275524375299413, 43.185087602409766 #500 m isol with 1 approx

    a0=a_a0*d0**2+b_a0*d0+c_a0
    # a0=c_try
    a1_front=0.0

    # a_a2, b_a2, c_a2= -7.37951417185169e-08,7.415276863337747e-06,0.0016999769569772587
    # a_a2, b_a2, c_a2= -4.8469380897279535e-08,4.735784181057871e-07,0.0020481248019126753 #500 m bis
    a_a2, b_a2, c_a2=  -6.016514714769361e-08,4.848420806082731e-06,0.001767031413096028 #500 m isol with 1 approx

    # a_a2, b_a2, c_a2= 0.0008982897531551649,1.6455812746636922, 0.0005404375020002251
    a2_front=a_a2*d0**2+b_a2*d0+c_a2


    def fnc(theta,phi):#определение функции с не заданными параметрами
            x_casc=((np.cos(theta)*np.cos(phi))*(x_front-x0_front)+(np.cos(theta)*np.sin(phi))*(y_front-y0_front)) #координаты x фотонов в системе ливня
            y_casc=(-np.sin(phi)*(x_front-x0_front)+np.cos(phi)*(y_front-y0_front)) #координаты y фотонов в системе ливня
            z_casc=((np.sin(theta)*np.cos(phi))*(x_front-x0_front)+(np.sin(theta)*np.sin(phi))*(y_front-y0_front)) #координаты z фотонов в системе ливня
            R=np.sqrt(x_casc**2+y_casc**2) # расстояние фотонов от центра оси в системе ливня
            tau=a0+a1_front*R+a2_front*R**2+z_casc/c_ns #аппроксимированое время
            s=(((t_fnc-tau)**2)) #функцию которую надо аппроксимировать
            Smin=np.sum(s)
            return Smin
        
    theta_multi_start=[]
    phi_multi_start=[]
    min_s=[]
    for th in range(0,len(theta_initt)):
            param_fix=Minuit(fnc, theta=theta_initt[th],phi=phi_initt[th])
            param_fix.limits['theta']=(0,30/180*np.pi)
            param_fix.limits['phi']=(0,2*np.pi)

            param_fix.migrad()
            theta_multi_start.append(param_fix.values[0]) 
            phi_multi_start.append(param_fix.values[1])
            min_s.append(param_fix.fval)
        
    theta=np.array(theta_multi_start)[np.where(np.array(min_s)==min(min_s))][0]
    phi=np.array(phi_multi_start)[np.where(np.array(min_s)==min(min_s))][0]

#==============

    #SAVE DATA
    theta_res.append(theta)
    phi_res.append(phi)
    d0_res.append(np.sqrt(x0_front**2+y0_front**2))
    sum_phot_res.append(sum(number_photons_i))
    files_res.append(name_files[files])
    sum_pmt.append(len(x_front))
    
    theta_fit=theta
    phi_fit=phi
    cx = np.sin(theta_real) * np.cos(phi_real) #theta - истинные углы
    cy = np.sin(theta_real) * np.sin(phi_real)
    cz = np.cos(theta_real)
    cx_fit = np.sin(theta_fit) * np.cos(phi_fit) #theta' -полученные углы моделированием
    cy_fit = np.sin(theta_fit) * np.sin(phi_fit)
    cz_fit = np.cos(theta_fit)
    delta=np.arccos(cx *cx_fit + cy * cy_fit + cz * cz_fit)*180/np.pi

    deltaa.append(delta)
    
    # monterr.append(mean_monter)
    
    
    
print(np.mean(deltaa))

plt.hist(deltaa)
plt.show()

# theta_delta=abs(np.array(theta_res)*(180/np.pi)-theta_real*180/np.pi)
plt.hist(np.array(theta_res)*(180/np.pi))
plt.axvline(theta_real*(180/np.pi),color='red')
plt.show()

plt.hist(np.array(phi_res)*(180/np.pi))
plt.axvline(phi_real*(180/np.pi),color='red')
plt.show()

plt.scatter(sum_pmt,deltaa)
plt.ylim(None,13)
plt.show()

