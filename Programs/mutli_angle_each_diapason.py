%%time
##=== IMPORT MODULES

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
import glob
import os


###=========================================================================

##=== !!! CHANGE ONLY HERE !!!

#Type of telescope
name_telesc='SPHERE-2' #SPHERE-2, SPHERE-3

#Altitude of the telescope
H=900 #500, 900, 1000 m

#Type of analyse
# type_analyse='only_sig'
# type_analyse='sig_&_bg'
type_analyse='electronic'

#Energy
En=30

part_of_int=0.5


###=========================================================================
###=========================================================================

##=== OPEN FILES

#===============
#== For all types of analyse

#PMT on mosaic: x (mm),y (mm), num
if name_telesc=='SPHERE-2':
    x_y_mos=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_background/Data/Data_init/mosaic_pmt_coords_df.csv')

if name_telesc=='SPHERE-3':
    pixel_data = pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/SPHERE3_pixel_data_A.dat',header=None, names=['x','y','z','nx','ny','nz','theta','phi'], sep='\s+') #number of PMT
    pixel_data['segment'] = pixel_data.index // 7 #== add num segment
    pixel_data['pixel'] =  pixel_data.index % 7 #== add num pixel

#===============

#Open the number of the PMT around a PMT

circle_of_pmt_around_pmt=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/circle_of_pmt_around_pmt.csv',)#===============

#===============

#Path of the directory of the events

#Only signal

if type_analyse!='electronic':
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/2_data_signal/only_sig_{}PeV_P_{}m'.format(En,H) #Signal only
    # path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/bad_res_500m'.format(En,H) #Signal only

#Signal with background  
 
    if type_analyse=='sig_&_bg':
        data_photons_bg = pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/for_BG-distribution_photons/{}m'.format(H))

#Electronics

if type_analyse=='electronic':
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/1_data_electronics/electronics_{}PeV_P_{}m'.format(En,H,En) #Electronic signal #044 #072
    # path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/1_data_electronics/electronics_30PeV_P_500m/mosaic_hits_m01_pro_30PeV_10-20_001_c031' #Electronic signal d0=147 m
    # path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/bad_res/mosaic_hits_m01_Fe_10PeV_10-20_038_c005' #Electronic signal #044 #072

name_files = glob.glob(os.path.join(path,'*'))

###=========================================================================
###=========================================================================

##=== FUNCTION

# __________________

# Return event with the type of analyse

def read_mosaic_hits_file(path, step = 12.5):

    a = pd.read_csv(path, header=None, skiprows=[0],sep='\s+')     
    a[4]=a[4]-min(a[4]) # номер столбца со временем - 4
       
    pd_event=pd.DataFrame({'pmt_ev':list(a[0]),'t_ev':list(a[4])}) #creation of DataFrame to facilite calculations
    q = np.zeros((109, 1020))  # empty array109 PMT * 1020 cells
    cells_here=np.arange(0,1021*step,step)
       
    shift = 500 #first bin of the signal
    
    for i in num_pmt:
        impulse = list((pd_event.query(f'pmt_ev=={i}'))['t_ev']) #impulse in the i pmt
        n_photons_ev,t_photon_evv=np.histogram(impulse,bins=cells_here) #number and time of photons in each cell of the i pmt
        n_phot_shift=[0]*shift+n_photons_ev[:-shift].tolist() #same (n) with shift
        q[i]=n_phot_shift

    return pd.DataFrame(q.T)

# ===

def read_mosaic_hits_file_background(path, step = 12.5):

    a = pd.read_csv(path, header=None, skiprows=[0],sep='\s+')     
    a[4]=a[4]-min(a[4]) # номер столбца со временем - 4
       
    pd_event=pd.DataFrame({'pmt_ev':list(a[0]),'t_ev':list(a[4])}) #creation of DataFrame to facilite calculations
    q = np.zeros((109, 1020))  # empty array109 PMT * 1020 cells
    cells_here=np.arange(0,1021*step,step)
       
    shift = 500 #first bin of the signal
    
    for i in num_pmt:
        impulse = list((pd_event.query(f'pmt_ev=={i}'))['t_ev']) #impulse in the i pmt
        n_photons_ev,t_photon_evv=np.histogram(impulse,bins=cells_here) #number and time of photons in each cell of the i pmt
        n_phot_shift=[0]*shift+n_photons_ev[:-shift].tolist() #same (n) with shift
        q[i]=n_phot_shift    
        q[i]=n_phot_shift+np.random.choice(data_photons_bg['{}'.format(i)].dropna(),size=1020)
            
    return pd.DataFrame(q.T)

# ===

def read_electonic_event_file(path):
    event = pd.read_csv(path, header=None, sep='\s+', skiprows=[0],on_bad_lines='skip')
    for pmt in range(event.shape[1]):
        event.iloc[0::2, pmt] -= event.iloc[0:400:2, pmt].mean()
        event.iloc[1::2, pmt] -= event.iloc[1:400:2, pmt].mean()
    return event

# __________________

# Return start values of multistart

def theta_initt(theta_init):
    return theta_init*len(theta_init)

##===================

#=== DIAPASON OF THE IMPULSE

def diapason_impulse(event):
    pulse = event.sum(axis=1) # total impulse in 1020 cells
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
    
    return event_diapason,cells_diapason

##===================

#=== TIME OF ARRIVAL

#Return PMT than can be analysed

def valuable_PMT(df_front_N_sup_10,df_front_N_inf_10,df_front_N):
    for i in range(0,20):
        isolated_num=[]
        new_PMT_ok=[]
        total_PMT_ok=list(df_front_N_sup_10['PMT'])
        total_t_ok=list(df_front_N_sup_10['t_max'])
    
        for num in range(0,len(df_front_N_inf_10)):

            if max(df_front_N_inf_10['N_max']) <= max_noise[num]:
                continue
            
            #PMT around the analysed PMT
            PMT_around=circle_of_pmt_around_pmt.loc[df_front_N_inf_10.iloc[num].name].dropna()[1:]
            
            #No PMT > thrs around the analysed PMT
            if df_front_N_sup_10['PMT'].isin(list(PMT_around)).any() == False and df_front_N_inf_10['PMT'].isin(list(PMT_around)).any() == False:
                continue
            if df_front_N_sup_10['PMT'].isin(list(PMT_around)).any() == False and df_front_N_inf_10['PMT'].isin(list(PMT_around)).any() == True:
                isolated_num.append(df_front_N_inf_10.iloc[num].name) #Save the PMT number without PMT around
          
                continue
            else:
                mean_t=df_front_N_sup_10['t_max'][df_front_N_sup_10['PMT'].isin(list(PMT_around))].mean()
          
            if df_front_N_inf_10.iloc[num]['t_max'] < mean_t - 100 or df_front_N_inf_10.iloc[num]['t_max'] > mean_t + 100:
                continue
            else:
                total_PMT_ok.append(df_front_N_inf_10.iloc[num]['PMT']) #new total PMT
                total_t_ok.append(df_front_N_inf_10.iloc[num]['t_max']) #new total time
                new_PMT_ok.append(df_front_N_inf_10.iloc[num]['PMT']) # new saved PMT
       
        df_front_N_sup_10=df_front_N[df_front_N['PMT'].isin(total_PMT_ok)]
        df_front_N_inf_10=df_front_N_inf_10.drop(new_PMT_ok)
        
        if len(isolated_num)==0:
            break
    
    x=np.array(df_front_N_sup_10['x'])
    y=np.array(df_front_N_sup_10['y'])
    t=np.array(df_front_N_sup_10['t_max'])
    N=np.array(df_front_N_sup_10['N_max'])
    PMT=np.array(df_front_N_sup_10['PMT'])
    index=np.array(df_front_N_sup_10['index'])
    
    return x, y, t, N, PMT, index


#Return time of arrival and N photons by part of integral

def t_part_of_int(PMT_total,index_total):
    t_front_reconstr=[]
    N_front_reconstr=[]
    x_front_reconstr=[]
    y_front_reconstr=[]
    PMT=[]
    for i in range(0,len(x_front)):
        # print(i)
        window=all_window[PMT_total[i]]
        event_pmt=window
    
        imax_pmt = int(index_total[i]) #index of the max impulse in the diapason
        start_pmt, stop_pmt = imax_pmt, imax_pmt
         
        
        threshold=mean_noise[PMT_total[i]]

        
        while event_pmt[start_pmt] > threshold and start_pmt>event_pmt.index[0]:
            start_pmt -= 1  
             
        while event_pmt[stop_pmt] > threshold and stop_pmt<event_pmt.index[-1]:
            stop_pmt += 1
             
        # #Impulse in the finded diapason
        impulse_filtrate=event_pmt[start_pmt+1-min(event_pmt.index):stop_pmt-min(event_pmt.index)]
        cells_impulse_filtrate=cells_diapason[start_pmt+1-min(event_pmt.index):stop_pmt-min(event_pmt.index)]
    
        if len(impulse_filtrate)<1:
            continue
    #___________________
    
    #Aoment of the arrival of shower in the PMT
        t_arrival=cells_impulse_filtrate[len(impulse_filtrate.cumsum()[impulse_filtrate.cumsum() < impulse_filtrate.sum()*part_of_int])-1]
    
        t_front_reconstr.append(t_arrival-t_path[int(PMT_total[i])])
        N_front_reconstr.append(impulse_filtrate.sum())
        x_front_reconstr.append(x_snow[PMT_total[i]])
        y_front_reconstr.append(y_snow[PMT_total[i]])
        PMT.append(PMT_total[i])

    return t_front_reconstr, N_front_reconstr, x_front_reconstr, y_front_reconstr, PMT

##===================

#=== ANGLES

def angles(x_front,y_front,t_fnc,N_photons):
    
    def fnc(theta,phi,a0,a1,a2):#определение функции с не заданными параметрами
            x_casc=((np.cos(theta)*np.cos(phi))*(x_front)+(np.cos(theta)*np.sin(phi))*(y_front)) #координаты x фотонов в системе ливня
            y_casc=(-np.sin(phi)*(x_front)+np.cos(phi)*(y_front)) #координаты y фотонов в системе ливня
            z_casc=((np.sin(theta)*np.cos(phi))*(x_front)+(np.sin(theta)*np.sin(phi))*(y_front)) #координаты z фотонов в системе ливня
            R=np.sqrt(x_casc**2+y_casc**2) # расстояние фотонов от центра оси в системе ливня
            tau=a0+a1*R+a2*R**2+z_casc/c_ns #аппроксимированое время
            s=(((t_fnc-tau)**2))*N_photons #функцию которую надо аппроксимировать
            Smin=np.sum(s)
            return Smin
        
    theta_multi_start=[]
    phi_multi_start=[]
    min_s=[]
    for th in range(0,len(theta_initt)):
            param_fix=Minuit(fnc, theta=theta_initt[th],phi=phi_initt[th],a0=0,a1=0,a2=0)
            param_fix.limits['theta']=(0,30/180*np.pi)
            param_fix.limits['phi']=(0,2*np.pi)
    
            param_fix.migrad()
            theta_multi_start.append(param_fix.values[0]) 
            phi_multi_start.append(param_fix.values[1])
            min_s.append(param_fix.fval)
        
    theta=np.array(theta_multi_start)[np.where(np.array(min_s)==min(min_s))][0]
    phi=np.array(phi_multi_start)[np.where(np.array(min_s)==min(min_s))][0]
    
    return theta, phi

###=========================================================================
###=========================================================================

##=== CONSTANT VALUES

cells=np.arange(0,1020*12.5,12.5)
num_pmt=np.arange(0,109)

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


# === VALUES FOR d0

if H==500:
    d0_center=180
    d0_before_c=140
    d_six=50
    d_double=100
if H==900:
    d0_center=330
    d0_before_c=260
    d_six=100
    d_double=150

#=====FNC

#a0 fnc
# a_a0=-0.630484 #more in the article "Алгоритм восстановления направления прихода ШАЛ для телескопа СФЕРА"
# b_a0=53.541113661

# a2=0.0005

#Paramatrization
theta_init=np.array([0.05,0,1,0.15, 0.2,0.25, 0.3,0.35]).tolist() #in rad
phi_init=np.arange(0.5,6.5,0.5) #in rad

theta_initt=theta_initt(theta_init)
phi_initt=phi_init.tolist()*len(theta_init)


#=====DELTA

real_theta_deg=0.1745*180/np.pi
real_phi_deg=3.0564*180/np.pi

theta_real=0.3277 #[0]
phi_real=3.0564 #[0]

# theta_real= 0.2718 #[files 008]
# phi_real=2.4855 #[files 008]

# theta_real=0.3332 #[14]
# phi_real=5.7728 #[14]

# theta_real=0.2175
# phi_real=2.982

# theta_real=0.1745
# phi_real=3.0564

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

###========================================================================
###========================================================================

delta_by_max=[]
delta_by_part_of_int=[]
files_saved=[]
total_nbr_pmt=[]
for files in range(0,len(name_files)):
    # print(name_files[files])

    #Open file
    # files=5
    if type_analyse=='only_sig':
        event = read_mosaic_hits_file(name_files[files])
        
    # Signal with the background
    if type_analyse=='sig_&_bg':
        event = read_mosaic_hits_file_background(name_files[files])

    #Signal with the electronic
    if type_analyse=='electronic':
        event = read_electonic_event_file(name_files[files])
    
    if len(event)<1020:
        continue

    ## START OF THE ALGORITHM
 
    #Return the DIAPASON of the impulse 
 
    found_diapason=diapason_impulse(event)
    event_diapason,cells_diapason = found_diapason[0],found_diapason[1]
    event_diapason.index=np.arange(len(event_diapason))
    
    noise_window=event[0:400].rolling(window=4).sum().dropna()
    mean_noise=noise_window[noise_window>0].mean()
    max_noise=noise_window.max()

    if type_analyse=='only_sig':
        big_trsh=20
    if type_analyse=='electronic':
        big_trsh=max_noise.max()+10

 
    #=========================
 
    #SLIDING WINDOW
 
    all_window=event_diapason.rolling(window=4).sum().dropna()
    all_window.index=np.arange(0,len(all_window))
 
    #Determine the position of the max
    x0_center=x_snow[all_window.max().idxmax()]
    y0_center=y_snow[all_window.max().idxmax()]
    d0_from_center=np.sqrt(x0_center**2+y0_center**2)
 
    if d0_from_center>d0_center:
        # print('not this')
        continue
 
    cells_window=cells_diapason[:-3]
 
 
    df_front_N=pd.DataFrame({'PMT':num_pmt,'x':x_snow,'y':y_snow,'N_max':list(all_window.max()),'t_max':cells_window[all_window.idxmax()],'index':all_window.idxmax()})
    df_front_N['N_max']=df_front_N['N_max'][df_front_N['N_max']>=big_trsh]
    df_front_N_sup_10=df_front_N.dropna()
 
 
    df_front_N=pd.DataFrame({'PMT':num_pmt,'x':x_snow,'y':y_snow,'N_max':list(all_window.max()),'t_max':cells_window[all_window.idxmax()],'index':all_window.idxmax()})
    df_front_N['N_max']=df_front_N['N_max'][df_front_N['N_max']<big_trsh]
    df_front_N_inf_10=df_front_N.dropna()
 
    df_front_N=pd.DataFrame({'PMT':num_pmt,'x':x_snow,'y':y_snow,'N_max':list(all_window.max()),'t_max':cells_window[all_window.idxmax()],'index':all_window.idxmax()})
 
    #Return PMT that can be analysed
 
    PMT_ok= valuable_PMT(df_front_N_sup_10,df_front_N_inf_10,df_front_N)
    x_front,y_front,t_max,N_max,PMT_max,index_max=PMT_ok[0],PMT_ok[1],PMT_ok[2],PMT_ok[3],PMT_ok[4],PMT_ok[5]
 
    if len(t_max)<=50:
        continue
    
    #Return angles by max of PMT
    new_pmt=t_max-t_path[PMT_max]

    both_angles=angles(x_front,y_front,t_max,N_max)
    theta=both_angles[0]
    phi=both_angles[1]
 
    theta_fit=theta
    phi_fit=phi
    cx = np.sin(theta_real) * np.cos(phi_real) #theta - истинные углы
    cy = np.sin(theta_real) * np.sin(phi_real)
    cz = np.cos(theta_real)
    cx_fit = np.sin(theta_fit) * np.cos(phi_fit) #theta' -полученные углы моделированием
    cy_fit = np.sin(theta_fit) * np.sin(phi_fit)
    cz_fit = np.cos(theta_fit)
    delta=np.arccos(cx *cx_fit + cy * cy_fit + cz * cz_fit)*180/np.pi
 
    delta_by_max.append(delta)

    #____________________________________
 
    #Return time of arrival by PART OF INTEGRAL

    time_part_of_int=t_part_of_int(PMT_max,index_max)  
    t_snow_int,N_front_reconstr,x_front_int,y_front_int, PMT_int=time_part_of_int[0],time_part_of_int[1] ,time_part_of_int[2] ,time_part_of_int[3], time_part_of_int[4]
        
    #Return ANGLES by part of integral

    both_angles=angles(np.array(x_front_int),np.array(y_front_int),np.array(t_snow_int)-min(t_snow_int),np.array(N_front_reconstr))
    
    theta=both_angles[0]
    phi=both_angles[1]
        
    # ======================================================================================
    # ======================================================================================
 
 
    # # #========================================== Delta
 
    # print('theta_fit=',theta*180/np.pi)
    # print('phi_fit=',phi*180/np.pi)
 
    # print('theta_real=',theta_real*180/np.pi)
    # print('phi_real=',phi_real*180/np.pi)
 
    theta_fit=theta
    phi_fit=phi
    cx = np.sin(theta_real) * np.cos(phi_real) #theta - истинные углы
    cy = np.sin(theta_real) * np.sin(phi_real)
    cz = np.cos(theta_real)
    cx_fit = np.sin(theta_fit) * np.cos(phi_fit) #theta' -полученные углы моделированием
    cy_fit = np.sin(theta_fit) * np.sin(phi_fit)
    cz_fit = np.cos(theta_fit)
    delta=np.arccos(cx *cx_fit + cy * cy_fit + cz * cz_fit)*180/np.pi
 


    files_saved.append(name_files[files])
    
    delta_by_part_of_int.append(delta)
    
    total_nbr_pmt.append(len(t_max))
 

    
###==================================================================
###==================================================================

print(np.mean(delta_by_max))
print(np.mean(delta_by_part_of_int))

# #========================================== Delta

# print('theta_fit=',theta*180/np.pi)
# print('phi_fit=',phi*180/np.pi)

# print('theta_real=',theta_real*180/np.pi)
# print('phi_real=',phi_real*180/np.pi)

# theta_fit=theta
# phi_fit=phi
# cx = np.sin(theta_real) * np.cos(phi_real) #theta - истинные углы
# cy = np.sin(theta_real) * np.sin(phi_real)
# cz = np.cos(theta_real)
# cx_fit = np.sin(theta_fit) * np.cos(phi_fit) #theta' -полученные углы моделированием
# cy_fit = np.sin(theta_fit) * np.sin(phi_fit)
# cz_fit = np.cos(theta_fit)
# delta=np.arccos(cx *cx_fit + cy * cy_fit + cz * cz_fit)*180/np.pi

   
# print(delta)


# #========================== Tau

# x_casc=((np.cos(theta)*np.cos(phi))*(x_front-x0_front)+(np.cos(theta)*np.sin(phi))*(y_front-y0_front)) #координаты x фотонов в системе ливня
# y_casc=(-np.sin(phi)*(x_front-x0_front)+np.cos(phi)*(y_front-y0_front)) #координаты y фотонов в системе ливня
# z_casc=((np.sin(theta)*np.cos(phi))*(x_front-x0_front)+(np.sin(theta)*np.sin(phi))*(y_front-y0_front)) #координаты z фотонов в системе ливня
# R=np.sqrt(x_casc**2+y_casc**2) # расстояние фотонов от центра оси в системе ливня
# tau=a0+a1_front*R+a2_front*R**2 +z_casc/c_ns


# ###==================================================================
# ###==================================================================

# #========= GRAPHICS

plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in"})

# ==================================

#========= APPROXIMATION OF T INCREASING

# All

# # plt.scatter(t_try_x,t_try,s=1)
# # plt.legend(title='Approximation \ntheta={:.2f}º \nE={} PeV \nH={} m'.format(18.78,En,H)) #18.78, or 10.4, 12,47
# # plt.ylim(-100,600)
# # plt.ylabel('t, ns')
# # plt.xlabel('d, m ; norm. on d0')
# # plt.show()

# plt.scatter(d_norm_same_all,t_fnc_reconstr_approx_same_all,s=1)
# plt.plot(x_approx_of_all_approx,approx_of_all_approx,color='red')
# plt.legend(title='Approximation \ntheta={:.2f}º \nE={} PeV \nH={} m \na={:.4f} ±0.002\nb={:.2f}±0.15 \nc=free'.format(18.78,En,H,a_mean,b_mean)) #18.78, or 10.4
# plt.ylim(-100,600)
# plt.ylabel('t, ns')
# plt.xlabel('d, m ; norm. on d0')
# plt.show()


# ##========= Parameters a0, a1, a2 when they aren't fixed

# plt.scatter(b_front,delta_total)
# # # plt.scatter(b_front,delta_total)
# # # plt.scatter(a_front,delta_total)
# plt.ylabel('delta, º')
# plt.xlabel('b')
# plt.show()

#DELTA

plt.hist(delta_by_max)
plt.show()

plt.hist(delta_by_part_of_int)
plt.show()

plt.scatter(total_nbr_pmt,delta_by_max)
plt.show()




