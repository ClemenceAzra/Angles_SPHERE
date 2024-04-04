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
H=500 #500, 900, 1000 m

#Type of analyse
# type_analyse='only_sig'
# type_analyse='sig_&_bg'
type_analyse='electronic'


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

real_theta=list(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_theta',header=None,sep='\s+')[0]) #number of PMT
real_phi=list(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_phi',header=None,sep='\s+')[0]) #number of PMT


#===============
#== For 1 type of analyse

#Path of the directory of the events

#Only signal

if type_analyse!='electronic':
    path='/Users/clemence/Documents/Scientific_research/Data/only_sig_10_PeV_P_900m/mosaic_hits_m01_pro_30PeV_10-20_001_c001' #Signal only

#Signal with background  
 
    if type_analyse=='sig_&_bg':
        data_photons_bg = pd.read_csv('/Users/clemence/Documents/Scientific_research/Data/for_BG-distribution_photons/{}m'.format(H))

#Electronics

if type_analyse=='electronic':
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/electronics/electronics_30PeV_P_{}m'.format(H) #Electronic signal

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

part_of_int=0.6 #part of integral 

#=====FNC

#a0 fnc
a_a0=-0.630484 #more in the article "Алгоритм восстановления направления прихода ШАЛ для телескопа СФЕРА"
b_a0=53.541113661

# a2=0.0006

#Multistart
theta_init=np.array([0.05,0,1,0.15, 0.2,0.25, 0.3,0.35]).tolist() #in rad
phi_init=np.arange(0.5,6.5,0.5) #in rad

def theta_initt(theta_init):
    return theta_init*len(theta_init)

theta_initt=theta_initt(theta_init)
phi_initt=phi_init.tolist()*len(theta_init)


#=====DELTA

# real_theta_deg=0.1745*180/np.pi
# real_phi_deg=3.0564*180/np.pi

# theta_real=0.3277
# phi_real=3.0564

# # theta_real=0.2175
# # phi_real=2.982

# # theta_real=0.1745
# # phi_real=3.0564

##=== PRELIMINARY CALCULATIONS

#==Translation x,y from mos to snow

if name_telesc=='SPHERE-2':
    x_snow=-b*x_y_mos['x']-a*x_y_mos['x']**2
    y_snow=-b*x_y_mos['y']-a*x_y_mos['y']**2

if name_telesc=='SPHERE-3':
    x_snow=-b*np.array(pixel_data.groupby('segment').mean()['x'])
    y_snow=-b*np.array(pixel_data.groupby('segment').mean()['y'])

#==Translation t from mos to snow (path from mos to snow)
  
t_path=((np.sqrt(H**2+(np.sqrt(np.array(x_snow)**2+np.array(y_snow)**2))**2))/c_ns) #path in ns


###=========================================================================

name_files = glob.glob(os.path.join(path,'*'))

##=== DEFINITION OF THE TYPE OF ANALYSE

#===Save to analyse final results
theta_res=[]
phi_res=[]
d0_res=[]
sum_phot_res=[]
files_res=[]

deltaa=[]

delta_choose=[]
a2_choose=[]
d0_choose=[]
mean_monter_choose=[]
sum_phot_choose=[]
 
#===
for files in range(0,len(name_files)):

    path=name_files[files]
    # print(path)

    #Only signal
    
    if type_analyse!='electronic':
    
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
                    q[i]=n_phot_shift+np.random.choice(data_photons_bg['{}'.format(i)].dropna(),size=1020)
                    
            return pd.DataFrame(q.T)
        event = read_mosaic_hits_file(path)
    
    #Electronic
    
    if type_analyse=='electronic':
    
        def read_electonic_event_file(path):
            event = pd.read_csv(path, header=None, sep='\s+', skiprows=[0],on_bad_lines='skip')
            for pmt in range(event.shape[1]):
                event.iloc[0::2, pmt] -= event.iloc[0:400:2, pmt].mean()
                event.iloc[1::2, pmt] -= event.iloc[1:400:2, pmt].mean()
            return event
        event = read_electonic_event_file(path)
    
    if len(event)<1020:
        continue
    
    
    ###========================================================================
    
    
    #Search the impulse in the отклик
    
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
    
    ###============================================
    
    #Save the PMT and bins in with no BG 
    
    saved_pmt=[]
    t_front_reconstr=[]
    x_front_reconstr=[]
    y_front_reconstr=[]
    number_photons=[]
    
    impulse_each_pmt=[]
    cells_each_pmt=[]
    
    for i in range(0,event_diapason.shape[1]):
        event_pmt=event_diapason[i]
        
        #Condition 1 : delete the PMT with max signal < max BG
        if max(event_pmt)<max(event[i][0:400]):
            continue
        
        imax_pmt = event_pmt.idxmax() #index of the max impulse
        start_pmt, stop_pmt = imax_pmt, imax_pmt

        bgd_pmt = max(event[i][0:400])
        while event_pmt[start_pmt] > bgd_pmt and start_pmt>event_pmt.index[0]:
            start_pmt -= 1  
            
        while event_pmt[stop_pmt] > bgd_pmt and stop_pmt<event_pmt.index[-1]:
            stop_pmt += 1
            
        # #Impulse in the finded diapason
        impulse=event_pmt[start_pmt-min(event_pmt.index):stop_pmt-min(event_pmt.index)]
        cells_impulse=cells_diapason[start_pmt-min(event_pmt.index):stop_pmt-min(event_pmt.index)]

        #Condition 2: delete bins < max BG
        impulse_filtrate=impulse[impulse>max(event[i][0:400])] 
        cells_impulse_filtrate=cells_impulse[impulse>max(event[i][0:400])] 

            
       #Condition 3 : delete the PMT with 0 bins
        if len(cells_impulse_filtrate)<1: #if number of bins in PMT less than 4 - low signal
            continue
        
        #Aoment of the arrival of shower in the PMT
        t_arrival=cells_impulse_filtrate[len(impulse_filtrate.cumsum()[impulse_filtrate.cumsum() < impulse_filtrate.sum()*part_of_int])-1]

        #SAVE DATA
        
        t_front_reconstr.append(t_arrival-t_path[i]) #t - obtain front
        x_front_reconstr.append(x_snow[i]) #x - obtain front
        y_front_reconstr.append(y_snow[i]) #y - obtain front
        
        saved_pmt.append(num_pmt[i]) # Nº of the saved PMT with signal
        number_photons.append(impulse_filtrate.sum()) #number of signal photons in each PMT

    ###============================================
    
    ###============================================
    
    #Axis of the shower
    
    x0_front=x_front_reconstr[number_photons.index(max(number_photons))] 
    y0_front=y_front_reconstr[number_photons.index(max(number_photons))]
    t0_front=t_front_reconstr[number_photons.index(max(number_photons))]
    
    d0=np.sqrt(x0_front**2+y0_front**2)
    
    t_fnc=t_front_reconstr-t0_front
    
    #Examination of the t
    
    
    #Delete the isolated dot from the front
    # #fit the 2D front : f(x)=t
    xx=np.arange(min(x_front_reconstr),max(x_front_reconstr),(max(x_front_reconstr)-min(x_front_reconstr))/len(x_front_reconstr)) #approximation - x
    x_y_t=pd.DataFrame({'x_fr':x_front_reconstr,'y_fr':y_front_reconstr,'t_fr':t_fnc,'nbr_phot':number_photons}).sort_values(by=['t_fr']) #x,y,t of the front sorted by t
    
    a_isol,b_isol,c_isol=np.polyfit(np.arange(0,len(x_y_t)),x_y_t['t_fr'],2)
    
    yy=a_isol*(np.arange(0,len(x_y_t)))**2+b_isol*np.arange(0,len(x_y_t))+c_isol #t approx
    
    t_discrep=abs(yy-x_y_t['t_fr'])
    
    x_y_t['discrep_t']=t_discrep
    
    x_front=list(x_y_t['x_fr'][x_y_t['discrep_t']<100])
    y_front=list(x_y_t['y_fr'][x_y_t['discrep_t']<100])
    t_fnc=list(x_y_t['t_fr'][x_y_t['discrep_t']<100])
    number_photons=list(x_y_t['nbr_phot'][x_y_t['discrep_t']<100])
    
    
    
    #Examination of the x
    
    x_t=pd.DataFrame({'x_fr':x_front,'y_fr':y_front,'t_fr':t_fnc,'nbr_phot':number_photons}).sort_values(by=['x_fr']) #x,y,t of the front sorted by t
    
    a_try,b_try,c_try=np.polyfit(x_t['x_fr'],x_t['t_fr'],2)
    yy_try=a_try*(np.array(x_t['x_fr']))**2+b_try*np.array(x_t['x_fr'])+c_try #t approx
    
    x_t['dicrepb']=abs(yy_try-x_t['t_fr'])
    
    x_front=list(x_t['x_fr'][x_t['dicrepb']<200])
    y_front=list(x_t['y_fr'][x_t['dicrepb']<200])
    t_fnc=list(x_t['t_fr'][x_t['dicrepb']<200])
    number_photons=list(x_t['nbr_phot'][x_t['dicrepb']<200])
    
    aaa,bbb,ccc=np.polyfit(x_front,t_fnc,2)
    
    
    #Examination of the y 
    
    x0=np.mean(y_front)
    y0=np.mean(t_fnc)
    
    aaaaa=pd.DataFrame({'x_fr':x_front,'y_fr':y_front,'t_fr':t_fnc})
    
    x_try_left=aaaaa[aaaaa['y_fr']<x0]
    x_left=x_try_left['y_fr'].mean()
    y_left=x_try_left['t_fr'].mean()
    
    x_try_right=aaaaa[aaaaa['y_fr']>x0]
    x_right=x_try_right['y_fr'].mean()
    y_right=x_try_right['t_fr'].mean()
    
    monter_left=abs(y_left-y0)
    monter_right=abs(y_right-y0)
    mean_monter=np.mean([monter_right,monter_left])
    
    # if mean_monter>10:
    #     continue
    #     y1_ok=y_left+mean_monter
        
    #     d1=np.sqrt((x_left-x0)**2+(y1_ok-y0)**2)
    #     d2=np.sqrt((x_left-x0)**2+(y_left-y0)**2)
    #     if d1<d2:
    #         alpha1= np.arccos( d1 / d2)
    #     else:
    #         alpha1= np.arccos( d2 / d1)
    #     x1_ok=np.cos(alpha1)*x_left
    #     y1_ok=np.sin(alpha1)*y_left
        
        
    #     y_front=np.cos(alpha1)*np.array(y_front)
    #     t_fnc=np.cos(alpha1)*np.array(t_fnc)
    
    # nbr_phot=pd.DataFrame({'x_fr':x_front,'y_fr':y_front,'t_fr':t_fnc,'nbr_phot':number_photons}).sort_values(by=['nbr_phot'])
    
    # perc=1#$
    
    # # x_frontt=list(nbr_phot['x_fr'])[len(nbr_phot)//3:]
    # # y_frontt=list(nbr_phot['y_fr'])[len(nbr_phot)//3:]
    # # t_frontt=list(nbr_phot['t_fr'])[len(nbr_phot)//3:]
    
    # x_front=list(nbr_phot['x_fr'])[int(np.floor(len(nbr_phot)*(1-perc))):]
    # y_front=list(nbr_phot['y_fr'])[int(np.floor(len(nbr_phot)*(1-perc))):]
    # t_fnc=list(nbr_phot['t_fr'])[int(np.floor(len(nbr_phot)*(1-perc))):]


    ###============================================
    
    # Algorithm for the retrieval of the direction
    
    ##=== Parameters of the front

    a0=a_a0*np.sqrt(x0_front**2+y0_front**2)+b_a0
    # a0=c_try
    a1_front=0
    # if mean_monter>10:
    #     a2_front=0.0001
    # else:
    #     a2_front=0.001
        

    a0=a_a0*np.sqrt(x0_front**2+y0_front**2)+b_a0
    # a1_front=0
    a2_front=np.arange(0.0001,0.0014,0.0001)
    
    delta_min=[]
    theta_min=[]
    phi_min=[]
    a2_min=[]
        
    for r in range(0,len(a2_front)):
        
        ##=== Algorithm
           
            
            def fnc(theta,phi):#определение функции с не заданными параметрами
                    x_casc=((np.cos(theta)*np.cos(phi))*(x_front-x0_front)+(np.cos(theta)*np.sin(phi))*(y_front-y0_front)) #координаты x фотонов в системе ливня
                    y_casc=(-np.sin(phi)*(x_front-x0_front)+np.cos(phi)*(y_front-y0_front)) #координаты y фотонов в системе ливня
                    z_casc=((np.sin(theta)*np.cos(phi))*(x_front-x0_front)+(np.sin(theta)*np.sin(phi))*(y_front-y0_front)) #координаты z фотонов в системе ливня
                    R=np.sqrt(x_casc**2+y_casc**2) # расстояние фотонов от центра оси в системе ливня
                    tau=a0+a1_front*R+a2_front[r]*R**2+z_casc/c_ns #аппроксимированое время
                    s=(((t_fnc-tau)**2)) #функцию которую надо аппроксимировать
                    Smin=np.sum(s)
                    return Smin
                
            theta_multi_start=[]
            phi_multi_start=[]
            min_s=[]
            #Multistart
            for th in range(0,len(theta_initt)):
                    param_fix=Minuit(fnc, theta=theta_initt[th],phi=phi_initt[th])
                    param_fix.limits['theta']=(0,30/180*np.pi)
                    param_fix.limits['phi']=(0,2*np.pi)
                    param_fix.migrad()
                    
                    theta_multi_start.append(param_fix.values[0]) #obtain theta with multistart
                    phi_multi_start.append(param_fix.values[1])  #obtain phi with multistart
                    min_s.append(param_fix.fval) #obtain min S with multistart
                
            theta=np.array(theta_multi_start)[np.where(np.array(min_s)==min(min_s))][0] #theta with the min S
            phi=np.array(phi_multi_start)[np.where(np.array(min_s)==min(min_s))][0]  #phi with the min S
            
            theta_fit=theta
            phi_fit=phi
            
            theta_real=real_theta[int(path[-8:-5])-1]
            phi_real=real_phi[int(path[-8:-5])-1]

            cx = np.sin(theta_real) * np.cos(phi_real) #theta - истинные углы
            cy = np.sin(theta_real) * np.sin(phi_real)
            cz = np.cos(theta_real)
            cx_fit = np.sin(theta_fit) * np.cos(phi_fit) #theta' -полученные углы моделированием
            cy_fit = np.sin(theta_fit) * np.sin(phi_fit)
            cz_fit = np.cos(theta_fit)
            delta=np.arccos(cx *cx_fit + cy * cy_fit + cz * cz_fit)*180/np.pi
            
            delta_min.append(delta)
            theta_min.append(theta)
            phi_min.append(phi)
            a2_min.append(a2_front[r])
            
#==============

    #SAVE DATA
         
    delta_choose.append(delta_min[delta_min.index(min(delta_min))])
    a2_choose.append(a2_min[delta_min.index(min(delta_min))])
    d0_choose.append(d0)
    mean_monter_choose.append(mean_monter)
    sum_phot_choose.append(sum(number_photons))
    
    files_res.append(path)
    
    
total=np.array([delta_choose,a2_choose,d0_choose,mean_monter_choose,sum_phot_choose]).T

# np.savetxt('choose_a2_500m_P_30PeV_res.csv',total)
# np.savetxt('choose_a2_500m_P_30PeV_files.csv',files_res,fmt='%s')



