%%time
##=== IMPORT MODULES

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit


###=========================================================================

##=== !!! CHANGE ONLY HERE !!!

#Type of telescope
name_telesc='SPHERE-2' #SPHERE-2, SPHERE-3

#Altitude of the telescope
H=500 #500, 900, 1000 m

#Type of analyse
type_analyse='only_sig'
# type_analyse='sig_&_bg'
# type_analyse='electronic'

#Energy
En=10

###=========================================================================
###=========================================================================

##=== OPEN FILES

#===============
#== For all types of analyse

#PMT on mosaic: x (mm),y (mm), num
if name_telesc=='SPHERE-2':
    x_y_mos=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_background/Data/Data_init/mosaic_pmt_coords_df.csv')
    total_number_of_pmt=109

if name_telesc=='SPHERE-3':
    pixel_data = pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/SPHERE3_pixel_data_A.dat',header=None, names=['x','y','z','nx','ny','nz','theta','phi'], sep='\s+') #number of PMT
    pixel_data['segment'] = pixel_data.index // 7 #== add num segment
    pixel_data['pixel'] =  pixel_data.index % 7 #== add num pixel
    total_number_of_pmt=379

#===============

#Open the number of the PMT around a PMT

circle_of_pmt_around_pmt=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/circle_of_pmt_around_pmt.csv',)#===============

#===============

#Path of the directory of the events

#SPHERE 2

#Only signal

if type_analyse!='electronic':
    # path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/2_data_signal/only_sig_{}PeV_P_{}m/mosaic_hits_m01_pro_{}PeV_10-20_001_c090'.format(En,H,En) #Signal only
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/bad_res_500m/mosaic_hits_m01_pro_10PeV_10-20_058_c030'.format(En,H,En) #Signal only

#Signal with background  
 
    if type_analyse=='sig_&_bg':
        data_photons_bg = pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/for_BG-distribution_photons/{}m'.format(H))

#Electronics

if type_analyse=='electronic':
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/1_data_electronics/electronics_{}PeV_P_{}m/mosaic_hits_m01_pro_{}PeV_10-20_001_c001'.format(En,H,En) #Electronic signal #044 #072

#SPHERE 3

if name_telesc=='SPHERE-3' and type_analyse=='only_sig':
    path='/Users/clemence/Documents/Scientific_research/Sphere_3/Data/{}m_{}PeV_Fe/moshits_Q2_atm01_5626_{}PeV_15_001_c006'.format(H,En,En)

real_theta=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_theta',header=None,sep='\s+')[0]) #number of PMT
real_phi=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_phi',header=None,sep='\s+')[0]) #number of PMT

theta_real=real_theta[57]
phi_real=real_theta[57]

###=========================================================================
###=========================================================================

##=== FUNCTION

# __________________

# Return event with the type of analyse
if name_telesc=='SPHERE-2':
    num_t_event=4
if name_telesc=='SPHERE-3': 
    num_t_event=5

def read_mosaic_hits_file(path, step = 12.5):

    a = pd.read_csv(path, header=None, skiprows=[0],sep='\s+') 
    a[num_t_event]=a[num_t_event]-min(a[num_t_event])  
    pd_event=pd.DataFrame({'pmt_ev':list(a[0]),'t_ev':list(a[num_t_event])}) #creation of DataFrame to facilite calculations
    q = np.zeros((total_number_of_pmt, 1020))  # empty array109 PMT * 1020 cells
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
    a[num_t_event]=a[num_t_event]-min(a[num_t_event]) # номер столбца со временем - 4
       
    pd_event=pd.DataFrame({'pmt_ev':list(a[0]),'t_ev':list(a[num_t_event])}) #creation of DataFrame to facilite calculations
    q = np.zeros((total_number_of_pmt, 1020))  # empty array109 PMT * 1020 cells
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

# __________________

#Return the diapason of the impulse 

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

# __________________

#Return the axis of the shower

def axis(x,y,t,N):
    x0_front=x_front_reconstr[number_photons_i.index(max(number_photons_i))] 
    y0_front=y_front_reconstr[number_photons_i.index(max(number_photons_i))]
    t0_front=t_front_reconstr[number_photons_i.index(max(number_photons_i))]
    
    d0_from_center=np.sqrt(x0_front**2+y0_front**2)
    
    d_x_to_x_max=np.sqrt( (np.array(x_front_reconstr)-x0_front)**2  + (np.array(y_front_reconstr)-y0_front)**2  )
    
    if d0_from_center>d0_before_c and d0_from_center<d0_center:
        max_d=d_six 
    
    elif d0_from_center<d0_before_c:
        max_d=d_double
        
    else:
        print('no this')
        
    x_circle_around_max=np.array(x_front_reconstr)[np.where(d_x_to_x_max<max_d)[0]]
    y_circle_around_max=np.array(y_front_reconstr)[np.where(d_x_to_x_max<max_d)[0]]
    t_circle_around_max=np.array(t_front_reconstr)[np.where(d_x_to_x_max<max_d)[0]]
    number_photons_around_max=np.array(number_photons_i)[np.where(d_x_to_x_max<max_d)[0]]
    
    x0=np.sum(x_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max)
    y0=np.sum(y_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max)
    t0=np.sum(t_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max)

    return x0, y0, t0,d0_from_center, x0_front, y0_front

# __________________

#Return the number of the PMT of isolated dots
def dots_t_increase(x,y,t,N,pmt):
     
    #Analyse t increase
    
    t_fnc_reconstr_d0=np.array(t)-t0


    df_analyse_t=pd.DataFrame({'x':x,'y':y,'t':t_fnc_reconstr_d0,'N_phot':N,'PMT':pmt}).sort_values(by=['t'])
       
    a_t,b_t,c_t=np.polyfit(np.arange(0,len(t_fnc_reconstr_d0)),df_analyse_t['t'],2,w=df_analyse_t['N_phot'])
    t_fnc_reconstr_approx=a_t*np.arange(0,len(t_fnc_reconstr_d0))**2+b_t*np.arange(0,len(t_fnc_reconstr_d0))+c_t

    df_analyse_t['delta_t']=abs(t_fnc_reconstr_approx-df_analyse_t['t']) #add delta t

    x_front_after_t=df_analyse_t['x'][df_analyse_t['delta_t']<df_analyse_t['delta_t'].mean()*3]
    y_front_after_t=df_analyse_t['y'][df_analyse_t['delta_t']<df_analyse_t['delta_t'].mean()*3]
    t_front_after_t=df_analyse_t['t'][df_analyse_t['delta_t']<df_analyse_t['delta_t'].mean()*3]

    saved_pmt_after_t=df_analyse_t['PMT'][df_analyse_t['delta_t']<df_analyse_t['delta_t'].mean()*3]

    pmt_isol_dots=list(df_analyse_t['PMT'][df_analyse_t['delta_t']>df_analyse_t['delta_t'].mean()*3])

    return saved_pmt_after_t, t_front_after_t, pmt_isol_dots

def dots_t_sinusoide(saved_pmt_after_t,t_front_after_t):
    new_t_sinusoide=t_front_after_t.sort_index()
    saved_pmt_after_t=list(saved_pmt_after_t.sort_index())

    circle_around_PMT_with_t=circle_of_pmt_around_pmt.loc[~circle_of_pmt_around_pmt['PMT'].isin(saved_pmt_after_t)== False].set_index(np.arange(0,len(saved_pmt_after_t)))
    circle_around_PMT_with_t['t_after_t']=list(new_t_sinusoide)
        
    def index(which):
        t_above=[]
        for i in range(0,len(circle_around_PMT_with_t)):
            index_above=list(circle_around_PMT_with_t[circle_around_PMT_with_t['PMT']==which[i]].index)
            
            if len(index_above)==0:
                t_above.append(np.nan)
            else:
                t_above.append(circle_around_PMT_with_t['t_after_t'][index_above[0]])
            
        return t_above
            
    circle_around_PMT_with_t['t_above']=index(circle_around_PMT_with_t['above'])
    circle_around_PMT_with_t['t_below']=index(circle_around_PMT_with_t['below'])
    circle_around_PMT_with_t['t_right']=index(circle_around_PMT_with_t['right'])
    circle_around_PMT_with_t['t_left']=index(circle_around_PMT_with_t['left'])
    circle_around_PMT_with_t['t_diag_left']=index(circle_around_PMT_with_t['diag_left'])
    circle_around_PMT_with_t['t_diag_right']=index(circle_around_PMT_with_t['diag_right'])


    circle_around_PMT_with_t['soustr_above']=abs(circle_around_PMT_with_t['t_after_t']-circle_around_PMT_with_t['t_above'])
    circle_around_PMT_with_t['soustr_below']=abs(circle_around_PMT_with_t['t_after_t']-circle_around_PMT_with_t['t_below'])
    circle_around_PMT_with_t['soustr_right']=abs(circle_around_PMT_with_t['t_after_t']-circle_around_PMT_with_t['t_right'])
    circle_around_PMT_with_t['soustr_left']=abs(circle_around_PMT_with_t['t_after_t']-circle_around_PMT_with_t['t_left'])
    circle_around_PMT_with_t['soustr_diag_left']=abs(circle_around_PMT_with_t['t_after_t']-circle_around_PMT_with_t['t_diag_left'])
    circle_around_PMT_with_t['soustr_diag_right']=abs(circle_around_PMT_with_t['t_after_t']-circle_around_PMT_with_t['t_diag_right'])

    circle_around_PMT_with_t['min_dist_neighb']=circle_around_PMT_with_t[['soustr_above','soustr_below','soustr_right','soustr_left','soustr_diag_left','soustr_diag_right']].min(axis=1)
    circle_around_PMT_with_t.dropna(subset=['min_dist_neighb'])

    pmt_isol_dots_2nd=list(circle_around_PMT_with_t['PMT'][circle_around_PMT_with_t['min_dist_neighb'] > circle_around_PMT_with_t['min_dist_neighb'].mean()*4])
    # pmt_dots=np.sort(pmt_isol_dots+pmt_isol_dots_2nd)

    return pmt_isol_dots_2nd, circle_around_PMT_with_t

# __________________

# Return the signal from the background (N, t)

def signal_in_pmt(event_pmt):
    imax_pmt = event_pmt.idxmax() #index of the max impulse
    start_pmt, stop_pmt = imax_pmt, imax_pmt
     
    if type_analyse=='only_sig':
        bgd_pmt=0
    else:
        bgd_pmt = np.mean(noise_widnow[noise_widnow>0])

     
    while event_pmt[start_pmt] > bgd_pmt and start_pmt>event_pmt.index[0]:
        start_pmt -= 1  
         
    while event_pmt[stop_pmt] > bgd_pmt and stop_pmt<event_pmt.index[-1]:
        stop_pmt += 1
         
    # #Impulse in the finded diapason
    impulse_filtrate=event_pmt[start_pmt+1-min(event_pmt.index):stop_pmt-min(event_pmt.index)]
    cells_impulse_filtrate=cells_diapason[start_pmt+1-min(event_pmt.index):stop_pmt-min(event_pmt.index)]
    return impulse_filtrate, cells_impulse_filtrate


def angles(x_front,y_front,t_fnc):
    
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
num_pmt=np.arange(0,total_number_of_pmt)

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

part_of_int=0.5


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

theta_real=0.3332 #[files 15]
phi_real=5.7728 #[files 15]

###=========================================================================
###=========================================================================


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


# Calculation of the diameter of the mosaic on snow

diameter=np.sqrt( (x_snow[104]-x_snow[95])**2+ (y_snow[104]-y_snow[95])**2  )
d_diameter=np.linspace(0,diameter,total_number_of_pmt)

###=========================================================================
###=========================================================================

##=== DEFINITION OF THE TYPE OF ANALYSE


# Only signal
if type_analyse=='only_sig':
    event = read_mosaic_hits_file(path)
    
# Signal with the background
if type_analyse=='sig_&_bg':
    event = read_mosaic_hits_file_background(path)

#Signal with the electronic
if type_analyse=='electronic':
    event = read_electonic_event_file(path)



###========================================================================
###========================================================================
###========================================================================


## START OF THE ALGORITHM

#Return the diapason of the impulse 

found_diapason=diapason_impulse(event)
event_diapason,cells_diapason = found_diapason[0],found_diapason[1]
event_diapason.index=np.arange(len(event_diapason))
###============================================

#Save the PMT and bins in with no BG 

saved_pmt=[]
t_front_reconstr=[]
x_front_reconstr=[]
y_front_reconstr=[]
number_photons_i=[]

impulse_each_pmt=[]
cells_each_pmt=[]
cells_impulse_filtrate_each_pmt=[]
impulse_filtrate_each_pmt=[]

event_window_each_pmt=[]
event_each_pmt=[]
noise_window_each_pmt=[]

index_sig_max=[]

# d_diameter_each_pmt=[]

event_diapason[0].index=np.arange(len(event_diapason))
event_pmt=event_diapason[0]

event_window=event_pmt.rolling(window=4).sum().dropna()
noise_widnow=event[0][0:400].rolling(window=4).sum().dropna()

index_max_0=event_window.idxmax()
 
for i in range(0,event_diapason.shape[1]):
    
    event_diapason[i].index=np.arange(len(event_diapason))
    event_pmt=event_diapason[i]
    
    event_window=event_pmt.rolling(window=4).sum().dropna()
    noise_widnow=event[i][0:400].rolling(window=4).sum().dropna()
    # event_pmt.append(num_pmt[i])
    
    #___________________

    #Condition 1 for PMT: delete the PMT with max signal < max BG

    event_pmt=event_window
    if max(event_pmt)<np.mean(noise_widnow[noise_widnow>0]):
        continue
    
    #___________________
    
    # Return the signal from the background (N, t)

    signal_in_pmt_t_N=signal_in_pmt(event_pmt)
    impulse_filtrate=signal_in_pmt_t_N[0]
    cells_impulse_filtrate=signal_in_pmt_t_N[1]
    
    #___________________
        
   #Condition 3 for bins: delete the PMT with 0 bins
   
    if len(cells_impulse_filtrate)<2: 
        continue
    if impulse_filtrate.max()<4:
        continue
    #___________________

    #Aoment of the arrival of shower in the PMT
    t_arrival=cells_impulse_filtrate[len(impulse_filtrate.cumsum()[impulse_filtrate.cumsum() < impulse_filtrate.sum()*part_of_int])-1]

    #_______________
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
    
    # index_max_0=new_index_max

    
    # d_diameter_each_pmt.append(d_diameter[i])
    
# ======================================================================================
# ======================================================================================

#============== Axis

x0_y0_t0=axis(x_front_reconstr,y_front_reconstr,t_front_reconstr,number_photons_i)
x0, y0, t0,d0_front= x0_y0_t0[0], x0_y0_t0[1], x0_y0_t0[2], x0_y0_t0[3]
x0_front, y0_front = x0_y0_t0[4], x0_y0_t0[5]

x0_front=x_front_reconstr[number_photons_i.index(max(number_photons_i))] 
y0_front=y_front_reconstr[number_photons_i.index(max(number_photons_i))]
t0_front=t_front_reconstr[number_photons_i.index(max(number_photons_i))]
    
d0=np.sqrt(x0**2+y0**2)

analyze_d0=pd.read_csv(path)
real_x0=float(analyze_d0.columns[0].split()[2])
real_y0=float(analyze_d0.columns[0].split()[3])

real_d0=np.sqrt(real_x0**2+real_y0**2)

# ============ Return number of PMT isolated dots

#Analysis method by ascent

res_t_by_increase=dots_t_increase(x_front_reconstr,y_front_reconstr,t_front_reconstr,number_photons_i,saved_pmt)
pmt_by_increase=res_t_by_increase[0]
t_by_increase=res_t_by_increase[1]
dots_t_by_increase=res_t_by_increase[2] 

#Analysis method by sinusoide

res_t_by_sinusoide=dots_t_sinusoide(pmt_by_increase,t_by_increase)
dots_t_by_sinusoide=res_t_by_sinusoide[0]
circle_around_PMT_with_t_sinusoide=res_t_by_sinusoide[1]

#Total of the isolated PMT

pmt_dots=dots_t_by_increase+dots_t_by_sinusoide

# ============ Return new x, y, t of the front without isolated dots

val=pd.DataFrame({'PMT':saved_pmt, 'x':x_front_reconstr, 'y':y_front_reconstr, 't':t_front_reconstr-t0, 'N_phot':number_photons_i})

x_front=list(val.loc[~val['PMT'].isin(pmt_dots)]['x'])
y_front=list(val.loc[~val['PMT'].isin(pmt_dots)]['y'])
N_photons=list(val.loc[~val['PMT'].isin(pmt_dots)]['N_phot'])
t_front=val.loc[~val['PMT'].isin(pmt_dots)]['t']
pmt_new=val.loc[~val['PMT'].isin(pmt_dots)]['PMT']

#============= Delete the isolated PMT on the mosaic (without any of the 6 neighbours around)


circle_around_PMT_with_t_2=dots_t_sinusoide(pmt_new,t_front)[1] #df total
circle_around_PMT_with_t_2['x']=x_front #add x front
circle_around_PMT_with_t_2['y']=y_front #add y front
circle_around_PMT_with_t_2['N']=N_photons #add y front

circle_around_PMT_with_t_2=circle_around_PMT_with_t_2.dropna(subset=['min_dist_neighb']) #delete the PMT with no neighbour

#t fromt

circle_around_PMT_with_t_2=circle_around_PMT_with_t_2.sort_values(by='N',ascending=False)

part_front=-1
x_front=np.array(circle_around_PMT_with_t_2['x'][0:part_front]) #save values
y_front=np.array(circle_around_PMT_with_t_2['y'][0:part_front]) 
t_front=np.array(circle_around_PMT_with_t_2['t_after_t'][0:part_front])
pmt_new=np.array(circle_around_PMT_with_t_2['PMT'][0:part_front])
N_photons=np.array(circle_around_PMT_with_t_2['N'][0:part_front])



#  # ============== Return the ANGLES

both_angles=angles(x_front,y_front,t_front)
theta=both_angles[0]
phi=both_angles[1]

# #========================================== Delta

print('theta_fit=',theta*180/np.pi)
print('phi_fit=',phi*180/np.pi)

print('theta_real=',theta_real*180/np.pi)
print('phi_real=',phi_real*180/np.pi)

theta_fit=theta
phi_fit=phi
cx = np.sin(theta_real) * np.cos(phi_real) #theta - истинные углы
cy = np.sin(theta_real) * np.sin(phi_real)
cz = np.cos(theta_real)
cx_fit = np.sin(theta_fit) * np.cos(phi_fit) #theta' -полученные углы моделированием
cy_fit = np.sin(theta_fit) * np.sin(phi_fit)
cz_fit = np.cos(theta_fit)
delta=np.arccos(cx *cx_fit + cy * cy_fit + cz * cz_fit)*180/np.pi

   
print(delta)


# #========================== Tau
# a0,a1,a2=param_fix.values[2], param_fix.values[3], param_fix.values[4]
# theta,phi=theta_real, phi_real
# x_casc=((np.cos(theta)*np.cos(phi))*(x_front-x0)+(np.cos(theta)*np.sin(phi))*(y_front-y0)) #координаты x фотонов в системе ливня
# y_casc=(-np.sin(phi)*(x_front-x0)+np.cos(phi)*(y_front-y0)) #координаты y фотонов в системе ливня
# z_casc=((np.sin(theta)*np.cos(phi))*(x_front-x0)+(np.sin(theta)*np.sin(phi))*(y_front-y0)) #координаты z фотонов в системе ливня
# R=np.sqrt(x_casc**2+y_casc**2) # расстояние фотонов от центра оси в системе ливня
# tau=a0+a1*R+a2*R**2 +z_casc/c_ns

# ###==================================================================
# ###==================================================================

# #========= GRAPHICS

plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in"})


# plt.bar(cells,event.sum(axis=1) ,width=50)
# # plt.axvline(x=cells[start], color='r')
# # plt.axvline(x=cells[stop],  color='r')
# plt.grid()
# plt.show()

num=60
plt.bar(cells_impulse_filtrate_each_pmt[num],impulse_filtrate_each_pmt[num],width=8)
plt.ylabel('n')
plt.xlabel('t, ns')
plt.xlim(cells_diapason[0],cells_diapason[-1])
plt.legend(title='num PMT:{}'.format(saved_pmt[num]))
plt.show()

plt.bar(cells_diapason[:-3],event_window_each_pmt[num],width=8)
plt.legend(title='sliding window \nnum PMT:{}'.format(saved_pmt[num]))
plt.show()

# # plt.bar(cells_diapason,event_each_pmt[num],width=8)
# # plt.legend(title='num PMT:{}'.format(saved_pmt[num]))
# # plt.show()

# plt.bar(cells_diapason[:-3],event_bad_0[0],width=8)
# plt.legend(title='num PMT:{} \n false impulse = 0'.format(saved_pmt[num]))
# plt.show()

# # plt.bar(cells[0:397],noise_window_each_pmt[num],width=10)
# # plt.legend(title='background \nnum PMT:{}'.format(saved_pmt[num]))
# # plt.show()

# #============================================================================= Front x. y t 
# #===== Initial 

# fig = plt.figure(figsize=(6,12))
# ax = fig.add_subplot(111, projection="3d")

# x_i=x_front_reconstr
# y_i=y_front_reconstr
# t_i=t_front_reconstr

# ax.scatter(x_i,y_i,t_i,s=10,label='front',color='red')

# ax.set_zlim(min(t_i),max(t_i))
# ax.view_init(0,270)  #gauche: droite gauche/droite: vision haut bas
# ax.set_xlabel("x, м")
# ax.set_ylabel("y, м")
# ax.set_zlabel("t, ns")
# # ax.legend(title='d0={:.2f} m / weight=1'.format(d0,len(x_front_reconstr)))
# ax.legend(title='initial')
# plt.show()

# # #===== after  after

# fig = plt.figure(figsize=(6,12))
# ax = fig.add_subplot(111, projection="3d")

# x=x_front
# y=y_front
# t=t_front


# ax.scatter(x,y,t,s=10,color='green', label='front')
# ax.scatter(x,y,tau,s=10,color='red', label='fitted front')
# # ax.scatter(x,y,tau,s=10,label='approx',color='red')
# # ax.set_zlim(min(t_i),max(t_i))
# ax.view_init(0,270)  #gauche: droite gauche/droite: vision haut bas
# ax.set_xlabel("x, м")
# ax.set_ylabel("y, м")
# ax.set_zlabel("t, ns")
# ax.legend(title='d0={:.2f} m \n{} PMT total \ndelta={:.2f}º'.format(d0,len(x_front),delta))
# # ax.set_zlim(-200,300)
# plt.show()

#==================================================================================

#Analyse d0
# plt.scatter(x_front,y_front)
# plt.scatter(real_x0,real_y0,color='red',label='real')
# plt.scatter(x0_front,y0_front,color='orange',label='by max')
# plt.scatter(x0,y0,color='blue',label='by gravity center')
# plt.axis('equal')
# plt.legend(fontsize=8,title='∆d={:.2f} m \nmax photons={}'.format(abs(real_d0-d0),max(number_photons_i)),title_fontsize=8)
# plt.show()

#Show the all the PMT of mosaic on snow

# plt.scatter(x_snow,y_snow)
# plt.axis('equal')
# plt.show()


#==================================================================================

#ANALYSE T increasing: first analyse to eliminate 

# plt.scatter(np.arange(0,len(df_analyse_t)),df_analyse_t['t'])
# plt.plot(np.arange(0,len(df_analyse_t)),t_fnc_reconstr_approx,c='red')
# plt.show()

#This one

# plt.scatter(np.arange(0,len(t_front)),np.sort(t_front))
# aaaa=t_N['d']
# plt.scatter(xx_approx,list(t_N['t']))
# plt.xlabel('d, m \nAxis of the front development')
# plt.ylabel('t, ns')
# # plt.xlim(-20,400)
# plt.ylim(-20,300)
# plt.legend(title='a={:.4f}, b={:.2f}, c={:.2f} \n Total number of PMT={}'.format(a_t,b_t,c_t, len(x_front)))
# plt.plot(xx_approx,t_fnc_reconstr_approx_same,c='red')
# plt.show()

# # #======================================

# # #SINUSOIDE

plt.scatter(saved_pmt, np.array(t_front_reconstr))
plt.xlabel('N PMT')
plt.ylim(None,4200)
plt.ylabel('t, ns')
plt.show()

plt.scatter(pmt_new,t_front)
plt.ylim(-200,1200)
plt.xlabel('N PMT')
plt.ylabel('t, ns')
plt.show()


