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
En=30

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
#== For 1 type of analyse

#Path of the directory of the events

#Only signal

if type_analyse!='electronic':
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/2_data_signal/only_sig_{}PeV_P_500m/mosaic_hits_m01_pro_{}PeV_10-20_001_c081'.format(En,En) #Signal only

#Signal with background  
 
    if type_analyse=='sig_&_bg':
        data_photons_bg = pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/for_BG-distribution_photons/{}m'.format(H))

#Electronics

if type_analyse=='electronic':
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/1_data_electronics/electronics_{}PeV_P_{}m/mosaic_hits_m01_pro_{}PeV_10-20_001_c001'.format(En,H,En) #Electronic signal #044 #072
    # path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/1_data_electronics/electronics_30PeV_P_500m/mosaic_hits_m01_pro_30PeV_10-20_001_c049' #Electronic signal #044 #072
    # path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/bad_res/mosaic_hits_m01_Fe_10PeV_10-20_038_c005' #Electronic signal #044 #072

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

part_of_int=0.5

#=====FNC

#a0 fnc
# a_a0=-0.630484 #more in the article "Алгоритм восстановления направления прихода ШАЛ для телескопа СФЕРА"
# b_a0=53.541113661

# a2=0.0005

#Paramatrization
theta_init=np.array([0.05,0,1,0.15, 0.2,0.25, 0.3,0.35]).tolist() #in rad
phi_init=np.arange(0.5,6.5,0.5) #in rad

def theta_initt(theta_init):
    return theta_init*len(theta_init)

theta_initt=theta_initt(theta_init)
phi_initt=phi_init.tolist()*len(theta_init)


#=====DELTA

real_theta_deg=0.1745*180/np.pi
real_phi_deg=3.0564*180/np.pi

theta_real=0.3277 #[0]
phi_real=3.0564 #[0]

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

#==Fonction of the approximation of the front

def tau_fnc(theta,phi,a0,a1,a2):
    x_casc=((np.cos(theta)*np.cos(phi))*(x_front-x0_front)+(np.cos(theta)*np.sin(phi))*(y_front-y0_front)) #координаты x фотонов в системе ливня
    y_casc=(-np.sin(phi)*(x_front-x0_front)+np.cos(phi)*(y_front-y0_front)) #координаты y фотонов в системе ливня
    z_casc=((np.sin(theta)*np.cos(phi))*(x_front-x0_front)+(np.sin(theta)*np.sin(phi))*(y_front-y0_front)) #координаты z фотонов в системе ливня
    R=np.sqrt(x_casc**2+y_casc**2) # расстояние фотонов от центра оси в системе ливня
    return a0+a1*R+a2*R**2+z_casc/c_ns


###=========================================================================

##=== DEFINITION OF THE TYPE OF ANALYSE

#Only signal

if type_analyse!='electronic':

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
number_photons_i=[]

impulse_each_pmt=[]
cells_each_pmt=[]
cells_impulse_filtrate_each_pmt=[]
impulse_filtrate_each_pmt=[]

for i in range(0,event_diapason.shape[1]):
    event_pmt=event_diapason[i]
    
    #Condition 1 for PMT: delete the PMT with max signal < max BG
    
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

    #Condition 2 for bins: delete bins < max BG
    
    impulse_filtrate=impulse[impulse>max(event[i][0:400])] 
    cells_impulse_filtrate=cells_impulse[impulse>max(event[i][0:400])] 

        
   #Condition 3 for bins: delete the PMT with 0 bins
   
    if len(cells_impulse_filtrate)<1: #if number of bins in PMT less than 4 - low signal
        continue
    
    #Aoment of the arrival of shower in the PMT
    t_arrival=cells_impulse_filtrate[len(impulse_filtrate.cumsum()[impulse_filtrate.cumsum() < impulse_filtrate.sum()*part_of_int])-1]

    #SAVE DATA
    
    t_front_reconstr.append(t_arrival-t_path[i]) #t - obtain front
    x_front_reconstr.append(x_snow[i]) #x - obtain front
    y_front_reconstr.append(y_snow[i]) #y - obtain front
    
    saved_pmt.append(num_pmt[i]) # Nº of the saved PMT with signal
    number_photons_i.append(impulse_filtrate.sum()) #number of signal photons in each PMT

    cells_impulse_filtrate_each_pmt.append(cells_impulse_filtrate) #to see image
    impulse_filtrate_each_pmt.append(impulse_filtrate) #to see image

###============================================

#Axis of the shower

# number_photons=np.array(number_photons)


x0_front=x_front_reconstr[number_photons_i.index(max(number_photons_i))] 
y0_front=y_front_reconstr[number_photons_i.index(max(number_photons_i))]
t0_front=t_front_reconstr[number_photons_i.index(max(number_photons_i))]

d0_from_center=np.sqrt(x0_front**2+y0_front**2)

d_x_to_x_max=np.sqrt( (np.array(x_front_reconstr)-x0_front)**2  + (np.array(y_front_reconstr)-y0_front)**2  )

if H==500:
    if d0_from_center>140 and d0_from_center<180:
        d_max=50
    if d0_from_center<140:
       d_max=100

       
if H==900:
    if d0_from_center>2 and d0_from_center<330:
        d_max=100
    if d0_from_center<260:
       d_max=150
    

x_circle_around_max=np.array(x_front_reconstr)[np.where(d_x_to_x_max<d_max)[0]]
y_circle_around_max=np.array(y_front_reconstr)[np.where(d_x_to_x_max<d_max)[0]]
number_photons_around_max=np.array(number_photons_i)[np.where(d_x_to_x_max<d_max)[0]]

x0_grav_center=np.sum(x_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max)
y0_grav_center=np.sum(y_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max)

x0_real=float(pd.read_csv(path).columns.tolist()[0].split()[2])
y0_real=float(pd.read_csv(path).columns.tolist()[0].split()[3])

all_window=event_diapason.rolling(window=4).sum().dropna()
all_window.index=np.arange(0,len(all_window))

x0_center=list(x_y_mos['x'])[all_window.max().idxmax()]
y0_center=list(x_y_mos['y'])[all_window.max().idxmax()]

print(x0_center,y0_center)

###==================================================================

#========= GRAPHICS

plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in"})

#Mosaic with last circle, before last circle, and interior

# plt.scatter(x_y_mos['x'],x_y_mos['y'])
# plt.scatter(-x0_center,-y0_center)
# plt.axis('square')
# plt.show()

#==========

#3 ZONES OF MOSAIC

exterior=180
before_last=140

d_x_to_center=np.sqrt( np.array(x_snow)**2  + np.array(y_snow)**2  )
 
x_last_circle=np.array(list(x_y_mos['x']))[np.where(d_x_to_center>exterior)[0]]
y_last_circle=np.array(list(x_y_mos['y']))[np.where(d_x_to_center>exterior)[0]]

x_before_last_circle=np.array(list(x_y_mos['x']))[np.where(d_x_to_center>before_last)[0]]
y_before_last_circle=np.array(list(x_y_mos['y']))[np.where(d_x_to_center>before_last)[0]]

plt.scatter(list(x_y_mos['x']),list(x_y_mos['y']),c='greenyellow',label='interior')
plt.scatter(x_before_last_circle,y_before_last_circle,label='before \nlast circle',color='blue')
plt.scatter(x_last_circle,y_last_circle,label='last circle',color='orange')
plt.axis('square')
plt.xlim(-350,350)
plt.ylim(-350,350)
plt.xlabel('x на мозаике, мм')
plt.ylabel('y на мозаике, мм')
plt.legend(fontsize=8)
plt.show()

#Circle around

#Modelisation of 3 max

#Snow
# x0_1, y0_1=-20.311293184751595, -35.180191474677365 #interior
# x0_2, y0_2=-79.79018458238589, 138.2006496359943 #before last
# x0_3, y0_3=39.903004803815, 204.21223889454336 #last

#mosaic
x0_1, y0_1=-22.186325, -115.76552 #interior
x0_2, y0_2=-88.154778,152.68855 #before last
x0_3, y0_3= 44.086131, 225.62029 #last


d_x_to_x_max_1=np.sqrt( (np.array(x_y_mos['x'])-x0_1)**2  + (np.array(x_y_mos['y'])-y0_1)**2  )
d_x_to_x_max_2=np.sqrt( (np.array(x_y_mos['x'])-x0_2)**2  + (np.array(x_y_mos['y'])-y0_2)**2  )

x_circle_around_max_1=np.array(x_y_mos['x'])[np.where(d_x_to_x_max_1<100)[0]] #interior
y_circle_around_max_1=np.array(x_y_mos['y'])[np.where(d_x_to_x_max_1<100)[0]] #interior

x_circle_around_max_2=np.array(x_y_mos['x'])[np.where(d_x_to_x_max_2<50)[0]] #before last
y_circle_around_max_2=np.array(x_y_mos['y'])[np.where(d_x_to_x_max_2<50)[0]] #before last

# 2 types of circle around 2 different max

plt.scatter(x_y_mos['x'],x_y_mos['y'])
plt.scatter(x_circle_around_max_1,y_circle_around_max_1,c='red')
plt.scatter(x_circle_around_max_2,y_circle_around_max_2,c='red')
plt.scatter(x0_1,y0_1,c='lime')
plt.scatter(x0_2,y0_2,c='darkblue')
plt.scatter(x0_3,y0_3,c='orange')
plt.axis('square')
plt.xlim(-350,350)
plt.ylim(-350,350)
plt.xlabel('x на снегу, мм')
plt.ylabel('y на снегу, мм')
plt.legend(fontsize=9)
plt.show()


exterior=180
before_last=140

d_x_to_center=np.sqrt( np.array(x_snow)**2  + np.array(y_snow)**2  )
 
x_last_circle=np.array(list(x_y_mos['x']))[np.where(d_x_to_center>exterior)[0]]
y_last_circle=np.array(list(x_y_mos['y']))[np.where(d_x_to_center>exterior)[0]]

x_before_last_circle=np.array(list(x_y_mos['x']))[np.where(d_x_to_center>before_last)[0]]
y_before_last_circle=np.array(list(x_y_mos['y']))[np.where(d_x_to_center>before_last)[0]]

plt.scatter(list(x_y_mos['x']),list(x_y_mos['y']),c='greenyellow',label='Внутренняя \nчасть')
plt.scatter(x_before_last_circle,y_before_last_circle,label='Предпоследнее \nкольцо',color='blue')
plt.scatter(x_last_circle,y_last_circle,label='Последнее \nкольцо',color='orange')
# plt.scatter(x_y_mos['x'],x_y_mos['y'])

plt.scatter(x_circle_around_max_1,y_circle_around_max_1,c='red')
plt.scatter(x_circle_around_max_2,y_circle_around_max_2,c='red')
plt.scatter(x0_1,y0_1,c='lime')
plt.scatter(x0_2,y0_2,c='darkblue')
plt.scatter(x0_3,y0_3,c='orange')
plt.axis('square')
plt.xlim(-400,400)
plt.ylim(-400,400)
plt.xlabel('x на мозаике, мм')
plt.ylabel('y на мозаике, мм')
plt.legend(fontsize=6)
# plt.savefig('foo.pdf')
plt.show()

