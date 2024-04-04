%%time
##=== IMPORT MODULES

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit


###==================================================================================

##=== !!! CHANGE ONLY HERE !!!

#Type of telescope
name_telesc='SPHERE-2' #SPHERE-2, SPHERE-3

#Altitude of the telescope
H=500 #500, 900, 1000 m

#Type of analyse
type_analyse='only_sig'
# type_analyse='sig_&_bg'
# type_analyse='electronic'

###==================================================================================

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
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/2_data_signal/only_sig_30_PeV_P_500m' #Signal only

#Signal with background  
 
    if type_analyse=='sig_&_bg':
        data_photons_bg = pd.read_csv('/Users/clemence/Documents/Scientific_research/Data/for_BG-distribution_photons/{}m'.format(H))

#Electronics

if type_analyse=='electronic':
    # path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/electronics/electronics_30PeV_P_{}m/mosaic_hits_m01_pro_30PeV_10-20_001_c003'.format(H) #Electronic signal #044 #072
    path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/RESULTS/bad_res/mosaic_hits_m01_pro_30PeV_10-20_015_c014' #Electronic signal #044 #072


###==================================================================================

###==================================================================================

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
    x_snow=-b*x_y_mos['x']-a*x_y_mos['x']**2
    y_snow=-b*x_y_mos['y']-a*x_y_mos['y']**2

if name_telesc=='SPHERE-3':
    x_snow=-b*np.array(pixel_data.groupby('segment').mean()['x'])
    y_snow=-b*np.array(pixel_data.groupby('segment').mean()['y'])

#==Translation t from mos to snow (path from mos to snow)
  
t_path=((np.sqrt(H**2+(np.sqrt(np.array(x_snow)**2+np.array(y_snow)**2))**2))/c_ns) #path in ns


###=================================================================================
name_files = glob.glob(os.path.join(path,'*'))

##=== DEFINITION OF THE TYPE OF ANALYSE

#===Save to analyse final results
sum_photons_event=[]
nbr_pmt_optimal=[]
min_err_d0=[]
d0_real_event=[]
#===
for files in range(0,len(name_files)):

    path=name_files[files]
    
    ##=== d0

    analyze_d0=pd.read_csv(path)

    real_x0=float(analyze_d0.columns[0].split()[2])
    real_y0=float(analyze_d0.columns[0].split()[3])

    real_d0=np.sqrt(real_x0**2+real_y0**2)
    
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
    
    # saved_PMT=list(event_diapason.max().where(event[0:400].max()<event_diapason.max()).dropna().index) #saved PMT where max impulse > max BG in the piedestal
    # event_not_BG=event_diapason[saved_PMT].set_index(np.arange(0,len(event_not_BG))) #delete PMT with the BG > max impulse
    
    # bins_not_BG_r=event_not_BG[3][event_not_BG[3].idxmax():].where(event_not_BG[3][event_not_BG[3].idxmax():]>event[3][0:400].max()).dropna()
    # bins_not_BG_l=event_not_BG[3][:event_not_BG[3].idxmax()].where(event_not_BG[3][:event_not_BG[3].idxmax()]>event[3][0:400].max()).dropna()
    
    # bins_not_BG=pd.concat([bins_not_BG_l, bins_not_BG_r])
    # aaaa=[]
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
        number_photons.append(impulse_filtrate.sum()) #number of signal photons in each PMT
    
    ###============================================
    
    #Axis of the shower
    
    x0_front=x_front_reconstr[number_photons.index(max(number_photons))] 
    y0_front=y_front_reconstr[number_photons.index(max(number_photons))]
    t0_front=t_front_reconstr[number_photons.index(max(number_photons))]
    
    d0=np.sqrt(x0_front**2+y0_front**2)
    
    x_y_t_phot=pd.DataFrame({'x_fr':x_front_reconstr,'y_fr':y_front_reconstr,'nbr_phot':number_photons}).sort_values(by=['nbr_phot']) #x,y,t of the front sorted by t
    
    test_nbr=np.arange(1,30,1)
    save_err=[]
    for test in range(0,len(test_nbr)):
        decoup=x_y_t_phot[-test_nbr[test]:]
        tryy_x,tryy_y=decoup['x_fr']*decoup['nbr_phot'],decoup['y_fr']*decoup['nbr_phot']
        x0_other=tryy_x.sum()/decoup['nbr_phot'].sum()
        y0_other=tryy_y.sum()/decoup['nbr_phot'].sum()
        d0_other=np.sqrt(x0_other**2+y0_other**2)
        save_err.append(abs(d0_other-real_d0))
        
    test_nbr_ok=test_nbr[save_err.index(min(save_err))]
    
    d0_real_event.append(real_d0)
    nbr_pmt_optimal.append(test_nbr_ok)
    min_err_d0.append(min(save_err))
    sum_photons_event.append(sum(number_photons))

# ###==================================================================

# #========= GRAPHICS

# plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in"})

# x=x_front_reconstr
# y=y_front_reconstr
# z=number_photons
# plt.scatter(x, y, z, cmap='Wistia')
# cb = plt.colorbar()
# cb.set_label('number of photons')
# plt.scatter(x0_front,y0_front)
# plt.scatter(real_x0,real_y0,color='yellow')
# plt.scatter(x0_other,y0_other,color='purple')
# plt.xlabel('x, m')
# plt.ylabel('y, m')
# plt.show()

# plt.scatter(test_nbr,save_err,label='real d0 = {:.2f} m \n ∑nbr phot={}'.format(real_d0,verifyy))
# plt.scatter(test_nbr_ok,min(save_err),color='red',label='min error of d0')
# plt.ylabel('err d0, m')
# plt.xlabel('N PMT')
# plt.legend()
# plt.show()
