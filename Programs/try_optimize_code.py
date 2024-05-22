%%time

###=========================================================================
###=========================================================================

### IMPORT MODULES

import pandas as pd
import numpy as np
from iminuit import Minuit
import matplotlib.pyplot as plt


###=========================================================================
###=========================================================================

##=== !!! CHANGE ONLY HERE !!!

telescope='SPHERE-2'        #Type of telescope
H=500                       #Altitude of the telescope, m
type_analyse='only_sig'   #Type of analyse --> only_sig / electronic / experiment
En=30                       #Energy
nuclei='pro'                #Type of nuclei

## !!! END CHANGE ONLY HERE !!!

###=========================================================================
###=========================================================================

### OPEN FILES

if telescope=='SPHERE-2':
    
    x_y_mos=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_background/Data/Data_init/mosaic_pmt_coords_df.csv')            #Coord. PMT on mosaic, mm
    circle_of_pmt_around_pmt=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/circle_of_pmt_around_pmt.csv')  #Neighbors of PMT
    coeff_amp=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/coeff_amp.csv',header=None,sep='\s+')          #Coeff. amplification

    #Path of the file
    
    if type_analyse=='only_sig':
        path=f'/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/2_data_signal/only_sig_{En}PeV_{nuclei}_{H}m/mosaic_hits_m01_pro_{En}PeV_10-20_001_c001' #Signal only --> #01 for diploma
    
    elif type_analyse=='electronic':
        path=f'/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/1_data_electronics/electronics_{En}PeV_{nuclei}_{H}m/mosaic_hits_m01_pro_{En}PeV_10-20_001_c001'
    
    elif type_analyse=='experiment':
        path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/3_experiment/01889.txt' #good 12524, 10677, too low 11846, / intense 13930/ for diploma low 11398

    #Real angles
    
    real_theta=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_theta',header=None,sep='\s+')[0])
    real_phi=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_phi',header=None,sep='\s+')[0]) 

if telescope=='SPHERE-3':
    
    pixel_data = pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/SPHERE3_pixel_data_A.dat',header=None, names=['x','y','z','nx','ny','nz','theta','phi'], sep='\s+')  #Coord. PMT on mosaic, mm
    pixel_data['segment'] = pixel_data.index // 7 #== add num segment
    pixel_data['pixel'] =  pixel_data.index % 7 #== add num pixel
    circle_of_pmt_around_pmt=pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/circle_of_pmt_around_pmt_sph3.csv')  #Neighbors of PMT

    #Path of the file
    
    if type_analyse=='only_sig':
        path=f'/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/sph3/{H}m_{En}PeV_{nuclei}/moshits_Q2_atm01_0014_{nuclei}PeV_15_001_c007' #07 for diploma
    
    elif type_analyse=='electronic':
        path=f'/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/1_data_electronics/electronics_{En}PeV_{nuclei}_{H}m/mosaic_hits_m01_pro_{En}PeV_10-20_001_c001' #Electronic signal #044 #072

    #Real angles
    
    real_theta=15/180*np.pi
    real_phi_sph3_Fe=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/angels_phi_Fe',header=None,sep='\s+')[0])
    real_phi_sph3_N=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/angels_phi_N',header=None,sep='\s+')[0])
    real_phi_sph3_P=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/angels_phi_P',header=None,sep='\s+')[0])


###=========================================================================
###=========================================================================

### CONSTANT VALUES

num_pmt=np.arange(0,total_number_of_pmt) #number of PMT for 0 to the last
c_ns=0.3                                 #speed of light m/ns
bins=12.5                                #ns
cells=np.arange(0,response*bins,bins)    #length of the electronic response, bin=12.5 ns *1200
pied=400
response=1020

#Start values of multistart
theta_initt=[0.05,0,1,0.15, 0.2,0.25, 0.3,0.35]*8 #in rad
phi_initt=np.arange(0.5,6.5,0.5).tolist()*8       #in rad

if telescope=='SPHERE-2':
    
    x_mos=np.array(x_y_mos['x'])    #x on mosaic, mm
    y_mos=np.array(x_y_mos['y'])    #y on mosaic, mm
    
    num_t_event=4                   #column of time in initial file
    total_number_of_pmt=109         #total number of PMT
    
    d0_center=192                   #last circle on mosaic
    d0_before_c=154                 #before last circle on mosaic
    d_six=50                        #distance PMT and 1st circle around
    d_double_six=100                #distance PMT and 2nd circles around
    
    lim_search=100                  #diapason for search the impulse by DFS, ns
    
    if type_analyse=='experiment':
        pied=395
        response=1015
    
if telescope=='SPHERE-3':
    
    x_mos=np.array(pixel_data.groupby('segment').mean()['x']) #x on mosaic
    y_mos=np.array(pixel_data.groupby('segment').mean()['y']) #y on mosaic
    
    num_t_event=5                                             #column of time in initial file
    total_number_of_pmt=379                                   #total number of PMT
    
    d0_center=292                                             #last circle on mosaic
    d0_before_c=262                                           #before last circle on mosaic
    d_six=32                                                  #distance PMT and 1st circle around
    d_double_six=70                                           #distance PMT and 2nd circles around
    
    lim_search=50                                             #diapason for search the impulse by DFS, ns

    
if telescope=='SPHERE-3': 
    
    x_mos=np.array(pixel_data.groupby('segment').mean()['x'])   #x on mosaic
    y_mos=np.array(pixel_data.groupby('segment').mean()['y'])   #y on mosaic
    
    num_t_event=5                   #column of time in initial file
    total_number_of_pmt=379         #total number of PMT
    
    d0_center=292                   #last circle on mosaic
    d0_before_c=262                 #before last circle on mosaic
    d_six=32                        #distance PMT and 1st circle around
    d_double_six=70                 #distance PMT and 2nd circles around
    
    lim_search=50                   #Diapason for search the impulse by DFS, ns
    
    pied=400                        #piedestal
    response=1020                   #lenght of the response

    
    if type_analyse=='electronic':  
        pied=100                    #piedestal
        response=500                #lenght of the respons



###=========================================================================
###=========================================================================

### ALGORITHM

class AngleRetrieval:
    
    def __init__(self,path):
        self.path=path
        
       
    if telescope=='SPHERE-2':
        
        if type_analyse=='only_sig':
        
            def read_mosaic_hits_file(self): #Open files
                step = 12.5
                shift = 500 #first bin of the signal
                q = np.zeros((total_number_of_pmt, 1020))  # empty array109 PMT * 1020 cells
                cells_here=np.arange(0,1021*step,step)
                
                a = pd.read_csv(self.path, header=None, skiprows=[0],sep='\s+')
                pd_event=np.array(pd.DataFrame({'pmt_ev':list(a[0]),'t_ev':list(a[num_t_event]-min(a[num_t_event]))})) #creation of DataFrame to facilite calculations
                
                for i in np.arange(num_pmt.shape[0]):
                  impulse = np.array(pd_event)[:,1][np.where(np.array(pd_event)[:,0]==i)]
                  n_photons_ev,t_photon_evv=np.histogram(impulse,bins=cells_here) #number and time of photons in each cell of the i pmt
                  n_phot_shift=[0]*shift+n_photons_ev[:-shift].tolist() #same (n) with shift
                  q[i]=n_phot_shift
                
                return pd.DataFrame(q.T)
        
        elif type_analyse=='electronic':
        
           def read_mosaic_hits_file(self):
               event = pd.read_csv(self.path, header=None, sep='\s+', skiprows=[0],on_bad_lines='skip').to_numpy()
               q = np.zeros((1020,total_number_of_pmt))  # empty array109 PMT * 1020 cells
               for pmt in range(event.shape[1]):
                  q[0::2, pmt] = event[0::2, pmt] - np.mean(event[0:pied:2, pmt])
                  q[1::2, pmt] = event[1::2, pmt] - np.mean(event[1:pied:2, pmt])
               for i in range(event.shape[1]):
                  q.T[i]=q.T[i]*coeff_amp[1].iloc[i]
               return pd.DataFrame(q)
           
        elif type_analyse=='experiment':
        
            def read_mosaic_hits_file(self):
                event = pd.read_csv(self.path,skiprows=40,skipfooter=2042,engine='python',sep='\s+',header=None)
                event=event.drop(columns=[0,110,111,112])
                event.columns-=1 #name of PMT from 0 to 108
                event=event.to_numpy()
                q = np.zeros((response,total_number_of_pmt))  # empty array109 PMT * 1020 cells
                for pmt in range(event.shape[1]):
                   q[0::2, pmt] = event[0::2, pmt] - np.mean(event[0:pied:2, pmt])
                   q[1::2, pmt] = event[1::2, pmt] - np.mean(event[1:pied:2, pmt])
                for i in range(event.shape[1]):
                   q.T[i]=q.T[i]*coeff_amp[1].iloc[i]
                return pd.DataFrame(q)
           
    elif telescope=='SPHERE-3':
        
        if type_analyse=='only_sig':
    
          def read_mosaic_hits_file(self):
              event = pd.read_csv(self.path, header=None, sep='\s+',on_bad_lines='skip').T
              event['pixel'],event['segment']=event.index%8,event.index//8 
              event=event[event['pixel']!=7].groupby('segment').mean().drop(['pixel'], axis=1).T  #delete each 8th pixel and average each segment + #drop column pixel 
              for pmt in range(event.shape[1]):
                  event[pmt] -= event[0:100][pmt].mean() #substraction of the piedestal
              return event
        
    
    def sum_of_impulse(self): #Summed impulse in each channel
        # event=final_results.read_mosaic_hits_file()
        event=self.read_mosaic_hits_file()
        return event.sum(axis=1) #summarize event
    
    def intensity(self): #Intensity of the event
                        
        event=self.read_mosaic_hits_file()
        summed_impulse=self.sum_of_impulse() #summed event
    
        max_signal=summed_impulse[1:-1].max() #amplitude max
        noise_pied=summed_impulse[1:pied][summed_impulse[1:pied]>0] #piedestal
        
        intensity_signal=(max_signal)/(noise_pied.mean()+3*noise_pied.std()) #ration max / noise
        
        if type_analyse!='only_sig' and intensity_signal<2:
            print('too low signal')
            return
        elif type_analyse=='only_sig' and max_signal<50:
            print('too low signal')
            return
        else:
            return max_signal,intensity_signal
        
    def diapason(self): #diapason of the event
        
        event=self.read_mosaic_hits_file()
        summed_impulse=self.sum_of_impulse() #summed event

        imax = summed_impulse[:-1].idxmax() #index of the max impulse --> delete the last cell -> > max of impulse
        start, stop = imax, imax
        
        maxpart = 1
        background=np.mean(summed_impulse[1:pied])+3*np.std(summed_impulse[1:pied])
        bgd = maxpart * background
        
        while summed_impulse[start] > bgd and start!=summed_impulse.index[0]:
            start -= 1
        while summed_impulse[stop] > bgd and stop!=summed_impulse.index[-1]:
            stop += 1
        
        margin_left,margin_right= 30, 30 
        
        stop  += margin_right
        start -= margin_left
        
        event_diapason=event[start:stop] #diapason of the impulse - n photons
        cells_diapason=cells[start:stop] #diapason of the impulse - time
        
        if cells_diapason[0]+30*12.5<pied*12.5:
            print('event in piedestal')
            return
        else:
            event_diapason.index=np.arange(len(event_diapason)) #reindex new finded diapason
            # self.save_results(event_diapason)
            return event_diapason, cells_diapason, background
        
    def translation_snow_mos(self): #Translation x,y,t from mosaic to snow
        
        if telescope=='SPHERE-2':
            if type_analyse!='experiment':
                High=H                                                      
                b=0.9051146902370498*(High/500)                             #coeff for translate x mos --> x snow
                x_snow=-b*x_mos                                             #x on snow
                y_snow=-b*y_mos                                             #y on snow

            else:
                High=float(list(pd.read_csv(self.path).iloc[4])[0])-456 #430 = Baikal altitude
                b=0.9051146902370498*(High/500) 
                x_snow=-b*x_mos #x on snow
                y_snow=-b*y_mos #y on snow
                t_path=((np.sqrt(H**2+(np.sqrt(np.array(x_snow)**2+np.array(y_snow)**2))**2))/c_ns) #path in ns

        if telescope=='SPHERE-3':
            High=H
            b=0.5657087664485793*(High/500)                                #coeff for translate x mos --> x snow
            x_snow=-b*x_mos                                             #x on snow
            y_snow=-b*y_mos                                             #y on snow
          
        t_path=((np.sqrt(High**2+(np.sqrt(np.array(x_snow)**2+np.array(y_snow)**2))**2))/c_ns) #path from mosaic to snow in ns

        return x_snow, y_snow, t_path, b
        
    def axis(self): #Axis of the event
        b=self.translation_snow_mos()[3]
        event_diapason=self.diapason()[0]
        
        x0_max,y0_max=x_mos[event_diapason.sum().idxmax()],y_mos[event_diapason.sum().idxmax()] #on mosaic, mm
        d0_from_center=np.sqrt(x0_max**2+y0_max**2)
        
        if d0_from_center>d0_center:
            print('axis too far')
            return
        else:
            d_pmt_to_d_max=np.sqrt( (np.array(x_mos)-x0_max)**2  + (np.array(y_mos)-y0_max)**2  ) #distance of each PMT from the max PMT, m
            
            if d0_from_center>d0_before_c and d0_from_center<d0_center: #if the max PMT is on the circle before the last one
                max_d=d_six #the minimum distance of the neighors around
            
            elif d0_from_center<d0_before_c: #if the max PMT is inside
                max_d=d_double_six #the minimum distance of the neighors around
        
            number_photons=list(event_diapason[event_diapason>0].sum()) #max of each PMT
          
            x_circle_around_max=x_mos[np.where(d_pmt_to_d_max<max_d)[0]] #x around the max PMT
            y_circle_around_max=y_mos[np.where(d_pmt_to_d_max<max_d )[0]] #y around the max PMT
            number_photons_around_max=np.array(number_photons)[np.where(d_pmt_to_d_max<max_d )[0]] #max PMT around the max PMT
            
            x0=np.sum(x_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max) #x by gravity center
            y0=np.sum(y_circle_around_max*number_photons_around_max)/np.sum(number_photons_around_max) #y by gravity center
            return x0*-b, y0*-b
    
    def lenght(self): #Lenght of the signal
        intensity_signal=self.intensity()[0]
        
        if type_analyse=='only_sig': #just the signal without background
            event=self.read_mosaic_hits_file()
            event_cells_sup_0=pd.DataFrame({'event':np.array(event.max(axis=1)),'cells':cells})
            event_cells_sup_0=event_cells_sup_0[event_cells_sup_0['event']>0]
            #Delete 90%
            cells_sup0=event_cells_sup_0['cells']
            lim_inf,lim_sup = np.percentile(np.array(cells_sup0), 5), np.percentile(np.array(cells_sup0), 95) #percent of the impulse saved
            idx_95=cells[cells_sup0[cells_sup0>lim_sup].index[0]-1]
            idx_5=cells[cells_sup0[cells_sup0<lim_inf].index[-1]+1]
            lenght_impulse=idx_95-idx_5 #lenght
            
        elif type_analyse!='only_sig': #signal with background
            cells_diapason=self.diapason()[1]
            lenght_impulse=cells_diapason[-1]-cells_diapason[0]-60*12.5  #lenght
        
        if lenght_impulse<100 or lenght_impulse>2500:
            print('lenght too short or long')
            return
        elif intensity_signal>10000 or lenght_impulse>2000:
            print('calibration event')
            return
        else:
            return lenght_impulse
    
    def amplification(self): #Amplification of the signal
        find_diapason=self.diapason()
        event_diapason,cells_diapason=find_diapason[0],find_diapason[1]
        
        amplif_signal=event_diapason.rolling(window=4).sum().dropna()#amplification of the signal by sliding window
        amplif_signal.index=np.arange(0,len(amplif_signal)) #reindex
        amplif_cells=cells_diapason[:-3] #t, ns
        return amplif_signal, amplif_cells
    

    def noise_pmt(self): #Noise in each PMT
        event=self.read_mosaic_hits_file()
        
        noise_window=event[1:pied].rolling(window=4).sum().dropna() #amplify the piedestal
        if type_analyse=='only_sig':
            mean_std=[0]*total_number_of_pmt #threshold <bg>+3 std
        else:
            mean_std=noise_window[noise_window>0].mean()+3*noise_window[noise_window>0].std() #threshold <bg>+3std
        return mean_std
    
    def front_up_noise(self): #Selection CRITERIA PMT: noise / Front x,y,t,N above the noise
        event_diapason=self.diapason()[0]
        N_t_amplification=self.amplification() #amplification of the new diapason
        all_window,cells_window=N_t_amplification[0],N_t_amplification[1] 
        noise_pmt=self.noise_pmt()
        translation_snow_mos_results=self.translation_snow_mos()
        x_snow,y_snow,t_path=translation_snow_mos_results[0],translation_snow_mos_results[1],translation_snow_mos_results[2]
        
        #Create DataFrame x,y,t,N
        x_y_t_N=pd.DataFrame({'PMT':num_pmt,'x':x_snow,'y':y_snow,'N_max':list(all_window.max()),'t_max':cells_window[all_window.idxmax()],'index':all_window.idxmax(),'noise_thrs':noise_pmt}).sort_values(by='N_max',ascending=False)
        x_y_t_N=x_y_t_N[x_y_t_N['N_max']>x_y_t_N['noise_thrs']] #save the PMT up to max noise
    
      
        x_y_t_N['t_max']=x_y_t_N['t_max']-t_path[np.array(x_y_t_N['PMT']).astype(int)]
        
        return x_y_t_N
        
    def DFS(self): #Selection CRITERIA PMT: diapason / DFS method: PMT with signal
    
        x_y_t_N=self.front_up_noise() #x,y,t,N
        
        x_y_t_N_good=x_y_t_N.iloc[0].to_frame().T #x,y,t,N: axis --> good PMT
        x_y_t_N_maybe=x_y_t_N.iloc[1:] #Others --> don't know
         
        PMT_around=circle_of_pmt_around_pmt.loc[x_y_t_N_good.index[0]].dropna()[1:]  #PMT around the max, delete the central PMT
        
        #Move them from maybe to good table
        x_y_t_N_good=pd.concat([x_y_t_N_good,x_y_t_N_maybe.loc[PMT_around[PMT_around.isin(x_y_t_N_maybe['PMT'])]]]) 
        x_y_t_N_maybe=x_y_t_N_maybe.drop(PMT_around[PMT_around.isin(x_y_t_N_maybe['PMT'])]) #Delete 6 neighbors around the max PMT
    
        for k in range(0,20):
            good_pmt=[]
            bad_pmt=[]
            for i in range(0,len(x_y_t_N_maybe)):
                looked=x_y_t_N_maybe.iloc[i] #number of the PMT in the table 2
                PMT_around=circle_of_pmt_around_pmt.loc[looked.name].dropna()[1:]  #PMT around the PMT, delete the central PMT
                
                PMT_around_in_table_1=PMT_around[PMT_around.isin(x_y_t_N_good['PMT'])]
                PMT_around_in_table_2=PMT_around[PMT_around.isin(x_y_t_N_maybe['PMT'])]
                
                if len(PMT_around_in_table_1)==0: #No neighbours of the PMT of table 2 in the table 1
                    continue
                
                else:
                    mean_time=x_y_t_N_good.loc[PMT_around_in_table_1]['t_max'].mean() #mean of the sure PMT around the examinated PMT
                    
                    if looked['t_max'] <= mean_time+lim_search and  looked['t_max'] >= mean_time-lim_search: 
                        good_pmt.append(looked['PMT'])
                    else:
                        bad_pmt.append(looked['PMT'])
                        continue
            
            x_y_t_N_good=pd.concat([x_y_t_N_good,x_y_t_N_maybe.loc[good_pmt]])
            x_y_t_N_maybe=x_y_t_N_maybe.drop(good_pmt+bad_pmt) #Delete sure PMT
         
        x_y_t_N_good['t_max_on_min']=np.array(x_y_t_N_good['t_max'])-min(np.array(x_y_t_N_good['t_max'])) #Readjust ti on t min
        x_y_t_N_good=x_y_t_N_good.sort_values(by='N_max')
        
        return x_y_t_N_good
         
    def neighbors(self): #Selection CRITERIA PMT: minimum 3 neighbors 
    
        x_y_t_N_good=self.DFS()
        
        bad_PMT=[]
        for i in range(0,len(x_y_t_N_good)):
            if len(x_y_t_N_good['PMT'][x_y_t_N_good['PMT'].isin(circle_of_pmt_around_pmt.loc[x_y_t_N_good.iloc[i][0]].dropna())].drop(x_y_t_N_good.iloc[i][0]))<2:
                bad_PMT.append(int(x_y_t_N_good.iloc[i][0]))
        
        x_y_t_N_neighbors=x_y_t_N_good.drop(bad_PMT)
        
        if len(x_y_t_N_neighbors)<20:
            print('not enough PMT')
            return
        else:
            # self.save_results(x_y_t_N_neighbors)
            return x_y_t_N_neighbors
    
    
    def angles(self):
        
        x_y_t_N_neighbors=self.neighbors()
        
        x_front,y_front,t_fnc,N_photons=np.array(x_y_t_N_neighbors['x']),np.array(x_y_t_N_neighbors['y']),np.array(x_y_t_N_neighbors['t_max_on_min']),np.array(x_y_t_N_neighbors['N_max']) #extract values of the x ,y, t, N of the front on snow
    
        
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
        
        a0_multi_start=[]
        a1_multi_start=[]
        a2_multi_start=[]
    
        
        min_s=[]
        for th in range(0,len(theta_initt)):
                param_fix=Minuit(fnc, theta=theta_initt[th],phi=phi_initt[th],a0=0,a1=0,a2=0)
                param_fix.limits['theta']=(0,60/180*np.pi)
                param_fix.limits['phi']=(0,2*np.pi)
    
        
                param_fix.migrad()
                theta_multi_start.append(param_fix.values[0]) 
                phi_multi_start.append(param_fix.values[1])
                
                a0_multi_start.append(param_fix.values[2]) 
                a1_multi_start.append(param_fix.values[3])
                a2_multi_start.append(param_fix.values[4])
    
                
                min_s.append(param_fix.fval)
            
        theta=np.array(theta_multi_start)[np.where(np.array(min_s)==min(min_s))][0]
        phi=np.array(phi_multi_start)[np.where(np.array(min_s)==min(min_s))][0]
        
        a0=np.array(a0_multi_start)[np.where(np.array(min_s)==min(min_s))][0]
        a1=np.array(a1_multi_start)[np.where(np.array(min_s)==min(min_s))][0]
        a2=np.array(a2_multi_start)[np.where(np.array(min_s)==min(min_s))][0]
        
        
        return theta, phi, a0, a1, a2, x_front, y_front, t_fnc
    
    if type_analyse!='experiment':
    
        def delta(self,theta_real,phi_real): #Delta (error)
        
            angles_res=self.angles()
            
            theta_fit=angles_res[0]
            phi_fit=angles_res[1]
            cx = np.sin(theta_real) * np.cos(phi_real) #theta - истинные углы
            cy = np.sin(theta_real) * np.sin(phi_real)
            cz = np.cos(theta_real)
            cx_fit = np.sin(theta_fit) * np.cos(phi_fit) #theta' -полученные углы моделированием
            cy_fit = np.sin(theta_fit) * np.sin(phi_fit)
            cz_fit = np.cos(theta_fit)
            delta=np.arccos(cx *cx_fit + cy * cy_fit + cz * cz_fit)*180/np.pi
            
            return delta

    def save_results(self):
        angles_and_param=self.angles()
        axis_values=self.axis()
        # lenght_value=self.lenght()

        table_final=np.array([angles_and_param[0],angles_and_param[1],angles_and_param[2],angles_and_param[3],angles_and_param[4],axis_values[0],axis_values[1]]).T
        # np.savetxt('/Users/clemence/Documents/test',table_final)


final_results=AngleRetrieval(path) 

final_results.angles()
# print(final_results.angles()[0])

class Figure:
    
    def __init__(self):
        pass
    
    def fig_sum_impulse(self):
        event=final_results.read_mosaic_hits_file()
        noise_mean_std=final_results.diapason()[2]
        
        plt.bar(cells,event.sum(axis=1),width=50)
        if type_analyse!='only_sig':
            plt.axhline(y=noise_mean_std,color='darkorange',linewidth=3)
        plt.grid()
        plt.ylabel('∑ Ni фот',fontsize=14)
        plt.xlabel('t, нс',fontsize=14)
        plt.xlim(0,None)
        # plt.ylim(-100,None)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.show()
        
    def front(self):
        
        front_3d=final_results.angles()
        
        fig = plt.figure(figsize=(8,10))
        ax = plt.axes(projection="3d")

        x_1=front_3d[5]
        y_1=front_3d[6]
        t_1=front_3d[7]

        ax.scatter(x_1,y_1,t_1,s=10,label='Experimental front',color='blue')
        # # ax.scatter(x_1,y_1,tau_fit,s=10,label='Аппроксимация',color='red')
        ax.view_init(10,140)  #For diploma: sph2 (10.80), sph3 (10,115)
        ax.set_xlabel("x, м",fontsize=14)
        ax.set_ylabel("y, м",fontsize=14)
        ax.set_zlabel("t, нс",fontsize=14)
        ax.xaxis.labelpad = 10
        ax.yaxis.labelpad = 10
        ax.zaxis.labelpad = 5
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

        plt.xticks(np.arange(-200,300,100))
        plt.yticks(np.arange(-200,300,100))
        ax.set_box_aspect(aspect=None, zoom=0.8)
        # # plt.savefig('/Users/clemence/Documents/Scientific_research/Diploma/Картинки/3d_front_sph2.pdf',bbox_inches='tight')
        plt.show() #3d_front_sph2

pictures=Figure()
# pictures.fig_sum_impulse()
# pictures.front()


