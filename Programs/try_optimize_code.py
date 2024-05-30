#%%time

###=========================================================================
###=========================================================================

### IMPORT MODULES

import pandas as pd
import numpy as np
from iminuit import Minuit
import matplotlib.pyplot as plt

from angleretrieval import *
from telescope import *
import figure 

        
##=========================================================================
###=========================================================================

##=== !!! CHANGE ONLY HERE !!!

telescope='SPHERE-2'        #Type of telescope
type_analyse='experiment'   #Type of analyse --> only_sig / electronic / experiment
H=500                       #Altitude of the telescope, m
En=30                       #Energy
nuclei='pro'                #Type of nuclei

## !!! END CHANGE ONLY HERE !!!

###=========================================================================
###=========================================================================

### path to FILES
maindir = "/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/"
maindir = "../data"

if telescope=='SPHERE-2':
    #Path of the file    
    if type_analyse=='only_sig':
        path=f'{maindir}/2_data_signal/only_sig_{En}PeV_{nuclei}_{H}m/mosaic_hits_m01_pro_{En}PeV_10-20_001_c001' #Signal only --> #01 for diploma   
    elif type_analyse=='electronic':
        path=f'{maindir}/1_data_electronics/electronics_{En}PeV_{nuclei}_{H}m/mosaic_hits_m01_pro_{En}PeV_10-20_001_c001'  
    elif type_analyse=='experiment':
        path = f'{maindir}/experiment/' #good 12524, 10677, too low 11846, / intense 13930/ for diploma low 11398
    #Real angles
    if type_analyse != 'experiment':
        real_theta=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_theta',header=None,sep='\s+')[0])
        real_phi=np.array(pd.read_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/Real_angles/angles_phi',header=None,sep='\s+')[0]) 

if telescope=='SPHERE-3':
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




#telescope = Telescope(name=telescope, type_analyse=type_analyse)        
telescope = Telescope()  

for filename in os.listdir(path):
    if not filename.endswith(".txt"):
        continue
    print(path + filename)
    
    angles_results = AngleRetrieval(path + filename, telescope) 
    angles_results.save_front_to_file()
    #angles_results.angles()
    # print(final_results.angles()[0])
        
    pictures = figure.Figure(angles_results)
    pictures.fig_sum_impulse()
#pictures.front()
