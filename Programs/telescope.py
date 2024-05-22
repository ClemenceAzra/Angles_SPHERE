import pandas as pd
import numpy as np


class Telescope:
    def __init__(self, name="SPHERE-2", type_analyse='experiment'):
        self.name = name
        self.bins = 12.5     #ns
        self.type_analyse = type_analyse
        print(f"Telescope {self.name} with \"{type_analyse}\" type of analyse")
        
        ## geometry
        if self.name == 'SPHERE-2':
            #Coord. PMT on mosaic, mm
            x_y_mos=pd.read_csv('../data/sphere-2/mosaic_pmt_coords_df.csv')  
            self.x_mos=np.array(x_y_mos['x'])    #x on mosaic, mm
            self.y_mos=np.array(x_y_mos['y'])    #y on mosaic, mm
            self.circle_of_pmt_around_pmt=pd.read_csv('../data/sphere-2/circle_of_pmt_around_pmt.csv')  #Neighbors of PMT
            self.coeff_amp = pd.read_csv('../data/sphere-2/coeff_amp.csv', header=None, sep='\s+')          #Coeff. amplification
        elif self.name =='SPHERE-3':
            pixel_data = pd.read_csv('../data/SPHERE3_pixel_data_A.dat', header=None, names=['x','y','z','nx','ny','nz','theta','phi'], sep='\s+')  #Coord. PMT on mosaic, mm
            pixel_data['segment'] = pixel_data.index // 7 #== add num segment
            #pixel_data['pixel']   = pixel_data.index % 7 #== add num pixel
            self.x_mos=np.array(pixel_data.groupby('segment').mean()['x']) #x on mosaic
            self.y_mos=np.array(pixel_data.groupby('segment').mean()['y']) #y on mosaic
            self.circle_of_pmt_around_pmt=pd.read_csv('../data/circle_of_pmt_around_pmt_sph3.csv')  #Neighbors of PMT            
        else:
            print("Erooor! Wrong name of telescope!!! SPHERE-4 ???")
            exit()
            
        if self.name == 'SPHERE-2':
            self.num_t_event=4                   ##  column of time in initial file
            self.total_number_of_pmt=109         ##  total number of PMT
            self.d0_center=192                   #last circle on mosaic
            self.d0_before_c=154                 #before last circle on mosaic
            self.d_six=50                        #distance PMT and 1st circle around
            self.d_double_six=100                #distance PMT and 2nd circles around
            self.lim_search=100                  #diapason for search the impulse by DFS, ns
        else:
            self.num_t_event=5                   #column of time in initial file
            self.total_number_of_pmt=379         #total number of PMT
            self.d0_center=292                   #last circle on mosaic
            self.d0_before_c=262                 #before last circle on mosaic
            self.d_six=32                        #distance PMT and 1st circle around
            self.d_double_six=70                 #distance PMT and 2nd circles around
            self.lim_search=50                   #diapason for search the impulse by DFS, ns

            
        ##  
        self.pied = 400
        self.response = 1020
        if self.name == 'SPHERE-2':    
            if self.type_analyse=='experiment':
                self.pied=395
                self.response=1015
    
        if self.name == 'SPHERE-3':    
            if self.type_analyse=='electronic':  
                self.pied=100                    #piedestal
                self.response=500                #lenght of the respons

        self.num_pmt = np.arange(0, self.total_number_of_pmt) #number of PMT for 0 to the last        
        ##  length of the electronic response, bin=12.5 ns *1200
        self.cells = np.arange(0, self.response * self.bins, self.bins)    
