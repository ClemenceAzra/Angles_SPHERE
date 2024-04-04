import pandas as pd
import os
import glob
import numpy as np

path='/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/1_data_electronics/electronics_{}PeV_{}_{}m'.format(En,nuclei,H) #Electronic signal #044 #072


name_files = glob.glob(os.path.join(path,'*'))

x0_real=[]
y0_real=[]

for files in range(0,len(name_files)):

    x0_real.append(float(pd.read_csv(name_files[files]).columns.tolist()[0].split()[2])) #real x of the axis
    y0_real.append(float(pd.read_csv(name_files[files]).columns.tolist()[0].split()[3])) #real y of the axis
