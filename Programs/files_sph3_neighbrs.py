import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

pixel_data = pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/SPHERE3_pixel_data_A.dat',header=None, names=['x','y','z','nx','ny','nz','theta','phi'], sep='\s+') #number of PMT
pixel_data['segment'] = pixel_data.index // 7 #== add num segment
pixel_data['pixel'] =  pixel_data.index % 7 #== add num pixel
total_number_of_pmt=379


H=500

if H==500:
    min_ngb=20
if H==1000:
    min_ngb=40
b=(1.1314175328971585)/(1000/H) #If H=500 - a --> a/2

x_snow=-b*np.array(pixel_data.groupby('segment').mean()['x'])
y_snow=-b*np.array(pixel_data.groupby('segment').mean()['y'])

PMT=[]
above=[]
below=[]
right=[]
left=[]
diag_left=[]
diag_right=[]
for i in range(0,len(x_snow)):
    coord=pd.DataFrame({'x':x_snow,'y':y_snow,'num':np.arange(0,379)})
    coord['d']=np.sqrt( (x_snow - coord['x'].iloc[i])**2+ (y_snow - coord['y'].iloc[i]) **2)
    
    coord_nump=coord['num'][ (coord['d']<min_ngb) & (coord['d']>0) ]
    
    PMT.append(coord['num'].loc[i])
    
    if len(coord_nump)==3:
        above.append(coord_nump.iloc[0])
        below.append(coord_nump.iloc[1])
        right.append(coord_nump.iloc[2])
        left.append(np.nan)
        diag_left.append(np.nan)
        diag_right.append(np.nan)

    elif len(coord_nump)==4:
        above.append(coord_nump.iloc[0])
        below.append(coord_nump.iloc[1])
        right.append(coord_nump.iloc[2])
        left.append(coord_nump.iloc[3])
        diag_left.append(np.nan)
        diag_right.append(np.nan)
    
    elif len(coord_nump)==5:
        above.append(coord_nump.iloc[0])
        below.append(coord_nump.iloc[1])
        right.append(coord_nump.iloc[2])
        left.append(coord_nump.iloc[3])
        diag_left.append(coord_nump.iloc[4])
        diag_right.append(np.nan)


    else:
        above.append(coord_nump.iloc[0])
        below.append(coord_nump.iloc[1])
        right.append(coord_nump.iloc[2])
        left.append(coord_nump.iloc[3])
        diag_left.append(coord_nump.iloc[4])
        diag_right.append(coord_nump.iloc[5])


total=pd.DataFrame({'PMT':PMT,'above':above,'below':below,'right':right,'left':left,'diag_left':diag_left,'diag_right':diag_right})

# total.to_csv('/Users/clemence/Documents/Scientific_research/Angle_retrieval/Data/circle_of_pmt_around_pmt_sph3.csv', index=False) 


d=np.sqrt(x_snow**2+y_snow**2)

coord_d=pd.DataFrame({'x':x_snow,'y':y_snow,'d':d})

if H==500:
    x_max=coord_d['x'][coord_d['d']>60]
    y_max=coord_d['y'][coord_d['d']>40]
if H==1000:
    x_max=coord_d['x'][coord_d['d']>310]
    y_max=coord_d['y'][coord_d['d']>310]

plt.scatter(x_snow,y_snow)
plt.scatter(x_max,y_max)
plt.axis('equal')
plt.show()




