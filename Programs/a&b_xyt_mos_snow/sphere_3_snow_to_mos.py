import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import numpy as np

###=========================================================================

##=== !!! CHANGE ONLY HERE !!!
#High

H=500

###=========================================================================
###=========================================================================

#Constant values

if H==500:
    vision_field=342/2 #m
if H==1000:
    vision_field=342
    
###=========================================================================
###=========================================================================

##===OPEN FILES

pixel_data = pd.read_csv('/Users/clemence/Documents/Scientific_research/Sphere_3/Data/SPHERE3_pixel_data_A.dat',header=None, names=['x','y','z','nx','ny','nz','theta','phi'], sep='\s+') #number of PMT

#== Open directory

directory_snow = '/Users/clemence/Documents/Scientific_research/Sphere_3/Data/phels2trace'
name_files_snow = glob.glob(os.path.join(directory_snow,'*'))

directory_mosaic = '/Users/clemence/Documents/Scientific_research/Sphere_3/Data/moshits'
name_files_mosaic = glob.glob(os.path.join(directory_mosaic,'*'))


###=========================================================================
###=========================================================================

initial_snow=[]
d_mosaic=[]

saved_snow=[]
saved_mos=[]

nbr_phot=[]

for files in range(0,len(name_files_snow)):
    
    #______________

    #Distance on SNOW
    
    xsnow_point=float(list(pd.read_csv(name_files_snow[files],sep='\s+'))[2])
    
    #SAVE DATA
    
    #Real d
    if H==1000:
        initial_snow.append(xsnow_point) # m 
    if H==500:
        initial_snow.append(xsnow_point/2) 
    
    #Number of the event
    saved_snow.append(int(name_files_snow[files][-3:]))

    #______________
    
    ##===OBTAIN distance on MOSAIC
        
    #DataFrame of segment, d
    coord_mos=pd.read_csv(name_files_mosaic[files],header=None,skiprows=[0], sep='\s+') #number of PMT
    
    #SAVE DATA
    
    d_mosaic.append(coord_mos[2].mean())

    #Number of the event
    saved_mos.append(int(name_files_mosaic[files][-3:]))
    
    #______________

    #SAVE DATA
    
    nbr_phot.append(len(coord_mos))

###=========================================================================
###=========================================================================

files_and_d_snow=pd.DataFrame({'files':saved_snow,'d_snow':initial_snow}).sort_values(by='files')
files_and_d_mosaic=pd.DataFrame({'files':saved_mos,'d_mosaic':d_mosaic,'N':nbr_phot}).sort_values(by='files')

total_distance=pd.DataFrame({'d_snow':list(files_and_d_snow['d_snow']),'d_mosaic':list(files_and_d_mosaic['d_mosaic']),'N':list(files_and_d_mosaic['N'])})

#To plot with the size 25 m
range_index=np.arange(0,100,5)
nbr_pmt_bins=[total_distance['N'][range_index[i]:range_index[i+1]].sum() for i in range(0,len(range_index)-1)]
nbr_pmt_bins=nbr_pmt_bins+[total_distance['N'][95:].sum()]

###=========================================================================
###=========================================================================

##=== APPROXIMATION

#r=aR
XX = np.vstack((np.array(total_distance[total_distance['d_snow']<vision_field]['d_snow']) *1, np.ones_like(np.array(total_distance[total_distance['d_snow']<vision_field]['d_snow'])))).T
p_no_offset = np.linalg.lstsq(XX[:, :-1], np.array(total_distance[total_distance['d_snow']<vision_field]['d_mosaic']))[0]
a=p_no_offset[0]

x_approx=np.arange(0,vision_field-1)
y_approx=a*x_approx

#R=ar

total_distance_new=total_distance[total_distance['d_snow']<vision_field]

XX = np.vstack((np.array(total_distance_new['d_mosaic']) *1, np.ones_like(total_distance_new['d_mosaic']))).T
p_no_offset = np.linalg.lstsq(XX[:, :-1], np.array(total_distance_new['d_snow']))[0]
a2=p_no_offset[0]

x_approx2=np.arange(0,31)
y_approx2=a2*x_approx2

###=========================================================================
###=========================================================================

# ##===GRAPHICS
params = {"xtick.direction": "in", "ytick.direction": "in"}
plt.rcParams.update(params)

#=======================

plt.scatter(list(total_distance['d_snow']),list(total_distance['d_mosaic']*100))
plt.plot(x_approx,y_approx*100,color='red',label='a={:.3f}'.format(a*100))
plt.axvline(x=vision_field,ymin=0, ymax=1,color='black')
plt.ylabel('r on mosaic, cm')
plt.xlabel('R on snow, m')
plt.grid(True)
plt.xlim(0,None)
plt.ylim(0,None)
plt.legend(title='{} m \n r=aR'.format(H))
plt.show()

#=======================

plt.scatter(list(total_distance_new['d_mosaic']*100),list(total_distance_new['d_snow']))
plt.plot(x_approx2,y_approx2/100,color='red',label='a={:.3f}'.format(a2/100))
# plt.axhline(y=vision_field,xmin=0, xmax=1,color='black')
plt.xlabel('r на мозаике, см')
plt.ylabel('R на снегу, м')
plt.grid(True)
plt.xlim(0,32)
plt.ylim(0,260)
plt.legend(title='{} м \n R=ar'.format(H),loc='upper left')
plt.show()

#=======================

#Number of photons in each R

plt.bar(np.arange(5,495,25)+8,nbr_pmt_bins,width=25,edgecolor='black',color='cadetblue')

plt.xlabel('R on snow, m')
plt.ylabel('N photons on mosaic')
plt.yscale('log')
plt.legend(title='{} m'.format(H))
plt.xlim(0,None)
plt.show()



#All SPHERE-2 and SPHERE-3

x_sph2=np.array([ 0.07679698,  1.07325244,  2.0941152 ,  3.17782316,  4.15476948,
        5.23068252,  6.27743836,  7.24842078,  8.26769575,  9.28504577,
       10.23617369, 11.17767333, 12.40028595, 13.25841488, 14.24109901,
       15.2293795 , 16.14962709, 17.05525302, 17.89341556, 18.84235461,
       19.71725262, 20.59167123, 21.38201666, 22.26667739, 23.13469526])
y_sph2=np.arange(0,250,10)

a,b,c=np.polyfit(x_sph2,y_sph2,2)
x_fit_sph2=np.arange(0,max(x_sph2),1)
y_fit_sph2=a*x_fit_sph2**2+b*x_fit_sph2


#== SPHERE 3

plt.scatter(x_sph2,y_sph2,color='green',label='СФРЕА-2')
plt.plot(x_fit_sph2,y_fit_sph2,color='red')


#== SPHERE 2

plt.scatter(list(total_distance_new['d_mosaic']*100),list(total_distance_new['d_snow']),color='blue',label='СФРЕА-3')
plt.plot(x_approx2,y_approx2/100,color='red')


# R=ar^2+br \na={:.3f}, b={:.1f}'.format(a,b)
# R=ar \na={:.2f}'.format(a2/100)

plt.xlabel('r на мозаике, см')
plt.ylabel('R на снегу, м')
plt.grid(True)
plt.xlim(0,32)
plt.ylim(0,260)
plt.legend(title='{} м'.format(H),loc='upper left')
plt.annotate('{:.2f}*r'.format(a2/100), xy=(20, 100), xytext=(24, 45),arrowprops={'facecolor':'blue', 'shrink':5} )
plt.annotate('{:.3f}*r2 + {:.1f}*r'.format(a,b), xy=(13, 150), xytext=(1, 130),arrowprops={'facecolor':'green', 'shrink':5} )
# plt.savefig('foo.pdf')
plt.show()




