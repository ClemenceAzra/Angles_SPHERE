import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class Figure:
    
    def __init__(self, angles_results):
        self.telescope = angles_results.telescope    
        self.cells = self.telescope.cells
        self.angles_results = angles_results
        self.results = {}
 

    def get_results(self):

        ## если нет результатов, то надо их получить
        if not self.results:
            if not self.angles_results.results:
                print("No results!!")
                angles_results.angles()
            self.results = self.angles_results.results
        self.x_front = self.angles_results.results['x_front']
        self.y_front = self.angles_results.results['y_front']
        self.t_fnc   = self.angles_results.results['t_fnc']
         
    
    def fig_sum_impulse(self):
        event = self.angles_results.event
        noise_mean_std=self.angles_results.diapason()[2]
        
        plt.bar(self.cells,event.sum(axis=1),width=50)
        if self.telescope.type_analyse!='only_sig':
            plt.axhline(y=noise_mean_std,color='darkorange',linewidth=3)
        plt.grid()
        plt.ylabel('∑ Ni фот',fontsize=14)
        plt.xlabel('t, нс',fontsize=14)
        plt.xlim(0,11000)
        plt.ylim(-50,250)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.show()

        
    def front(self):      
        front_3d = self.angles_results.angles()
        
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
