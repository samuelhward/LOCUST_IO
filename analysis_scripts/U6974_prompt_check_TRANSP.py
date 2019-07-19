#PLOT PROMPT LOSS MAPS FOR EACH RUN AND CODE

#prompt loss maps one row, each for each value of limiter size with ascot, transp and locust prompt losses on each

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm #get colourmaps
import copy
import context
from classes.output_classes.particle_list import Final_Particle_List as fpl
from classes.input_classes.wall import Wall

extensions=['U69','U70','U71','U72','U73','U74']
locust_files=[]
ascot_files=[]
radii=['1.05','1.10','1.20','1.30','1.40','1.50']
colors=['r','g','b','y','m','k']
file_numbers=[1] 
times=['0.100']
wall_file='LCFS_DIII-D.dat_'
colourmaps=[matplotlib.cm.get_cmap('Greens'),matplotlib.cm.get_cmap('Reds'),matplotlib.cm.get_cmap('Blues'),matplotlib.cm.get_cmap('Oranges'),matplotlib.cm.get_cmap('Purples'),matplotlib.cm.get_cmap('Greys')]


#separate plots
for time,file_number in zip(times,file_numbers):    
    fig,(ax1)=plt.subplots(1)

    for extension,radius,colourmap in zip(extensions,radii,colourmaps): #add locust_file and ascot_file to this eventually
    
        ax.set_xlim(1.,2.7)
        ax.set_ylim(-1.8,.6)
        
        #TRANSP runs
        loss_list=fpl(ID=extension+'   radius='+radius+'   time='+time+'s',data_format='TRANSP',filename='157418'+extension+'_fi_'+str(file_number)+'_gc.cdf')
        wall=Wall(radius,'LOCUST_2D',wall_file+radius)

        for quantity in ['time','R','Z','weight','status_flag','E','V_pitch']: #first values are junk
            loss_list[quantity]=np.delete(loss_list[quantity],0)
        
        loss_list['time']=(loss_list['time']-np.min(loss_list['time']))/(np.max(loss_list['time']-np.min(loss_list['time']))) #normalise along time dimension
        scatter=loss_list.plot(some_equilibrium=wall,limiters=True,style='scatter',axes=['R','Z'],fig=fig,ax=ax1,colmap=colourmap,status_flags=['all_losses'],colfield='time')
    

plt.colorbar(scatter)
ax1.set_title('TRANSP loss map t = 100ms')    
ax1.legend(radii)
legend = ax1.get_legend()
for i,colourmap in enumerate(colourmaps):
    legend.legendHandles[i].set_color(colourmap(.8))
    
plt.show()