#PLOT PROMPT LOSS MAPS FOR EACH RUN AND CODE

#prompt loss maps one row, each for each value of limiter size with ascot, transp and locust prompt losses on each

from scipy.io import netcdf as ncdf
import numpy as np
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib import cm #get colourmaps
cmap_default=matplotlib.cm.get_cmap('jet') #set default colourmap
import copy
import context
from classes.output_classes.particle_list import Final_Particle_List as fpl
from processing import plot_output as plou 

extensions=['U61','U62','U63','U64','U65','U66','U67','U68']
locust_files=[]
ascot_files=[]
radii=[1.01,1.05,1.1,1.15,1.2,1.3,1.4,1.5]
file_numbers=[1] 
times=['0.001']

#seperate plots
for time,file_number in zip(times,file_numbers):    
    fig,((ax1,ax2,ax3,ax4),(ax5,ax6,ax7,ax8))=plt.subplots(2,4)
    axes=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
    for extension,ax,radius in zip(extensions,axes,radii): #add locust_file and ascot_file to this eventually
        legend=[]
        ax.set_title('limiter radius = {}'.format(radius))
        #TRANSP runs
        loss_list=fpl(ID=extension+'   XBMBND='+str(radius)+'   time='+time+'s',data_format='TRANSP',filename='157418'+extension+'_fi_'+str(file_number)+'_gc.cdf')
        plou.plot_final_particle_list(loss_list,type='scatter',axes=['R','Z'],fig=fig,ax=ax,colmap='r',status_flags=['all_losses'],real_scale=True)
        legend+=['TRANSP']
        
        #LOCUST runs
        
        #ASCOT runs
    
    plt.legend(legend)
    plt.show()