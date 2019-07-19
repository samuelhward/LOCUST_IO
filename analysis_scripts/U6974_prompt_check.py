#PLOT PROMPT LOSS MAPS FOR EACH RUN AND CODE

#prompt loss maps one row, each for each value of limiter size with ascot, transp and locust prompt losses on each

from scipy.io import netcdf as ncdf
import numpy as np
import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib import cm #get colourmaps
from mpl_toolkits.axes_grid1 import make_axes_locatable #to sort size of colorbar
import copy
import context
from classes.output_classes.particle_list import Final_Particle_List as fpl
from classes.input_classes.equilibrium import Equilibrium
from classes.input_classes.wall import Wall

cmap_default=matplotlib.cm.get_cmap('jet') #set default colourmap

filename_eq='g157418.03000'
equi=Equilibrium(filename_eq,'GEQDSK',filename_eq)

cmap_TRANSP=matplotlib.cm.get_cmap('Blues_r')
cmap_LOCUST=matplotlib.cm.get_cmap('Greens_r')
cmap_ASCOT=matplotlib.cm.get_cmap('Reds_r')


extensions=['U69','U70','U71','U72','U73','U74']
locust_files=[]
ascot_files=[]
radii=['1.05','1.10','1.20','1.30','1.40','1.50']
file_numbers=[1] 
times=['0.100']
wall_file='input.wall_2d_'

ASCOT_files=['ascot_freia_9842185.h5','ascot_freia_9842187.h5','ascot_freia_9842189.h5','ascot_freia_9842190.h5','ascot_freia_9842193.h5','ascot_freia_9842195.h5']
ASCOT_run='ascot/run_1/' #this is with old kinetic profiles which are not extrapolated, ORBITMETHOD=1

ASCOT_files=['ascot_freia_9842186.h5','ascot_freia_9842188.h5','','ascot_freia_9842191.h5','ascot_freia_9842194.h5','ascot_freia_9842196.h5']
ASCOT_run='ascot/run_2/' #changed ORBITMETHOD to 4, added extrapolated kinetic profiles

ASCOT_files=['ascot_freia_9836979.h5','','ascot_freia_9837273.h5','ascot_freia_9837277.h5','ascot_freia_9837280.h5','ascot_freia_9837284.h5']
ASCOT_run='ascot/run_3/' #changed ORBITMETHOD back to 1, keep extrapolated kinetic profiles

ASCOT_files=['ascot_freia_9836980.h5','ascot_freia_9837054.h5','ascot_freia_9837274.h5','ascot_freia_9837278.h5','ascot_freia_9837281.h5','ascot_freia_9837285.h5']
ASCOT_run='ascot/run_4/' #ORBITMETHOD 4 and non-extrapolated kinetic profiles

#separate plots
for time,file_number in zip(times,file_numbers):    
    fig,((ax1,ax2,ax3),(ax4,ax5,ax6))=plt.subplots(2,3)
    axes=[ax1,ax2,ax3,ax4,ax5,ax6]
    for extension,ax,radius,ASCOT_file in zip(extensions,axes,radii,ASCOT_files): #add locust_file and ascot_file to this eventually
        
        wall=Wall(ID=radius,data_format='ASCOT_2D_input',filename=wall_file+radius)

        #TRANSP runs
        loss_list=fpl(ID='TRANSP - '+radius,data_format='TRANSP',filename='157418'+extension+'_fi_'+str(file_number)+'_gc.cdf')
        
        for quantity in ['time','R','Z','weight','status_flag','E','V_pitch']: #first values are junk
            loss_list[quantity]=np.delete(loss_list[quantity],0)

        loss_list['status_flag']+=0.5 #add some number (0<1) to give TRANSP PFC intercept status flag its own value so it has it's own colour when plotting
        loss_list['status_flags']['all_losses']+=0.5 #also need to change definition of status_flags accordingly
        
        loss_list['time']=(loss_list['time']-np.min(loss_list['time']))/(np.max(loss_list['time']-np.min(loss_list['time']))) #normalise along time dimension
        loss_list.plot(limiters=wall,LCFS=equi,style='scatter',axes=['R','Z'],fig=fig,ax=ax,colmap=cmap_TRANSP,status_flags=['all_losses'],colfield='time')
        
        #LOCUST runs
        
        #ASCOT runs 
        loss_list=fpl(ID='ASCOT - '+radius,data_format='ASCOT',filename=ASCOT_run+ASCOT_file)
        scatter=loss_list.plot(limiters=wall,LCFS=equi,style='scatter',axes=['R','Z'],fig=fig,ax=ax,colmap=cmap_ASCOT,status_flags=['wall_collision'],colfield='time')            

        ax.set_title('{}'.format(radius))
        ax.set_xlim(1.,2.7)
        ax.set_ylim(-1.8,.6)
        ax.set_aspect('equal')


    #legends and colourbars
    ax1.legend(['TRANSP','ASCOT'])#XXX add LOCUST to this
    legend=ax1.get_legend()
    for i,colourmap in enumerate([cmap_TRANSP,cmap_ASCOT]): #XXX add LOCUST to this
        legend.legendHandles[i].set_color(colourmap(.2))

    divider = make_axes_locatable(ax6) #make colourbar the right size
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(scatter,cax=cax)
    plt.show()