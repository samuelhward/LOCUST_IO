#perturbation_animation.py
"""
Samuel Ward
08/07/2019
----
quick plot of ITER FILD study MARS-F field_types
---
usage:
    see README.md for usage
 
notes:         
---
"""

import context
from classes.input_classes.perturbation import Perturbation as pert  
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import pathlib
import support

quantity_to_plot='B_field_R_imag' #magnitude of perturbation

def draw_perturbation(perturbation,ax1,ax2,fig):

    ax1.collections=[]
    ax2.cla()
    mesh=perturbation.plot(key=quantity_to_plot,fig=fig,ax=ax1,number_bins=100)#,vminmax=[0,0.01]) 
    cbar=fig.colorbar(mesh,cax=ax2,orientation='vertical')
    ax2.set_xlim([0,20])
    ax2.tick_params(axis="y",direction="in", pad=-150,length=0)

response=True
vacuum=True
rob_data=True
ascot_data=False

data_format='LOCUST'

filepaths=[]
field_types=[]

if ascot_data:

    if response:
        filepaths.append((support.dir_input_files / 'ITER_data' / 'B_field_Amp_90').glob('*'))
        field_types+=['response']
    if vacuum:
        filepaths.append((support.dir_input_files / 'ITER_data' / 'B_field_Amp_90_vac').glob('*'))
        field_types+=['vacuum']

    for filepath_group,field in zip(filepaths,field_types): #list of lists, where each list in the list is all files belonging to a type of perturbation field e.g. vacuum   

        fig,(ax1,ax2)=plt.subplots(1,2)
        ax1.set_xscale('linear')
        ax1.set_yscale('linear')

        perturbations=[] #pre-read and dump all perturbations before animating
        for filepath in filepath_group: 
            perturbation=pert(ID='ASCOT - {filename} - {field} - {quantity}'.format(filename=filepath.parts[-1],field=field,quantity=quantity_to_plot),data_format=data_format,filename=filepath)
            perturbations.append(perturbation)
            #perturb.dump_data(data_format='LOCUST',filename='LOCUST/coil_scans/BPLASMA_{harmonic}_{phase}_resist_resp.dat'.format(harmonic=harmonic,phase=phase))

        animation=FuncAnimation(fig,draw_perturbation,frames=perturbations,fargs=[ax1,ax2,fig],repeat=True,interval=1500) #cycle through phases and make animation
        plt.show()
        #animation.save('MARSF_{field}_field_types_{quantity}.gif'.format(field=field,quantity=quantity_to_plot),writer='pillow')


if rob_data:

    fig,(ax1,ax2)=plt.subplots(1,2)
    ax1.set_xscale('linear')
    ax1.set_yscale('linear')

    perturbations=[]
    filepaths=list((support.dir_input_files / 'rob_data').glob('*'))
    for filepath in filepaths:
        perturbation=pert(ID='LOCUST - {filename} - {quantity}'.format(filename=filepath.parts[-1],quantity=quantity_to_plot),data_format=data_format,filename=filepath)
        perturbations.append(perturbation)
    animation=FuncAnimation(fig,draw_perturbation,frames=perturbations,fargs=[ax1,ax2,fig],repeat=True,interval=1500) #cycle through phases and make animation
    plt.show()
    animation.save('LOCUST_field_{quantity}.gif'.format(field=filepath.parts[-1],quantity=quantity_to_plot),writer='pillow')
