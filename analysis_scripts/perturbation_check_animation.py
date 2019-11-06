#perturbation_check_animation.py
"""
Samuel Ward
16/08/2019
----
check perturbation with upper/lower coils rotated
---
usage:
    see README.md for usage
 
notes:         
    assumes that n2 and n6 perturbations will rotate together
---
"""

import context
from classes.input_classes.perturbation import Perturbation as pert 
from classes.input_classes.equilibrium import Equilibrium 
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import pathlib



quantity_to_plot='dB_field_R_real' #magnitude of perturbation
phase_min=0 #in degrees
phase_max=360
number_phases=100
phases=np.linspace(phase_min,phase_max,number_phases) #for each harmonic vary the phase between the upper and lower coils
n2=True #read, plot and dump n=2 harmonic?
n6=False #True #read, plot and dump n=6 harmonic?
response=True #include plasma response?
response_tag = 'response' if response else 'vacuum'  
ideal=False #ideal or resistive?
ideal_tag = 'ideal' if ideal else 'resistive'  
data_format_input='MARSF_bplas' #data format of data source
data_format_output='LOCUST' #data format to dump data to
dump_data=False #toggle whether to dump harmonics data
vminmax=[-0.0005,0.0005]

filepaths=[]
harmonics=[]
mode_numbers=[]

eq=Equilibrium(ID='Antti Snicker ITPA GEQDSK equilibrium',data_format='GEQDSK',filename=pathlib.Path('ASCOT')/'asnicker_data'/'g033143.02730')


def draw_perturbation(perturbation,fig,ax_array):
    perturbation.plot_components(R=1.7,Z=0.,phi=0.0,phase=0,i3dr=-1,LCFS=eq,limiters=False,vminmax=vminmax,absolute=False,ax_array=ax_array,fig=fig)


if n2:
    filepaths.append(pathlib.Path('ASCOT')/'dryan_data'/'33143_2730_mars_data')
    harmonics.append('n2')
    mode_numbers.append(2)
if n6:
    filepaths.append(pathlib.Path('ASCOT')/'dryan_data'/'33143_2730_n_6')
    harmonics.append('n6')
    mode_numbers.append(6)

for filepath,harmonic,mode_number in zip(filepaths,harmonics,mode_numbers): #cycle through harmonics   

    fig,((ax1,ax2),(ax3,ax4))=plt.subplots(2,2)
    ax_array=list((ax1,ax2,ax3,ax4))

    perturbations=[] #pre-read and dump all perturbations before animating
    for phase in phases: 
        perturbation=pert(ID='harmonic - {harmonic}, phase - {phase}'.format(harmonic=harmonic,phase=phase),data_format=data_format_input,filename=filepath,response=response,ideal=ideal,phase_shift=phase,bcentr=1.75660107,rmaxis=1.70210874,mode_number=mode_number)
        perturbations.append(perturbation)

    animation=FuncAnimation(fig,draw_perturbation,frames=perturbations,fargs=[fig,ax_array],repeat=True) #cycle through phases and make animation
    #plt.show()
    animation.save('a_beautiful_perturbation_check_{}.gif'.format(harmonic),writer='pillow')
