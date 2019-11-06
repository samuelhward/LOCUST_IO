#LOCUST_field_data_gen_animation.py
"""
Samuel Ward
07/04/2019
----
generate multiple LOCUST fields with upper/lower coils rotated
---
usage:
    see README.md for usage
 
notes:         
    assumes that n2 and n6 perturbations will rotate together
---
"""


import context
from classes.input_classes.perturbation import Perturbation as pert  
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import pathlib



quantity_to_plot='dB_field_R_real' #magnitude of perturbation
phase_min=0 #in degrees
phase_max=180
number_phases=10
n2=False#True #read, plot and dump n=2 harmonic?
n6=True #read, plot and dump n=6 harmonic?
response=False #include plasma response?
response_tag = 'response' if response else 'vacuum'  
ideal=True #ideal or resistive?
ideal_tag = 'ideal' if ideal else 'resistive'  
data_format_input='MARSF_bplas' #data format of data source
data_format_output='LOCUST' #data format to dump data to
dump_data=False #toggle whether to dump harmonics data

def draw_perturbation(perturbation,ax,fig):    
    ax.collections=[]
    perturbation.plot(key=quantity_to_plot,fig=fig,ax=ax)

phases=np.linspace(phase_min,phase_max,number_phases) #for each harmonic vary the phase between the upper and lower coils

filepaths=[]
harmonics=[]

if n2:
    filepaths.append(pathlib.Path('ASCOT')/'dryan_data'/'33143_2730_mars_data')
    harmonics.append('n2')
if n6:
    filepaths.append(pathlib.Path('ASCOT')/'dryan_data'/'33143_2730_n_6')
    harmonics.append('n6')

for filepath,harmonic in zip(filepaths,harmonics): #cycle through harmonics   

    fig,ax=plt.subplots(1)
    ax.set_xscale('linear')
    ax.set_yscale('linear')

    perturbations=[] #pre-read and dump all perturbations before animating
    for phase in phases: 
        perturbation=pert(ID='harmonic - {harmonic}, phase - {phase}'.format(harmonic=harmonic,phase=phase),data_format=data_format_input,filename=filepath,response=response,ideal=ideal,phase_shift=phase,bcentr=1.75660107,rmaxis=1.70210874)
        perturbations.append(perturbation)
        perturbation.dump_data(data_format=data_format_output,filename=pathlib.Path('LOCUST') / '3D_field_scan' / 'BPLASMA_{harmonic}_{phase}_{ideal}_{response}.dat'.format(harmonic=harmonic,phase=phase,ideal=ideal_tag,response=response_tag))

    animation=FuncAnimation(fig,draw_perturbation,frames=perturbations,fargs=[ax,fig],repeat=True) #cycle through phases and make animation
    plt.show()
    animation.save('a_beautiful_perturbation_{}.gif'.format(harmonic),writer='pillow')
