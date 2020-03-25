#script to animate R Z distribution function for TRANSP runs

import context
from classes.input_classes.equilibrium import Equilibrium as eq
from run_scripts import utils
import matplotlib.pyplot as plt

run_1='U46'
run_2='U47'

filename_eq='g157418.03000'
equi=eq('','GEQDSK',filename_eq)

number_files=13 #number files per run
FI_array=[]

for file in range(number_files): #add files
    file+=1 #starts from file 1
    filename='157418{}_fi_{}_gc.cdf'.format(run_1,file)
    FI_array.append(utils.TRANSP_output_FI(filename))

for file in range(number_files): #add second set of files
    file+=1 #starts from file 1
    filename='157418{}_fi_{}_gc.cdf'.format(run_2,file)
    FI_array.append(utils.TRANSP_output_FI(filename)) 


fig,(ax)=plt.subplots(1,1) #initialise plot

for FI in FI_array:
    mesh=FI.plot(axes=['R','Z'],some_equilibrium=equi,LCFS=True,limiters=True,real_scale=True,ax=ax,fig=fig)
    cbar=fig.colorbar(mesh,ax=ax)
    plt.draw()
    plt.pause(0.1)
    cbar.remove()
    plt.cla()