#script to plot 1D energy distribution of U46 and U47 TRANSP runs
import sys
import numpy as np
import context
import matplotlib.pyplot as plt
from classes.input_classes.equilibrium import Equilibrium
from classes.output_classes.distribution_function import Distribution_Function
from processing import plot_input
from processing import plot_output
from processing import process_input
from processing import process_output
import run_scripts.utils


#read equilibrium here
filename_eq='g157418.03000'
equi=Equilibrium(filename_eq,'GEQDSK',filename_eq)

#read and integrate the TRANSP fast ion files here
number_files=13 #number files per run
TRANSP_FI_array=[]
for file in range(number_files): #add the U46 files
    file+=1 #starts from file 1
    filename='157418U46_fi_{}_gc.cdf'.format(file)
    TRANSP_FI_array.append(run_scripts.utils.TRANSP_output_FI('TRANSP fast ion density',filename))

for file in range(number_files): #add the U47 files
    file+=1 #starts from file 1
    filename='157418U47_fi_{}_gc.cdf'.format(file)
    TRANSP_FI_array.append(run_scripts.utils.TRANSP_output_FI('TRANSP fast ion density',filename))

time_offset=3.0
for FI in TRANSP_FI_array:
    FI['TIME']-=time_offset #subtract shot time offset

#read LOCUST fast ion diffusion here
filename_LOCUST='FINT.dat' #default
LOCUST_dfn=run_scripts.utils.FINT_LOCUST('LOCUST fast ion density',filename_LOCUST)




#start by creating some axis objects
fig,(ax1,ax2)=plt.subplots(1,2)

#plot TRANSP energy diffusion
TRANSP_mesh=TRANSP_FI_array[0].plot(axes=['E','time'],TRANSP_output_FI_list=[FI for FI in TRANSP_FI_array],ax=ax1,fig=fig)

#plot LOCUST energy diffusion
LOCUST_mesh=LOCUST_dfn.plot(some_equilibrium=equi,axes=['E','time'],ax=ax2,fig=fig)

#fig.colorbar(TRANSP_mesh,ax=ax1)
#fig.colorbar(LOCUST_mesh,ax=ax2)

ax1.set_ylim([0,100000])
ax1.set_xlim([0,.26])
ax2.set_ylim([0,100000])
ax2.set_xlim([0,.26])

ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Energy [eV]')
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('Energy [eV]')


plt.show()