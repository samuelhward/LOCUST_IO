
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm #get colourmaps

import context
from classes.input_classes.number_density import Number_Density as nd
from classes.input_classes.temperature import Temperature as t
from classes.input_classes.equilibrium import Equilibrium
import run_scripts.utils
import copy

input_Te_file='profile_Te.dat'
input_Ti_file='profile_Ti.dat'
input_ne_file='profile_ne.dat'
output_file='LOCUST_15-11-2018_01-00-56.342.h5'

Te_before=t(ID='electron temperature',data_format='LOCUST',filename=input_Te_file,species='electrons')
Ti_before=t(ID='ion density',data_format='LOCUST',filename=input_Ti_file,species='ions')
ne_before=nd(ID='electron density',data_format='LOCUST',filename=input_ne_file,species='electrons')

Te_after=t(ID='electron temperature',data_format='LOCUST_h5',filename=output_file,species='electrons')
Ti_after=t(ID='ion density',data_format='LOCUST_h5',filename=output_file,species='ions')
ne_after=nd(ID='electron density',data_format='LOCUST_h5',filename=output_file,species='electrons')

fig,(ax1,ax2,ax3)=plt.subplots(1,3)

Te_after.plot(ax=ax1,fig=fig,colmap='red')
Te_before.plot(ax=ax1,fig=fig,colmap='blue')

Ti_after.plot(ax=ax2,fig=fig,colmap='red')
Ti_before.plot(ax=ax2,fig=fig,colmap='blue')

ne_after.plot(ax=ax3,fig=fig,colmap='red')
ne_before.plot(ax=ax3,fig=fig,colmap='blue')

ax1.legend(['after simulation','before simulation'])

for ax in [ax1,ax2,ax3]:
    ax.set_ylabel('')
plt.show()



#optional writing to file

ni_after=copy.deepcopy(ne_after)
file_eq='g157418.03000'
eq=Equilibrium(ID='equilibrium',data_format='GEQDSK',filename=file_eq)
rot=np.zeros(len(Ti_after['T']))
run_scripts.utils.dump_profiles_ASCOT(filename='input.plasma_1d_extrap',temperature_i=Ti_after,temperature_e=Te_after,density_i=ni_after,density_e=ne_after,rotation_toroidal=rot)

