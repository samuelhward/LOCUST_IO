import context
from classes.input_classes.beam_deposition import Beam_Deposition
from classes.input_classes.equilibrium import Equilibrium
import processing.utils
import run_scripts.utils
import processing.process_input as prin
import constants
import numpy as np
import matplotlib.pyplot as plt

filename_eq='g157418.03000'
equi_0=Equilibrium('GEQDSKFIX = 0','GEQDSK',filename_eq,GEQDSKFIX=0)
equi_1=Equilibrium('GEQDSKFIX = 1','GEQDSK',filename_eq,GEQDSKFIX=1)


for eq in [equi_0,equi_1]:
    eq.set(B_field=prin.B_calc(eq))
    eq.set(B_field_R=eq['B_field'][:,:,0])
    eq.set(B_field_tor=eq['B_field'][:,:,1])
    eq.set(B_field_Z=eq['B_field'][:,:,2])


fig,(ax1,ax2)=plt.subplots(1,2)
key='psirz'
equi_0.plot(key=key,ax=ax1,LCFS=True)
equi_1.plot(key=key,ax=ax2,LCFS=True)
plt.show()
