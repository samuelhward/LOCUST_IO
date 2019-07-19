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

filename_eq='g157418.03000'
equi=Equilibrium(filename_eq,'GEQDSK',filename_eq)

#Tmax=100ms
LOCUST_file='F_24-09-2018_11-55-41.801_TOTL.dfn'
TRANSP_file='157418U46_fi_10_gc.cdf'
ASCOT_file='ascot_freia_6496371.h5'

LOCUST_dfn=Distribution_Function('LOCUST fast ion density (0.1s)','LOCUST',filename=LOCUST_file)
TRANSP_dfn=run_scripts.utils.TRANSP_output_FI('TRANSP fast ion density (0.1s)',filename=TRANSP_file)
ASCOT_dfn=run_scripts.utils.ASCOT_output('ASCOT fast ion density (0.1s)',filename=ASCOT_file,datatype='distribution_function')
ASCOT_dfn['V_pitch']*=-1.



#start by creating some axis objects
fig,(ax1,ax2,ax3)=plt.subplots(1,3)




#R Z
'''
axes=['R','Z']
TRANSP_mesh=TRANSP_dfn.dfn_plot(axes=axes,some_equilibrium=equi,LCFS=True,limiters=True,ax=ax1,fig=fig)
ASCOT_mesh=ASCOT_dfn.dfn_plot(axes=axes,some_equilibrium=equi,LCFS=True,limiters=True,ax=ax2,fig=fig)
LOCUST_mesh=plot_output.plot_distribution_function(LOCUST_dfn,some_equilibrium=equi,axes=axes,LCFS=True,limiters=True,ax=ax3,fig=fig)

ax1.set_xlim([np.amin(equi['R_1D']),np.amax(equi['R_1D'])])
ax2.set_xlim([np.amin(equi['R_1D']),np.amax(equi['R_1D'])])
ax3.set_xlim([np.amin(equi['R_1D']),np.amax(equi['R_1D'])])
ax1.set_ylim([np.amin(equi['Z_1D']),np.amax(equi['Z_1D'])])
ax2.set_ylim([np.amin(equi['Z_1D']),np.amax(equi['Z_1D'])])
ax3.set_ylim([np.amin(equi['Z_1D']),np.amax(equi['Z_1D'])])

ax1.set_xlabel('R [m]')
ax1.set_ylabel('Z [m]')
ax2.set_xlabel('R [m]')
ax2.set_ylabel('Z [m]')
ax3.set_xlabel('R [m]')
ax3.set_ylabel('Z [m]')
'''


#ENERGY PITCH 
axes=['E','V_pitch']
TRANSP_mesh=TRANSP_dfn.dfn_plot(axes=axes,ax=ax1,fig=fig)
ASCOT_mesh=ASCOT_dfn.dfn_plot(axes=axes,ax=ax2,fig=fig)
LOCUST_mesh=plot_output.plot_distribution_function(LOCUST_dfn,axes=axes,ax=ax3,fig=fig)

#fig.colorbar(TRANSP_mesh,ax=ax1)
#fig.colorbar(ASCOT_mesh,ax=ax2)
#fig.colorbar(LOCUST_mesh,ax=ax3)

ax1.set_xlabel('E [eV]')
ax1.set_ylabel('V$_{par}$/V')
ax2.set_xlabel('E [eV]')
ax2.set_ylabel('V$_{par}$/V')
ax3.set_xlabel('E [eV]')
ax3.set_ylabel('V$_{par}$/V')

ax1.set_xlim([0,90000])
ax1.set_ylim([-1.,1.])
ax2.set_xlim([0,90000])
ax2.set_ylim([-1.,1.])
ax3.set_xlim([0,90000])
ax3.set_ylim([-1.,1.])

plt.show()