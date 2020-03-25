import sys
import numpy as np
import context
import matplotlib.pyplot as plt
from classes.input_classes.equilibrium import Equilibrium
from classes.output_classes.distribution_function import Distribution_Function
from processing import process_input
from processing import process_output
import run_scripts.utils

filename_eq='g157418.03000'
equi=Equilibrium(filename_eq,'GEQDSK',filename_eq)

#Tmax=100ms
LOCUST_file='F_24-09-2018_11-55-41.801_TOTL.dfn'
TRANSP_file='157418U46_fi_10_gc.cdf'
ASCOT_file='ascot_freia_6496371.h5'

#Tmax=260ms
LOCUST_file='F_24-09-2018_15-24-51.759_TOTL.dfn'
TRANSP_file='157418U47_fi_13_gc.cdf'
ASCOT_file='ascot_freia_6434831.h5'

LOCUST_dfn=Distribution_Function('LOCUST fast ion density (0.1s)','LOCUST',filename=LOCUST_file)
TRANSP_dfn=run_scripts.utils.TRANSP_output_FI('TRANSP fast ion density (0.1s)',filename=TRANSP_file)
ASCOT_dfn=run_scripts.utils.ASCOT_output('ASCOT fast ion density (0.1s)',filename=ASCOT_file,datatype='distribution_function')
ASCOT_dfn['V_pitch']*=-1.



#start by creating some axis objects
fig,((ax1,ax2,ax3),(ax4,ax5,ax6))=plt.subplots(2,3)

r_values=np.linspace(2,2.3,10)
z=0.

for r in r_values:

    #TOTAL R Z

    axes=['R','Z']
    TRANSP_mesh=TRANSP_dfn.plot(axes=axes,some_equilibrium=equi,LCFS=True,limiters=True,ax=ax1,fig=fig)
    ASCOT_mesh=ASCOT_dfn.plot(axes=axes,some_equilibrium=equi,LCFS=True,limiters=True,ax=ax2,fig=fig)
    LOCUST_mesh=LOCUST_dfn.plot(some_equilibrium=equi,axes=axes,LCFS=True,limiters=True,ax=ax3,fig=fig)

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

    #for point we have chosen, find closest dfn bin in real space
    diff_r=np.abs(TRANSP_dfn['R2D']-r)**2 #irregular grid so need to minimise distance to every point
    diff_z=np.abs(TRANSP_dfn['Z2D']-z)**2
    index_rz_transp=(diff_r+diff_z).argmin()
    r_transp=TRANSP_dfn['R2D'][index_rz_transp]
    z_transp=TRANSP_dfn['Z2D'][index_rz_transp]

    index_r_ascot=np.abs(ASCOT_dfn['R']-r).argmin() #regular grids so can use this method
    index_z_ascot=np.abs(ASCOT_dfn['Z']-z).argmin()-1 #required to match grids since lack of numerical precision
    r_ascot=ASCOT_dfn['R'][index_r_ascot]
    z_ascot=ASCOT_dfn['Z'][index_z_ascot]

    index_r_locust=np.abs(LOCUST_dfn['R']-r).argmin()
    index_z_locust=np.abs(LOCUST_dfn['Z']-z).argmin()
    r_locust=LOCUST_dfn['R'][index_r_locust]
    z_locust=LOCUST_dfn['Z'][index_z_locust]

    ax1.scatter(r_transp,z_transp,color='m',s=0.5)
    ax1.scatter(r_ascot,z_ascot,color='g',s=0.5)
    ax1.scatter(r_locust,z_locust,color='b',s=0.5)


    #LOCALISED ENERGY PITCH 
    axes=[index_rz_transp,slice(None),slice(None)]
    TRANSP_mesh=TRANSP_dfn.plot(axes=axes,ax=ax4,fig=fig) #plotting with slices and specifying the point plots from the non-integrated distribution function

    axes=[slice(None),slice(None),index_r_ascot,index_z_ascot]
    ASCOT_mesh=ASCOT_dfn.plot(axes=axes,ax=ax5,fig=fig)

    axes=[0,slice(None),slice(None),index_r_locust,index_z_locust]
    LOCUST_mesh=LOCUST_dfn.plot(axes=axes,ax=ax6,fig=fig)

    #fig.colorbar(TRANSP_mesh,ax=ax1)
    #fig.colorbar(ASCOT_mesh,ax=ax2)
    #fig.colorbar(LOCUST_mesh,ax=ax3)

    ax4.set_xlabel('E [eV]')
    ax4.set_ylabel('V$_{par}$/V')
    ax5.set_xlabel('E [eV]')
    ax5.set_ylabel('V$_{par}$/V')
    ax6.set_xlabel('E [eV]')
    ax6.set_ylabel('V$_{par}$/V')

    ax4.set_xlim([0,90000])
    ax4.set_ylim([-1.,1.])
    ax5.set_xlim([0,90000])
    ax5.set_ylim([-1.,1.])
    ax6.set_xlim([0,90000])
    ax6.set_ylim([-1.,1.])

    plt.draw()
    plt.pause(0.1)
    plt.cla()
    #plt.show()