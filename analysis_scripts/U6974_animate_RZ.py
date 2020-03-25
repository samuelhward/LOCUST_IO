import sys
import numpy as np
import context
import matplotlib.pyplot as plt
from classes.input_classes.equilibrium import Equilibrium
from classes.output_classes.distribution_function import Distribution_Function
from classes.input_classes.wall import Wall
from processing import process_input
from processing import process_output
import run_scripts.utils

filename_eq='g157418.03000'
equi=Equilibrium(filename_eq,'GEQDSK',filename_eq)

wall_files='LCFS_DIII-D.dat'
radii=['1.05','1.10','1.20','1.30','1.40','1.50'] #radii for limiter profiles
LOCUST_moments=['LOCUST_15-11-2018_01-00-56.342.h5','LOCUST_15-11-2018_00-43-29.482.h5','LOCUST_15-11-2018_00-48-33.477.h5','LOCUST_15-11-2018_00-53-28.875.h5','LOCUST_15-11-2018_11-38-05.889.h5','LOCUST_15-11-2018_11-38-00.187.h5']

TRANSP_files=['157418U69_fi_1_gc.cdf','157418U70_fi_1_gc.cdf','157418U71_fi_1_gc.cdf','157418U72_fi_1_gc.cdf','157418U73_fi_1_gc.cdf','157418U74_fi_1_gc.cdf']

ASCOT_files=['ascot_freia_8849968.h5','ascot_freia_8849970.h5','ascot_freia_8857826.h5','ascot_freia_8857830.h5','ascot_freia_8850120.h5','ascot_freia_8850121.h5']
ASCOT_run='ascot/run_1/'

ASCOT_files=['ascot_freia_9242197.h5','ascot_freia_9242198.h5','ascot_freia_9242199.h5','ascot_freia_9242200.h5','ascot_freia_9242201.h5','ascot_freia_9242202.h5']
ASCOT_run='ascot/run_2/' #changed ORBITMETHOD to 4, added new kinetic profiles

ASCOT_files=['ascot_freia_9305131.h5','ascot_freia_9305176.h5','ascot_freia_9305207.h5','ascot_freia_9305216.h5','ascot_freia_9305222.h5','ascot_freia_9305226.h5']
ASCOT_run='ascot/run_3/' #changed ORBITMETHOD back to 1

LOCUST_files=['F_15-11-2018_01-00-56.342_TOTL.dfn','F_15-11-2018_00-43-29.482_TOTL.dfn','F_15-11-2018_00-48-33.477_TOTL.dfn','F_15-11-2018_00-53-28.875_TOTL.dfn','F_15-11-2018_11-38-05.889_TOTL.dfn','F_15-11-2018_11-38-00.187_TOTL.dfn']
LOCUST_moments=['LOCUST_15-11-2018_01-00-56.342.h5','LOCUST_15-11-2018_00-43-29.482.h5','LOCUST_15-11-2018_00-48-33.477.h5','LOCUST_15-11-2018_00-53-28.875.h5','LOCUST_15-11-2018_11-38-05.889.h5','LOCUST_15-11-2018_11-38-00.187.h5']
LOCUST_run='locust/run_1/'

#LOCUST_files=['F_19-11-2018_17-20-33.749_TOTL.dfn','F_19-11-2018_14-08-59.547_TOTL.dfn','F_19-11-2018_14-09-09.941_TOTL.dfn','F_19-11-2018_14-10-19.388_TOTL.dfn','F_19-11-2018_14-10-29.187_TOTL.dfn','F_19-11-2018_21-47-44.702_TOTL.dfn']
#LOCUST_moments=['LOCUST_19-11-2018_17-20-33.749.h5','LOCUST_19-11-2018_14-08-59.547.h5','LOCUST_19-11-2018_14-09-09.941.h5','LOCUST_19-11-2018_14-10-19.388.h5','LOCUST_19-11-2018_14-10-29.187.h5','LOCUST_19-11-2018_21-47-44.702.h5']
#LOCUST_run='locust/run_2/' #added -DSOLCOL

#start by creating some axis objects
fig,((ax1,ax2,ax3),(ax4,ax5,ax6))=plt.subplots(2,3)

for LOCUST_file,TRANSP_file,ASCOT_file,radius in zip(LOCUST_files,TRANSP_files,ASCOT_files,radii):

    wall=Wall('limiter - '+radius,data_format='LOCUST_2D',filename=wall_files+'_'+radius)

    LOCUST_dfn=Distribution_Function(ID='LOCUST - '+radius,data_format='LOCUST',filename=LOCUST_run+LOCUST_file)
    TRANSP_dfn=run_scripts.utils.TRANSP_output_FI(ID='TRANSP - '+radius,filename=TRANSP_file)
    ASCOT_dfn=run_scripts.utils.ASCOT_output(ID='ASCOT - '+radius,filename=ASCOT_run+ASCOT_file,datatype='distribution_function')
    ASCOT_dfn['V_pitch']*=-1.

    r_values=np.linspace(1.5,2.3,1)
    z=0.

    for r in r_values:

        #TOTAL R Z

        axes=['R','Z']
        TRANSP_mesh=TRANSP_dfn.plot(axes=axes,limiters=wall,LCFS=equi,ax=ax1,fig=fig)
        ASCOT_mesh=ASCOT_dfn.plot(axes=axes,limiters=wall,LCFS=equi,ax=ax2,fig=fig)
        LOCUST_mesh=LOCUST_dfn.plot(axes=axes,limiters=wall,LCFS=equi,ax=ax3,fig=fig)


        for ax in [ax1,ax2,ax3]:

            ax.set_xlim([np.amin(equi['R_1D']),np.amax(equi['R_1D'])])
            ax.set_ylim([np.amin(equi['Z_1D']),np.amax(equi['Z_1D'])])
            ax.set_xlabel('R [m]')
            ax.set_ylabel('Z [m]')

        #for point we have chosen, find closest dfn bin in real space
        diff_r=np.abs(TRANSP_dfn['R2D']-r) #irregular grid so need to minimise distance to every point
        diff_z=np.abs(TRANSP_dfn['Z2D']-z)
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

        ax1.scatter(r_transp,z_transp,color='m',s=1.5)
        ax1.scatter(r_ascot,z_ascot,color='g',s=1.5)
        ax1.scatter(r_locust,z_locust,color='b',s=1.5)


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

        for ax in [ax4,ax5,ax5]:
            ax.set_xlabel('E [eV]')
            ax.set_ylabel('V$_{par}$/V')
            ax.set_xlim([0,90000])
            ax.set_ylim([-1.,1.])

        plt.draw()
        plt.pause(0.1)

    for ax in [ax1,ax2,ax3,ax4,ax5,ax6]:
        ax.cla()
        #plt.show()
