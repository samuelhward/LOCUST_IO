#quick test to check normalisation between distribution functions produced by ASCOT and LOCUST for U74 inputs
#for the first time, these ASCOT runs use scaled beam deposition weights such that power deposited into the plasma adds up to Pdep(LOCUST)=1W 


import sys
import numpy as np
import context
import matplotlib.pyplot as plt
from classes.input_classes.equilibrium import Equilibrium
from classes.input_classes.beam_deposition import Beam_Deposition
from classes.input_classes.wall import Wall
from classes.output_classes.distribution_function import Distribution_Function
from classes.output_classes.particle_list import Final_Particle_List
from processing import process_input
from processing import process_output
import run_scripts.utils
import constants

filename_eq='g157418.03000'
equi=Equilibrium(filename_eq,'GEQDSK',filename_eq)

#no DSOLCOL, DSOLCOL
ascot_files=['ascot/test_2/ascot_freia_7022.h5','ascot/test_4/ascot_freia_7023.h5'] 
locust_files=['locust/run_1/F_04-12-2018_14-11-36.625_TOTL.dfn','locust/run_2/F_05-12-2018_01-15-34.057_TOTL.dfn'] 
locust_FINTs=['locust/run_1/FINT.dat','locust/run_2/FINT.dat']
run_types=['no -DSOLCOL','-DSOLCOL']

fig,((ax1,ax2,ax3),(ax4,ax5,ax6))=plt.subplots(2,3)

axes=[[ax1,ax2,ax3],[ax4,ax5,ax6]]

for ASCOT_file,LOCUST_file,ax_group,LOCUST_FINT,run_type in zip(ascot_files,locust_files,axes,locust_FINTs,run_types):

	#distribution functions
	LOCUST_dfn=Distribution_Function('LOCUST fast ion density (0.1s) - limiter radius = 1.5','LOCUST',filename=LOCUST_file)
	ASCOT_dfn=run_scripts.utils.ASCOT_output('ASCOT fast ion density (0.1s) - limiter radius = 1.5',filename=ASCOT_file,datatype='distribution_function')

	axes=['R','Z']
    LOCUST_mesh=LOCUST_dfn.plot(axes=axes,ax=ax_group[0],LCFS=equi,fig=fig)
    ASCOT_mesh=ASCOT_dfn.dfn_plot(axes=axes,ax=ax_group[1],LCFS=equi,fig=fig)

	for ax,mesh in zip([ax_group[0],ax_group[1]],[LOCUST_mesh,ASCOT_mesh]):
	    cbar=fig.colorbar(mesh,ax=ax,orientation='vertical')

	#particle lists
    mom_locust=run_scripts.utils.FINT_LOCUST(ID='LOCUST moments - radius = 1.5',filename=LOCUST_FINT)
    mom_ascot=Final_Particle_List(ID='ASCOT final particle list - radius = 1.5',data_format='ASCOT',filename=ASCOT_file)
	i=np.where(mom_ascot['status_flag']==3)[0] #separate ASCOT markers which have hit PFC components
    mom_ascot['PFC_power']=np.full(len(mom_locust['time']),np.sum(mom_ascot['E'][i]*constants.species_charge*mom_ascot['weight'][i])) #place against same time base as LOCUST

    ax_group[2].plot(mom_locust['time'],mom_locust['PFC_power'],'g-') #ASCOT and LOCUST on same time base
    ax_group[2].plot(mom_locust['time'],mom_ascot['PFC_power'],'b-') 
	ax_group[2].legend(('LOCUST','ASCOT'))
	ax_group[2].set_title('PFC power loss (watts)')

	ax_group[0].set_ylabel(run_type)

plt.show()


