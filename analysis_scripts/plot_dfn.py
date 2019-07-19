#plot_dfn.py
 
"""
Samuel Ward
18/04/2019
----
script to plot time slices for MAST V and V
---
usage:
 
notes:         
    uses get_filenames.py to retrieve the filenames first via globs
---
"""

import numpy as np
import context 
from classes.output_classes.distribution_function import Distribution_Function as dfn
from classes.input_classes.equilibrium import Equilibrium as equi
from classes.output_classes.particle_list import Final_Particle_List as fpl
from classes.input_classes.beam_deposition import Beam_Deposition as bd
from classes.input_classes.wall import Wall as wall
from run_scripts.utils import  TRANSP_output as output_T
from run_scripts.utils import TRANSP_output_FI as dfn_t
import support
import constants
import settings
import pathlib

shot='29034'
run='W03'
dimension_LOCUST='2D'
dimension_ASCOT='2D'
extra_info_LOCUST='NUBEAM_birth_GC'
extra_info_ASCOT='NUBEAM_birth_GC'
ASCOT_beam_depo_data_format='ASCOT_GC'
ASCOT_beam_depo_filename=pathlib.Path('ASCOT') / (shot+run) / 'input.particles_BBNBI_FO_reduced'
ASCOT_beam_depo_filename=pathlib.Path('ASCOT') / (shot+run) / 'input.particles_NUBEAM_birth_GC'
time_slices=np.linspace(0.37,0.46,10)
tinit=0.36 #start time
cmap_LOCUST=settings.cmap_g
cmap_ASCOT=settings.cmap_b
cmap_TRANSP=settings.cmap_r
number_bins=5

from get_filepaths import get_filepaths #grab all the filenames 
files_locust,files_ascot,files_transp=get_filepaths(shot=shot,run=run,dimension_LOCUST=dimension_LOCUST,dimension_ASCOT=dimension_ASCOT,extra_info_LOCUST=extra_info_LOCUST,extra_info_ASCOT=extra_info_ASCOT)

file_eq=pathlib.Path('TRANSP') / (shot+run) / ('g'+shot)
file_wall=pathlib.Path('TRANSP') / (shot+run) / ('OMF'+shot+'.LIM')

eq=equi('',data_format='GEQDSK',filename=file_eq)
wall_MAST=wall('',data_format='UFILE',filename=file_wall)

#grab the ASCOT beam deposition froom the inistate to renormalise the ASCOT dfn to the correct power deposition
beam_depo=bd(ID='ASCOT beam depo',data_format=ASCOT_beam_depo_data_format,filename=ASCOT_beam_depo_filename)
beam_power=1. #desired beam power
if 'E' not in beam_depo.data: 
    beam_depo['E']=0.5*constants.mass_deuteron*(beam_depo['V_R']**2+beam_depo['V_tor']**2+beam_depo['V_Z']**2)/constants.species_charge
Pdep_ASCOT=np.sum(beam_depo['E']*constants.species_charge*beam_depo['weight'])

for file_locust,file_ascot,file_transp,time_slice in zip(files_locust,files_ascot,files_transp,time_slices):

    #read DFNs
    dfn_LOCUST=dfn('LOCUST time = {}'.format(time_slice-tinit),data_format='LOCUST',filename=file_locust)
    dfn_TRANSP=dfn_t('TRANSP time = {}'.format(time_slice-tinit),filename=file_transp)
    dfn_ASCOT=dfn('ASCOT time = {}'.format(time_slice-tinit),data_format='ASCOT',filename=file_ascot)

    #scale DFNs to match 1W beam deposited power
    N_LOCUST=dfn_LOCUST.transform(axes=['N'])['dfn']
    N_TRANSP=dfn_TRANSP.dfn_integrate()['dfn']
    N_ASCOT=dfn_ASCOT.transform(axes=['N'])['dfn']
    file_transp_output=pathlib.Path(file_transp).parents[0] / (shot+run+'.CDF')

    output_TRANSP=output_T('',filename=file_transp_output)
    BPCAP_index=np.abs(output_TRANSP['TIME']-time_slice).argmin()
    BPCAP=output_TRANSP['BPCAP'][BPCAP_index]
    dfn_TRANSP['dfn']/=BPCAP
    dfn_ASCOT['dfn']*=beam_power/(Pdep_ASCOT) #these are calculated above outside the loop
    dfn_ASCOT['V_pitch']*=-1
    #dfn_TRANSP['dfn']/=N_TRANSP #normalise distribution functions against eachother
    #dfn_LOCUST['dfn']/=N_LOCUST
    #dfn_ASCOT['dfn']/=dfn_ASCOT

    import matplotlib.pyplot as plt

    fig,((ax1,ax2,ax3),(ax4,ax5,ax6))=plt.subplots(2,3)

    axes=['E','V_pitch']
    for ax in [ax1,ax2,ax3]:
        ax.set_xlim(0,80000)
        ax.set_ylim(-1.1,1.1)
    real_scale=False
    limiters=False
    LCFS=False
    vminmax=[1.2e6,1.2e8]
    LOCUST_mesh=dfn_LOCUST.plot(axes=axes,fig=fig,ax=ax1,LCFS=LCFS,limiters=limiters,real_scale=real_scale,vminmax=vminmax,fill=False,colmap=cmap_LOCUST,number_bins=number_bins)
    TRANSP_mesh=dfn_TRANSP.dfn_plot(axes=axes,fig=fig,ax=ax1,LCFS=LCFS,limiters=limiters,real_scale=real_scale,vminmax=vminmax,fill=False,colmap=cmap_TRANSP,number_bins=number_bins)
    ASCOT_mesh=dfn_ASCOT.plot(axes=axes,fig=fig,ax=ax1,LCFS=LCFS,limiters=limiters,real_scale=real_scale,vminmax=vminmax,fill=False,colmap=cmap_ASCOT,number_bins=number_bins)

    if True:
        for ax,mesh in zip([ax1,ax2,ax3],[LOCUST_mesh,TRANSP_mesh,ASCOT_mesh]):
            cbar=fig.colorbar(mesh,ax=ax,orientation='vertical')

    #'''
    axes=['R','Z']
    real_scale=True
    limiters=wall_MAST
    LCFS=eq
    vminmax=[3e10,3e12]
    LOCUST_mesh=dfn_LOCUST.plot(axes=axes,fig=fig,ax=ax4,LCFS=LCFS,limiters=limiters,real_scale=real_scale,vminmax=vminmax,fill=False,colmap=cmap_LOCUST,number_bins=number_bins)
    TRANSP_mesh=dfn_TRANSP.dfn_plot(axes=axes,fig=fig,ax=ax4,LCFS=LCFS,limiters=limiters,real_scale=real_scale,vminmax=vminmax,fill=False,colmap=cmap_TRANSP,number_bins=number_bins)
    ASCOT_mesh=dfn_ASCOT.plot(axes=axes,fig=fig,ax=ax4,LCFS=LCFS,limiters=limiters,real_scale=real_scale,vminmax=vminmax,fill=False,colmap=cmap_ASCOT,number_bins=number_bins)
    #'''

    #'''
    axes=['V_pitch']
    dfn_LOCUST.plot(axes=axes,fig=fig,ax=ax5,colmap=cmap_LOCUST)
    dfn_TRANSP.dfn_plot(axes=axes,fig=fig,ax=ax5,colmap=cmap_TRANSP)
    dfn_ASCOT.plot(axes=axes,fig=fig,ax=ax5,colmap=cmap_ASCOT)
    #'''
   
    #'''
    axes=['E']
    dfn_LOCUST.plot(axes=axes,fig=fig,ax=ax6,colmap=cmap_LOCUST)
    dfn_TRANSP.dfn_plot(axes=axes,fig=fig,ax=ax6,colmap=cmap_TRANSP)
    dfn_ASCOT.plot(axes=axes,fig=fig,ax=ax6,colmap=cmap_ASCOT)
    #'''

    #'''
    if True:
        for ax,mesh in zip([ax4,ax5,ax6],[LOCUST_mesh,TRANSP_mesh,ASCOT_mesh]):
            cbar=fig.colorbar(mesh,ax=ax,orientation='vertical')
    #'''

    plt.show()

#################################
 
##################################################################
 
###################################################################################################
