#plot_dfn_single.py
 
"""
Samuel Ward
18/04/2019
----
script to plot time slices for MAST V and V - single plot
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
import matplotlib.pyplot as plt

cmap_LOCUST=settings.cmap_g
cmap_ASCOT=settings.cmap_b
cmap_TRANSP=settings.cmap_r

shot='29034'
run='W03'

dimensions=['2D','2D','2D','2D']
codes=['ASCOT','ASCOT','LOCUST','TRANSP']
extra_infos=['BBNBI_birth_FO','BBNBI_birth_GC','BBNBI_birth_GC_noDCTDC2','NUBEAM_birth_GC']
data_formats=['ASCOT','ASCOT','LOCUST',None]
cmaps=[cmap_ASCOT,cmap_ASCOT,cmap_LOCUST,cmap_TRANSP]

ASCOT_beam_depo_data_formats=['ASCOT_FO','ASCOT_FO',None]#,]'ASCOT_FO'] #empty spaces indicate other code is currently being plotted
ASCOT_beam_depo_filenames=[
                            pathlib.Path('ASCOT') / (shot+run) / 'input.particles_BBNBI_FO_reduced',
                            pathlib.Path('ASCOT') / (shot+run) / 'input.particles_BBNBI_FO_reduced',
                            None,
                            None] #path in input_files

time_slices=np.linspace(0.37,0.46,10)
tinit=0.36 #start time

file_eq=pathlib.Path('TRANSP') / (shot+run) / ('g'+shot)
file_wall=pathlib.Path('TRANSP') / (shot+run) / ('OMF'+shot+'.LIM')
eq=equi('',data_format='GEQDSK',filename=file_eq)
wall_MAST=wall('',data_format='UFILE',filename=file_wall)

from get_filepaths_single import get_filepaths_single #grab all the filenames 
files=[] #list of list, each high list belongs to a run, with each sublist a time slice
for dimension,extra_info,code in zip(dimensions,extra_infos,codes):
    files.append(get_filepaths_single(shot=shot,run=run,dimension=dimension,extra_info=extra_info,code=code))

for counter_time_slice,time_slice in enumerate([time_slices[-1]]): #XXX hack to plot final time slice

    fig,(ax1)=plt.subplots(1)

    for counter_code,(time_slices_for_this_run,code,data_format,cmap) in enumerate(zip(files,codes,data_formats,cmaps)):
        file_this_time_slice=time_slices_for_this_run[len(time_slices)-1]#counter_time_slice] XXX hack to plot final time slice 

        #read DFNs
        if code=='TRANSP':
            DFN=dfn_t('TRANSP time = {}ms'.format(int(1000*(time_slice-tinit))),filename=file_this_time_slice)
            N=DFN.dfn_integrate()['dfn']
            #scale DFNs to match 1W beam deposited power
            file_transp_output=pathlib.Path(file_this_time_slice).parents[0] / (shot+run+'.CDF')
            output_TRANSP=output_T('',filename=file_transp_output)
            BPCAP_index=np.abs(output_TRANSP['TIME']-time_slice).argmin()
            BPCAP=output_TRANSP['BPCAP'][BPCAP_index]
            DFN['dfn']/=BPCAP
        else:
            DFN=dfn('time = {}ms'.format(int(1000*(time_slice-tinit))),data_format=data_format,filename=file_this_time_slice)
            N=DFN.transform(axes=['N'])['dfn']
            
            if code=='ASCOT': #grab the ASCOT beam deposition froom the inistate to renormalise the ASCOT dfn to the correct power deposition
                beam_power=1. #desired beam power
                ASCOT_beam_depo_data_format=ASCOT_beam_depo_data_formats[counter_code]
                ASCOT_beam_depo_filename=ASCOT_beam_depo_filenames[counter_code]
                ASCOT_beam_depo=bd(ID='ASCOT beam depo',data_format=ASCOT_beam_depo_data_format,filename=ASCOT_beam_depo_filename)
                if 'E' not in ASCOT_beam_depo.data: 
                    ASCOT_beam_depo['E']=0.5*constants.mass_deuteron*(ASCOT_beam_depo['V_R']**2+ASCOT_beam_depo['V_tor']**2+ASCOT_beam_depo['V_Z']**2)/constants.species_charge
                Pdep_ASCOT=np.sum(ASCOT_beam_depo['E']*constants.species_charge*ASCOT_beam_depo['weight'])
                DFN['dfn']*=beam_power/Pdep_ASCOT 
                DFN['V_pitch']*=-1

        DFN['E']/=1000. #get things in KeV
        DFN['dE']*=1000.

        axes=['E']
        #ax1.set_xlim(0,100)
        #ax1.set_ylim(-1.1,1.1)
        real_scale=False
        limiters=False#wall_MAST
        LCFS=False#eq
        fill=False
        number_bins=10
        vminmax=[0,1.e8]
        mesh=DFN.plot(axes=axes,fig=fig,ax=ax1,LCFS=LCFS,limiters=limiters,real_scale=real_scale,fill=fill,vminmax=vminmax,colmap=cmap,number_bins=number_bins)
        #TRANSP_mesh=dfn_TRANSP.dfn_plot(axes=axes,fig=fig,ax=ax1,LCFS=LCFS,limiters=limiters,real_scale=real_scale,vminmax=vminmax,fill=fill,colmap=cmap_TRANSP,number_bins=number_bins)
        #ASCOT_mesh=dfn_ASCOT.plot(axes=axes,fig=fig,ax=ax1,LCFS=LCFS,limiters=limiters,real_scale=real_scale,vminmax=vminmax,fill=fill,colmap=cmap_ASCOT,number_bins=number_bins)

        if True and len(axes)>1:
            cbar=fig.colorbar(mesh,ax=ax1,orientation='vertical')

    plt.show()

#################################
 
##################################################################
 
###################################################################################################
