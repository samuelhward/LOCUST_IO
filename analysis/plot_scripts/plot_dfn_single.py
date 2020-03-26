#plot_dfn_single.py
 
"""
Samuel Ward
18/04/2019
----
script to plot time slices for MAST V and V - single plot
---
usage:
 
notes:         
    can plot arbitrary number of DFNs using this script
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
import matplotlib

def cmap_custom(from_rgb,to_rgb):
    """
    generate custom colormaps
    args:
        from_rgb - list of starting r g b values  
        to_rgb - list of final r g b values
    notes:
    """
    from matplotlib.colors import LinearSegmentedColormap
    r1,g1,b1=from_rgb
    r2,g2,b2=to_rgb
    cdict={'red':((0,r1,r1),(1,r2,r2)),
            'green':((0,g1,g1),(1,g1,g1)),
            'blue':((0,b1,b1),(1,b1,b1))}
    cmap=LinearSegmentedColormap('custom cmap - {from_rgb}/{to_rgb}'.format(from_rgb=str(from_rgb),to_rgb=str(to_rgb)),cdict)
    return cmap
cmap_default=matplotlib.cm.get_cmap('jet') #set default colourmap
cmap_r=cmap_custom([1,0,0],[1,0,0]) #red
cmap_g=cmap_custom([0,1,0],[0,1,0]) #green
cmap_b=cmap_custom([0,0,1],[0,0,1]) #blue
cmap_y=cmap_custom([1,1,0],[1,1,0]) #yellow
cmap_m=cmap_custom([1,0,1],[1,0,1]) #magenta
cmap_c=cmap_custom([0,1,1],[0,1,1]) #cyan
cmap_w=cmap_custom([1,1,1],[1,1,1]) #white
cmap_k=cmap_custom([0,0,0],[0,0,0]) #black

cmap_LOCUST=settings.cmap_g
cmap_ASCOT=settings.cmap_b
cmap_TRANSP=settings.cmap_r

shot='29034'
run='W04'
tinit=0.36 #start time
axes=['E']
real_scale=False
limiters=False#MAST_wall
LCFS=False#MAST_eq
fill=False
number_bins=10
vminmax=False#[0.1e15,6.e15]
legend=True
title='1W BBNBI injection F(E) at 100ms - MAST #29034'
xlabel='energy [keV]'
ylabel='density [#/eV]'

file_eq=pathlib.Path('TRANSP') / (shot+run) / ('g'+shot) #grab the wall and equilibrium used in these simulations
file_wall=pathlib.Path('TRANSP') / (shot+run) / ('OMF'+shot+'.LIM')
MAST_eq=equi('',data_format='GEQDSK',filename=file_eq)
MAST_wall=wall('',data_format='UFILE',filename=file_wall)

dimensions=[] #intialise lists holding run metadata
codes=[]
extra_infos=[]
data_formats=[]
cmaps=[]
ASCOT_beam_depo_data_formats=[]
ASCOT_beam_depo_filenames=[]
IDs=[] #some more information to put as plot legend or misc

################################################################## W03

time_slices=np.linspace(0.37,0.46,10)

################################# BBNBI depositions
'''
dimensions.append('2D') #append a run's metadata
codes.append('ASCOT')
extra_infos.append('BBNBI_birth_FO')
data_formats.append('ASCOT')
cmaps.append(settings.cmap_b)
ASCOT_beam_depo_data_formats.append('ASCOT_FO') #if this is an ASCOT run, append the beam deposition data format so we can rescale according to true power deposited
ASCOT_beam_depo_filenames.append(pathlib.Path('ASCOT') / (shot+run) / 'input.particles_BBNBI_FO_reduced')
IDs.append('ASCOT FO')

dimensions.append('2D') #append a run's metadata
codes.append('ASCOT')
extra_infos.append('BBNBI_birth_GC')
data_formats.append('ASCOT')
cmaps.append(settings.cmap_r)
ASCOT_beam_depo_data_formats.append('ASCOT_FO')
ASCOT_beam_depo_filenames.append(pathlib.Path('ASCOT') / (shot+run) / 'input.particles_BBNBI_FO_reduced')
IDs.append('ASCOT GC')


dimensions.append('2D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('BBNBI_birth_GC_noDCTDC2')
data_formats.append('LOCUST')
cmaps.append(cmap_custom([0,0.7,0],[0,0.7,0]))
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST GC')


dimensions.append('2D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('BBNBI_birth_FO')
data_formats.append('LOCUST')
cmaps.append(settings.cmap_k)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST FO truncated')
'''

################################# NUBEAM depositions
'''
dimensions.append('2D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('NUBEAM_birth_GC')
data_formats.append('LOCUST')
cmaps.append(settings.cmap_k)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST GC truncated') 

dimensions.append('2D') #append a run's metadata
codes.append('ASCOT')
extra_infos.append('NUBEAM_birth_GC')
data_formats.append('ASCOT')
cmaps.append(cmap_ASCOT)
ASCOT_beam_depo_data_formats.append('ASCOT_GC')
ASCOT_beam_depo_filenames.append(pathlib.Path('ASCOT') / (shot+run) / 'input.particles_NUBEAM_birth_GC')
IDs.append('ASCOT GC')

dimensions.append('2D') #append a run's metadata
codes.append('TRANSP')
extra_infos.append('NUBEAM_birth_GC')
data_formats.append(None) #TRANSP distribution functions are unique here
cmaps.append(cmap_TRANSP)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('TRANSP GC')

dimensions.append('2D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('NUBEAM_birth_GC_noDCTDC2')
data_formats.append('LOCUST')
cmaps.append(cmap_custom([0,0.7,0],[0,0.7,0]))
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST GC')

'''

################################################################## W04
#'''

time_slices=np.array([1,2,3,4,5,10,20,40,70,100])

dimensions.append('3D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('NUBEAM_birth_GC')
data_formats.append('LOCUST')
cmaps.append(settings.cmap_k)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST GC truncated') 

dimensions.append('2D') #append a run's metadata
codes.append('ASCOT')
extra_infos.append('NUBEAM_birth_GC')
data_formats.append('ASCOT')
cmaps.append(cmap_ASCOT)
ASCOT_beam_depo_data_formats.append('ASCOT_GC')
ASCOT_beam_depo_filenames.append(pathlib.Path('ASCOT') / (shot+run) / 'input.particles_NUBEAM_birth_GC')
IDs.append('ASCOT GC')

#'''
dimensions.append('2D') #append a run's metadata
codes.append('TRANSP')
extra_infos.append('NUBEAM_birth_GC')
data_formats.append(None) #TRANSP distribution functions are unique here
cmaps.append(cmap_TRANSP)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('TRANSP GC')
#'''

#'''


'''
dimensions.append() #append a run's metadata
codes.append()
extra_infos.append()
data_formats.append()
cmaps.append()
ASCOT_beam_depo_data_formats.append()
ASCOT_beam_depo_filenames.append()
IDs.append()
'''
        



from get_filepaths_single import get_filepaths_single #use external script for globbing all the filenames
files=[] #list of lists, each list is a set of timeslices for a given run
for dimension,extra_info,code in zip(dimensions,extra_infos,codes):
    files.append(get_filepaths_single(shot=shot,run=run,dimension=dimension,extra_info=extra_info,code=code))


print(files)

#iterate over time slices in top level
for counter_time_slice,time_slice in enumerate([time_slices[-1]]): #XXX hack to plot final time slice
    fig,(ax1)=plt.subplots(1)
    #for this time slice, iterate over each code - for each code we need to pick out the simulation corresponding to the current time slice
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
        #ax1.set_xlim(0,100)
        #ax1.set_ylim(-1.1,1.1)
        mesh=DFN.plot(axes=axes,fig=fig,ax=ax1,LCFS=LCFS,limiters=limiters,real_scale=real_scale,fill=fill,colmap=cmap,number_bins=number_bins,vminmax=vminmax)

        if True and len(axes)>1:
            cbar=fig.colorbar(mesh,ax=ax1,orientation='vertical')
    
    if legend:
        plt.legend(([ID for ID in IDs]))
    
    if title:
        ax1.set_title(title)
    if xlabel:
        ax1.set_xlabel(xlabel)
    if ylabel:
        ax1.set_ylabel(ylabel)
    plt.show()

#################################
 
##################################################################
 
###################################################################################################
