#plot_dfn_single.py
 
"""
Samuel Ward
18/04/2019
----
script to plot single time slice for MAST V and V
---
usage:
 
notes:         
    can plot arbitrary number of DFNs using this script

    rules for line_styles:
        FO = solid cmap_G
        add GC --> add dashes
        lnLA --> make cmap_G_
        trunc --> add dots
---
"""

import numpy as np
import context 
import copy
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
from get_filepaths_single import get_filepaths_single #use external script for globbing all the filenames
import processing.utils

#define some colourmaps
cmap_r=settings.colour_custom([194,24,91,1])
cmap_g=settings.colour_custom([76,175,80,1])
cmap_g_=settings.colour_custom([205,220,57,1])
cmap_b=settings.colour_custom([33,150,243,1])
cmap_grey=settings.colour_custom([97,97,97,1])

tinit=0.36 #start time
axes=['E']
real_scale=False
limiters=False#MAST_wall
LCFS=False#MAST_eq
fill=False
number_bins=10
vminmax=False#[0.1e15,6.e15]
legend=True
legend_fontsize=15
linewidth=3
title=''
xlabel=''
ylabel=''
residual=False #toggle plotting f(E) or f(E)-mean{f(E)}
normalise_inject=False
normalise=False

def plot_dfns(axes=['E'],normalise_inject=False,normalise=True,E_inject=62000,time_slice_index=-1,residual=False,return_dfns=True,same_base=True):
    """
    plots the dfns but uses lots of global variables as free variables e.g. figs,ax etc. defined outside this function scope

    args:
        axes - 
        normalise_inject - toggle whether to normalise all dfns to injection energy
        normalise - toggle whether to normalise all dfns compared to average
        E_inject - specify injection energy
        time_slice_index - index of time slice to plot
        residual - toggle plotting DFNs as residuals of average
        return_dfns - toggle whether to return dict of integrated and plotted dfns
        same_base - interpolate/plot DFNs on same axis
    notes:
    """

    if residual or normalise or same_base: average_DFN,average_axes=dfn_find_average(axes=axes,normalise_inject=normalise_inject)

    if return_dfns: DFNs={}

    #plot single time slice in top level
    for time_slice in [time_slices[time_slice_index]]: #XXX hack to plot final time slice
        
        #for this time slice, iterate over each code - for each code we need to pick out the simulation corresponding to the current time slice
        for counter_code,(time_slices_for_this_run,code,data_format,ID,line_style,cmap) in enumerate(zip(files,codes,data_formats,IDs,line_styles,cmaps)):
            file_this_time_slice=time_slices_for_this_run[time_slice_index]
            if file_this_time_slice is None: continue

            #read DFNs
            if code=='TRANSP':
                DFN=dfn_t(ID=ID,filename=file_this_time_slice)
                #N=DFN.dfn_integrate()['dfn']
                #scale DFNs to match 1W beam deposited power
                file_transp_output=pathlib.Path(file_this_time_slice).parents[0] / (shot+run+'.CDF')
                output_TRANSP=output_T('',filename=file_transp_output)
                BPCAP_index=np.abs(output_TRANSP['TIME']-time_slice).argmin()
                BPCAP=output_TRANSP['BPCAP'][BPCAP_index]
                DFN['dfn']/=BPCAP

                if normalise_inject:
                    inject_index=np.abs(DFN['E']-E_inject).argmin()
                    DFN['dfn']/=DFN.dfn_integrate(energy=False)['dfn'][inject_index]
            else:
                DFN=dfn(ID=ID,time_string='time = {}ms'.format(int(1000*(time_slice-tinit))),data_format=data_format,filename=file_this_time_slice)
                #N=DFN.transform(axes=['N'])['dfn']
                
                if code=='ASCOT': #grab the ASCOT beam deposition froom the inistate to renormalise_inject the ASCOT dfn to the correct power deposition

                    beam_power=1. #desired beam power
                    ASCOT_beam_depo_data_format=ASCOT_beam_depo_data_formats[counter_code]
                    ASCOT_beam_depo_filename=ASCOT_beam_depo_filenames[counter_code]
                    ASCOT_beam_depo=bd(ID='ASCOT beam depo',data_format=ASCOT_beam_depo_data_format,filename=ASCOT_beam_depo_filename)
                    if 'E' not in ASCOT_beam_depo.data: 
                        ASCOT_beam_depo['E']=0.5*constants.mass_deuteron*(ASCOT_beam_depo['V_R']**2+ASCOT_beam_depo['V_phi']**2+ASCOT_beam_depo['V_Z']**2)/constants.species_charge
                    Pdep_ASCOT=np.sum(ASCOT_beam_depo['E']*constants.species_charge*ASCOT_beam_depo['weight'])
                    DFN['dfn']*=beam_power/Pdep_ASCOT 
                    DFN['V_pitch']*=-1

                    #calculate % power loss
                    ffppll=fpl('time = {}ms'.format(int(1000*(time_slice-tinit))),data_format='ASCOT',filename=file_this_time_slice)
                    PFC_hits=np.where(ffppll['status_flag']==ffppll['status_flags']['wall_collision'])[0]
                    print(f"PFC flux {file_this_time_slice}={100.*np.sum(ffppll['E'][PFC_hits]*ffppll['weight'][PFC_hits])*constants.charge_e*beam_power/(Pdep_ASCOT*time_slice)}%")

                if normalise_inject:
                    inject_index=np.abs(DFN['E']-E_inject).argmin()
                    DFN['dfn']/=DFN.transform(axes=['E'])['dfn'][inject_index]

            transform=True #have we pre-integrated the distribution function? not yet

            if residual or same_base: 
                if code=='TRANSP':
                    DFN=DFN.dfn_integrate(energy=False)
                else: 
                    DFN=DFN.transform(axes=axes)
                transform=False
    
                if len(axes)==1: interpolator=processing.utils.interpolate_1D(*[DFN[axis] for axis in axes],DFN['dfn'],type='RBF',function='linear')
                elif len(axes)==2: interpolator=processing.utils.interpolate_2D(*[DFN[axis] for axis in axes],DFN['dfn'],type='RBF',function='linear')
                else:
                    print("ERROR: axes must be 1D or 2D!\nreturning\n")
                    return
                DFN['dfn']=interpolator(*average_axes)
                for axis,average_axis in zip(axes,average_axes):
                    DFN[axis]=copy.deepcopy(average_axis)
                if residual:
                    DFN['dfn']=DFN['dfn']-average_DFN
        
            DFN['E']/=1000. #get axes in keV - .plot integrates the DFNs using .transform with 'dE' so should not affect plotting
            if normalise: DFN['dfn']/=average_DFN #measure DFN against un-normalised average
            mesh=DFN.plot(axes=axes,fig=fig,ax=ax,LCFS=LCFS,limiters=limiters,real_scale=real_scale,fill=fill,colmap=cmap,number_bins=number_bins,vminmax=vminmax,transform=transform,label=DFN.ID,line_style=line_style)

            if True and len(axes)>1:
                cbar=fig.colorbar(mesh,ax=ax,orientation='vertical')

            if return_dfns: DFNs[ID]=DFN['dfn']

        if residual and len(axes)==1: #plot 0 mark and dfn average
            ax.axhline(0.,color='k',linestyle='-',label='$f_{\mathrm{mean}}$')
            ax_av=ax.twinx()
            if normalise: average_DFN/=np.max(average_DFN)
            ax_av.plot(*[average_ax/1000 for average_ax in average_axes],average_DFN,'k-',label='average $f$')
            ax_av.tick_params(axis='y', labelcolor='k')

            label_mean='$f_{\mathrm{mean}}$'
            if normalise: label_mean+=' [a.u.]'
            ax_av.set_ylabel(label_mean,color='k') 
            ax_av.spines['right'].set_position(('outward', 0))
            ax.set_ylim(-0.15,0.15)

        ax.legend(fontsize=legend_fontsize,loc='lower center',ncol=2)
        ax.set_xlim(0,67.)
        
        if title:
            ax.set_title(title,fontsize=25)
        if xlabel:
            ax.set_xlabel(xlabel,fontsize=25)
        if ylabel:
            ax.set_ylabel(ylabel,fontsize=25)

    if return_dfns: return DFNs

def dfn_find_average(axes=['E'],normalise_inject=False,E_inject=62000,time_slice_index=-1):
    
    DFNs=[]
    DFNs_axes=[]

    for time_slice in [time_slices[time_slice_index]]: #XXX hack to plot final time slice

        #for this time slice, iterate over each code - for each code we need to pick out the simulation corresponding to the current time slice
        for counter_code,(time_slices_for_this_run,code,data_format,ID,cmap) in enumerate(zip(files,codes,data_formats,IDs,cmaps)):
            file_this_time_slice=time_slices_for_this_run[time_slice_index]
            if file_this_time_slice is None: continue

            #read DFNs
            if code=='TRANSP':
                DFN=dfn_t(ID=ID,filename=file_this_time_slice)
                #scale DFNs to match 1W beam deposited power
                file_transp_output=pathlib.Path(file_this_time_slice).parents[0] / (shot+run+'.CDF')
                output_TRANSP=output_T('',filename=file_transp_output)
                BPCAP_index=np.abs(output_TRANSP['TIME']-time_slice).argmin()
                BPCAP=output_TRANSP['BPCAP'][BPCAP_index]
                DFN['dfn']/=BPCAP

                if normalise_inject:
                    inject_index=np.abs(DFN['E']-E_inject).argmin()
                    DFN['dfn']/=DFN.dfn_integrate(energy=False)['dfn'][inject_index]
                
                DFN=DFN.dfn_integrate(energy=False)

            else:
                DFN=dfn(ID=ID,time_string='time = {}ms'.format(int(1000*(time_slice-tinit))),data_format=data_format,filename=file_this_time_slice)
                
                if code=='ASCOT': #grab the ASCOT beam deposition froom the inistate to renormalise_inject the ASCOT dfn to the correct power deposition

                    beam_power=1. #desired beam power
                    ASCOT_beam_depo_data_format=ASCOT_beam_depo_data_formats[counter_code]
                    ASCOT_beam_depo_filename=ASCOT_beam_depo_filenames[counter_code]
                    ASCOT_beam_depo=bd(ID='ASCOT beam depo',data_format=ASCOT_beam_depo_data_format,filename=ASCOT_beam_depo_filename)
                    if 'E' not in ASCOT_beam_depo.data: 
                        ASCOT_beam_depo['E']=0.5*constants.mass_deuteron*(ASCOT_beam_depo['V_R']**2+ASCOT_beam_depo['V_phi']**2+ASCOT_beam_depo['V_Z']**2)/constants.species_charge
                    Pdep_ASCOT=np.sum(ASCOT_beam_depo['E']*constants.species_charge*ASCOT_beam_depo['weight'])
                    DFN['dfn']*=beam_power/Pdep_ASCOT 
                    DFN['V_pitch']*=-1

                if normalise_inject:
                    inject_index=np.abs(DFN['E']-E_inject).argmin()
                    DFN['dfn']/=DFN.transform(axes=['E'])['dfn'][inject_index]

                DFN=DFN.transform(axes=axes)

            DFNs.append(DFN)
            DFNs_axes.append([DFN[axis] for axis in axes])
 
        #find narrowest range to interpolate over out of all axes
        axes_maxes=[]
        axes_mins=[]
        for DFNs_axis in DFNs_axes:
            axes_maxes.append([np.max(axis) for axis in DFNs_axis])
            axes_mins.append([np.min(axis) for axis in DFNs_axis])
        maxes=np.min(np.array(axes_maxes,ndmin=len(axes)),axis=0)
        mins=np.max(np.array(axes_mins,ndmin=len(axes)),axis=0)
        number_bins=100
        average_axes=[np.linspace(min_,max_,number_bins) for min_,max_ in zip(mins,maxes)]
        DFN_average=np.zeros(shape=[number_bins]*len(axes)) #store n-dimensional dfn average

        for DFN in DFNs:
            if len(axes)==1:
                DFN_interpolator=processing.utils.interpolate_1D(*[DFN[axis] for axis in axes],DFN['dfn'],type='RBF',function='linear')
                DFN_average+=DFN_interpolator(*[average_axis for average_axis in average_axes])/len(DFNs)   
            elif len(axes)==2: 
                DFN_interpolator=processing.utils.interpolate_2D(*[DFN[axis] for axis in axes],DFN['dfn'],type='RBF',function='linear')
                DFN_average+=DFN_interpolator(np.meshgrid(*[average_axis for average_axis in average_axes]))/len(DFNs)   
            else:
                print("ERROR: axes must be 1D or 2D!\nreturning\n")
                return
    
    return DFN_average,average_axes

################################################################## W03

dimensions=[] #intialise lists holding run metadata
codes=[]
extra_infos=[]
data_formats=[]
cmaps=[]
ASCOT_beam_depo_data_formats=[]
ASCOT_beam_depo_filenames=[]
IDs=[] #some more information to put as plot legend or misc
line_styles=[]
time_slices=np.linspace(0.37,0.46,10)
shot='29034'
run='W03'

################################# BBNBI depositions

#'''

dimensions.append('2D') #append a run's metadata
codes.append('ASCOT')
extra_infos.append('BBNBI_birth_GC')
data_formats.append('ASCOT')
cmaps.append(cmap_b)
ASCOT_beam_depo_data_formats.append('ASCOT_FO')
ASCOT_beam_depo_filenames.append(pathlib.Path('ASCOT') / (shot+run) / 'input.particles_BBNBI_FO_reduced')
IDs.append('ASCOT GC')
line_styles.append('dashed')

dimensions.append('2D') #append a run's metadata
codes.append('ASCOT')
extra_infos.append('BBNBI_birth_FO')
data_formats.append('ASCOT')
cmaps.append(cmap_b)
ASCOT_beam_depo_data_formats.append('ASCOT_FO') #if this is an ASCOT run, append the beam deposition data format so we can rescale according to true power deposited
ASCOT_beam_depo_filenames.append(pathlib.Path('ASCOT') / (shot+run) / 'input.particles_BBNBI_FO_reduced')
IDs.append('ASCOT FO')
line_styles.append('solid')

dimensions.append('2D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('BBNBI_birth_GC')
data_formats.append('LOCUST')
cmaps.append(cmap_g)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST GC')
line_styles.append('dashed')

dimensions.append('2D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('BBNBI_birth_GC_trunc')
data_formats.append('LOCUST')
cmaps.append(cmap_g)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST GC truncated')
line_styles.append('dashdot')

dimensions.append('3D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('BBNBI_birth_GC_ASCOTlnL')
data_formats.append('LOCUST')
cmaps.append(cmap_g_)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST GC ASCOT ln($\Lambda$)')
line_styles.append('dashed')

dimensions.append('3D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('BBNBI_birth_GC_ASCOTlnL_trunc')
data_formats.append('LOCUST')
cmaps.append(cmap_g_)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST GC ASCOT ln($\Lambda$) truncated')
line_styles.append('dashdot')

dimensions.append('3D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('BBNBI_birth_FO')
data_formats.append('LOCUST')
cmaps.append(cmap_g)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST FO')
line_styles.append('solid')

dimensions.append('2D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('BBNBI_birth_FO_trunc')
data_formats.append('LOCUST')
cmaps.append(cmap_g)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST FO truncated')
line_styles.append('dotted')

dimensions.append('3D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('BBNBI_birth_FO_ASCOTlnL')
data_formats.append('LOCUST')
cmaps.append(cmap_g_)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST FO ASCOT ln($\Lambda$)')
line_styles.append('solid')

dimensions.append('3D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('BBNBI_birth_FO_ASCOTlnL_trunc')
data_formats.append('LOCUST')
cmaps.append(cmap_g_)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST FO ASCOT ln($\Lambda$) truncated')
line_styles.append('dotted')

#'''

################################# NUBEAM depositions (monoenergetic)
'''
dimensions.append('2D') #append a run's metadata
codes.append('ASCOT')
extra_infos.append('NUBEAM_birth_GC')
data_formats.append('ASCOT')
cmaps.append(cmap_b)
ASCOT_beam_depo_data_formats.append('ASCOT_GC')
ASCOT_beam_depo_filenames.append(pathlib.Path('ASCOT') / (shot+run) / 'input.particles_NUBEAM_birth_GC')
IDs.append('ASCOT GC')
line_styles.append('dashed')

dimensions.append('2D') #append a run's metadata
codes.append('TRANSP')
extra_infos.append('NUBEAM_birth_GC_trunc')
data_formats.append(None) #TRANSP distribution functions are unique here
cmaps.append(cmap_r)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('TRANSP GC truncated')
line_styles.append('dashdot')

dimensions.append('2D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('NUBEAM_birth_GC_trunc')
data_formats.append('LOCUST')
cmaps.append(cmap_g)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST GC truncated') 
line_styles.append('dashdot')

dimensions.append('2D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('NUBEAM_birth_GC')
data_formats.append('LOCUST')
cmaps.append(cmap_g)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST GC')
line_styles.append('dashed')
'''

file_eq=pathlib.Path('TRANSP') / (shot+run) / ('g'+shot) #grab the wall and equilibrium used in these simulations
file_wall=pathlib.Path('TRANSP') / (shot+run) / ('OMF'+shot+'.LIM')
MAST_eq=equi('',data_format='GEQDSK',filename=file_eq)
MAST_wall=wall('',data_format='UFILE',filename=file_wall)

files=[] #list of lists, each list is a set of timeslices for a given run
for dimension,extra_info,code in zip(dimensions,extra_infos,codes):
    files.append(get_filepaths_single(shot=shot,run=run,dimension=dimension,extra_info=extra_info,code=code))

fig,(ax)=plt.subplots(1)
title=f'BBNBI injection $f(E)$ at 100ms - MAST #{shot}{run}' if normalise else f'1W BBNBI injection $F(E)$ at 100ms - MAST #{shot}{run}'
xlabel='Energy [keV]'
if residual: ylabel+=' residual'
if normalise: ylabel+= ' [a.u.]' 
else: ylabel += ' [#/eV]'
W03_DFNS=plot_dfns(axes=axes,normalise_inject=normalise_inject,normalise=normalise,residual=residual)
plt.show()

################################################################## W04
#'''

dimensions=[] #intialise lists holding run metadata
codes=[]
extra_infos=[]
data_formats=[]
cmaps=[]
ASCOT_beam_depo_data_formats=[]
ASCOT_beam_depo_filenames=[]
IDs=[] #some more information to put as plot legend or misc
line_styles=[]
time_slices=np.array([1,2,3,4,5,10,20,40,70,100])
shot='29034'
run='W04'

dimensions.append('2D') #append a run's metadata
codes.append('ASCOT')
extra_infos.append('NUBEAM_birth_GC')
data_formats.append('ASCOT')
cmaps.append(cmap_b)
ASCOT_beam_depo_data_formats.append('ASCOT_GC')
ASCOT_beam_depo_filenames.append(pathlib.Path('ASCOT') / (shot+run) / 'input.particles_NUBEAM_birth_GC')
IDs.append('ASCOT GC')
line_styles.append('dashed')

dimensions.append('2D') #append a run's metadata
codes.append('TRANSP')
extra_infos.append('NUBEAM_birth_GC')
data_formats.append(None) #TRANSP distribution functions are unique here
cmaps.append(cmap_r)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('TRANSP GC')
line_styles.append('dashed')

dimensions.append('3D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('NUBEAM_birth_GC')
data_formats.append('LOCUST')
cmaps.append(cmap_g)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST GC')
line_styles.append('dashed')

dimensions.append('3D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('NUBEAM_birth_GC_trunc')
data_formats.append('LOCUST')
cmaps.append(cmap_g)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST GC truncated') 
line_styles.append('dashdot')

dimensions.append('3D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('NUBEAM_birth_GC_ASCOTlnL')
data_formats.append('LOCUST')
cmaps.append(cmap_g_)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST GC ASCOT ln($\Lambda$)')
line_styles.append('dashed')

dimensions.append('3D') #append a run's metadata
codes.append('LOCUST')
extra_infos.append('NUBEAM_birth_GC_ASCOTlnL_trunc')
data_formats.append('LOCUST')
cmaps.append(cmap_g_)
ASCOT_beam_depo_data_formats.append(None)
ASCOT_beam_depo_filenames.append(None)
IDs.append('LOCUST GC ASCOT ln($\Lambda$) truncated')
line_styles.append('dashdot')

#'''

file_eq=pathlib.Path('TRANSP') / (shot+run) / ('g'+shot) #grab the wall and equilibrium used in these simulations
file_wall=pathlib.Path('TRANSP') / (shot+run) / ('OMF'+shot+'.LIM')
MAST_eq=equi('',data_format='GEQDSK',filename=file_eq)
MAST_wall=wall('',data_format='UFILE',filename=file_wall)

files=[] #list of lists, each list is a set of timeslices for a given run
for dimension,extra_info,code in zip(dimensions,extra_infos,codes):
    files.append(get_filepaths_single(shot=shot,run=run,dimension=dimension,extra_info=extra_info,code=code))

fig,(ax)=plt.subplots(1)
title=f'NUBEAM injection $F(E)$ at 100ms - MAST #{shot}{run}' if normalise else f'1W NUBEAM injection $F(E)$ at 100ms - MAST #{shot}{run}'
xlabel='Energy [keV]'
ylabel='Density'

if residual: ylabel+=' residual'
if normalise: ylabel+= ' [a.u.]' 
else: ylabel += ' [#/eV]'

W04_DFNS=plot_dfns(axes=axes,normalise_inject=normalise_inject,normalise=normalise,residual=residual)
plt.show()

##################################################################

'''
dimensions.append() #append a run's metadata
codes.append()
extra_infos.append()
data_formats.append()
cmaps.append()
ASCOT_beam_depo_data_formats.append()
ASCOT_beam_depo_filenames.append()
IDs.append()
line_styles.append()
'''

##################################################################
#perform KS tests

#requires residual=False,normalise=False,normalise_inject=False

KS_pairs_W03=[
            ['LOCUST GC','LOCUST FO'],
            ['LOCUST GC','ASCOT GC'],
            ['LOCUST FO','ASCOT FO'],
            ['LOCUST FO ASCOT ln($\Lambda$)','ASCOT FO'],
            ['LOCUST GC ASCOT ln($\Lambda$)','ASCOT GC'],
            ['LOCUST FO ASCOT ln($\Lambda$) truncated','ASCOT FO'], #check KS test with this
            ['LOCUST GC ASCOT ln($\Lambda$) truncated','ASCOT GC'], #check KS test with this
            ]

KS_pairs_W04=[
            ['LOCUST GC','ASCOT GC'],
            ['LOCUST GC ASCOT ln($\Lambda$)','ASCOT GC'],
            ['LOCUST GC','TRANSP GC'],
            ['LOCUST GC ASCOT ln($\Lambda$) truncated','ASCOT GC'], #check KS test with this
            ]

'''
number_markers={}
number_markers['W03']={}
number_markers['W04']={}
number_markers['W03']['ASCOT']=20000    
number_markers['W04']['ASCOT']=70507

W03 

#do not care about NUBEAM_birth* since single energy 
BBNBI_birth_FO=32768 
BBNBI_birth_FO_ASCOTlnL=32768 
BBNBI_birth_FO_ASCOTlnL_trunc=65536 
BBNBI_birth_FO_trunc=4194304 
BBNBI_birth_GC=524288 
BBNBI_birth_GC_ASCOTlnL=32768 
BBNBI_birth_GC_ASCOTlnL_trunc=65536 
BBNBI_birth_GC_trunc=1048576  

W04

NUBEAM_birth_GC=65536 
NUBEAM_birth_GC_ASCOTlnL=32768 
NUBEAM_birth_GC_ASCOTlnL_trunc=65536 
NUBEAM_birth_GC_trunc=65536 
'''

for alpha in [0.1,0.05,0.01,0.001]:
    
    print(f'********alpha={alpha}********\n')
    
    for run_dfns,run_pairs,run_name, in zip([W03_DFNS,W04_DFNS],[KS_pairs_W03,KS_pairs_W04],['W03','W04']):
    
        print(f'****{run_name} KS tests****')
        for pair in run_pairs:
            print(pair)
            D,P,reject=processing.utils.KS_test(70507,70507,run_dfns[pair[0]],run_dfns[pair[1]],alpha=alpha)
            print(f'D={D}')
            print(f'P={P}')
            print(f'reject? = {reject}')
            print('\n')

#################################

##################################################################

###################################################################################################
