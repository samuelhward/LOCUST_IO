#plot_scan_DNB_bounds.py
 
"""
Samuel Ward
06/11/21
----
script for plotting scan_DNB
---
 
notes:
remember to comment out any continue blocks in the launch script
---
"""

###################################################################################################
#Preamble
 
import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway

try:
    import numpy as np
    import pathlib
    import copy
    import matplotlib.pyplot as plt
    import os
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

if __name__=='__main__':

    try:
        import context
    except:
        raise ImportError("ERROR: context.py could not be imported!\nreturning\n")
        sys.exit(1)
try:
    from classes.input_classes.equilibrium import Equilibrium
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/equilibrium.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.beam_deposition import Beam_Deposition
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/perturbation.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/processing/utils.py could not be imported!\nreturning\n") 
    sys.exit(1)

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/src/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import constants
except:
    raise ImportError("ERROR: LOCUST_IO/src/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import settings
except:
    raise ImportError("ERROR: LOCUST_IO/src/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)

try:
    cwd=pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    sys.path.append(str(cwd.parents[1]))
    import templates.plot_mod
except:
    raise ImportError("ERROR: templates/template_mod.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## 
#Main 

cmap_r=settings.colour_custom([194,24,91,1])
cmap_g=settings.colour_custom([76,175,80,1])
cmap_b=settings.colour_custom([33,150,243,1])

colours=plt.rcParams['axes.prop_cycle'].by_key()['color']

import scan_DNB_launch as batch_data

output_dfns=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'dfn',batch_data,processes=16,chunksize=1)).reshape(
    len(batch_data.parameters__databases),
    len(batch_data.parameters__toroidal_mode_numbers),
    len(batch_data.parameters__phases_uppers_case7[0]),
    len(batch_data.parameters__currents_upper),
    )

output_moms=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'mom',batch_data,processes=32,chunksize=1)).reshape(
    len(batch_data.parameters__databases),
    len(batch_data.parameters__toroidal_mode_numbers),
    len(batch_data.parameters__phases_uppers_case7[0]),
    len(batch_data.parameters__currents_upper),
    )

fig=plt.figure(figsize=(10, 15),constrained_layout=False)
fig_mom=plt.figure(figsize=(10, 15),constrained_layout=False)
fig_2D_n3=plt.figure(figsize=(10, 15),constrained_layout=False)
fig_2D_n4=plt.figure(figsize=(10, 15),constrained_layout=False)
dimensions=(len(batch_data.parameters__databases),len(batch_data.parameters__currents_uppers[0])-1)
ax_total=fig.add_subplot(1,1,1,frameon=False)
ax_total.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_total.spines['top'].set_color('none')
ax_total.spines['bottom'].set_color('none')
ax_total.spines['left'].set_color('none')
ax_total.spines['right'].set_color('none')
ax_total_mom=fig_mom.add_subplot(1,1,1,frameon=False)
ax_total_mom.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_total_mom.spines['top'].set_color('none')
ax_total_mom.spines['bottom'].set_color('none')
ax_total_mom.spines['left'].set_color('none')
ax_total_mom.spines['right'].set_color('none')
ax_total_2D_n3=fig_2D_n3.add_subplot(1,1,1,frameon=False)
ax_total_2D_n3.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_total_2D_n3.spines['top'].set_color('none')
ax_total_2D_n3.spines['bottom'].set_color('none')
ax_total_2D_n3.spines['left'].set_color('none')
ax_total_2D_n3.spines['right'].set_color('none')
ax_total_2D_n4=fig_2D_n4.add_subplot(1,1,1,frameon=False)
ax_total_2D_n4.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_total_2D_n4.spines['top'].set_color('none')
ax_total_2D_n4.spines['bottom'].set_color('none')
ax_total_2D_n4.spines['left'].set_color('none')
ax_total_2D_n4.spines['right'].set_color('none')

axs=[]
axs_mom=[]
axs_2D_n3=[]
axs_2D_n4=[]
_=1
for dim in dimensions: _*=dim
for _ in range(_):
    axs.append(fig.add_subplot(*dimensions,_+1))
    axs_mom.append(fig_mom.add_subplot(*dimensions,_+1))
    axs_2D_n3.append(fig_2D_n3.add_subplot(*dimensions,_+1))
    axs_2D_n4.append(fig_2D_n4.add_subplot(*dimensions,_+1))

axs=np.array(axs).T.reshape(dimensions)
axs_mom=np.array(axs_mom).T.reshape(dimensions)
axs_2D_n3=np.array(axs_2D_n3).T.reshape(dimensions)
axs_2D_n4=np.array(axs_2D_n4).T.reshape(dimensions)

plot_axes=['E','V_pitch']
plot_axes=['R','Z']
plot_axes=['R']
plot_axes=['E']
plot_axes=['V_pitch']
plot_axes_2D=['E','V_pitch']
moment_quantity='NBI-heating-power(TOT)'
moment_quantity='energy_para'
moment_quantity='density'
moment_axis='flux_pol_norm'
plot_mode=2

for database_counter,(
        parameters__database,parameters__sheet_name_kinetic_prof,
        parameters__phases_uppers,parameters__phases_middles,parameters__phases_lowers,
        parameters__currents_upper,parameters__currents_middle,parameters__currents_lower 
        ) in enumerate(zip(
        batch_data.parameters__databases,batch_data.parameters__sheet_names_kinetic_prof,
        batch_data.parameters__phases_uppers_cases_all,batch_data.parameters__phases_middles_cases_all,batch_data.parameters__phases_lowers_cases_all,
        batch_data.parameters__currents_uppers,batch_data.parameters__currents_middles,batch_data.parameters__currents_lowers
        )): 
    #try:
    DFN_2D=output_dfns[database_counter,0,0,0].transform(axes=plot_axes)
    DFN_2D['E']/=1.e3 #plot in keV
    mom_2D=output_moms[database_counter,0,0,0]
    for parameters__kinetic_prof_tF_tE,parameters__kinetic_prof_Pr in zip(batch_data.parameters__kinetic_profs_tF_tE,batch_data.parameters__kinetic_profs_Pr):
        for mode_counter,(parameters__toroidal_mode_number,parameters__phases_upper,parameters__phases_middle,parameters__phases_lower) in enumerate(zip(batch_data.parameters__toroidal_mode_numbers,parameters__phases_uppers,parameters__phases_middles,parameters__phases_lowers)):
            for parameters__rotation_upper,parameters__rotation_middle,parameters__rotation_lower in zip(batch_data.parameters__rotations_upper,batch_data.parameters__rotations_middle,batch_data.parameters__rotations_lower): 
                for current_counter,(parameters__current_upper,parameters__current_middle,parameters__current_lower) in enumerate(zip(
                    parameters__currents_upper,
                    parameters__currents_middle,
                    parameters__currents_lower)):                    
                    if parameters__current_upper==0: continue
                    axs[database_counter,current_counter-1].set_xlabel('')
                    axs[database_counter,current_counter-1].set_ylabel('')
                    #V_pitch
                    #axs[database_counter,current_counter-1].set_xticks([-1.,-0.5,0.,0.5,1.])
                    #axs[database_counter,current_counter-1].set_xlim([-1,1])
                    #axs[database_counter,current_counter-1].set_ylim([-0.5,0.7])
                    #E
                    #axs[database_counter,current_counter-1].set_yticks([-1.,-0.5,0.])
                    #axs[database_counter,current_counter-1].set_xlim([0,1.5e2]) 
                    #axs[database_counter,current_counter-1].set_ylim([-1.,0.3])
                    #R
                    #axs[database_counter,current_counter-1].set_yticks([-1.,-0.5,0.])
                    #axs[database_counter,current_counter-1].set_xlim([4.,9.5])
                    #axs[database_counter,current_counter-1].set_ylim([-1.,0.3])
                    DFN_3Ds=[]
                    mom_3Ds=[]
                    for phase_counter,(parameters__phase_upper,parameters__phase_middle,parameters__phase_lower) in enumerate(zip(parameters__phases_upper,parameters__phases_middle,parameters__phases_lower)): #nest at same level == offset them together rigidly 
                        if output_dfns[database_counter,mode_counter,phase_counter,current_counter]:
                            DFN_3Ds.append(output_dfns[database_counter,mode_counter,phase_counter,current_counter].transform(axes=plot_axes)['dfn'])            
                            mom_3Ds.append(output_moms[database_counter,mode_counter,phase_counter,current_counter][moment_quantity])            
                            #either plot each DFN here + 2D axisymmetric 
                            marker='+' if -3 in parameters__toroidal_mode_number else '.'
                            if plot_mode==1:
                                #plot distribution
                                axs[database_counter,current_counter-1].plot(
                                DFN_2D[plot_axes[0]],
                                DFN_3Ds[-1],
                                #label=r'$\Phi_{\mathrm{m}}={}$'.format(str(parameters__phase_middle)),
                                label='{}'.format(str(parameters__phase_middle)),
                                color=settings.cmap_inferno(phase_counter/len(parameters__phases_upper)),
                                marker=marker,
                                )
                                axs_mom[database_counter,current_counter-1].plot(
                                mom_2D[moment_axis],
                                mom_3Ds[-1],
                                #label=r'$\Phi_{\mathrm{m}}={}$'.format(str(parameters__phase_middle)),
                                label='{}'.format(str(parameters__phase_middle)),
                                color=settings.cmap_inferno(phase_counter/len(parameters__phases_upper)),
                                marker=marker,
                                )
                            if plot_mode==2:
                                label=r'$n$={}'.format(np.abs(parameters__toroidal_mode_number[0])) if phase_counter==0 else ''
                                colours=plt.rcParams['axes.prop_cycle'].by_key()['color']
                                colour=colours[0] if -3 in parameters__toroidal_mode_number else colours[1]
                                #plot single point against absolute phase over absolute phase
                                axs[database_counter,current_counter-1].scatter(
                                parameters__phase_middle,
                                (DFN_3Ds[-1][int(len(DFN_3Ds[-1])/2)]-DFN_2D['dfn'][int(len(DFN_3Ds[-1])/2)]),#/DFN_2D['dfn'][int(len(DFN_3Ds[-1])/2)],
                                color=colour,#settings.cmap_viridis(phase_counter/len(parameters__phases_upper)),
                                marker=marker,
                                label=label
                                )
                        if plot_mode==1:
                            axs[database_counter,current_counter-1].plot(
                            DFN_2D[plot_axes[0]],
                            DFN_2D['dfn'],
                            linestyle='dashed',
                            color='black',
                            )                                
                            axs_mom[database_counter,current_counter-1].plot(
                            mom_2D[moment_axis],
                            mom_2D[moment_quantity],
                            linestyle='dashed',
                            color='black',
                            )   
                    #or plot range between DFN differences here
                    if plot_mode==3:
                        label=f'{parameters__database.split("_")[-1]} '
                        label=f'$n$ = {"+".join(str(abs(int(mode))) for mode in parameters__toroidal_mode_number)}'
                        DFN_3Ds=np.array(DFN_3Ds)                    
                        DFN_min=np.min(DFN_3Ds,axis=0)
                        DFN_max=np.max(DFN_3Ds,axis=0)
                        axs[database_counter,current_counter-1].fill_between(
                            DFN_2D[plot_axes[0]],
                            np.nan_to_num((DFN_max-DFN_2D['dfn'])),#/DFN_2D['dfn']),
                            np.nan_to_num((DFN_min-DFN_2D['dfn'])),#/DFN_2D['dfn']),
                            label=label,
                            alpha=0.45,
                            #color=colours[2*mode_counter],
                            )        
                        mom_3Ds=np.array(mom_3Ds)                    
                        mom_min=np.min(mom_3Ds,axis=0)
                        mom_max=np.max(mom_3Ds,axis=0)
                        axs_mom[database_counter,current_counter-1].fill_between(
                            mom_2D[moment_axis],
                            np.nan_to_num((mom_max-mom_2D[moment_quantity])/mom_2D[moment_quantity]),
                            np.nan_to_num((mom_min-mom_2D[moment_quantity])/mom_2D[moment_quantity]),
                            label=label,
                            alpha=0.45,
                            #color=colours[2*mode_counter],
                            )  
                    #or plot range between DFNs here
                    if plot_mode==4:
                        label=f'{parameters__database.split("_")[-1]} '
                        label=f'$n$ = {"+".join(str(abs(int(mode))) for mode in parameters__toroidal_mode_number)}'
                        DFN_3Ds=np.array(DFN_3Ds)                    
                        DFN_min=np.min(DFN_3Ds,axis=0)
                        DFN_max=np.max(DFN_3Ds,axis=0)
                        axs[database_counter,current_counter-1].fill_between(
                            DFN_2D[plot_axes[0]],
                            np.nan_to_num(DFN_max),
                            np.nan_to_num(DFN_min),
                            label=label,
                            alpha=0.45,
                            #color=colours[2*mode_counter],
                            )       
                        label_axisymmetric='axisymmetric' if -3 in parameters__toroidal_mode_number else None
                        axs[database_counter,current_counter-1].plot(
                            DFN_2D[plot_axes[0]],
                            DFN_2D['dfn'],
                            label=label_axisymmetric,
                            linestyle='dashed',
                            color='black'
                            )        
                        mom_3Ds=np.array(mom_3Ds)                    
                        mom_min=np.min(mom_3Ds,axis=0)
                        mom_max=np.max(mom_3Ds,axis=0)
                        axs_mom[database_counter,current_counter-1].fill_between(
                            mom_2D[moment_axis],
                            mom_max,
                            mom_min,
                            label=label,
                            alpha=0.45,
                            #color=colours[2*mode_counter],
                            )  
                        axs_mom[database_counter,current_counter-1].plot(
                            mom_2D[moment_axis],
                            mom_2D[moment_quantity],
                            label=label_axisymmetric,
                            linestyle='dashed',
                            color='black'
                            )  
                    axs[database_counter,0].get_shared_y_axes().join(axs[database_counter,0],axs[database_counter,current_counter-1])
                    axs_mom[database_counter,0].get_shared_y_axes().join(axs_mom[database_counter,0],axs_mom[database_counter,current_counter-1])
                    #axs[database_counter,current_counter-1].set_xticks(axs[database_counter,0].get_xticks())
                    #axs[database_counter,current_counter-1].set_yticks(axs[database_counter,0].get_yticks())
                    #axs_mom[database_counter,current_counter-1].set_xticks(axs_mom[database_counter,0].get_xticks())
                    #axs_mom[database_counter,current_counter-1].set_yticks(axs_mom[database_counter,0].get_yticks())
                    axs[database_counter,current_counter-1].text(
                        #f'{parameters__database.split("_")[-1]}'+', $I_{\mathrm{c}}$'+f'={int(parameters__current_upper/1000)}kAt',
                        #loc='right',y=1.,
                        fontsize=10,
                        color=settings.colour_custom(rgba=[100,100,100,1])(0.),
                        x=0.05,
                        y=0.05,
                        s=f'{parameters__database.split("_")[-1]}'+', $I_{\mathrm{c}}$'+f'={int(parameters__current_upper/1000)}kAt',
                        ha='left',
                        va='bottom',
                        rotation='horizontal',
                        transform=axs[database_counter,current_counter-1].transAxes,
                    )
                    axs_mom[database_counter,current_counter-1].text(
                        #f'{parameters__database.split("_")[-1]}'+', $I_{\mathrm{c}}$'+f'={int(parameters__current_upper/1000)}kAt',
                        #loc='right',y=1.,
                        fontsize=10,
                        color=settings.colour_custom(rgba=[100,100,100,1])(0.),
                        x=0.05,
                        y=0.05,
                        s=f'{parameters__database.split("_")[-1]}'+', $I_{\mathrm{c}}$'+f'={int(parameters__current_upper/1000)}kAt',
                        ha='left',
                        va='bottom',
                        rotation='horizontal',
                        transform=axs_mom[database_counter,current_counter-1].transAxes,
                    )
                    axs[database_counter,current_counter-1].tick_params(
                        axis='both',          
                        which='both',      
                        labelbottom=False,
                        labelleft=False,
                        left=False,
                        )
    #except:
    #    pass
    
for axis in axs[-1,:]:
    axis.tick_params(
        axis='both',          
        which='both',      
        labelbottom=True,
        labelsize=10,
        )

for axis in axs[:,0]:
    axis.tick_params(
        axis='both',          
        which='both',      
        labelleft=True,
        labelsize=10,
        left=True,
        )

ax_total.set_xlabel('R [m]',fontsize=15)
ax_total.set_xlabel('Energy [keV]',fontsize=15)
ax_total.set_xlabel(r'Pitch ($v_{\parallel}/v$)',fontsize=15)
ax_total.set_xlabel('Absolute phase shift of RMP ($\Phi_{\mathrm{m}}$) [deg]',fontsize=15) #\Phi for absolute

ax_total.set_ylabel(r'Fast-ion density [$\mathrm{m}^{-1}$]',fontsize=15)
ax_total.set_ylabel('Fractional change in fast-ion density',fontsize=15)
ax_total.set_ylabel(r'Fast-ion density [$\mathrm{eV}^{-1}$]',fontsize=15,labelpad=-0)
ax_total.set_ylabel(r'Fast-ion density [$\lambda^{-1}$]',fontsize=15)
ax_total.set_ylabel(r'Change in fast-ion density at $\lambda=0$ [$\lambda^{-1}$]',fontsize=15)

axs[0,0].legend(bbox_to_anchor=(0.5,1.05),ncol=2,loc='center',bbox_transform=ax_total.transAxes,fontsize=10)

#plt.savefig(
#    fname=f'/home/ITER/wards2/pics/paper_3_DNB_loss_bounds_E.pdf',
#    dpi=800,
#)


DFNs_for_2D_plots=[]
for dfn in output_dfns[...,1:].flatten():
    DFNs_for_2D_plots.append(dfn.transform(axes=plot_axes_2D)['dfn'])

DFNs_for_2D_plots=np.array(DFNs_for_2D_plots).reshape(
    len(batch_data.parameters__databases),
    len(batch_data.parameters__toroidal_mode_numbers),
    len(batch_data.parameters__phases_uppers_case7[0]),
    len(batch_data.parameters__currents_upper)-1, #got rid of zero current case
    len(dfn[plot_axes_2D[0]]),
    len(dfn[plot_axes_2D[1]]),
    )

for database_counter,(
        parameters__database,parameters__sheet_name_kinetic_prof,
        parameters__phases_uppers,parameters__phases_middles,parameters__phases_lowers,
        parameters__currents_upper,parameters__currents_middle,parameters__currents_lower 
        ) in enumerate(zip(
        batch_data.parameters__databases,batch_data.parameters__sheet_names_kinetic_prof,
        batch_data.parameters__phases_uppers_cases_all,batch_data.parameters__phases_middles_cases_all,batch_data.parameters__phases_lowers_cases_all,
        batch_data.parameters__currents_uppers,batch_data.parameters__currents_middles,batch_data.parameters__currents_lowers
        )): 
    DFN_2D=output_dfns[database_counter,0,0,0].transform(axes=plot_axes_2D)
    DFN_2D['E']/=1.e3 #plot in keV
    #set ticks in first ax of the row
    for parameters__kinetic_prof_tF_tE,parameters__kinetic_prof_Pr in zip(batch_data.parameters__kinetic_profs_tF_tE,batch_data.parameters__kinetic_profs_Pr):
        for mode_counter,(ax_group,fig_2D,parameters__toroidal_mode_number,parameters__phases_upper,parameters__phases_middle,parameters__phases_lower) in enumerate(zip([axs_2D_n3,axs_2D_n4],[fig_2D_n3,fig_2D_n4],batch_data.parameters__toroidal_mode_numbers,parameters__phases_uppers,parameters__phases_middles,parameters__phases_lowers)):
            ax_group[database_counter,0].set_xlim([0,1.6e2]) 
            ax_group[database_counter,0].set_ylim([-1,1])
            ax_group[database_counter,0].set_xticks(np.linspace(0,150,11))
            ax_group[database_counter,0].set_yticks([-1.,-0.5,0.,0.5,1.])
            for parameters__rotation_upper,parameters__rotation_middle,parameters__rotation_lower in zip(batch_data.parameters__rotations_upper,batch_data.parameters__rotations_middle,batch_data.parameters__rotations_lower): 
                for current_counter,(parameters__current_upper,parameters__current_middle,parameters__current_lower) in enumerate(zip(
                    parameters__currents_upper,
                    parameters__currents_middle,
                    parameters__currents_lower)): 
                    if parameters__current_upper==0: continue
                    ax_group[database_counter,0].get_shared_y_axes().join(ax_group[database_counter,0],ax_group[database_counter,current_counter-1])
                    ax_group[database_counter,current_counter-1].set_xticks(ax_group[database_counter,0].get_xticks())
                    ax_group[database_counter,current_counter-1].set_yticks(ax_group[database_counter,0].get_yticks())
                    DFN_2D_plot_min=np.min(DFNs_for_2D_plots[database_counter,mode_counter,:,current_counter-1,...],axis=0)
                    DFN_2D_plot_max=np.max(DFNs_for_2D_plots[database_counter,mode_counter,:,current_counter-1,...],axis=0)
                    DFN_2_plot=copy.deepcopy(output_dfns[database_counter,0,0,0])
                    DFN_2_plot['E']/=1000.
                    DFN_diff=np.nan_to_num(DFN_2D_plot_max-DFN_2D_plot_min)
                    DFN_2_plot.set(dfn=(DFN_diff-np.min(DFN_diff))/(np.max(DFN_diff)-np.min(DFN_diff))),#/DFN_2_plot.transform(axes=plot_axes_2D)['dfn']))
                    mesh=DFN_2_plot.plot(
                        fig=fig_2D,
                        ax=ax_group[database_counter,current_counter-1],
                        transform=False,
                        axes=plot_axes_2D,
                        colmap=settings.cmap_inferno,
                        number_bins=11,
                    )
                    ax_group[database_counter,current_counter-1].text(x=0.95,y=0.05,s=f'{parameters__database.split("_")[-1]}'+', $I_{\mathrm{c}}$'+f'={int(parameters__current_upper/1000)}kAt',fontsize=10,
                    ha='right',
                    va='bottom',
                    rotation='horizontal',
                    transform=ax_group[database_counter,current_counter-1].transAxes,
                    color=settings.colour_custom(rgba=[242,242,242,1])(0.),
                    )
                    ax_group[database_counter,current_counter-1].set_xlabel('')
                    ax_group[database_counter,current_counter-1].set_ylabel('')
                    ax_group[database_counter,current_counter-1].set_xlim([0,150])
                    ax_group[database_counter,current_counter-1].set_ylim([-1.,1.])
                    ax_group[database_counter,current_counter-1].set_title('')
                    ax_group[database_counter,current_counter-1].tick_params(
                        axis='both',          
                        which='both',      
                        labelbottom=False,
                        labelleft=False,    
                        )

            for axis in ax_group[-1,:]:
                axis.tick_params(
                    axis='both',          
                    which='both',      
                    labelbottom=True,
                    labelsize=10,
                    )
            for axis in ax_group[:,0]:
                axis.tick_params(
                    axis='both',          
                    which='both',      
                    labelleft=True,
                    labelsize=10,
                    )

title_n3=r'Variation in DNB density due to $n=3$ rotating RMP [a.u.]'
title_n4=r'Variation in DNB density due to $n=4$ rotating RMP [a.u.]'

ax_total_2D_n3.set_title(title_n3,fontsize=15)
ax_total_2D_n4.set_title(title_n4,fontsize=15)
ax_total_2D_n3.set_xlabel('Energy [keV]',fontsize=15)
ax_total_2D_n3.set_ylabel(r'Pitch ($v_{\parallel}/v$)',fontsize=15)
ax_total_2D_n4.set_xlabel('Energy [keV]',fontsize=15)
ax_total_2D_n4.set_ylabel(r'Pitch ($v_{\parallel}/v$)',fontsize=15)

cbar_ax = fig_2D_n3.add_axes([0.91, 0.1, 0.025, 0.7])
cbar=fig_2D_n3.colorbar(mesh,cax=cbar_ax,orientation='vertical',)#,pad=0.2)
cbar.ax.tick_params(labelsize=10)
cbar_ax = fig_2D_n4.add_axes([0.91, 0.1, 0.025, 0.7])
cbar=fig_2D_n4.colorbar(mesh,cax=cbar_ax,orientation='vertical')#,pad=0.2)
cbar.ax.tick_params(labelsize=10)


if __name__ == '__main__':
    plt.show()





