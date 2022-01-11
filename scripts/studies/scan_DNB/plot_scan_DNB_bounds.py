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

import scan_DNB_launch as batch_data

output_dfns=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'dfn',batch_data,processes=16,chunksize=1)).reshape(
    len(batch_data.parameters__databases),
    len(batch_data.parameters__toroidal_mode_numbers),
    len(batch_data.parameters__phases_uppers_case7[0]),
    len(batch_data.parameters__currents_upper),
    )

fig=plt.figure(figsize=(10, 15),constrained_layout=False)
dimensions=(len(batch_data.parameters__databases),len(batch_data.parameters__currents_uppers[0])-1)
ax_total=fig.add_subplot(1,1,1,frameon=False)
ax_total.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_total.spines['top'].set_color('none')
ax_total.spines['bottom'].set_color('none')
ax_total.spines['left'].set_color('none')
ax_total.spines['right'].set_color('none')
axs=[]
_=1
for dim in dimensions: _*=dim
for _ in range(_):
    axs.append(fig.add_subplot(*dimensions,_+1))

axs=np.array(axs).T.reshape(dimensions)

plot_axes=['E','V_pitch']
plot_axes=['R','Z']
plot_axes=['V_pitch']
plot_axes=['E']

for database_counter,(
        parameters__database,parameters__sheet_name_kinetic_prof,
        parameters__phases_uppers,parameters__phases_middles,parameters__phases_lowers,
        parameters__currents_upper,parameters__currents_middle,parameters__currents_lower 
        ) in enumerate(zip(
        batch_data.parameters__databases,batch_data.parameters__sheet_names_kinetic_prof,
        batch_data.parameters__phases_uppers_cases_all,batch_data.parameters__phases_middles_cases_all,batch_data.parameters__phases_lowers_cases_all,
        batch_data.parameters__currents_uppers,batch_data.parameters__currents_middles,batch_data.parameters__currents_lowers
        )): 
    try:
        DFN_2D=output_dfns[database_counter,0,0,0].transform(axes=plot_axes)
        DFN_2D['E']/=1.e3 #plot in keV
        for parameters__kinetic_prof_tF_tE,parameters__kinetic_prof_Pr in zip(batch_data.parameters__kinetic_profs_tF_tE,batch_data.parameters__kinetic_profs_Pr):
            for mode_counter,(parameters__toroidal_mode_number,parameters__phases_upper,parameters__phases_middle,parameters__phases_lower) in enumerate(zip(batch_data.parameters__toroidal_mode_numbers,parameters__phases_uppers,parameters__phases_middles,parameters__phases_lowers)):
                for parameters__rotation_upper,parameters__rotation_middle,parameters__rotation_lower in zip(batch_data.parameters__rotations_upper,batch_data.parameters__rotations_middle,batch_data.parameters__rotations_lower): 
                    for current_counter,(parameters__current_upper,parameters__current_middle,parameters__current_lower) in enumerate(zip(
                        parameters__currents_upper,
                        parameters__currents_middle,
                        parameters__currents_lower)):                    
                        if parameters__current_upper==0: continue
                        DFN_3Ds=[]
                        for phase_counter,(parameters__phase_upper,parameters__phase_middle,parameters__phase_lower) in enumerate(zip(parameters__phases_upper,parameters__phases_middle,parameters__phases_lower)): #nest at same level == offset them together rigidly 
                            if output_dfns[database_counter,mode_counter,phase_counter,current_counter]:
                                DFN_3Ds.append(output_dfns[database_counter,mode_counter,phase_counter,current_counter].transform(axes=plot_axes)['dfn'])            
                                #either plot each DFN here + 2D axisymmetric 
                            """ 
                                marker='+' if -3 in parameters__toroidal_mode_number else '.'
                                axs[database_counter,current_counter-1].plot(
                                DFN_2D[plot_axes[0]],
                                DFN_3Ds[-1],
                                #label=r'$\Phi_{\mathrm{m}}={}$'.format(str(parameters__phase_middle)),
                                label='{}'.format(str(parameters__phase_middle)),
                                color=settings.cmap_inferno(phase_counter/len(parameters__phases_upper)),
                                marker=marker,
                                )
                            axs[database_counter,current_counter-1].plot(
                            DFN_2D[plot_axes[0]],
                            DFN_2D['dfn'],
                            linestyle='dashed',
                            color='black',
                            )                                
                            """
                        #or plot range between DFNs here
                        #"""
                        label=f'{parameters__database.split("_")[-1]} '
                        label=f'$n$ = {"+".join(str(abs(int(mode))) for mode in parameters__toroidal_mode_number)}'
                        DFN_3Ds=np.array(DFN_3Ds)                    
                        DFN_min=np.min(DFN_3Ds,axis=0)
                        DFN_max=np.max(DFN_3Ds,axis=0)
                        axs[database_counter,current_counter-1].fill_between(
                            DFN_2D[plot_axes[0]],
                            np.nan_to_num((DFN_max-DFN_2D['dfn'])/DFN_2D['dfn']),
                            np.nan_to_num((DFN_min-DFN_2D['dfn'])/DFN_2D['dfn']),
                            label=label,
                            alpha=0.45,
                            )        
                        #"""
                        axs[database_counter,current_counter-1].text(x=0.1,y=0.05,s=f'{parameters__database.split("_")[-1]}'+', $I_{\mathrm{p}}$'+f'={int(parameters__current_upper/1000)}kAt',fontsize=10)
                        axs[database_counter,current_counter-1].set_xlabel('')
                        axs[database_counter,current_counter-1].set_ylabel('')
                        #axs[database_counter,current_counter-1].set_yticks([-1.,-0.5,0.])
                        ##axs[database_counter,current_counter-1].set_xlim([-1,1])
                        #axs[database_counter,current_counter-1].set_xlim([0,1e2]) #beam energy is 100keV
                        #axs[database_counter,current_counter-1].set_ylim([-1.,0.])
                        axs[database_counter,current_counter-1].tick_params(
                            axis='both',          
                            which='both',      
                            labelbottom=False,
                            labelleft=False,
                            )
    except:
        pass
    
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
        )

ax_total.set_xlabel('Pitch')
ax_total.set_xlabel('Energy [keV]',fontsize=15)
ax_total.set_ylabel('Fractional change in fast-ion density',fontsize=15)
axs[0,0].legend(loc='lower right',fontsize=10)

#plt.savefig(
#    fname=f'/home/ITER/wards2/pics/paper_3_DNB_loss_bounds_E.pdf',
#    dpi=800,
#)

if __name__ == '__main__':
    plt.show()