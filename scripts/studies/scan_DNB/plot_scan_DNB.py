#plot_scan_DNB.py
 
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

Pinj=.13e6

import scan_DNB_launch as batch_data

PFC_power=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.calc_PFC_power,'rund',batch_data,processes=32,chunksize=1)).reshape(
    len(batch_data.parameters__databases),
    len(batch_data.parameters__toroidal_mode_numbers),
    len(batch_data.parameters__phases_uppers_case7[0]),
    len(batch_data.parameters__currents_upper),
    )

fig,axs=plt.subplots(2,constrained_layout=False)
for scenario_counter,scenario in enumerate(batch_data.parameters__databases):
    for mode_number_counter,mode_number in enumerate(batch_data.parameters__toroidal_mode_numbers):
        for current_counter,current in enumerate(batch_data.parameters__currents_uppers[scenario_counter][1:]):
            
            label=f'{scenario.split("_")[-1]}'
            label+=f', $n$ = {"+".join(str(abs(int(mode))) for mode in mode_number)}'
            axs[current_counter].plot(batch_data.parameters__phases_middles_cases_all[scenario_counter][mode_number_counter],100.*PFC_power[scenario_counter,mode_number_counter,:,current_counter+1]/Pinj,
            label=label,
            )
fig.subplots_adjust(right=0.85)
axs[0].set_xlabel('$\Delta\Phi$ [deg]') #\Phi for absolute
axs[1].set_xlabel('$\Delta\Phi$ [deg]') #\Phi for absolute
axs[0].set_ylabel('Deposited power lost [%]')
axs[1].legend(loc='center',bbox_to_anchor=(1.1,0.5),ncol=1)
plt.show()


output_dfns=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'dfn',batch_data,processes=8,chunksize=1)).reshape(
    len(batch_data.parameters__databases),
    len(batch_data.parameters__toroidal_mode_numbers),
    len(batch_data.parameters__phases_uppers_case7[0]),
    len(batch_data.parameters__currents_upper),
    )

fig=plt.figure(constrained_layout=False)
dimensions=(len(batch_data.parameters__toroidal_mode_numbers),len(batch_data.parameters__databases))
ax_total=fig.add_subplot(1,1,1,frameon=False)
ax_total.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_total.spines['top'].set_color('none')
ax_total.spines['bottom'].set_color('none')
ax_total.spines['left'].set_color('none')
ax_total.spines['right'].set_color('none')
axs=[]

plot_axes=['E','V_pitch']
plot_axes=['R','Z']
plot_axes=['V_pitch']
plot_axes=['E']
ax_counter=0
for database_counter, (parameters__database,parameters__sheet_name_kinetic_prof,plasma_species,
    parameters__phases_upper,parameters__phases_middle,parameters__phases_lower) in enumerate(zip(
    batch_data.parameters__databases,
    batch_data.parameters__sheet_names_kinetic_prof,
    batch_data.plasmas_species,
    batch_data.parameters__phases_uppers_cases_all,
    batch_data.parameters__phases_middles_cases_all,
    batch_data.parameters__phases_lowers_cases_all)):
    DFN_2D=output_dfns[database_counter,0,0,0].transform(axes=plot_axes)
    for mode_number_counter,(parameters__toroidal_mode_number,parameters__phase_upper,parameters__phase_middle,parameters__phase_lower) in enumerate(zip(batch_data.parameters__toroidal_mode_numbers,batch_data.parameters__phases_upper,batch_data.parameters__phases_middle,batch_data.parameters__phases_lower)):
        ax_counter+=1    
        axs.append(fig.add_subplot(*dimensions,ax_counter))
        #remove labels
        axs[-1].tick_params(
            axis='both',          
            which='both',      
            labelbottom=False,
            labelleft=False,
            )
        
        print(f'total difference for case={parameters__database}, n={parameters__toroidal_mode_number}: {output_dfns[database_counter,mode_number_counter,0,-1].transform(axes=["N"])["dfn"]-output_dfns[database_counter,0,0,0].transform(axes=["N"])["dfn"]}')
        DFN_diff=output_dfns[database_counter,mode_number_counter,0,-1].transform(axes=plot_axes)
        DFN_diff['dfn']=np.nan_to_num((DFN_diff['dfn']-DFN_2D['dfn'])/DFN_2D['dfn'],nan=0.)
        DFN_diff['dfn'][np.abs(DFN_diff['dfn'])>1e10]=0
        """
        output_dfns[database_counter,mode_number_counter,0,-1].plot(
            axes=plot_axes,
            ax=axs[-1],
            fig=fig,
            colmap=settings.cmap_r,
        )
        output_dfns[database_counter,0,0,0].plot(
            axes=plot_axes,
            ax=axs[-1],
            fig=fig,
            colmap=settings.cmap_g,
        )
        """ 
        DFN_diff['E']/=1.e3 #convert to keV
        mesh=DFN_diff.plot(
            axes=plot_axes,
            transform=False,
            ax=axs[-1],
            fig=fig,
            colmap=settings.cmap_k,
            #vminmax=[-5,5],
            #real_scale=True
        )
        #plt.colorbar(mesh,ax=axs[-1])
        axs[-1].set_title('')
        axs[-1].set_xlabel('')
        axs[-1].set_ylabel('')
        axs[-1].set_xlim([-1,1])
        axs[-1].set_xlim([0,1e2]) #beam energy is 100keV
        axs[-1].set_ylim([-1.,0.])

axs=np.array(axs).reshape(dimensions)
#re-add labels for edge plots
for axis in axs[-1,:]:
    axis.tick_params(
        axis='both',          
        which='both',      
        labelbottom=True,
        )
for axis in axs[:,0]:
    axis.tick_params(
        axis='both',          
        which='both',      
        labelleft=True,
        )

ax_total.set_xlabel('Energy [keV]')
ax_total.set_ylabel('Pitch')

if __name__ == '__main__':
    plt.show()