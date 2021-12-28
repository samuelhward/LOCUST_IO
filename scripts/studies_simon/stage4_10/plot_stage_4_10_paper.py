#plot_stage_4_10_paper.py
 
"""
Samuel Ward
29/06/21
----
script for plotting stage 4.10 of Simon's studies 
used to create plot for paper 3
make sure the 1_10 launch script you're importing only the 60kAt runs
---
 
notes:         
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
    import scipy.constants
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

import stage_4_10_launch as batch_data

method='rund'

Pinj=33.e6
PFC_power=templates.plot_mod.apply_func_parallel(templates.plot_mod.calc_PFC_power,method,batch_data,processes=16,chunksize=1)
PFC_power=np.array(PFC_power).reshape(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers),len(batch_data.parameters__phases_upper))

#account for possible difference in mass of injected isotope
if method=='fpl':
    PFC_power/=constants.mass_deuteron #different masses in some cases, so divide by mass here for re-multiplication later
    for scenario_counter,(scenario,beam_species) in enumerate(zip(batch_data.parameters__databases,batch_data.configs_beam_species)):
        if beam_species=='hydrogen':
            PFC_power[scenario_counter]*=scipy.constants.physical_constants['proton mass'][0]
        elif beam_species=='deuterium':
            PFC_power[scenario_counter]*=scipy.constants.physical_constants['deuteron mass'][0]

fig,ax=plt.subplots(1,constrained_layout=True)
for scenario_counter,scenario in enumerate(batch_data.parameters__databases):
    for mode_number_counter,mode_number in enumerate(batch_data.parameters__toroidal_mode_numbers):
        #find relative phase between coils for this toroidal mode number
        relative_phases_upper_middle=batch_data.parameters__phases_uppers_cases_all[scenario_counter][mode_number_counter]-batch_data.parameters__phases_middles_cases_all[scenario_counter][mode_number_counter]
        relative_phases_lower_middle=batch_data.parameters__phases_lowers_cases_all[scenario_counter][mode_number_counter]-batch_data.parameters__phases_middles_cases_all[scenario_counter][mode_number_counter]
        label=f'{scenario.split("_")[-1]} ({batch_data.configs_beam_species[scenario_counter][0].capitalize()})'
        label+=f', $n$ = {"+".join(str(abs(int(mode))) for mode in mode_number)}'
        #label+=', '+r'$n\Delta\Phi_{\mathrm{u,m}}$'
        #label+=f'={int(np.abs(mode_number[0])*relative_phases_upper_middle[0])}'
        #label+=', '+r'$n\Delta\Phi_{\mathrm{l,m}}$='
        #label+=f'{int(np.abs(mode_number[0])*relative_phases_lower_middle[0])}'
        #plot relative phase
        ax.plot(batch_data.parameters__phases_middles_cases_all[scenario_counter][mode_number_counter],100.*PFC_power[scenario_counter,mode_number_counter]/Pinj,label=label
        #label=f'{scenario} ({batch_data.configs_beam_species[scenario_counter]}), $n$ = {"+".join(str(abs(int(mode))) for mode in mode_number)}, U:M={int(np.abs(mode_number[0])*relative_phases_upper_middle[0])}, L:M={int(np.abs(mode_number[0])*relative_phases_lower_middle[0])}'
        )

#now add 1_10 data to plot (or could re-run 1_10 to fill in missing simulations...it doesn't matter - but 1_10 data are huge! so think about editing the launch file...)
#'''
import studies_simon.stage1_10.stage_1_10_launch as batch_data_1_10

inds_to_del=[counter for counter,val in enumerate(batch_data_1_10.args_batch['LOCUST_run__dir_output']) if '4_10' in str(val)]
inds_to_del.reverse()
for key,value in batch_data_1_10.args_batch.items():
    if value:
        for ind in inds_to_del:
            del(batch_data_1_10.args_batch[key][ind])
        
PFC_power_1_10=templates.plot_mod.apply_func_parallel(templates.plot_mod.calc_PFC_power,method,batch_data_1_10,processes=4,chunksize=1)

PFC_power_1_10=np.array(PFC_power_1_10).reshape(
                    len(batch_data_1_10.parameters__kinetic_profs_Pr),
                    len(batch_data_1_10.parameters__currents_upper),
                    len(batch_data_1_10.parameters__toroidal_mode_numbers),
                    len(batch_data_1_10.parameters__phases_upper)
                    )

for scenario_counter,scenario in enumerate(batch_data_1_10.parameters__databases):
    for plasma_state_counter,(Pr,tftE) in enumerate(zip(batch_data_1_10.parameters__kinetic_profs_Pr,batch_data_1_10.parameters__kinetic_profs_tF_tE)):
        for current_counter,(current,linestyle) in enumerate(zip(batch_data_1_10.parameters__currents_upper,['solid','solid'])):
            for mode_number_counter,mode_number in enumerate(batch_data_1_10.parameters__toroidal_mode_numbers):

                relative_phases_upper_middle=batch_data_1_10.parameters__phases_uppers[mode_number_counter]-batch_data_1_10.parameters__phases_middles[mode_number_counter]
                relative_phases_lower_middle=batch_data_1_10.parameters__phases_lowers[mode_number_counter]-batch_data_1_10.parameters__phases_middles[mode_number_counter]
                label=f'{scenario.split("_")[-1]} ({batch_data_1_10.configs_beam_species[scenario_counter][0].capitalize()})'
                label+=f', $n$ = {"+".join(str(abs(int(mode))) for mode in mode_number)}'
                #label+=', '+r'$n\Delta\Phi_{\mathrm{u,m}}$'
                #label+=f'={int(np.abs(mode_number[0])*relative_phases_upper_middle[0])}'
                #label+=', '+r'$n\Delta\Phi_{\mathrm{l,m}}$='
                #label+=f'{int(np.abs(mode_number[0])*relative_phases_lower_middle[0])}'
                ax.plot(batch_data_1_10.parameters__phases_middles[mode_number_counter],100*PFC_power_1_10[plasma_state_counter,current_counter,mode_number_counter]/Pinj,label=label)

ax.set_xlabel('Absolute phase shift of RMP ($\Phi_{\mathrm{M}}$) [deg]') #\Phi for absolute
ax.set_ylabel('NBI power loss [%]')
ax.set_ylim([0,8.1])
#fig.subplots_adjust(right=0.8)
ax.legend()#loc='center',bbox_to_anchor=(1.15,0.5),ncol=1)

#'''
plt.show()

#################################
 
##################################################################
 
###################################################################################################