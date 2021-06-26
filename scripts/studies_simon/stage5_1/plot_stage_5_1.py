#plot_stage_5_1.py
 
"""
Samuel Ward
02/06/21
----
script for plotting stage 5 of Simon's studies 
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
    from random_scripts.plot_X_point_displacement import plot_X_point_displacement
except:
    raise ImportError("ERROR: templates/template_mod.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## 
#Main 

import stage_5_1_launch as batch_data

outputs=templates.plot_mod.get_output_files(batch_data,'fpl')

Pinj=33.e6
PFC_power=[]
for output in outputs:
    if output:
        i=np.where(output['status_flag']=='PFC_intercept_3D')[0]
        PFC_power.append([1.e6*output['f']*np.sum((output['V_R'][i]**2+output['V_phi'][i]**2+output['V_Z'][i]**2)*output['FG'][i])*0.5])
    else:
        PFC_power.append([-10.])

'''
outputs=templates.plot_mod.get_output_files(batch_data,'rund')

Pinj=33.e6
PFC_power=[]
for output in outputs:
    if output is not None and output['run_status']=='completed':
        PFC_power.append([output['PFC_power']['total']])
    else:
        PFC_power.append([-10.])

'''

PFC_power=np.array(PFC_power).reshape(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers),len(batch_data.parameters__phases_upper))

for scenario_counter,(scenario,beam_species) in enumerate(zip(batch_data.parameters__databases,batch_data.configs_beam_species)):
    if beam_species=='hydrogen':
        PFC_power[scenario_counter]*=scipy.constants.physical_constants['proton mass'][0]
    elif beam_species=='deuterium':
        PFC_power[scenario_counter]*=scipy.constants.physical_constants['deuteron mass'][0]
#determine unique scenarios we are studying
scenarios={scenario for scenario in batch_data.parameters__databases}

# one figure per plasma
fig,axs=plt.subplots(len(scenarios)*len(batch_data.parameters__toroidal_mode_numbers),2,constrained_layout=True)
axs=np.array([axs]).reshape(len(scenarios),len(batch_data.parameters__toroidal_mode_numbers),2)
axes={scenario:ax for scenario,ax in zip(scenarios,axs)}
for plasma_state_counter,(database,beam_species) in enumerate(zip(batch_data.parameters__databases,batch_data.configs_beam_species)):
    for mode_number_counter,(mode_number) in enumerate(batch_data.parameters__toroidal_mode_numbers):
        if beam_species == 'hydrogen':
            colour='r'
        elif beam_species == 'deuterium':
            colour='b'
        else:
            colour='k'
        case=int(batch_data.parameters__databases[plasma_state_counter][-1])
        n=int(np.abs(mode_number[0]))
        ax_left,ax_right=axes[database][mode_number_counter]        
        plot_X_point_displacement(case=case,n=n,LVV=90,fig=fig,ax=ax_left)
        ax_right.scatter(np.arange(len(PFC_power[plasma_state_counter][mode_number_counter])),100*PFC_power[plasma_state_counter][mode_number_counter]/Pinj,label=f'({beam_species} plasma)',color=colour)
        ax_right.set_xlabel('Point',fontsize=15)
        ax_right.set_ylabel(f'Deposited power lost [%]',fontsize=15)
        ax_right.set_ylabel('Ploss [%]',fontsize=15)

#axes[batch_data.parameters__databases[0]][0][1].legend() 

# remove ticks from total ax
#ax_total = fig.add_subplot(111,frameon=False)
#ax_total.tick_params(axis='both',which='both',bottom=False,labelbottom=False,left=False,labelleft=False)

plt.show()

#################################
 
##################################################################
 
###################################################################################################