#plot_stage_4_10.py
 
"""
Samuel Ward
29/06/21
----
script for plotting stage 4.10 of Simon's studies 
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

def calc_PFC_power(filename,classtype='fpl'):
    fpl=templates.plot_mod.read_locust_io_obj(filename=filename,classtype=classtype)
    if fpl:
        i=np.where(fpl['status_flag']=='PFC_intercept_3D')[0]
        PFC_power=1.e6*fpl['f']*np.sum((fpl['V_R'][i]**2+fpl['V_phi'][i]**2+fpl['V_Z'][i]**2)*fpl['FG'][i])*0.5 #leave out species mass for now, since applied later
    else:
        PFC_power=-1
    return PFC_power 

Pinj=33.e6
PFC_power=templates.plot_mod.apply_func_parallel(calc_PFC_power,'fpl',batch_data,processes=8,chunksize=1)
PFC_power=np.array(PFC_power).reshape(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers),len(batch_data.parameters__phases_upper))

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
        #plot relative phase
        ax.scatter(np.abs(mode_number[0])*(batch_data.parameters__phases_uppers_cases_all[scenario_counter][mode_number_counter]),100.*PFC_power[scenario_counter,mode_number_counter]/Pinj,
        label=f'{scenario} ({batch_data.configs_beam_species[scenario_counter]}), $n$ = {"+".join(str(abs(int(mode))) for mode in mode_number)}, U:M={int(np.abs(mode_number[0])*relative_phases_upper_middle[0])}, L:M={int(np.abs(mode_number[0])*relative_phases_lower_middle[0])}')

ax.set_xlabel('Relative rigid phase $\Delta\phi$') #\Phi for absolute
ax.set_ylabel('Deposited power lost [%]')
ax.legend(bbox_to_anchor=(.5,.3),ncol=1,loc='best')
plt.show()

#################################
 
##################################################################
 
###################################################################################################