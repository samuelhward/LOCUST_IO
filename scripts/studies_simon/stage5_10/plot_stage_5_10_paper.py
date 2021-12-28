#plot_stage_5_10_paper.py
 
"""
Samuel Ward
29/06/21
----
script for plotting stage 5.10 of Simon's studies 
used to create plot for paper 3
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

import stage_5_10_launch as batch_data

method='rund'

Pinj=33.e6
PFC_power=templates.plot_mod.apply_func_parallel(templates.plot_mod.calc_PFC_power,method,batch_data,processes=32,chunksize=1)
PFC_power=np.array(PFC_power).reshape(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers),len(batch_data.parameters__phases_upper))

if method=='fpl':
    PFC_power/=constants.mass_deuteron #different masses in some cases, so divide by mass here for re-multiplication later
    for scenario_counter,(scenario,beam_species) in enumerate(zip(batch_data.parameters__databases,batch_data.configs_beam_species)):
        if beam_species=='hydrogen':
            PFC_power[scenario_counter]*=scipy.constants.physical_constants['proton mass'][0]
        elif beam_species=='deuterium':
            PFC_power[scenario_counter]*=scipy.constants.physical_constants['deuteron mass'][0]

fig,ax=plt.subplots(1)
for plasma_state_counter,(scenario,beam_species) in enumerate(zip(batch_data.parameters__databases,batch_data.configs_beam_species)):
    for mode_number_counter,(mode_number) in enumerate(batch_data.parameters__toroidal_mode_numbers):
        label=f'{scenario.split("_")[-1]}'+f' ({batch_data.configs_beam_species[plasma_state_counter][0].capitalize()})'
        label+=f', $n$ = {"+".join(str(abs(int(mode))) for mode in mode_number)}'
        ax.plot(np.arange(len(PFC_power[plasma_state_counter][mode_number_counter])),100*PFC_power[plasma_state_counter][mode_number_counter]/Pinj,label=label)


#now add 1_10_1 data to plot (or could re-run 1_10 to fill in missing simulations...it doesn't matter - but 1_10_1 data are huge! so think about editing the launch file...)
del(batch_data)
import studies_simon.stage3_10.stage_3_10_launch as batch_data

PFC_power_3_10=templates.plot_mod.apply_func_parallel(templates.plot_mod.calc_PFC_power,method,batch_data,processes=32,chunksize=1)
PFC_power_3_10=np.array(PFC_power_3_10).reshape(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers),len(batch_data.parameters__phases_upper))
for plasma_state_counter,(scenario,beam_species) in enumerate(zip(batch_data.parameters__databases,batch_data.configs_beam_species)):
    for mode_number_counter,mode_number in enumerate(batch_data.parameters__toroidal_mode_numbers):
        label=f'{scenario.split("_")[-1]}'+f' ({batch_data.configs_beam_species[plasma_state_counter][0].capitalize()})'
        label+=f', $n$ = {"+".join(str(abs(int(mode))) for mode in mode_number)}'
        ax.plot(np.arange(len(PFC_power_3_10[plasma_state_counter][mode_number_counter])),100*PFC_power_3_10[plasma_state_counter][mode_number_counter]/Pinj,label=label)

ax.set_xlabel('XPD contour point')
ax.set_ylabel('NBI power loss [%]')
ax.set_ylim([0,8.1])
ax.legend(loc='center',bbox_to_anchor=(1.1,0.5),ncol=1)
plt.show()

#################################
 
##################################################################
 
###################################################################################################