#plot_stage_4.py
 
"""
Samuel Ward
25/05/21
----
script for plotting stage 4 of Simon's studies 
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

import stage_4_launch as batch_data

outputs=templates.plot_mod.get_output_files(batch_data,'fpl')

Pinj=33.e6
PFC_power=[]
for output in outputs:
    if output:
        i=np.where(output['status_flag']=='PFC_intercept_3D')[0]
        PFC_power.append([1.e6*output['f']*np.sum((output['V_R'][i]**2+output['V_phi'][i]**2+output['V_Z'][i]**2)*output['FG'][i])*0.5*constants.mass_deuteron])
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

fig,ax=plt.subplots(1)
for scenario_counter,scenario in enumerate(batch_data.parameters__databases):
    for mode_number_counter,mode_number in enumerate(batch_data.parameters__toroidal_mode_numbers):
        #find relative phase between coils for this toroidal mode number
        relative_phases_upper_middle=batch_data.parameters__phases_uppers[mode_number_counter]-batch_data.parameters__phases_middles[mode_number_counter]
        relative_phases_lower_middle=batch_data.parameters__phases_lowers[mode_number_counter]-batch_data.parameters__phases_middles[mode_number_counter]
        #plot relative phase
        ax.scatter(np.abs(mode_number[0])*(batch_data.parameters__phases_middles[mode_number_counter]-26.7),100.*PFC_power[scenario_counter,mode_number_counter]/33.e6,label=f'{scenario}, $n$ = {mode_number}')

ax.set_xlabel('Absolute rigid phase $\Delta\Phi$',fontsize=20)
ax.set_ylabel('PFC power flux [%P_dep]',fontsize=20)
ax.legend(fontsize=10)
plt.show()

#################################
 
##################################################################
 
###################################################################################################