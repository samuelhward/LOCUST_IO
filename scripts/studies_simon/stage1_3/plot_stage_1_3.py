#plot_stage_1_3.py
 
"""
Samuel Ward
22/06/21
----
script for plotting stage 1.3 of Simon's studies 
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

import stage_1_3_launch as batch_data

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

PFC_power=np.array(PFC_power).reshape(len(batch_data.parameters__kinetic_profs_Pr),len(batch_data.parameters__toroidal_mode_numbers),len(batch_data.parameters__phases_upper))

# one figure per plasma
fig,axs=plt.subplots(len(batch_data.parameters__kinetic_profs_Pr),constrained_layout=True)
colours=['r','b'] #one colour per mode number

for plasma_state_counter,(ax,power,Pr,tftE) in enumerate(zip([axs],PFC_power,batch_data.parameters__kinetic_profs_Pr,batch_data.parameters__kinetic_profs_tF_tE)):
    for mode_number_counter,(mode_number,colour) in enumerate(zip(batch_data.parameters__toroidal_mode_numbers,colours)):
        #find relative phase between coils for this toroidal mode number
        relative_phases_upper_middle=batch_data.parameters__phases_uppers[mode_number_counter]-batch_data.parameters__phases_middles[mode_number_counter]
        relative_phases_lower_middle=batch_data.parameters__phases_lowers[mode_number_counter]-batch_data.parameters__phases_middles[mode_number_counter]
        ax.scatter(np.abs(mode_number[0])*(batch_data.parameters__phases_middles[mode_number_counter]),100*power[mode_number_counter]/Pinj,label=f'$n$ = {mode_number}, U:M={int(np.abs(mode_number[0])*relative_phases_upper_middle[0])}, L:M={int(np.abs(mode_number[0])*relative_phases_lower_middle[0])}',color=colour)

ax.legend()
ax.set_xlabel('Relative rigid phase $\Delta\phi$',fontsize=20) #\Phi for absolute
ax.set_ylabel('Deposited power lost [%]',fontsize=20)
# remove ticks from total ax
#ax_total = fig.add_subplot(111,frameon=False)
#ax_total.tick_params(axis='both',which='both',bottom=False,labelbottom=False,left=False,labelleft=False)

plt.show()

#################################
 
##################################################################
 
###################################################################################################