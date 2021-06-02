#plot_stage_3.py
 
"""
Samuel Ward
02/06/21
----
script for plotting stage 3 of Simon's studies 
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
    from random_scripts.plot_X_point_displacement import plot_X_point_displacement
except:
    raise ImportError("ERROR: templates/template_mod.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## 
#Main 

import stage_3_launch as batch_data

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

# one figure per plasma
fig,axs=plt.subplots(len(batch_data.parameters__databases)*len(batch_data.parameters__toroidal_mode_numbers),2,constrained_layout=True)
axs=np.array([axs]).reshape(len(batch_data.parameters__databases),len(batch_data.parameters__toroidal_mode_numbers),2)
colours=['r','b'] #one colour per mode number

for plasma_state_counter,(row,database) in enumerate(zip(axs,batch_data.parameters__databases)):
    for mode_number_counter,(subrow,mode_number,colour) in enumerate(zip(row,batch_data.parameters__toroidal_mode_numbers,colours)):
        case=int(batch_data.parameters__databases[plasma_state_counter][-1])
        n=int(np.abs(mode_number[0]))
        ax_left,ax_right=subrow        
        plot_X_point_displacement(case=case,n=n,LVV=90,fig=fig,ax=ax_left)
        ax_right.scatter(np.arange(len(PFC_power[plasma_state_counter][mode_number_counter])),100*PFC_power[plasma_state_counter][mode_number_counter]/Pinj)
        ax_right.set_xlabel('Point',fontsize=15)
        ax_right.set_ylabel('Deposited power lost [%]',fontsize=15)

# remove ticks from total ax
#ax_total = fig.add_subplot(111,frameon=False)
#ax_total.tick_params(axis='both',which='both',bottom=False,labelbottom=False,left=False,labelleft=False)

plt.show()

#################################
 
##################################################################
 
###################################################################################################