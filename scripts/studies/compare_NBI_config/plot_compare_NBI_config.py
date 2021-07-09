#plot_compare_NBI_config.py
 
"""
Samuel Ward
08/07/20
----
plot script to compare NBI geometry
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
    from classes.input_classes.equilibrium import Equilibrium
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/equilibrium.py could not be imported!\nreturning\n") 
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

import compare_NBI_config_launch as batch_data
outputs=templates.plot_mod.get_output_files(batch_data,'fpl')
fig1,ax1=plt.subplots(1)
fig2,ax2=plt.subplots(1)
fig3,ax3=plt.subplots(1)
Pinj=33.e6
PFC_power=[]

run_number=0
for parameters__database,parameters__sheet_name_kinetic_prof in zip(
        batch_data.parameters__databases,batch_data.parameters__sheet_names_kinetic_prof): 
    for parameters__kinetic_prof_tF_tE,parameters__kinetic_prof_Pr in zip(batch_data.parameters__kinetic_profs_tF_tE,batch_data.parameters__kinetic_profs_Pr):
        for parameters__toroidal_mode_number,parameters__phases_upper,parameters__phases_middle,parameters__phases_lower in zip(batch_data.parameters__toroidal_mode_numbers,batch_data.parameters__phases_uppers,batch_data.parameters__phases_middles,batch_data.parameters__phases_lowers):
            for parameters__phase_upper,parameters__phase_middle,parameters__phase_lower in zip(batch_data.parameters__phases_upper,batch_data.parameters__phases_middle,batch_data.parameters__phases_lower): #nest at same level == offset them together rigidly 
                for parameters__rotation_upper,parameters__rotation_middle,parameters__rotation_lower in zip(batch_data.parameters__rotations_upper,batch_data.parameters__rotations_middle,batch_data.parameters__rotations_lower): #nest at same level == rotating them together rigidly
                    for parameters__current_upper,parameters__current_middle,parameters__current_lower in zip(batch_data.parameters__currents_upper,batch_data.parameters__currents_middle,batch_data.parameters__currents_lower):

                        output=next(outputs)
                        if output: 
                            if 'B3D_EX' not in batch_data.args_batch['LOCUST_run__flags'][run_number]: #this is the 2D case
                                line_style='dashed'
                                label='axisymmetric'
                            else:
                                line_style='solid'                                
                                label=parameters__toroidal_mode_number
                            output.plot(fig=fig1,ax=ax1,axes=['time'],fill=False,label=label,number_bins=200,weight=False,line_style=line_style)
                            output['E']/=1000. #convert to keV
                            output.plot(fig=fig2,ax=ax2,axes=['E'],fill=False,label=label,number_bins=200,weight=False,line_style=line_style)
                            #calculate losses from particle list
                            i=np.where(output['status_flag']=='PFC_intercept_3D')[0]
                            power_loss=1.e6*output['f']*np.sum((output['V_R'][i]**2+output['V_phi'][i]**2+output['V_Z'][i]**2)*output['FG'][i])*0.5*constants.mass_deuteron
                            PFC_power.append([power_loss/1.e6])
                        else:
                            PFC_power.append([-10.])

                        run_number+=1

ax1.legend()
ax2.legend()
ax1.set_xlabel('Time [s]')
ax2.set_xlabel('Energy [keV]')
ax1.set_ylabel('losses')
ax2.set_ylabel('losses')
ax1.set_title('')
ax2.set_title('')
plt.show()
        
#################################
 
##################################################################
 
###################################################################################################