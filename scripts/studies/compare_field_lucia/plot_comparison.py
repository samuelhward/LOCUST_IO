#plot_comparison.py
 
"""
Samuel Ward
25/07/21
----
compare input fields to Lucia/ASCOT fields
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
    from random_scripts.plot_X_point_displacement import plot_X_point_displacement
except:
    raise ImportError("ERROR: templates/template_mod.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## 
#Main 

import compare_field_lucia_launch as batch_data

perturbations=templates.plot_mod.get_input_files(batch_data,'pert')
poincare_maps=templates.plot_mod.get_output_files(batch_data,'poinc')


run_number=0
for parameters__database,parameters__sheet_name_kinetic_prof in zip(
        batch_data.parameters__databases,batch_data.parameters__sheet_names_kinetic_prof): 
    for parameters__kinetic_prof_tF_tE,parameters__kinetic_prof_Pr in zip(batch_data.parameters__kinetic_profs_tF_tE,batch_data.parameters__kinetic_profs_Pr):
        #note - different sets of toroidal mode numbers require scanning different absolute phase ranges
        for parameters__plasma_response in batch_data.parameters__plasma_responses:
            for parameters__toroidal_mode_number,parameters__i3dr in zip(batch_data.parameters__toroidal_mode_numbers,batch_data.parameters__i3drs):
                for parameters__phase_upper,parameters__phase_middle,parameters__phase_lower in zip(batch_data.parameters__phases_uppers,batch_data.parameters__phases_middles,batch_data.parameters__phases_lowers):
                    for parameters__rotation_upper,parameters__rotation_middle,parameters__rotation_lower in zip(batch_data.parameters__rotations_upper,batch_data.parameters__rotations_middle,batch_data.parameters__rotations_lower): #nest at same level == rotating them together rigidly
                        for parameters__current_upper,parameters__current_middle,parameters__current_lower in zip(batch_data.parameters__currents_upper,batch_data.parameters__currents_middle,batch_data.parameters__currents_lower):

                            perturbation=next(perturbations)
                            poincare_map=next(poincare_maps)

                            if perturbation and poincare_map:
                                
                                equilibrium=Equilibrium(ID='',data_format='GEQDSK',filename=batch_data.args_batch['RMP_study__filepath_equilibrium'][run_number],GEQDSKFIX1=True,GEQDSKFIX2=True)
                                
                                fig,axs=plt.subplots(1,2)
                                
                                perturbation.mode_number=parameters__toroidal_mode_number[0] 
                                perturbation.plot(key='dB_field_R',ax=axs[0],fig=fig,i3dr=parameters__i3dr)

                                #poincare_map.plot(ax=axs[1],fig=fig,real_scale=True)
                                templates.plot_mod.plot_poincare_psi_theta(poincare_map,equilibrium,phi=0.,ax=axs[1],fig=fig)

                                plt.title(f'i3dr={parameters__i3dr}, $n$={perturbation.mode_number}')
                                plt.show()

                            run_number+=1