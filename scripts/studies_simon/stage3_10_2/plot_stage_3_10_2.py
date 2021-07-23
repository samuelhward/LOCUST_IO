#plot_stage_3_10_2.py
 
"""
Samuel Ward
16/07/21
----
script for plotting stage 3_10_2 of Simon's studies 
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
    from matplotlib.animation import FuncAnimation
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
    from random_scripts.plot_X_point_displacement import plot_X_point_displacement
except:
    raise ImportError("ERROR: templates/template_mod.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## 
#Main 

import stage_3_10_2_launch as batch_data
total_runs=len(batch_data.args_batch['LOCUST_run__flags'])





Pinj=33.e6

#fig=plt.figure()
#ax=fig.add_subplot(111,polar=False)
outputs=templates.plot_mod.get_output_files(batch_data,'fpl')
for parameters__database,parameters__sheet_name_kinetic_prof in zip(
        batch_data.parameters__databases,batch_data.parameters__sheet_names_kinetic_prof): 
    for parameters__kinetic_prof_tF_tE,parameters__kinetic_prof_Pr in zip(batch_data.parameters__kinetic_profs_tF_tE,batch_data.parameters__kinetic_profs_Pr):
        #note - different sets of toroidal mode numbers require scanning different absolute phase ranges
        for parameters__toroidal_mode_number,parameters__phases_upper,parameters__phases_middle,parameters__phases_lower in zip(batch_data.parameters__toroidal_mode_numbers,batch_data.parameters__phases_uppers,batch_data.parameters__phases_middles,batch_data.parameters__phases_lowers):

            point=0
            number_points=len(parameters__phases_lower)    
            #fig,ax=plt.subplots(number_points)
            #ax=np.array(ax,ndmin=1)
            for parameters__phase_upper,parameters__phase_middle,parameters__phase_lower in zip(parameters__phases_upper,parameters__phases_middle,parameters__phases_lower): #nest at same level == offset them together rigidly 
                for parameters__rotation_upper,parameters__rotation_middle,parameters__rotation_lower in zip(batch_data.parameters__rotations_upper,batch_data.parameters__rotations_middle,batch_data.parameters__rotations_lower): #nest at same level == rotating them together rigidly
                    for parameters__current_upper,parameters__current_middle,parameters__current_lower in zip(batch_data.parameters__currents_upper,batch_data.parameters__currents_middle,batch_data.parameters__currents_lower):

                        output=next(outputs)
                        if output:
                            output.set(theta=processing.utils.angle_pol(R_major=0.639277243e1,R=output['R'],Z=output['Z'],Z_major=0.597384943))
                            output['weight']*=0.5*constants.mass_deuteron*(output['V_R']**2+output['V_phi']**2+output['V_Z']**2)/Pinj
                            output['phi']*=180/np.pi
                            output['theta']*=180/np.pi
                            output.plot(axes=['phi','Z'],style='scatter',fill=False,label=str(point),colmap=settings.cmap_default,colfield='V_phi',number_bins=200,weight=True)#XXX ,fig=fig,ax=ax[point])
                        ax[point].set_xlim(-180,180)
                        #ax[point].set_ylim(240,270) #theta
                        ax[point].set_title(f'point {point}')
                        point+=1
            #XXX plt.show()

#ax.set_xlabel('R [m]',fontsize=25)  
#ax.set_xlabel('Energy [keV]',fontsize=25)  
#ax.set_ylabel('marker loss fraction',fontsize=25)  
#ax.set_ylabel('marker loss fraction',fontsize=25)  
#ax.set_title('')
#ax.set_title('')
#ax.legend()
#ax.legend()






outputs=templates.plot_mod.get_output_files(batch_data,'dfn')
for parameters__database,parameters__sheet_name_kinetic_prof in zip(
        batch_data.parameters__databases,batch_data.parameters__sheet_names_kinetic_prof): 
    for parameters__kinetic_prof_tF_tE,parameters__kinetic_prof_Pr in zip(batch_data.parameters__kinetic_profs_tF_tE,batch_data.parameters__kinetic_profs_Pr):
        #note - different sets of toroidal mode numbers require scanning different absolute phase ranges
        for parameters__toroidal_mode_number,parameters__phases_upper,parameters__phases_middle,parameters__phases_lower in zip(batch_data.parameters__toroidal_mode_numbers,batch_data.parameters__phases_uppers,batch_data.parameters__phases_middles,batch_data.parameters__phases_lowers):

            point=0
            number_points=len(parameters__phases_lower)            
            fig,ax=plt.subplots(1)
            
            for parameters__phase_upper,parameters__phase_middle,parameters__phase_lower in zip(parameters__phases_upper,parameters__phases_middle,parameters__phases_lower): #nest at same level == offset them together rigidly 
                for parameters__rotation_upper,parameters__rotation_middle,parameters__rotation_lower in zip(batch_data.parameters__rotations_upper,batch_data.parameters__rotations_middle,batch_data.parameters__rotations_lower): #nest at same level == rotating them together rigidly
                    for parameters__current_upper,parameters__current_middle,parameters__current_lower in zip(batch_data.parameters__currents_upper,batch_data.parameters__currents_middle,batch_data.parameters__currents_lower):

                        point+=1
                        output=next(outputs)
                        if output:
                            colmap=settings.colour_custom(list(settings.cmap_inferno(point/number_points,bytes=True))[0:3]+[1])
                            output.plot(axes=['V_pitch'],fill=False,label=str(point),colmap=colmap,fig=fig,ax=ax)

            plt.show()





outputs=templates.plot_mod.get_output_files(batch_data,'fpl')
def plot(frame,fig,ax): 
    output,run_number=frame
    run_number-=1

    n=batch_data.args_batch['MARS_read__flags'][run_number]['N0']
    #phi_U=float(batch_data.args_batch['MARS_read__flags'][run_number]['UPHASE']).replace('D0')
    #phi_M=float(batch_data.args_batch['MARS_read__flags'][run_number]['MPHASE']).replace('D0')
    #phi_L=float(batch_data.args_batch['MARS_read__flags'][run_number]['LPHASE']).replace('D0')

    try:
        if output:   
            ax.cla()
            output.set(theta=processing.utils.angle_pol(R_major=0.639277243e1,R=output['R'],Z=output['Z'],Z_major=0.597384943))
            output.plot(axes=['phi','R'],style='scatter',fill=False,label=str(run_number-2),colmap=settings.cmap_default,colfield='V_phi',number_bins=200,weight=True,fig=fig,ax=ax)
            ax.title(f'$n$ = {n}, Point {run_number-2}')
            if 'B3D_EX' not in batch_data.args_batch['LOCUST_run__flags'][run_number]: #axisymmetric run
                ax.title(f'2D field')
        else:
            ax.cla()
    except:
        pass

fig=plt.figure()
ax=fig.add_subplot(111,polar=True)
animation=FuncAnimation(fig,plot,frames=zip(outputs,list(range(total_runs))),fargs=[fig,ax],repeat=True,interval=1)
plt.show()



#################################
 
##################################################################
 
###################################################################################################