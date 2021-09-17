#plot_stage_1_10_2.py
 
"""
Samuel Ward
23/07/21
----
script for plotting stage 1.10.1 of Simon's studies 
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
    from classes.input_classes.beam_deposition import Beam_Deposition
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/perturbation.py could not be imported!\nreturning\n") 
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
except:
    raise ImportError("ERROR: templates/template_mod.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## 
#Main 

cmap_r=settings.colour_custom([194,24,91,1])
cmap_g=settings.colour_custom([76,175,80,1])
cmap_b=settings.colour_custom([33,150,243,1])

import stage_1_10_2_launch as batch_data

#outputs=templates.plot_mod.get_output_files(batch_data,'fpl')

Pinj=33.e6


#read input data
input_dBs=np.array(list(templates.plot_mod.get_output_files(batch_data,'pert'))).reshape(len(batch_data.parameters__phases_upper),len(batch_data.parameters__toroidal_mode_numbers[0]))
equilibrium=Equilibrium(ID='',data_format='GEQDSK',filename=batch_data.args_batch['RMP_study__filepath_equilibrium'][0],GEQDSKFIX1=True,GEQDSKFIX2=True)
beam_deposition=Beam_Deposition(ID='',data_format='LOCUST_FO_weighted',filename=pathlib.Path(batch_data.args_batch['LOCUST_run__dir_cache'][0].strip('\''))/'ptcles.dat')
number_markers=8*batch_data.LOCUST_run__settings_prec_mod['threadsPerBlock']*batch_data.LOCUST_run__settings_prec_mod['blocksPerGrid']
for key in beam_deposition.data.keys():
    try:
        beam_deposition[key]=beam_deposition[key][0:number_markers]
    except:
        pass

equilibrium['q_rz']=processing.utils.flux_func_to_RZ(psi=equilibrium['flux_pol'],quantity=equilibrium['qpsi'],equilibrium=equilibrium)

#read outputs
output_fpls=list(templates.plot_mod.get_output_files(batch_data,'fpl'))

#calculate some data
for output_fpl,input_dB in zip(output_fpls,input_dBs):
    if output_fpl and np.all(input_dB):
        #set toroidal mode number
        for mode in input_dB:
            mode.mode_number=-int(str(mode.filename)[-1])    
        output_fpl.set(theta_initial=processing.utils.angle_pol(R_major=0.639277243e1,R=output_fpl['R_initial'],Z=output_fpl['Z_initial'],Z_major=0.597384943))
        output_fpl.set(theta_final=processing.utils.angle_pol(R_major=0.639277243e1,R=output_fpl['R'],Z=output_fpl['Z'],Z_major=0.597384943))
        output_fpl.set(dTheta=output_fpl['theta_final']-output_fpl['theta_initial'])
        output_fpl.set(dPhi=output_fpl['phi']-output_fpl['phi_initial'])
        output_fpl['power']=output_fpl['weight']*0.5*constants.mass_deuteron*(output_fpl['V_R']**2+output_fpl['V_phi']**2+output_fpl['V_Z']**2)
        #output_fpl['weight']=output_fpl['weight']*0.5*constants.mass_deuteron*(output_fpl['V_R']**2+output_fpl['V_phi']**2+output_fpl['V_Z']**2)
        output_fpl.set(V_pitch_initial_2D=processing.utils.pitch_calc({key.replace('_initial',''):value for key,value in output_fpl.data.items() if 'initial' in key},equilibria=[equilibrium],perturbations=None,i3dr=-1,phase=0.))
        output_fpl.set(V_pitch_initial_3D=processing.utils.pitch_calc({key.replace('_initial',''):value for key,value in output_fpl.data.items() if 'initial' in key},equilibria=[equilibrium],perturbations=[mode for mode in input_dB],i3dr=-1,phase=0.))
        output_fpl.set(V_pitch_final_2D=processing.utils.pitch_calc(output_fpl,equilibria=[equilibrium],perturbations=None,i3dr=-1,phase=0.))
        output_fpl.set(V_pitch_final_3D=processing.utils.pitch_calc(output_fpl,equilibria=[equilibrium],perturbations=[mode for mode in input_dB],i3dr=-1,phase=0.))
        output_fpl.set(dV_pitch_initial=output_fpl['V_pitch_initial_3D']-output_fpl['V_pitch_initial_2D'])
        output_fpl.set(dV_pitch_final=output_fpl['V_pitch_final_3D']-output_fpl['V_pitch_final_2D'])
        output_fpl.set(dV_pitch=output_fpl['V_pitch_final_3D']-output_fpl['V_pitch_initial_3D'])
        output_fpl.set(ddV_pitch=output_fpl['dV_pitch_final']-output_fpl['dV_pitch_initial'])
        output_fpl.set(psi_norm_sqrt_initial=np.sqrt(output_fpl['psi_norm_initial']))
        output_fpl.set(q_initial=processing.utils.value_at_RZ(R=output_fpl['R_initial'],Z=output_fpl['Z_initial'],quantity=equilibrium['q_rz'],grid=equilibrium))
        #to avoid discontinuities in the colour map since lots of markers near to phi=0 and theta=0
        for key in output_fpl.data:
            if np.any([angle in key for angle in ['theta','phi']]):
                output_fpl[key][output_fpl[key]>np.pi]-=2.*np.pi

if __name__=='__main__':

    run_number=0
    for parameters__database,parameters__sheet_name_kinetic_prof in zip(
            parameters__databases,parameters__sheet_names_kinetic_prof): 
        for parameters__kinetic_prof_tF_tE,parameters__kinetic_prof_Pr in zip(parameters__kinetic_profs_tF_tE,parameters__kinetic_profs_Pr):
            
            '''
            equilibrium=Equilibrium(ID='',data_format='GEQDSK',filename=batch_data.args_batch['RMP_study__filepath_equilibrium'][run_number],GEQDSKFIX1=True,GEQDSKFIX2=True)
            beam_deposition=Beam_Deposition(ID='',data_format='LOCUST_FO_weighted',filename=pathlib.Path(batch_data.args_batch['LOCUST_run__dir_cache'][run_number].strip('\''))/'ptcles.dat')
            #always generate more markers that used, so crop the beam deposition down
            number_markers=8*batch_data.LOCUST_run__settings_prec_mod['threadsPerBlock']*batch_data.LOCUST_run__settings_prec_mod['blocksPerGrid']
            for key in beam_deposition.data.keys():
                try:
                    beam_deposition[key]=beam_deposition[key][0:number_markers]
                except:
                    pass
            '''

            for parameters__rotation_upper,parameters__rotation_middle,parameters__rotation_lower in zip(batch_data.parameters__rotations_upper,batch_data.parameters__rotations_middle,batch_data.parameters__rotations_lower): #nest at same level == rotating them together rigidly
                for parameters__current_upper,parameters__current_middle,parameters__current_lower in zip(batch_data.parameters__currents_upper,batch_data.parameters__currents_middle,batch_data.parameters__currents_lower):
                    for parameters__toroidal_mode_number,parameters__phases_upper,parameters__phases_middle,parameters__phases_lower in zip(batch_data.parameters__toroidal_mode_numbers,batch_data.parameters__phases_uppers,parameters__phases_middles,batch_data.parameters__phases_lowers):

                        number_points=len(batch_data.parameters__phases_lower)    
                        #fig,ax=plt.subplots(number_points)
                        #ax=np.array(ax,ndmin=1)

                        for parameters__phase_upper,parameters__phase_middle,parameters__phase_lower in zip(batch_data.parameters__phases_upper,batch_data.parameters__phases_middle,batch_data.parameters__phases_lower): #nest at same level == offset them together rigidly 

                            output_fpls[run_number].plot(axes=['phi','Z'],style='scatter',fill=False,colmap=settings.cmap_default,colfield='V_phi',number_bins=200,weight=True)
                        
                            run_number+=1

                        plt.show()

    #find markers quickly lost
    for output_fpl in output_fpls:
        markers_prompt_loss,=np.where((output_fpl['time']<=0.005) & (output_fpl['status_flag']=='PFC_intercept_3D'))
        output_fpl['status_flag'][list(markers_prompt_loss)]='quickly_lost'

    #find markers lost to wall panels

    #find markers that started inside/outside psi_n=0.9
    for output_fpl in output_fpls:
        markers_start_core,=np.where((output_fpl['psi_initial_norm']<=0.9) & (output_fpl['status_flag']=='PFC_intercept_3D'))
        markers_start_edge,=np.where((output_fpl['psi_initial_norm']>0.9) & (output_fpl['status_flag']=='PFC_intercept_3D'))
        output_fpl['status_flag'][list(markers_start_edge)]='edge_born'

    #find markers which lost at all phases
    markers_lost_always = set(list(np.where(np.all(np.array([output_fpl['status_flag'] for output_fpl in output_fpls])=='PFC_intercept_3D',axis=0))[0]))
    markers_lost_sometimes = set(list(np.where(np.any(np.array([output_fpl['status_flag'] for output_fpl in output_fpls])=='PFC_intercept_3D',axis=0))[0]))-markers_lost_always
    #reset flags of markers always lost
    #can infer those only sometimes lost as retaining the PFC_intercept_3D flag
    #do not want to overwrite markers_lost_sometimes since we lose information as to which phase they're lost at
    for output_fpl in output_fpls:
        output_fpl['status_flag'][list(markers_lost_always)]='always_lost'

    #find max-min variations throughout phases
    for quantity in ['time','R','Z','phi','theta_final','V_pitch_final_3D']:
        for output_fpl in output_fpls:
            output_fpl[f'd{quantity}']=np.max(np.array([output_fpl[quantity] for output_fpl in output_fpls],),axis=0)-np.min(np.array([output_fpl[quantity] for output_fpl in output_fpls],),axis=0)


    fig,axs = plt.subplots(3,2,subplot_kw={'projection': 'polar'})

    for counter,(ax,status_flag) in enumerate(zip(axs,['always_lost','PFC_intercept_3D'])):

        output_fpl.plot(axes=['phi','R'],style='scatter',status_flags=['always_lost'],fill=False,colmap=settings.cmap_inferno,colmap_val=counter/len(output_fpls),colfield='time',number_bins=200,weight=True,fig=fig,ax=ax)


    fig,axs = plt.subplots(2)

    for counter,(ax,status_flag) in enumerate(zip(axs,['always_lost','PFC_intercept_3D'])):

        output_fpls[0].plot(axes=['phi','theta_final'],style='histogram',status_flags=[status_flag],fill=False,colmap=settings.cmap_inferno,colmap_val=counter/len(output_fpls),colfield='time',number_bins=200,weight=True,fig=fig,ax=ax)


    fig,axs = plt.subplots(3,2)

    for counter,(output_fpl,ax) in enumerate(zip(output_fpls,[ax for row in axs for ax in row])):

        output_fpl.plot(axes=['theta_final'],style='histogram',status_flags=['always_lost'],fill=False,colmap=settings.cmap_inferno,colmap_val=0.8,colfield='time',number_bins=200,weight=True,fig=fig,ax=ax)
        output_fpl.plot(axes=['theta_final'],style='histogram',status_flags=['PFC_intercept_3D'],fill=False,colmap=settings.cmap_inferno,colmap_val=0.2,colfield='time',number_bins=200,weight=True,fig=fig,ax=ax)

        output_fpl.plot(axes=['phi','theta_final'],style='scatter',status_flags=['always_lost'],fill=False,colmap=settings.cmap_inferno,colmap_val=counter/len(output_fpls),colfield='time',number_bins=200,weight=True,fig=fig,ax=ax)

        output_fpl.plot(axes=['E','V_pitch_final_3D'],style='scatter',status_flags=['always_lost'],fill=False,colmap=settings.cmap_inferno,colmap_val=counter/len(output_fpls),colfield='theta_final',number_bins=200,weight=True,fig=fig,ax=ax)




    fig,ax = plt.subplots(1)

    for counter,output_fpl in enumerate(output_fpls):
        output_fpl.plot(axes=['time'],style='histogram',status_flags=['always_lost'],fill=False,colmap=settings.cmap_inferno,colmap_val=counter/len(output_fpls),colfield='V_phi',number_bins=1000,weight=True,fig=fig,ax=ax)

        output_fpl.plot(axes=['theta_final'],style='histogram',status_flags=['PFC_intercept_3D'],fill=False,colmap=settings.cmap_inferno,colmap_val=counter/len(output_fpls),colfield='time',number_bins=200,weight=True,fig=fig,ax=ax)


#################################

##################################################################

###################################################################################################