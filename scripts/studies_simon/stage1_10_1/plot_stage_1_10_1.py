#plot_stage_1_10_1.py
 
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

import stage_1_10_1_launch as batch_data


Pinj=33.e6

#read input data
#cannot read input_dBs in parallel, since sometimes variable number of input filenames per simulation
input_dBs=np.array(list(templates.plot_mod.get_io_files(batch_data,'pert',return_lists=True)),dtype='object').reshape(
    len(batch_data.parameters__currents_upper),
    len(batch_data.parameters__toroidal_mode_numbers),
    len(batch_data.parameters__phases_uppers[0]),
    )
#input_dBs=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'pert',batch_data,processes=16,chunksize=1)).reshape(
#    len(batch_data.parameters__currents_upper),
#    len(batch_data.parameters__toroidal_mode_numbers),
#    len(batch_data.parameters__phases_uppers[0]),
#    )

equilibrium=Equilibrium(ID='',data_format='GEQDSK',filename=batch_data.args_batch['RMP_study__filepath_equilibrium'][0],GEQDSKFIX1=True,GEQDSKFIX2=True)
equilibrium['q_rz']=processing.utils.flux_func_to_RZ(psi=equilibrium['flux_pol'],quantity=equilibrium['qpsi'],equilibrium=equilibrium)
beam_deposition=Beam_Deposition(ID='',data_format='LOCUST_FO_weighted',filename=pathlib.Path(batch_data.args_batch['LOCUST_run__dir_cache'][0].strip('\''))/'ptcles.dat')
number_markers=8*batch_data.LOCUST_run__settings_prec_mod['threadsPerBlock']*batch_data.LOCUST_run__settings_prec_mod['blocksPerGrid']
for key in beam_deposition.data.keys():
    try:
        beam_deposition[key]=beam_deposition[key][0:number_markers]
    except:
        pass
beam_deposition.set(theta=processing.utils.angle_pol(R_major=0.639277243e1,R=beam_deposition['R'],Z=beam_deposition['Z'],Z_major=0.597384943))
beam_deposition['power']=beam_deposition['weight']*0.5*constants.mass_deuteron*(beam_deposition['V_R']**2+beam_deposition['V_phi']**2+beam_deposition['V_Z']**2)
beam_deposition.set(V_pitch_2D=processing.utils.pitch_calc(beam_deposition,equilibria=[equilibrium],perturbations=None,i3dr=-1,phase=0.))
beam_deposition.set(q=processing.utils.value_at_RZ(R=beam_deposition['R'],Z=beam_deposition['Z'],quantity=equilibrium['q_rz'],grid=equilibrium))
for key in beam_deposition.data:
    if np.any([angle in key for angle in ['theta','phi']]):
        beam_deposition[key][beam_deposition[key]>np.pi]-=2.*np.pi
        beam_deposition[key]*=180./np.pi

#read outputs
#output_fpls=np.array(templates.plot_mod.get_io_files(batch_data,'fpl')).reshape(
#    len(batch_data.parameters__currents_upper),
#    len(batch_data.parameters__toroidal_mode_numbers),
#    len(batch_data.parameters__phases_uppers[0]),
#    )
output_fpls=np.array(templates.plot_mod.apply_func_parallel(templates.plot_mod.read_locust_io_obj,'fpl',batch_data,processes=4,chunksize=1)).reshape(
    len(batch_data.parameters__currents_upper),
    len(batch_data.parameters__toroidal_mode_numbers),
    len(batch_data.parameters__phases_uppers[0]),
    )
#output_orbits_3D=templates.plot_mod.get_io_files(batch_data,'orbit3D')
#output_orbits_2D=templates.plot_mod.get_io_files(batch_data,'orbit2D')

#calculate some data
for current_counter,(parameters__current_upper,parameters__current_middle,parameters__current_lower) in enumerate(zip(batch_data.parameters__currents_upper,batch_data.parameters__currents_middle,batch_data.parameters__currents_lower)):
    for mode_counter,(parameters__toroidal_mode_number,parameters__phases_upper,parameters__phases_middle,parameters__phases_lower) in enumerate(zip(batch_data.parameters__toroidal_mode_numbers,batch_data.parameters__phases_uppers,batch_data.parameters__phases_middles,batch_data.parameters__phases_lowers)):
        for phase_counter,(parameters__phase_upper,parameters__phase_middle,parameters__phase_lower) in enumerate(zip(parameters__phases_upper,parameters__phases_middle,parameters__phases_lower)): #nest at same level == offset them together rigidly 

            if output_fpls[current_counter,mode_counter,phase_counter] and np.all(input_dBs[current_counter,mode_counter,phase_counter]):

                for mode in input_dBs[current_counter,mode_counter,phase_counter]:
                    mode.mode_number=-int(str(mode.filename)[-1])

                output_fpls[current_counter,mode_counter,phase_counter].set(theta_initial=processing.utils.angle_pol(R_major=0.639277243e1,R=output_fpls[current_counter,mode_counter,phase_counter]['R_initial'],Z=output_fpls[current_counter,mode_counter,phase_counter]['Z_initial'],Z_major=0.597384943))
                output_fpls[current_counter,mode_counter,phase_counter].set(theta_final=processing.utils.angle_pol(R_major=0.639277243e1,R=output_fpls[current_counter,mode_counter,phase_counter]['R'],Z=output_fpls[current_counter,mode_counter,phase_counter]['Z'],Z_major=0.597384943))
                output_fpls[current_counter,mode_counter,phase_counter].set(dTheta=output_fpls[current_counter,mode_counter,phase_counter]['theta_final']-output_fpls[current_counter,mode_counter,phase_counter]['theta_initial'])
                output_fpls[current_counter,mode_counter,phase_counter].set(dPhi=output_fpls[current_counter,mode_counter,phase_counter]['phi']-output_fpls[current_counter,mode_counter,phase_counter]['phi_initial'])
                output_fpls[current_counter,mode_counter,phase_counter]['power']=output_fpls[current_counter,mode_counter,phase_counter]['weight']*0.5*constants.mass_deuteron*(output_fpls[current_counter,mode_counter,phase_counter]['V_R']**2+output_fpls[current_counter,mode_counter,phase_counter]['V_phi']**2+output_fpls[current_counter,mode_counter,phase_counter]['V_Z']**2)
                #output_fpls[current_counter,mode_counter,phase_counter]['weight']=output_fpls[current_counter,mode_counter,phase_counter]['weight']*0.5*constants.mass_deuteron*(output_fpls[current_counter,mode_counter,phase_counter]['V_R']**2+output_fpls[current_counter,mode_counter,phase_counter]['V_phi']**2+output_fpls[current_counter,mode_counter,phase_counter]['V_Z']**2)
                output_fpls[current_counter,mode_counter,phase_counter].set(V_pitch_initial_2D=processing.utils.pitch_calc({key.replace('_initial',''):value for key,value in output_fpls[current_counter,mode_counter,phase_counter].data.items() if 'initial' in key},equilibria=[equilibrium],perturbations=None,i3dr=-1,phase=0.))
                output_fpls[current_counter,mode_counter,phase_counter].set(V_pitch_initial_3D=processing.utils.pitch_calc({key.replace('_initial',''):value for key,value in output_fpls[current_counter,mode_counter,phase_counter].data.items() if 'initial' in key},equilibria=[equilibrium],perturbations=[mode for mode in input_dBs[current_counter,mode_counter,phase_counter]],i3dr=-1,phase=0.))
                output_fpls[current_counter,mode_counter,phase_counter].set(V_pitch_final_2D=processing.utils.pitch_calc(output_fpls[current_counter,mode_counter,phase_counter],equilibria=[equilibrium],perturbations=None,i3dr=-1,phase=0.))
                output_fpls[current_counter,mode_counter,phase_counter].set(V_pitch_final_3D=processing.utils.pitch_calc(output_fpls[current_counter,mode_counter,phase_counter],equilibria=[equilibrium],perturbations=[mode for mode in input_dBs[current_counter,mode_counter,phase_counter]],i3dr=-1,phase=0.))
                output_fpls[current_counter,mode_counter,phase_counter].set(dV_pitch_initial=output_fpls[current_counter,mode_counter,phase_counter]['V_pitch_initial_3D']-output_fpls[current_counter,mode_counter,phase_counter]['V_pitch_initial_2D'])
                output_fpls[current_counter,mode_counter,phase_counter].set(dV_pitch_final=output_fpls[current_counter,mode_counter,phase_counter]['V_pitch_final_3D']-output_fpls[current_counter,mode_counter,phase_counter]['V_pitch_final_2D'])
                output_fpls[current_counter,mode_counter,phase_counter].set(dV_pitch=output_fpls[current_counter,mode_counter,phase_counter]['V_pitch_final_3D']-output_fpls[current_counter,mode_counter,phase_counter]['V_pitch_initial_3D'])
                output_fpls[current_counter,mode_counter,phase_counter].set(ddV_pitch=output_fpls[current_counter,mode_counter,phase_counter]['dV_pitch_final']-output_fpls[current_counter,mode_counter,phase_counter]['dV_pitch_initial'])
                output_fpls[current_counter,mode_counter,phase_counter].set(psi_norm_sqrt_initial=np.sqrt(output_fpls[current_counter,mode_counter,phase_counter]['psi_norm_initial']))
                output_fpls[current_counter,mode_counter,phase_counter].set(q_initial=processing.utils.value_at_RZ(R=output_fpls[current_counter,mode_counter,phase_counter]['R_initial'],Z=output_fpls[current_counter,mode_counter,phase_counter]['Z_initial'],quantity=equilibrium['q_rz'],grid=equilibrium))
                output_fpls[current_counter,mode_counter,phase_counter].set(q_final=processing.utils.value_at_RZ(R=output_fpls[current_counter,mode_counter,phase_counter]['R'],Z=output_fpls[current_counter,mode_counter,phase_counter]['Z'],quantity=equilibrium['q_rz'],grid=equilibrium))
                #to avoid discontinuities in the colour map since lots of markers near to phi=0 and theta=0
                for key in output_fpls[current_counter,mode_counter,phase_counter].data:
                    if np.any([angle in key for angle in ['theta','phi']]):
                        output_fpls[current_counter,mode_counter,phase_counter][key][output_fpls[current_counter,mode_counter,phase_counter][key]>np.pi]-=2.*np.pi
                        output_fpls[current_counter,mode_counter,phase_counter][key]*=180./np.pi

if __name__=='__main__':



    #find markers which lost at all phases
    markers_always_lost = set(list(np.where(np.all(np.array([output_fpl['status_flag'] for output_fpl in output_fpls])=='PFC_intercept_3D',axis=0))[0]))
    markers_sometimes_lost = set(list(np.where(np.any(np.array([output_fpl['status_flag'] for output_fpl in output_fpls])=='PFC_intercept_3D',axis=0))[0]))-markers_always_lost

    for output_fpl in output_fpls:
        output_fpl['status_flag'][list(markers_always_lost)]='always_lost'


    fig,axs = plt.subplots(3,2,subplot_kw={'projection': 'polar'})


    fig,axs = plt.subplots(3,2)

    for counter,(output_fpl,ax) in enumerate(zip(output_fpls,[ax for row in axs for ax in row])):
        output_fpl.plot(axes=['phi','theta_final'],style='scatter',status_flags=['always_lost'],fill=False,colmap=settings.cmap_inferno,colmap_val=counter/len(output_fpls),colfield='time',number_bins=200,weight=True,fig=fig,ax=ax)
        output_fpl.plot(axes=['E','V_pitch_final_3D'],style='scatter',status_flags=['always_lost'],fill=False,colmap=settings.cmap_inferno,colmap_val=counter/len(output_fpls),colfield='theta_final',number_bins=200,weight=True,fig=fig,ax=ax)


    fig,ax = plt.subplots(1)

    for counter,output_fpl in enumerate(output_fpls):
        output_fpl.plot(axes=['time'],style='histogram',fill=False,colmap=settings.cmap_inferno,colmap_val=counter/len(output_fpls),colfield='V_phi',number_bins=500,weight=True,fig=fig,ax=ax)


#################################

##################################################################

###################################################################################################