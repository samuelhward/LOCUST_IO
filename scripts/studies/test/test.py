#test.py
 
"""
Samuel Ward
18/05/20
----
test project
---
 
notes:         
    look at template_mod for additional variables!
    launches the workflow defined in template_run.py
    since only flat arrays are passed to Batch class, __ labels are designed to help distiniguish what variables are needed for 

    directory structure goes as:
        LOCUST_IO
            data
                input_files
                    path_to_master_data - unchanging
                    path_to_additional_data - unchanging
                    database_folder_name
                        RMP_study
                            metadata_folder_name - stores all the inputs just before a run is performed, deleted after a run
                output_files
                    database_folder_name
                        RMP_study
                            parameter_folder_name - stores all the outputs permanently
                cache_files
                    database_folder_name - store parameter-agnostic cache here
                        parameter_folder_name - store parameter-specific cache here

        /tmp/<uname>/ (LOCUST root dir)
            InputFiles
            OutputFiles
            CacheFiles

    assumes PureVac fields are renamed to:
        mv BPLASMA_MARSF_n3_LV BPLASMA_MARSF_n3_cL_V  
        mv BPLASMA_MARSF_n3_UV BPLASMA_MARSF_n3_cU_V  
        mv BPLASMA_MARSF_n6_MV BPLASMA_MARSF_n6_cM_V
        mv BPLASMA_MARSF_n3_MV BPLASMA_MARSF_n3_cM_V  
        mv BPLASMA_MARSF_n6_LV BPLASMA_MARSF_n6_cL_V
        mv BPLASMA_MARSF_n6_UV BPLASMA_MARSF_n6_cU_V
---
"""

###################################################################################################
#Preamble
 
import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway

try:
    import numpy as np
    import pathlib
    import copy
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
    from run_scripts.batch import Batch
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/batch.py could not be imported!\nreturning\n")
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
    from templates.template_mod import *
except:
    raise ImportError("ERROR: templates/template_mod.py could not be imported!\nreturning\n") 
    sys.exit(1)


################################################################## 
#Main 

#################################
#define study name 

RMP_study__name='test'

#################################
#define options and dispatch tables for helping choosing settings

parameters__toroidal_mode_numbers__options={} #XXX needs studying - toroidal mode number combinations 
parameters__toroidal_mode_numbers__options['n=3']=[-3,-6]
parameters__toroidal_mode_numbers__options['n=4']=[-4,-5]

##################################################################
#choose the scenarios we will want to examine

parameters__databases=['ITER_15MAQ10_case5'] #these all zipped at same level in main loop because source data is not consistent enough
parameters__sheet_names_kinetic_prof=["'Flat n'"]

##################################################################
#define the parameter space for a given scenario

#kinetic profile parameters which vary independently
parameters__kinetic_profs_Pr=[1.] #highest rotation
parameters__kinetic_profs_tF_tE=[0.5]

#3D field parameters which vary independently - if you want to vary these together then put them into the same loop nesting below
#2D arrays, each element has length = number of modes
parameters__toroidal_mode_numbers=[[-3,-6]]
parameters__phases_upper=np.linspace(0,120,9)[:-1]+30. #86,0,34 = default for maximmum stochasticity in coil coordinate system
parameters__phases_middle=np.linspace(0,120,9)[:-1]+26.7
parameters__phases_lower=np.linspace(0,120,9)[:-1]+30.
parameters__rotations_upper=np.array([0.])
parameters__rotations_middle=np.array([0.])
parameters__rotations_lower=np.array([0.])
parameters__currents_upper=np.array([90.])*1000.
parameters__currents_middle=np.array([90.])*1000.
parameters__currents_lower=np.array([90.])*1000.

config_beam_1='on'
config_beam_2='on'

##################################################################
#define the workflow commands in order we want to execute them

RMP_study__workflow_commands="\"['mkdir','save_args','kin_get','3D_get','3D_calc','input_get','IDS_create','kin_extrap','run_BBNBI','depo_get','run_LOCUST','clean_input']\""

##################################################################
#create every valid combination of parameter, returned in flat lists
#use zip and nest levels to define specific combinations which cannot be varied

run_number=0
parameter_strings=[]
#first level are the data which remain constant for a parameter scan
for parameters__database,parameters__sheet_name_kinetic_prof in zip(
        parameters__databases,parameters__sheet_names_kinetic_prof): 
    for parameters__kinetic_prof_tF_tE,parameters__kinetic_prof_Pr in zip(parameters__kinetic_profs_tF_tE,parameters__kinetic_profs_Pr):
        for parameters__toroidal_mode_number in parameters__toroidal_mode_numbers:
            for parameters__phase_upper,parameters__phase_middle,parameters__phase_lower in zip(parameters__phases_upper,parameters__phases_middle,parameters__phases_lower): #nest at same level == offset them together rigidly 
                for parameters__rotation_upper,parameters__rotation_middle,parameters__rotation_lower in zip(parameters__rotations_upper,parameters__rotations_middle,parameters__rotations_lower): #nest at same level == rotating them together rigidly
                    for parameters__current_upper,parameters__current_middle,parameters__current_lower in zip(parameters__currents_upper,parameters__currents_middle,parameters__currents_lower):

                        run_number+=1 #increment run counter              
                                                   
                        #create a string of variables identifying this run
                        parameters__kinetic_prof_tF_tE_string=parameters__kinetic_profs_tF_tE__dispatch[parameters__kinetic_prof_tF_tE] #generate some variable string equivalents for later
                        parameters__kinetic_prof_Pr_string=parameters__kinetic_profs_Pr__dispatch[parameters__kinetic_prof_Pr]

                        parameters__parameter_string=''
                        parameters__parameter_string+='_'.join(['{}_{}'.format(parameter,str(value)) for parameter,value in zip([
                                        'tFtE',
                                        'Pr'
                                        ],[
                                        parameters__kinetic_prof_tF_tE,
                                        parameters__kinetic_prof_Pr])])

                        parameters__parameter_string+='_ntor_'
                        for mode in parameters__toroidal_mode_number:
                            parameters__parameter_string+='{}_'.format(str(mode)) #add toroidal mode information    
                        parameters__parameter_string+='_'.join(['{}_{}'.format(parameter,str(value)) for parameter,value in zip([
                                'phaseu',
                                'phasem',
                                'phasel',
                                'rotu',
                                'rotm',
                                'rotl',
                                'ikatu',
                                'ikatm',
                                'ikatl'],[
                                parameters__phase_upper,
                                parameters__phase_middle,
                                parameters__phase_lower,
                                parameters__rotation_upper,
                                parameters__rotation_middle,
                                parameters__rotation_lower,
                                parameters__current_upper,
                                parameters__current_middle,
                                parameters__current_lower])])

                        parameters__parameter_string+=f'_beams_{str(config_beam_1)}_{str(config_beam_2)}'

                        parameter_strings.append(parameters__parameter_string)

                        #################################
                        #define corresponding workflow args passed to batch (denoted wth __batch)

                        #run-specific settings

                        
                        LOCUST_run__flags=LOCUST_run__flags_default
                        LOCUST_run__flags['TIMAX']='2.0D0' #increase max tracking time
                        #XXX CURRENTLY WAITING FOR FIX LOCUST_run__flags['I3DR']=-1 
                        LOCUST_run__settings_prec_mod={}
                        LOCUST_run__settings_prec_mod['nmde']=len(parameters__toroidal_mode_number) #number of total toroidal harmonics = number of modes
                        LOCUST_run__settings_prec_mod['Ab']='AD' 
                        LOCUST_run__settings_prec_mod['Zb']='+1.0_gpu' 
                        LOCUST_run__settings_prec_mod['file_tet']="'locust_wall'" 
                        LOCUST_run__settings_prec_mod['file_eqm']="'locust_eqm'" 
                        LOCUST_run__settings_prec_mod['threadsPerBlock']=32
                        LOCUST_run__settings_prec_mod['blocksPerGrid']=32
                        LOCUST_run__settings_prec_mod['root']="'/tmp/{username}/{study}/{params}'".format(username=settings.username,study=RMP_study__name,params=parameters__parameter_string)
                        LOCUST_run__settings_prec_mod['i3dr']=-1 #XXX WHILST I3DR FLAG IS BROKE
                        LOCUST_run__settings_prec_mod['niter']=1

                        LOCUST_run__settings_prec_mod['incl']=np.zeros(32) #opt to run more markers on fewer GPUs in parallel
                        LOCUST_run__settings_prec_mod['incl'][run_number]=1 
                        LOCUST_run__settings_prec_mod['incl']=f"{LOCUST_run__settings_prec_mod['incl']}"

                        MARS_read__flags={}
                        MARS_read__flags['TOKAMAK']=1
                        MARS_read__flags['N0']=parameters__toroidal_mode_number[0]
                        MARS_read__flags['PLS']=True
                        MARS_read__flags['UPHASE']=f'{parameters__phase_upper}D0' #XXX does this account for counter-rotating harmonics?
                        MARS_read__flags['MPHASE']=f'{parameters__phase_middle}D0'
                        MARS_read__flags['LPHASE']=f'{parameters__phase_lower}D0'
                        MARS_read__settings={}
                        MARS_read__settings['TAIL']="{}".format(MARS_read__tails)
                        MARS_read__settings['IKATN']=f'[{parameters__current_upper/1000.}_gpu,{parameters__current_middle/1000.}_gpu,{parameters__current_lower/1000.}_gpu]'
                        MARS_read__settings['dXR']=f'{0.010}_gpu'
                        MARS_read__settings['dXZ']=f'{0.010}_gpu'
                        
                        NEMO_run__xml_settings={}
                        NEMO_run__xml_settings['nmarker']=LOCUST_run__settings_prec_mod['threadsPerBlock']*LOCUST_run__settings_prec_mod['blocksPerGrid']*16
                        NEMO_run__xml_settings['fokker_flag']=0

                        BBNBI_run__xml_settings={}
                        BBNBI_run__number_particles=LOCUST_run__settings_prec_mod['threadsPerBlock']*LOCUST_run__settings_prec_mod['blocksPerGrid']*16
                        BBNBI_run__dir_BBNBI=support.dir_bbnbi
                        
                        #if all coilsets do not rotate together we must split them up individually!
                        if all(rotation==parameters__rotation_upper for rotation in [parameters__rotation_upper,parameters__rotation_middle,parameters__rotation_lower]): 
                            #if coils rotate together but we still want one row offset with others then define relative phase for mars_read
                            nnum_string=str(['{}'.format(mode) for mode in parameters__toroidal_mode_number]).replace('\'','')
                            phase_string='['
                            omega_string='['
                            for mode in parameters__toroidal_mode_number:
                                for phase,omega in zip([parameters__phase_upper],
                                                        [parameters__rotation_upper]):
                                    phase_string+='{}_gpu,'.format(0.0) #setting phase to 0 since this now controlled in mars_read preprocessor code 
                                    omega_string+='{}_gpu,'.format(omega)
                            phase_string=phase_string[:-1]+']'
                            omega_string=omega_string[:-1]+']'

                        else:
                            LOCUST_run__flags['NCOILS']=3 #if we have separate coils then multiply up quantities for each coilset
                            LOCUST_run__settings_prec_mod['nmde']*=LOCUST_run__flags['NCOILS'] #number of total toroidal harmonics = number of modes * number of coilsets
                            nnum_string=str(['{}'.format(mode) for mode in parameters__toroidal_mode_number for coil_row in range(LOCUST_run__flags['NCOILS'])]).replace('\'','')
                            phase_string='['
                            omega_string='['
                            for mode in parameters__toroidal_mode_number:
                                for phase,omega in zip([parameters__phase_upper,parameters__phase_middle,parameters__phase_lower],
                                                        [parameters__rotation_upper,parameters__rotation_middle,parameters__rotation_lower]):
                                    phase_string+='{}_gpu,'.format(0.0) #setting phase to 0 since this now controlled in mars_read preprocessor code
                                    omega_string+='{}_gpu,'.format(omega)
                            phase_string=phase_string[:-1]+']'
                            omega_string=omega_string[:-1]+']'

                        LOCUST_run__settings_prec_mod['phase']=phase_string
                        LOCUST_run__settings_prec_mod['omega']=omega_string
                        LOCUST_run__settings_prec_mod['nnum']=nnum_string

                        args_batch['parameters__sheet_name_kinetic_prof'].append(copy.deepcopy(parameters__sheet_name_kinetic_prof))
                        args_batch['parameters__sheet_name_rotation'].append(copy.deepcopy('"{}"'.format(target_IDS_dispatch[parameters__database][parameters__kinetic_prof_Pr_string][parameters__kinetic_prof_tF_tE_string]['sheet_name_rotation'])))
                        args_batch['parameters__var_name_rotation'].append(copy.deepcopy('"{}"'.format(target_IDS_dispatch[parameters__database][parameters__kinetic_prof_Pr_string][parameters__kinetic_prof_tF_tE_string]['var_name_rotation'])))
                        args_batch['parameters__toroidal_mode_numbers'].append(copy.deepcopy("'{}'".format(parameters__toroidal_mode_number)))
                        args_batch['LOCUST_run__dir_LOCUST'].append(copy.deepcopy("'{}'".format(str(support.dir_locust / parameters__database / RMP_study__name / parameters__parameter_string))))
                        args_batch['LOCUST_run__dir_LOCUST_source'].append(copy.deepcopy("'{}'".format(str(support.dir_locust / 'source'))))
                        args_batch['LOCUST_run__dir_input'].append(copy.deepcopy("'{}'".format(str(support.dir_input_files / parameters__database / RMP_study__name / parameters__parameter_string))))
                        args_batch['LOCUST_run__dir_output'].append(copy.deepcopy("'{}'".format(str(support.dir_output_files / parameters__database / RMP_study__name / parameters__parameter_string))))
                        args_batch['LOCUST_run__dir_cache'].append(copy.deepcopy("'{}'".format(str(support.dir_cache_files / parameters__database / RMP_study__name)))) #one level less to pool cache files into same directory across simulations
                        args_batch['LOCUST_run__environment_name'].append(copy.deepcopy(LOCUST_run__environment_name))
                        args_batch['LOCUST_run__repo_URL'].append(copy.deepcopy(LOCUST_run__repo_URL))
                        args_batch['LOCUST_run__commit_hash'].append(copy.deepcopy(LOCUST_run__commit_hash))
                        args_batch['LOCUST_run__settings_prec_mod'].append(copy.deepcopy(LOCUST_run__settings_prec_mod))
                        args_batch['LOCUST_run__flags'].append(copy.deepcopy(LOCUST_run__flags))
                        args_batch['NEMO_run__dir_NEMO'].append(copy.deepcopy(support.dir_nemo))
                        args_batch['NEMO_run__xml_settings'].append(copy.deepcopy(NEMO_run__xml_settings))
                        args_batch['BBNBI_run__dir_BBNBI'].append(copy.deepcopy(BBNBI_run__dir_BBNBI))
                        args_batch['BBNBI_run__xml_settings'].append(copy.deepcopy(BBNBI_run__xml_settings))
                        args_batch['BBNBI_run__number_particles'].append(copy.deepcopy(BBNBI_run__number_particles))
                        args_batch['MARS_read__tail_U'].append(copy.deepcopy(MARS_read__tail_U))
                        args_batch['MARS_read__tail_M'].append(copy.deepcopy(MARS_read__tail_M))
                        args_batch['MARS_read__tail_L'].append(copy.deepcopy(MARS_read__tail_L))
                        args_batch['MARS_read__settings'].append(copy.deepcopy(MARS_read__settings))
                        args_batch['MARS_read__flags'].append(copy.deepcopy(MARS_read__flags))
                        args_batch['MARS_read__dir_MARS_builder'].append(copy.deepcopy("'{}'".format(str(support.dir_cache_files / parameters__database / RMP_study__name / (parameters__parameter_string + '_MARS_builder') ))))
                        args_batch['RMP_study__name'].append(copy.deepcopy(RMP_study__name))
                        args_batch['RMP_study__filepath_kinetic_profiles'].append(copy.deepcopy(list((RMP_study__dir_input_database / parameters__database / folder_name_DataEq).glob('*.xlsx'))[0])) #determine path to current kinetic profiles
                        args_batch['RMP_study__filepath_equilibrium'].append(copy.deepcopy(list((RMP_study__dir_input_database / parameters__database / folder_name_DataEq).glob('*eqdsk*'))[0])) #determine path to current equilibrium
                        args_batch['RMP_study__filepath_additional_data'].append(copy.deepcopy(RMP_study__filepaths_additional_data))
                        args_batch['RMP_study__workflow_commands'].append(RMP_study__workflow_commands)
                        
                        #find paths to 3D fields corresponding to desired parameters depending on requested field type

                        RMP_study__field_type='plasma_response' #response or vacuum

                        for coil_row in ['cU','cM','cL']:

                            if RMP_study__field_type is 'vacuum':
                                RMP_study__filepaths_3D_field_head=RMP_study__dir_input_database / 'ITER_15MAQ10_case5'/ 'DataVac' / 'PureVac'
                                field_filepath_string='"{}"'.format([str(RMP_study__filepaths_3D_field_head/f'BPLASMA_MARSF_n{np.abs(mode)}_{coil_row}_V') for mode in parameters__toroidal_mode_number])

                            elif RMP_study__field_type is 'vacuum_resistive_wall':
                                RMP_study__filepaths_3D_field_head=RMP_study__dir_input_database / 'ITER_15MAQ10_case5'/ 'DataVac' / 'WithRW'
                                field_filepath_string='"{}"'.format([str(RMP_study__filepaths_3D_field_head/f'BPLASMA_MARSF_n{np.abs(mode)}_{coil_row}_Pr{parameters__kinetic_prof_Pr_string}_tfte{parameters__kinetic_prof_tF_tE_string}.IN') for mode in parameters__toroidal_mode_number])
                            
                            elif RMP_study__field_type is 'plasma_response':
                                RMP_study__filepaths_3D_field_head=RMP_study__dir_input_database / parameters__database / folder_name_DataMarsf
                                field_filepath_string='"{}"'.format([str(RMP_study__filepaths_3D_field_head/f'BPLASMA_MARSF_n{np.abs(mode)}_{coil_row}_Pr{parameters__kinetic_prof_Pr_string}_tfte{parameters__kinetic_prof_tF_tE_string}.IN') for mode in parameters__toroidal_mode_number])

                            args_batch[f'RMP_study__filepaths_3D_fields_{coil_row[-1]}'].append(field_filepath_string)

                        args_batch['IDS__shot'].append(copy.deepcopy(IDS__shot))
                        args_batch['IDS__run'].append(copy.deepcopy(IDS__run))
                        args_batch['IDS__username'].append(copy.deepcopy(IDS__username))
                        args_batch['IDS__imasdb'].append(copy.deepcopy(IDS__imasdb))
                        args_batch['IDS__target_IDS_shot'].append(copy.deepcopy(target_IDS_dispatch[parameters__database][parameters__kinetic_prof_Pr_string][parameters__kinetic_prof_tF_tE_string]['shot']))
                        args_batch['IDS__target_IDS_run'].append(copy.deepcopy(target_IDS_dispatch[parameters__database][parameters__kinetic_prof_Pr_string][parameters__kinetic_prof_tF_tE_string]['run']))

                        args_batch['IDS__NBI_shot'].append(copy.deepcopy(config_beam_dispatch[config_beam_1][config_beam_2]['shot']))
                        args_batch['IDS__NBI_run'].append(copy.deepcopy(config_beam_dispatch[config_beam_1][config_beam_2]['run']))
                        args_batch['IDS__NBI_imasdb'].append(copy.deepcopy(config_beam_dispatch[config_beam_1][config_beam_2]['imasdb']))
                        args_batch['IDS__NBI_username'].append(copy.deepcopy(config_beam_dispatch[config_beam_1][config_beam_2]['user']))

##################################################################
#define and launch the batch scripts

if __name__=='__main__':

    def Batch_wrapper(args_batch,**kwargs):
        batch_instance=Batch(**args_batch)
        batch_instance.launch(**kwargs) 

    split_node=False
    if split_node:

        import multiprocessing,time
        ngpus=8
        ngpus_per_group=1
        ngpu_groups=int(ngpus/ngpus_per_group)
        n_sims_per_gpu_group=int(len(list(args_batch.values())[0])/ngpu_groups)

        for gpu_group in range(ngpu_groups):
            for sim in np.arange(n_sims_per_gpu_group*gpu_group,n_sims_per_gpu_group*(gpu_group+1)):
                args_batch['LOCUST_run__settings_prec_mod'][sim]['incl']=np.zeros(32,dtype=int) #opt to run more markers on fewer GPUs in parallel
                args_batch['LOCUST_run__settings_prec_mod'][sim]['incl'][gpu_group*ngpus_per_group:(gpu_group+1)*ngpus_per_group]=1 
                args_batch['LOCUST_run__settings_prec_mod'][sim]['incl']=np.array2string(args_batch['LOCUST_run__settings_prec_mod'][sim]['incl'],separator=',')

        processes = [multiprocessing.Process(target=Batch_wrapper, args=({key:value[slice(gpu*n_sims_per_gpu_group,(gpu+1)*n_sims_per_gpu_group)] for key,value in args_batch.items()},), kwargs={'workflow_filepath':path_template_run,'environment_name_batch':LOCUST_run__environment_name,'environment_name_workflow':LOCUST_run__environment_name,'interactive':True}) for gpu in range(ngpus)]

        for p in processes:
            p.start()
            time.sleep(300)
        for p in processes:
            p.join()

    else:

        Batch_wrapper(args_batch,
            workflow_filepath=path_template_run,
            environment_name_batch=LOCUST_run__environment_name,
            environment_name_workflow=LOCUST_run__environment_name,
            interactive=False)  

#################################
 
##################################################################
 
###################################################################################################