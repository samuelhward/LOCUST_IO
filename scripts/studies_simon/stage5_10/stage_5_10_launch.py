#stage_5_10_launch.py
 
"""
Samuel Ward
29/06/21
----
stage 5.10

    - note that max X-point displacement is optimal
    - fix Pr=0.3 and τφ/τE=2.
    - fix current at 45kAt
    - fix central coil phase from stage 4.10
    - vary upper/lower coil phase to trace out contour on 2/3 max X-point displacement in 8 steps
    - vary n=3,n=4 + sideband
    - vary 7.5MA/2.65T (case 2) hydrogen/Hbeam + deuterium-deuterium/Dbeam and 7.5MA/5.3T (case 3) deuterium-deuterium/Dbeam and 7.5MA/4.5T (case 7) deuterium-deuterium/Dbeam scenarios
    - 64 runs
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

RMP_study__name='stage_5_10'

#################################
#define options and dispatch tables for helping choosing settings

parameters__toroidal_mode_numbers__options={} #XXX needs studying - toroidal mode number combinations 
parameters__toroidal_mode_numbers__options['n=3']=[-3,-6]
parameters__toroidal_mode_numbers__options['n=4']=[-4,-5]

##################################################################
#choose the scenarios we will want to examine

parameters__databases=['ITER_7d5MAHalfB_case2','ITER_7d5MAHalfB_case2','ITER_7d5MAFullB_case3','ITER_7d5MA4d5T_case7'] #these all zipped at same level in main loop because source data is not consistent enough
parameters__sheet_names_kinetic_prof=["'nT=0.76ne'","'nT=0.76ne'","'iterDD.iterFSBMI'","'Pr=0.3,tF=2tE'"]
configs_beam_species=['hydrogen','deuterium','deuterium','deuterium']
plasmas_species=[['hydrogen'],['deuterium'],['deuterium'],['deuterium']] #"\"['']\""

##################################################################
#define the parameter space for a given scenario

#kinetic profile parameters which vary independently
parameters__kinetic_profs_Pr=[0.3]
parameters__kinetic_profs_tF_tE=[2.]


#3D field parameters which vary independently - if you want to vary these together then put them into the same loop nesting below
#2D arrays, each element has length = number of modes

#n=3 and 4 harmonic settings - both first side band and corresponding rigid phases to scan
parameters__toroidal_mode_numbers=[[-3,-6],[-4,-5]]

# since optimal phases now depend on both case and toroidal mode number need to introduce another iteration level

# coordinates along contour of 2/3 max X-point displacement using plot_X_point_displacement.py
# angles are relative to middle coil position
# divide by mode number here to transform to real space
# assuming that XPD is same for H and D plasmas
contour_relative_phase_upper_case2_n3 = -np.array([74.63723481476995, 94.63311692520949, 157.18125169115623, 245.6339460009049, 304.83859661955285, 284.85559201054417, 222.31397935590854, 133.86030211070397])/3.
contour_relative_phase_lower_case2_n3 = -np.array([90.0, 180.43614789702443, 251.14892495783212, 270.51724339742515, 202.10359194497602, 111.67821291052805, 40.964071279354286, 21.57431145425649])/3.
contour_relative_phase_upper_case2_n4 = -np.array([-235.25428228245076, -216.41234418780357, -159.42596117030766, -74.09198985324687, -14.34647749802351, -33.180026710639886, -90.16255166751492, -175.50484225141875])/4.
contour_relative_phase_lower_case2_n4 = -np.array([20.0, 110.68656027750981, 185.14699648712272, 214.1580849547286, 149.0387661808968, 58.34688734132995, -16.122178881430294, -45.10064442788318])/4.
contour_relative_phase_upper_case3_n3 = -np.array([-145.2716393483703, -120.29546165786665, -58.86720067794775, 26.453077446119536, 85.44889915416834, 60.479195149170856, -0.9504587861743696, -86.27183216038702])/3.
contour_relative_phase_lower_case3_n3 = -np.array([-40.0, 47.38213631307841, 115.75256832617907, 145.10254993460026, 85.3990492991743, -1.9832026288900206, -70.35180144888817, -99.7339733621876])/3.
contour_relative_phase_upper_case3_n4 = -np.array([16.551006721206214, 101.63738829604354, 186.16856300536085, 272.622759500644, 206.2821589689221, 101.73275455393512, 360.0, 283.9907025144561])/4.
contour_relative_phase_lower_case3_n4 = -np.array([0.0, 86.70346757605165, 173.95152266509925, 259.290976102509, 360.0, 258.0480590049009, 153.2833530850786, 75.78994993180221])/4.
contour_relative_phase_upper_case7_n3 = -np.array([-180.8982266, -163.6965917 , -106.09998978, -24.07638456, 46.97770407, 29.81883198, -27.77752547, -109.7981419])/3.
contour_relative_phase_lower_case7_n3 = -np.array([-25., 64.24894848, 136.30523227, 176.08312545, 131.3915865, 42.15514932, -29.90019883, -69.69885907])/3.
contour_relative_phase_upper_case7_n4 = -np.array([-125.07069558, -103.38580963, -44.16881469, 38.68319062, 126.51926018, 104.94660554, 45.71372546, -37.14116784])/4.
contour_relative_phase_lower_case7_n4 = -np.array([-120., -23.67811381, 57.0005977, 112.43866172, 96.51266678, 0.20322827, -80.46864141, -135.90423969])/4.


# leave as negative values since easier to translate to positive than translate whole X-displacement pattern
#contour_relative_phase_upper_case2_n3[contour_relative_phase_upper_case2_n3<0]+=360.
#contour_relative_phase_lower_case2_n3[contour_relative_phase_lower_case2_n3<0]+=360.
#contour_relative_phase_upper_case2_n4[contour_relative_phase_upper_case2_n4<0]+=360.
#contour_relative_phase_lower_case2_n4[contour_relative_phase_lower_case2_n4<0]+=360.
#contour_relative_phase_upper_case3_n3[contour_relative_phase_upper_case3_n3<0]+=360.
#contour_relative_phase_lower_case3_n3[contour_relative_phase_lower_case3_n3<0]+=360.
#contour_relative_phase_upper_case3_n4[contour_relative_phase_upper_case3_n4<0]+=360.
#contour_relative_phase_lower_case3_n4[contour_relative_phase_lower_case3_n4<0]+=360.

# index of optimal absolute phase (lowest losses) from stage_4_10 are:
# case 2, H n=3: 1st simulation 
# case 2, H n=4: 2nd simulation
# case 2, D n=3: 1st simulation 
# case 2, D n=4: 1st simulation
# case 3, D n=3: 2nd simulation 
# case 3, D n=4: 5th simulation 
# case 7, D n=3: 4th simulation 
# case 7, D n=4: 1st simulation 

# to check these are correct, verify that U-M and L-M give the above relative phases and parameters__phases_middles* are the optimal values from stage_4
# assuming that XPD is same for H and D plasmas
parameters__phases_uppers_case2_H=[
    contour_relative_phase_upper_case2_n3+np.linspace(0,120,6)[:-1][0]+30.-3.3,
    contour_relative_phase_upper_case2_n4+np.linspace(0,90,6)[:-1][1]+30.-3.3]
parameters__phases_middles_case2_H=[
    np.full(len(contour_relative_phase_upper_case2_n3),np.linspace(0,120,6)[:-1][0])+26.7,
    np.full(len(contour_relative_phase_upper_case2_n4),np.linspace(0,90,6)[:-1][1])+26.7]
parameters__phases_lowers_case2_H=[
    contour_relative_phase_lower_case2_n3+np.linspace(0,120,6)[:-1][0]+30.-3.3,
    contour_relative_phase_lower_case2_n4+np.linspace(0,90,6)[:-1][1]+30.-3.3]

parameters__phases_uppers_case2_D=[
    contour_relative_phase_upper_case2_n3+np.linspace(0,120,6)[:-1][0]+30.-3.3,
    contour_relative_phase_upper_case2_n4+np.linspace(0,90,6)[:-1][0]+30.-3.3]
parameters__phases_middles_case2_D=[
    np.full(len(contour_relative_phase_upper_case2_n3),np.linspace(0,120,6)[:-1][0])+26.7,
    np.full(len(contour_relative_phase_upper_case2_n4),np.linspace(0,90,6)[:-1][0])+26.7]
parameters__phases_lowers_case2_D=[
    contour_relative_phase_lower_case2_n3+np.linspace(0,120,6)[:-1][0]+30.-3.3,
    contour_relative_phase_lower_case2_n4+np.linspace(0,90,6)[:-1][0]+30.-3.3]

parameters__phases_uppers_case3_D=[
    contour_relative_phase_upper_case3_n3+np.linspace(0,120,6)[:-1][1]+30.-3.3,
    contour_relative_phase_upper_case3_n4+np.linspace(0,90,6)[:-1][4]+30.-3.3]
parameters__phases_middles_case3_D=[
    np.full(len(contour_relative_phase_upper_case3_n3),np.linspace(0,120,6)[:-1][1])+26.7,
    np.full(len(contour_relative_phase_upper_case3_n4),np.linspace(0,90,6)[:-1][4])+26.7]
parameters__phases_lowers_case3_D=[
    contour_relative_phase_lower_case3_n3+np.linspace(0,120,6)[:-1][1]+30.-3.3,
    contour_relative_phase_lower_case3_n4+np.linspace(0,90,6)[:-1][4]+30.-3.3]

parameters__phases_uppers_case7_D=[
    contour_relative_phase_upper_case3_n3+np.linspace(0,120,6)[:-1][3]+30.-3.3,
    contour_relative_phase_upper_case3_n4+np.linspace(0,90,6)[:-1][0]+30.-3.3]
parameters__phases_middles_case7_D=[
    np.full(len(contour_relative_phase_upper_case3_n3),np.linspace(0,120,6)[:-1][3])+26.7,
    np.full(len(contour_relative_phase_upper_case3_n4),np.linspace(0,90,6)[:-1][0])+26.7]
parameters__phases_lowers_case7_D=[
    contour_relative_phase_lower_case3_n3+np.linspace(0,120,6)[:-1][3]+30.-3.3,
    contour_relative_phase_lower_case3_n4+np.linspace(0,90,6)[:-1][0]+30.-3.3]


parameters__phases_uppers_cases_all=[parameters__phases_uppers_case2_H,parameters__phases_uppers_case2_D,parameters__phases_uppers_case3_D,parameters__phases_uppers_case7_D]
parameters__phases_middles_cases_all=[parameters__phases_middles_case2_H,parameters__phases_middles_case2_D,parameters__phases_middles_case3_D,parameters__phases_middles_case7_D]
parameters__phases_lowers_cases_all=[parameters__phases_lowers_case2_H,parameters__phases_lowers_case2_D,parameters__phases_lowers_case3_D,parameters__phases_lowers_case7_D]

parameters__rotations_upper=np.array([0.])
parameters__rotations_middle=np.array([0.])
parameters__rotations_lower=np.array([0.])
parameters__currents_upper=np.array([45.])*1000.
parameters__currents_middle=np.array([45.])*1000.
parameters__currents_lower=np.array([45.])*1000.

config_beam_1='off'
config_beam_2='on'

##################################################################
#define the workflow commands in order we want to execute them

RMP_study__workflow_commands="\"['mkdir','save_args','kin_get','3D_get','3D_calc','input_get','IDS_create','run_BBNBI','depo_get','calc_poinc','run_LOCUST','clean_input']\""

##################################################################
#create every valid combination of parameter, returned in flat lists
#use zip and nest levels to define specific combinations which cannot be varied

run_number=0
parameter_strings=[]
#first level are the data which remain constant for a parameter scan
for (parameters__database,parameters__sheet_name_kinetic_prof,config_beam_species,plasma_species,
        parameters__phases_uppers,parameters__phases_middles,parameters__phases_lowers) in zip(
        parameters__databases,parameters__sheet_names_kinetic_prof,configs_beam_species,plasmas_species,
        parameters__phases_uppers_cases_all,parameters__phases_middles_cases_all,parameters__phases_lowers_cases_all): 
    for parameters__kinetic_prof_tF_tE,parameters__kinetic_prof_Pr in zip(parameters__kinetic_profs_tF_tE,parameters__kinetic_profs_Pr):
        #note - different sets of toroidal mode numbers require scanning different absolute phase ranges
        for parameters__toroidal_mode_number,parameters__phases_upper,parameters__phases_middle,parameters__phases_lower in zip(parameters__toroidal_mode_numbers,parameters__phases_uppers,parameters__phases_middles,parameters__phases_lowers):
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

                        parameters__parameter_string+=f'_beams_{str(config_beam_1)}_{str(config_beam_2)}_{config_beam_species[0]}'

                        parameter_strings.append(parameters__parameter_string)

                        #################################
                        #define corresponding workflow args passed to batch (denoted wth __batch)

                        #run-specific settings

                        
                        LOCUST_run__flags=copy.deepcopy(LOCUST_run__flags_default)
                        LOCUST_run__flags['TIMAX']='0.25D0'
                        LOCUST_run__flags['PSIMGRLX']=True #equilibria are relatively poor so relax magnetic axis location
                        #XXX CURRENTLY WAITING FOR FIX LOCUST_run__flags['I3DR']=-1 
                        LOCUST_run__settings_prec_mod={}
                        LOCUST_run__settings_prec_mod['nmde']=len(parameters__toroidal_mode_number) #number of total toroidal harmonics = number of modes
                        LOCUST_run__settings_prec_mod['Ab']=table_species_LOCUST[config_beam_species] 
                        LOCUST_run__settings_prec_mod['Zb']=f'+{table_species_AZ[config_beam_species][-1]}_gpu'  
                        LOCUST_run__settings_prec_mod['file_tet']="'locust_wall'" 
                        LOCUST_run__settings_prec_mod['file_eqm']="'locust_eqm'" 
                        LOCUST_run__settings_prec_mod['threadsPerBlock']=64
                        LOCUST_run__settings_prec_mod['blocksPerGrid']=64
                        LOCUST_run__settings_prec_mod['root']="'/tmp/{username}/{study}/{params}'".format(username=settings.username,study=RMP_study__name,params=parameters__parameter_string)
                        LOCUST_run__settings_prec_mod['i3dr']=-1 #XXX WHILST I3DR FLAG IS BROKE
                        LOCUST_run__settings_prec_mod['niter']=1
                        MARS_read__flags={}
                        MARS_read__flags['TOKAMAK']=1
                        MARS_read__flags['N0']=np.abs(parameters__toroidal_mode_number[0])
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

                        args_batch['parameters__plasma_species'].append('"{}"'.format(plasma_species))
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

                        args_batch['IDS__NBI_shot'].append(copy.deepcopy(config_beam_dispatch[config_beam_1][config_beam_2][config_beam_species]['shot']))
                        args_batch['IDS__NBI_run'].append(copy.deepcopy(config_beam_dispatch[config_beam_1][config_beam_2][config_beam_species]['run']))
                        args_batch['IDS__NBI_imasdb'].append(copy.deepcopy(config_beam_dispatch[config_beam_1][config_beam_2][config_beam_species]['imasdb']))
                        args_batch['IDS__NBI_username'].append(copy.deepcopy(config_beam_dispatch[config_beam_1][config_beam_2][config_beam_species]['user']))

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