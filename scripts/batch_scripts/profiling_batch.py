#profiling_batch.py
 
"""
Samuel Ward
22/01/2019
----
script for controlling and launching LOCUST profiling workflow
--- 
notes:         
    single species
---
"""

###################################################################################################
#Preamble
 
import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway

try:
    import numpy as np
    import pathlib
    import copy
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

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

################################################################## 
#Main 

#define parameters which are fixed throughout a parameter scan - if we want to vary then add as a layer in the for loops
study_name='scan_threads_blocks'

LOCUST_run__environment_name='TITAN'
LOCUST_run__repo_URL=f"'{settings.repo_URL_LOCUST}'"
LOCUST_run__commit_hash=f"'{settings.commit_hash_default_LOCUST}'"

GPU_card_dispatch={} #all possible options listed here
GPU_card_dispatch['GPU6']='GTX_TITAN'
GPU_card_dispatch['GPU9']='K80'
GPU_card_dispatch['TITAN']='P100'
GPU_card_dispatch['VIKING']='V100'

#################################
#initialise args needed for Batch class

args_batch={}
args_batch['LOCUST_run__dir_LOCUST']=[]
args_batch['LOCUST_run__dir_LOCUST_source']=[]
args_batch['LOCUST_run__dir_input']=[]
args_batch['LOCUST_run__dir_output']=[]
args_batch['LOCUST_run__dir_cache']=[]
args_batch['LOCUST_run__environment_name']=[]
args_batch['LOCUST_run__repo_URL']=[]
args_batch['LOCUST_run__commit_hash']=[]
args_batch['LOCUST_run__settings_prec_mod']=[]
args_batch['LOCUST_run__flags']=[]

##################################################################
#define variable parameters

parameter__threads=np.array([2**n for n in [4,5,6]])
parameter__blocks=np.array([2**n for n in [4,5,6,7,8,9,10,11]])
parameter__timax=np.array([0.001,0.01,0.1,1.0])

##################################################################
#create every valid combination of parameter, returned in flat lists
#use zip and nest levels to define specific combinations which cannot be varied

for thread in parameter__threads:
    for block in parameter__blocks:
        for timax in parameter__timax:
        
            #create string holding varying parameters                                       
            parameter_string=''
            parameter_string+='_'.join(['{}_{}'.format(parameter,str(value)) for parameter,value in zip([
                            'threads',
                            'blocks',
                            'timax'
                            ],[
                            thread,
                            block,
                            timax])])

            LOCUST_run__flags={}
            LOCUST_run__flags['TOKAMAK']=1
            LOCUST_run__flags['LEIID']=8
            LOCUST_run__flags['WLIST']=True
            LOCUST_run__flags['BRELAX']=True
            LOCUST_run__flags['UNBOR']=100
            LOCUST_run__flags['OPENMESH']=True
            LOCUST_run__flags['OPENTRACK']=True
            LOCUST_run__flags['PFCMOD']=True
            LOCUST_run__flags['NOPFC']=True
            LOCUST_run__flags['TOKHEAD']=True
            LOCUST_run__flags['JXB2']=True
            LOCUST_run__flags['PROV']=True
            LOCUST_run__flags['PITCHCUR']=True
            LOCUST_run__flags['EBASE']=True
            LOCUST_run__flags['UHST']=True
            LOCUST_run__flags['LNLBT']=True
            LOCUST_run__flags['GEQDSKFIX1']=True
            LOCUST_run__flags['BP']=True
            LOCUST_run__flags['TIMAX']='{}D0'.format(timax)
            LOCUST_run__flags['B3D']=True
            LOCUST_run__flags['B3D_EX']=True
            LOCUST_run__flags['SMALLEQ']=True #XXX test whether we need this when using mesh
            LOCUST_run__flags['NCOILS']=3
            LOCUST_run__flags['STATIC']=True
            LOCUST_run__settings_prec_mod={}
            LOCUST_run__settings_prec_mod['Ab']='AD' 
            LOCUST_run__settings_prec_mod['Zb']='+1.0_gpu' 
            LOCUST_run__settings_prec_mod['file_tet']="'locust_wall'" 
            LOCUST_run__settings_prec_mod['file_eqm']="'locust_eqm'" 
            LOCUST_run__settings_prec_mod['threadsPerBlock']=thread
            LOCUST_run__settings_prec_mod['blocksPerGrid']=block
            LOCUST_run__settings_prec_mod['root']="'/tmp/{username}/{study}/{params}'".format(username=settings.username,study=study_name,params=parameter_string)
            LOCUST_run__settings_prec_mod['nnum']='[-3,-3,-3,-6,-6,-6]'
            LOCUST_run__settings_prec_mod['nmde']=6
            LOCUST_run__settings_prec_mod['phase']='[0.0_gpu,0.0_gpu,0.0_gpu,0.0_gpu,0.0_gpu,0.0_gpu]'
            LOCUST_run__settings_prec_mod['omega']='[0.0_gpu,0.0_gpu,0.0_gpu,0.0_gpu,0.0_gpu,0.0_gpu]'
            LOCUST_run__settings_prec_mod['fi']='[1.0_gpu]'
            LOCUST_run__settings_prec_mod['Ai']='[AD]'
            LOCUST_run__settings_prec_mod['Zi']='[+1.0_gpu]'
            LOCUST_run__settings_prec_mod['nion']=1

            args_batch['LOCUST_run__dir_LOCUST'].append(copy.deepcopy("'{}'".format(str(support.dir_locust / study_name / parameter_string))))
            args_batch['LOCUST_run__dir_LOCUST_source'].append(copy.deepcopy("'{}'".format(str(support.dir_locust / 'source'))))
            args_batch['LOCUST_run__dir_input'].append(copy.deepcopy("'{}'".format(str(support.dir_input_files / study_name))))
            args_batch['LOCUST_run__dir_output'].append(copy.deepcopy("'{}'".format(str(support.dir_output_files / study_name / GPU_card_dispatch[LOCUST_run__environment_name] / parameter_string))))
            args_batch['LOCUST_run__dir_cache'].append(copy.deepcopy("'{}'".format(str(support.dir_cache_files / study_name)))) #one level less to pool cache files into same directory across simulations
            args_batch['LOCUST_run__environment_name'].append(copy.deepcopy(LOCUST_run__environment_name))
            args_batch['LOCUST_run__repo_URL'].append(copy.deepcopy(LOCUST_run__repo_URL))
            args_batch['LOCUST_run__commit_hash'].append(copy.deepcopy(LOCUST_run__commit_hash))
            args_batch['LOCUST_run__settings_prec_mod'].append(copy.deepcopy(LOCUST_run__settings_prec_mod))
            args_batch['LOCUST_run__flags'].append(copy.deepcopy(LOCUST_run__flags))

##################################################################
#define and launch the batch scripts

if __name__=='__main__':

    profiling_batch=Batch(**args_batch)
    profiling_batch.launch(workflow_filepath='profiling_workflow.py',environment_name_batch=LOCUST_run__environment_name,environment_name_workflow=LOCUST_run__environment_name,interactive=False)   

#################################
 
##################################################################
 
###################################################################################################