#profiling_workflow.py

'''
Samuel Ward
22/01/2019
----
workflow designed for profiling LOCUST
---
notes: 

    assumes all input files already located in locust dir_input
---
'''


##################################################################
#Preamble

try:
    import sys
    import subprocess
    import pathlib
    import shlex
    import copy
    import numpy as np
    import ast
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
    import run_scripts.workflow
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/workflow.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import run_scripts.environment
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/environment.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    from run_scripts.LOCUST_run import LOCUST_run
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/LOCUST_run.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import run_scripts.NEMO_run
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/NEMO_run.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    from run_scripts.MARS_builder_run import MARS_builder_run
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/MARS_builder_run.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import run_scripts.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/utils.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)

try:
    from classes.input_classes.temperature import Temperature
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/temperature.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.number_density import Number_Density
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/number_density.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.rotation import Rotation
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/rotation.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.equilibrium import Equilibrium
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/equilibrium.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.perturbation import Perturbation
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/perturbation.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.beam_deposition import Beam_Deposition
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/perturbation.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.wall import Wall
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/wall.py could not be imported!\nreturning\n") 
    sys.exit(1)

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/src/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import settings
except:
    raise ImportError("ERROR: LOCUST_IO/src/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)

##################################################################
#Main

class profiling_workflow(run_scripts.workflow.Workflow):
    """
    defines a workflow which runs a performance parameter scan

    notes:
    usage:
    """ 

    def __init__(self,
                *args,
                **kwargs):
        """
        notes:
        """

        super().__init__()

        self.args={**kwargs}

        #################################
        #derive some information
        
        #################################
        #add workflow stages

        self.add_command(command_name='mkdir',command_function=self.setup_study_dirs,position=1) #add new workflow stages
        self.add_command(command_name='run_LOCUST',command_function=self.run_LOCUST,position=2)
        self.add_command(command_name='delete_outputs',command_function=self.delete_outputs,position=3)

    def setup_study_dirs(self,*args,**kwargs):
        """
        notes:
        """

        for direc in [
                        self.args['LOCUST_run__dir_input'],
                        self.args['LOCUST_run__dir_output'],
                        self.args['LOCUST_run__dir_cache']
                        ]:
            if not pathlib.Path(direc).is_dir(): pathlib.Path(direc).mkdir(parents=True)

    def run_LOCUST(self,*args,**kwargs):
        """
        notes:
        """

        LOCUST_workflow=LOCUST_run(environment_name=self.args['LOCUST_run__environment_name'],
            repo_URL=self.args['LOCUST_run__repo_URL'],
            commit_hash=self.args['LOCUST_run__commit_hash'],
            dir_LOCUST=self.args['LOCUST_run__dir_LOCUST'],
            dir_input=self.args['LOCUST_run__dir_input'],
            dir_output=self.args['LOCUST_run__dir_output'],
            dir_cache=self.args['LOCUST_run__dir_cache'],
            settings_prec_mod=self.args['LOCUST_run__settings_prec_mod'],
            flags=self.args['LOCUST_run__flags'])

        LOCUST_workflow.run()

    def delete_outputs(self,*args,**kwargs):
        """
        notes:
        """

        for file in pathlib.Path(self.args['LOCUST_run__dir_output']).glob('*'):
            if 'rundata' not in str(file):
                subprocess.run(shlex.split('rm {file}'.format(file=str(file))),shell=False) #delete all main outputs except rundata messages


if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser(description='run profiling_workflow from command line in Python!')
    
    parser.add_argument('--LOCUST_run__dir_LOCUST',type=str,action='store',dest='LOCUST_run__dir_LOCUST',help="",default=support.dir_locust)
    parser.add_argument('--LOCUST_run__dir_input',type=str,action='store',dest='LOCUST_run__dir_input',help="",default=support.dir_input_files)
    parser.add_argument('--LOCUST_run__dir_output',type=str,action='store',dest='LOCUST_run__dir_output',help="",default=support.dir_output_files)
    parser.add_argument('--LOCUST_run__dir_cache',type=str,action='store',dest='LOCUST_run__dir_cache',help="",default=support.dir_cache_files)
    parser.add_argument('--LOCUST_run__environment_name',type=str,action='store',dest='LOCUST_run__environment_name',help="",default='TITAN')
    parser.add_argument('--LOCUST_run__repo_URL',type=str,action='store',dest='LOCUST_run__repo_URL',help="",default=settings.repo_URL_LOCUST)
    parser.add_argument('--LOCUST_run__commit_hash',type=str,action='store',dest='LOCUST_run__commit_hash',help="",default=None)
    parser.add_argument('--LOCUST_run__settings_prec_mod',nargs='+',type=str,action='store',dest='LOCUST_run__settings_prec_mod',help="",default={})
    parser.add_argument('--LOCUST_run__flags',nargs='+',type=str,action='store',dest='LOCUST_run__flags',help="",default={})

    args=parser.parse_args()

    #provide some extra parsing steps to dict-like and array-like input arguments
    args.LOCUST_run__settings_prec_mod=run_scripts.utils.command_line_arg_parse_dict(args.LOCUST_run__settings_prec_mod)
    args.LOCUST_run__flags=run_scripts.utils.command_line_arg_parse_dict(args.LOCUST_run__flags)

    workflow=profiling_workflow(**{key:arg for key,arg in args._get_kwargs()})
    workflow.run()

#################################

##################################################################

###################################################################################################