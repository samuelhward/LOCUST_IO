#batch.py

'''
Samuel Ward
03/10/2019
----
tools to batch run workflows 
---
notes:
    essentially translates given settings into command line args for hard-coded batch scripts to parse
---
'''


##################################################################
#Preamble

try:
    import sys
    import subprocess
    import pathlib
    import shlex
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import run_scripts.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/utils.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import run_scripts.environment
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/environment.py could not be imported!\nreturning\n")
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

class Batch:
    """
    class to batch run many LOCUST_runs

    notes:

    usage:
    """ 

    batch_system={}

    batch_system['TITAN']={}
    batch_system['TITAN']['shebang']='#!/bin/bash'
    batch_system['TITAN']['command_execute_workflow']='python3'
    batch_system['TITAN']['command_launch_batch']='sbatch'
    batch_system['TITAN']['flag_command']='#SBATCH'
    batch_system['TITAN']['flags']={}
    batch_system['TITAN']['flags']['job-name']='LOCUST'
    batch_system['TITAN']['flags']['gres']='gpu:8'
    batch_system['TITAN']['flags']['partition']='gpu_p100_titan'
    batch_system['TITAN']['flags']['cpus-per-task']=1
    batch_system['TITAN']['flags']['exclusive']=True
    batch_system['TITAN']['flags']['mail-user']='samuel.ward@york.ac.uk'
    batch_system['TITAN']['flags']['mail-type']='END,FAIL,TIME_LIMIT'
    batch_system['TITAN']['flags']['o']='LOCUST_SLURM.out'
    batch_system['TITAN']['flags']['e']='LOCUST_SLURM.err'

    batch_system['CUMULUS']={}
    batch_system['CUMULUS']['shebang']='#!/bin/bash'
    batch_system['CUMULUS']['command_execute_workflow']='python3'
    batch_system['CUMULUS']['command_launch_batch']='qsub'
    batch_system['CUMULUS']['flag_command']='#PBS'
    batch_system['CUMULUS']['flags']={}
    batch_system['CUMULUS']['flags']['V']=True
    batch_system['CUMULUS']['flags']['l']='select=1:ncpus=8:place=excl:walltime=2400:00:00'
    batch_system['CUMULUS']['flags']['q']='gpu'
    batch_system['CUMULUS']['flags']['N']='LOCUST'
    batch_system['CUMULUS']['flags']['o']='LOCUST_PBS.out'
    batch_system['CUMULUS']['flags']['e']='LOCUST_PBS.err'

    batch_system['GPU5']={}
    batch_system['GPU5']['command_execute_workflow']='python'
    batch_system['GPU6']={}
    batch_system['GPU6']['command_execute_workflow']='python'
    batch_system['GPU7']={}
    batch_system['GPU7']['command_execute_workflow']='python'
    batch_system['GPU8']={}
    batch_system['GPU8']['command_execute_workflow']='python'
    batch_system['GPU9']={}
    batch_system['GPU9']['command_execute_workflow']='python'
    batch_system['GPU10']={}
    batch_system['GPU10']['command_execute_workflow']='python'

    batch_system['VIKING']={}
    batch_system['VIKING']['shebang']='#!/bin/bash'
    batch_system['VIKING']['command_execute_workflow']='python'
    batch_system['VIKING']['command_launch_batch']='sbatch'
    batch_system['VIKING']['flag_command']='#SBATCH'
    batch_system['VIKING']['flags']={}
    batch_system['VIKING']['flags']['job-name']='LOCUST'
    batch_system['VIKING']['flags']['gres']='gpu:2'
    batch_system['VIKING']['flags']['partition']='gpu'
    batch_system['VIKING']['flags']['cpus-per-task']=1
    batch_system['VIKING']['flags']['mail-user']='samuel.ward@york.ac.uk'
    batch_system['VIKING']['flags']['mail-type']='END,FAIL,TIME_LIMIT'
    batch_system['VIKING']['flags']['output']='LOCUST_SLURM.out'
    batch_system['VIKING']['flags']['error']='LOCUST_SLURM.err'
    batch_system['VIKING']['flags']['ntasks']=1
    batch_system['VIKING']['flags']['mem']='12gb'
    batch_system['VIKING']['flags']['time']='02:00:00'


    def __init__(self,**batch_settings):
        """
        notes:
            every list MUST be the same length
            args containing multiple settings should be set using list of dicts e.g. flags=[{'TOKAMAK':1,'STDOUT':True}] or settings_prec_mod=[{'root':"'some/file/path'"}] - this equates to a single workflow run with --flags TOKAMAK=1 STDOUT  
        args:
            batch_settings - kwargs storing lists which describe command line args for chosen workflow e.g. filepaths_in=[/some/path/1,some/path/2] if you want to launch a workflow twice, once with --filepaths_in=/some/path/1 and once with --filepaths_in=/some/path/2 
        """

        self.batch_settings=batch_settings
        for setting in self.batch_settings:
            self.number_runs=len(batch_settings[setting]) #infer how many runs the user wants to launch by the lenght of the lists stored in batch_settings

    def launch(self,workflow_filepath,environment_name_batch=None,environment_name_workflow='TITAN',interactive=False):
        """
        notes:
            the goal is to move through all the lists stored in batch_settings and generate a corresponding batch run script which calls some workflow with args that match those defined in batch_settings
        args:
            workflow_filepath - filepath to file storing workflow to run e.g. support.dir_run_scripts / 'LOCUST_run.py'
            environment_name_batch - string identifier to choose batch script environment 
            environment_name_workflow - string identifier to choose workflow environment 
            interactive - launch workflow interactively without batch system
        """

        if not environment_name_workflow: 
            print("ERROR: batch.launch() must know environment_name_workflow!\nreturning\n!")
            return

        if environment_name_batch:
            try:
                environment_batch=run_scripts.environment.Environment(environment_name_batch) #attach an environment
            except:
                print(f"WARNING: Batch.launch() could not find valid environment_name_batch ({environment_name_batch}) - running with default!")
                environment_batch=None

        ################################# loop over each run in this set of batch runs

        for run in range(self.number_runs): #equivalent to the length of lists held in batch_settings

            ################################# create arg string to pass to workflow script

            workflow_args=run_scripts.utils.command_line_arg_parse_generate_string(command_number_=run,**self.batch_settings) #string to hold args passed to run script called from batch script

            if interactive is False:

                ################################# create string holding batch script file settings

                system_batch_flags=''
                for batch_flag,value in Batch.batch_system[environment_name_batch]['flags'].items(): #add '-' or '--' to batch system flags e.g. num_processors --> --num_processorss                    
                    batch_flag=''.join(['-',batch_flag]) if len(batch_flag)==1 else ''.join(['--',batch_flag]) 
                    batch_flag='{} {}'.format(batch_flag,str(value)) if value is not True else '{}'.format(batch_flag) #create strings for batch system flags                
                    system_batch_flags=''.join([system_batch_flags,Batch.batch_system[environment_name_batch]['flag_command'],' ',batch_flag,'\n'])

                ################################# combine to generate final batch script


                #create a dummy batch script file
                batch_file="""\
{shebang}
{system_batch_flags}

{environment_batch}

{command_execute_workflow} {workflow_filepath} {workflow_args}
                """.format(
                shebang=Batch.batch_system[environment_name_batch]['shebang'],
                system_batch_flags=system_batch_flags,
                environment_batch=''.join(f'{command}\n' for command in environment_batch.create_command_string(listed=True)) if environment_batch else '',
                command_execute_workflow=Batch.batch_system[environment_name_batch]['command_execute_workflow'],
                workflow_filepath=str(workflow_filepath),
                workflow_args=workflow_args
                ).lstrip(' ')

                batch_file_path=support.dir_run_scripts / '.Batch_submission_script_{}'.format(run)
                with open(batch_file_path,'w') as file:
                    file.write(batch_file.lstrip(' '))

                #launch batch scripts here   
                try:
                    subprocess.run([str(Batch.batch_system[environment_name_batch]['command_launch_batch']),str(batch_file_path)])
                except subprocess.CalledProcessError as err:
                    raise(err)

            else:

                ################################# if not running with batch just execute each workflow individually

                try:
                    subprocess.run(
                        environment_batch.create_command_string()+'; '+
                        str(Batch.batch_system[environment_name_batch]['command_execute_workflow'])+' '+
                        str(workflow_filepath)+
                        workflow_args, 
                        shell=True)
                except subprocess.CalledProcessError as err:
                    raise(err)

#################################

##################################################################

###################################################################################################