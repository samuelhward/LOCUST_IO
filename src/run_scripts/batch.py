#batch.py

'''
Samuel Ward
03/10/2019
----
tools to batch run workflows 
---
notes:
    essentially translates given settings into command line args for hard-coded batch scripts to call
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
    batch_system['TITAN']['flags']['exclusive']=''
    batch_system['TITAN']['flags']['mail-user']='samuel.ward@iter.org'
    batch_system['TITAN']['flags']['mail-type']='END,FAIL,TIME_LIMIT'
    batch_system['TITAN']['flags']['o']='LOCUST_SLURM.out'
    batch_system['TITAN']['flags']['e']='LOCUST_SLURM.err'

    batch_system['CUMULUS']={}
    batch_system['CUMULUS']['shebang']='#!/bin/bash'
    batch_system['CUMULUS']['command_execute_workflow']='python3'
    batch_system['CUMULUS']['command_launch_batch']='qsub'
    batch_system['CUMULUS']['flag_command']='#PBS'
    batch_system['CUMULUS']['flags']={}
    batch_system['CUMULUS']['flags']['V']=None
    batch_system['CUMULUS']['flags']['l']='select=1:ncpus=8:place=excl:walltime=2400:00:00'
    batch_system['CUMULUS']['flags']['q']='gpu'
    batch_system['CUMULUS']['flags']['N']='LOCUST'
    batch_system['CUMULUS']['flags']['o']='LOCUST_PBS.out'
    batch_system['CUMULUS']['flags']['e']='LOCUST_PBS.err'

    def __init__(self,**batch_settings):
        """
        notes:
            every list MUST be the same length
            args containing multiple settings should be set using list of dicts e.g. flags=[{'TOKAMAK':1,'STDOUT':None}] or settings_prec_mod=[{'root':"'some/file/path'"}] - this equates to a single workflow run with --flags TOKAMAK=1 STDOUT  
        args:
            batch_settings - kwargs storing lists which describe command line args for chosen workflow e.g. filepaths_in=[/some/path/1,some/path/2] if you want to launch a workflow twice, once with --filepaths_in=/some/path/1 and once with --filepaths_in=/some/path/2 
        """

        self.batch_settings=batch_settings
        for setting in self.batch_settings:
            self.number_runs=len(batch_settings[setting])

    def launch(self,workflow_filepath,system_name='TITAN'):
        """
        notes:
            the goal is to move through all the lists stored in batch_settings and generate a corresponding batch run script which calls some workflow with args defined in batch_settings
        args:
            workflow_filepath - filepath to file storing workflow to run e.g. support.dir_run_scripts / 'LOCUST_run.py'
            system_name - string identifier to choose from selection of environments stored as class attributes 
        """

        ################################# loop over each run in this set of batch runs

        for run in range(self.number_runs): #equivalent to the length of lists held in batch_settings

            ################################# create string holding batch script file settings

            system_batch_flags=''
            for batch_flag,value in Batch.batch_system[system_name]['flags'].items(): #add '-' or '--' to batch system flags e.g. num_processors --> --num_processorss                    
                batch_flag=''.join(['-',batch_flag]) if len(batch_flag)==1 else ''.join(['--',batch_flag]) 
                batch_flag='{} {}'.format(batch_flag,str(value)) if value else '{}'.format(batch_flag) #create strings for batch system flags                
                system_batch_flags=''.join([system_batch_flags,Batch.batch_system[system_name]['flag_command'],' ',batch_flag,'\n'])

            ################################# create arg string to pass to workflow script

            workflow_args='' #string to hold args passed to run script called from batch script
            for counter,(workflow_arg,workflow_settings) in enumerate(self.batch_settings.items()): #workflow_arg is the name of the arguement supplied to workflow when run from command line (e.g. filepath_input), workflow_settings are corresponding settings (e.g. /some/file/path)
                setting_values_this_run=workflow_settings[run] #at this point we have picked single element of list describing a particular setting over multiple runs e.g. compile flags 

                if setting_values_this_run:
                    if type(setting_values_this_run)==type({}): #if type is dict then these workflow settings are passed at command line differently - in the form --settings setting1=value1 setting2=value2
                        for setting,value in setting_values_this_run.items(): #if value1 above is a string, will need extra set of quotes to maintain continuity
                            if type(value)==type(''): #add some extra formatting to insert quotes to maintain continuity when calling workflows
                                setting_values_this_run[setting]='"{value}"'.format(value=value)

                        #parse args in form --settings setting1=value1 setting2=value2 setting3
                        workflow_args_this_setting=' '.join(['{}={}'.format(setting,value) if value else '{}'.format(setting) for setting,value in setting_values_this_run.items()])                               
                        workflow_args=' '.join([workflow_args,'--'+str(workflow_arg),workflow_args_this_setting])

                    else:

                        #parse args in form --settings setting1
                        workflow_args=' '.join([workflow_args,'--'+str(workflow_arg),str(setting_values_this_run)]) 

            ################################# combine to generate final batch script

            #create a dummy batch script file
            batch_file="""\
{shebang}
{system_batch_flags}

{command_execute_workflow} {workflow_filepath} {workflow_args}
            """.format(
            shebang=Batch.batch_system[system_name]['shebang'],
            system_batch_flags=system_batch_flags,
            command_execute_workflow=Batch.batch_system[system_name]['command_execute_workflow'],
            workflow_filepath=str(workflow_filepath),
            workflow_args=workflow_args
            ).lstrip(' ')

            batch_file_path=support.dir_run_scripts / '.Batch_submission_script_{}'.format(run)
            with open(batch_file_path,'w') as file:
                file.write(batch_file.lstrip(' '))

            #launch batch scripts here   
            try:
                subprocess.run([str(Batch.batch_system[system_name]['command_launch_batch']),str(batch_file_path)])
            except subprocess.CalledProcessError as err:
                raise(err)

#################################

##################################################################

###################################################################################################