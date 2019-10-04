#LOCUST_batch.py

'''
Samuel Ward
03/10/2019
----
tools to batch submit many LOCUST_runs 
---
notes: 
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
    import run_scripts.LOCUST_run
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/LOCUST_run.py could not be imported!\nreturning\n")
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

class LOCUST_batch:
    """
    class to batch run many LOCUST_runs

    notes:

    usage:
    """ 

    def __init__(self,dirs_LOCUST=[support.dir_locust],dirs_input=[support.dir_input_files],dirs_output=[support.dir_output_files],dirs_cache=[support.dir_cache_files],system_names=['TITAN'],repo_URLs=[settings.repo_URL_LOCUST],commit_hashes=[None],settings_prec_mod=[{}],flags=[{}]):
        """
        notes:
            args list exactly same as LOCUST_run except all provided as lists
            every list MUST be the same length 
        args:
            dirs_LOCUST - list of directories for temporarily storing source code
            dirs_input - list of directories to read input data from 
            dirs_output - list of directories to write results to 
            dirs_cache - list of directories to read cache data from 
            system_names - list of string identifiers to choose from selection of environments stored as class attributes 
            repo_URLs - list of SSH addresses of remote target git repo 
            commit_hashes - list of commit hashes identifying this build - defaults to latest commit on settings.branch_default_LOCUST branch
            settings_prec_mod - list of dicts denoting variable names and values to set them to within prec_mod.f90
            flags - list of dicts denoting compile flags (no '-D' please e.g. STDOUT TOKAMAK=3)
        """

        self.dirs_LOCUST=dirs_LOCUST
        self.dirs_input=dirs_input
        self.dirs_output=dirs_output
        self.dirs_cache=dirs_cache
        self.system_names=system_names
        self.repo_URLs=repo_URLs
        self.commit_hashes=commit_hashes
        self.settings_prec_mod=settings_prec_mod
        self.flags=flags
        self.system_batch_names_avail=['SLURM','PBS',None]

    def launch(self,system_batch_name='SLURM',**system_batch_flags):
        """
        notes:
        args:
            system_batch_name - type of batch system installed on running machine e.g. SLURM or PBS...None specifies interactive runs (useful for many short runs)
            system_batch_flags - batch system runtime flags (given to sbatch)
        """

        supplied_system_batch_flags=True if system_batch_flags else False

        for (counter,
            (dir_LOCUST,
            dir_input,
            dir_output,
            dir_cache,
            system_name,
            repo_URL,
            commit_hash,
            set_of_prec_mod_settings,
            set_of_flags)) in enumerate(zip(
            self.dirs_LOCUST,
            self.dirs_input,
            self.dirs_output,
            self.dirs_cache,
            self.system_names,
            self.repo_URLs,
            self.commit_hashes,
            self.settings_prec_mod,
            self.flags)):

            ##################################################################
            if system_batch_name=='SLURM':

                shebang='#!/bin/bash'

                if not supplied_system_batch_flags: #add default batch system flags
                    system_batch_flags={}
                    system_batch_flags['job-name']='LOCUST_{}'.format(counter)
                    system_batch_flags['gres']='gpu:8'
                    system_batch_flags['partition']='gpu_p100_titan'
                    system_batch_flags['cpus-per-task']=1
                    system_batch_flags['mail-user']='samuel.ward@iter.org'
                    system_batch_flags['mail-type']='END,FAIL,TIME_LIMIT'
                    system_batch_flags['o']='LOCUST_SLURM_{}.out'.format(counter)
                    system_batch_flags['e']='LOCUST_SLURM_{}.err'.format(counter)
 
                SLURM_flags='' #turn batch system flags into strings
                for batch_flag,value in system_batch_flags.items(): #add '-' or '--' to batch system flags e.g. num_processors --> --num_processorss                    
                    batch_flag=''.join(['-',batch_flag]) if len(batch_flag)==1 else ''.join(['--',batch_flag]) 
                    batch_flag='{}={}'.format(batch_flag,str(value)) if value else '{}'.format(batch_flag) #create strings for batch system flags                
                    SLURM_flags=''.join([SLURM_flags,'#SBATCH ',batch_flag,'\n'])


                #need to provide extra formatting since we are calling python from command line
                for setting,value in set_of_prec_mod_settings.items(): #strings will need extra set of quotes to maintain continuity
                    if type(value)==type(''): #since we need to match the input arguement formatting of LOCUST_run.py - which requires multiple quotes to format correctly
                        set_of_prec_mod_settings[setting]='"{value}"'.format(value=value)
                    else:
                        set_of_prec_mod_settings[setting]="'{value}'".format(value=value)

                LOCUST_run_flags="""\
--dir_LOCUST {dir_LOCUST} 
--dir_in {dir_input} 
--dir_out {dir_output} 
--dir_cache {dir_cache} 
--sys_name {system_name} 
--repo_URL {repo_URL} 
{commit} {commit_hash} 
{settings_prec_mod} {set_of_prec_mod_settings} 
{flags} {set_of_flags} 
                """.format(
                    dir_LOCUST=dir_LOCUST,
                    dir_input=dir_input,
                    dir_output=dir_output,
                    dir_cache=dir_cache,
                    system_name=system_name,
                    repo_URL=repo_URL,
                    commit='--commit' if commit_hash else '', #do not print if empty/not supplied
                    commit_hash=commit_hash if commit_hash else '',
                    settings_prec_mod='--settings_prec_mod ' if set_of_prec_mod_settings else '', #do not print if empty/not supplied
                    set_of_prec_mod_settings=' '.join(['{} {}'.format(setting,value) for setting,value in set_of_prec_mod_settings.items()]),
                    flags='--flags ' if set_of_flags else '', #do not print if empty/not supplied
                    set_of_flags=' '.join(['{}={}'.format(flag,value) if value else '{}'.format(flag) for flag,value in set_of_flags.items()])
                    ).replace('\n','')
                        
                #create a dummy batch script file
                batch_file="""\
{shebang}
{SLURM_flags}

python3 {LOCUST_run} {LOCUST_run_flags}
                """.format(
                shebang=shebang,
                SLURM_flags=SLURM_flags,
                LOCUST_run=support.dir_run_scripts / 'LOCUST_run.py',
                LOCUST_run_flags=LOCUST_run_flags
                ).lstrip(' ')
 
                batch_file_path=support.dir_run_scripts / '.LOCUST_batch_submission_script_{}'.format(counter)
                with open(batch_file_path,'w') as file:
                    file.write(batch_file.lstrip(' '))

                #launch batch scripts here   
                try:
                    subprocess.run(['sbatch',str(batch_file_path)])
                except subprocess.CalledProcessError as err:
                    raise(err)

            ##################################################################
            elif system_batch_name=='PBS':

                shebang='#!/bin/bash'

                if not system_batch_flags: #add default batch system flags
                    system_batch_flags={}
                    system_batch_flags['V']=None
                    system_batch_flags['l']='nodes=1:ppn=16'
                    system_batch_flags['q']='gpu'
                    system_batch_flags['N']='LOCUST_{}'.format(counter)
                    system_batch_flags['cpus-per-task']=1
                    system_batch_flags['o']='LOCUST_PBS_{}.out'.format(counter)
                    system_batch_flags['e']='LOCUST_PBS_{}.err'.format(counter)
 
                PBS_flags='' #turn batch system flags into strings
                for batch_flag,value in system_batch_flags.items(): #add '-' or '--' to batch system flags e.g. num_processors --> --num_processorss                    
                    batch_flag=''.join(['-',batch_flag]) if len(batch_flag)==1 else ''.join(['--',batch_flag]) 
                    batch_flag='{}={}'.format(batch_flag,str(value)) if value else '{}'.format(batch_flag) #create strings for batch system flags                
                    PBS_flags=''.join([PBS_flags,'#PBS ',batch_flag,'\n'])


                #need to provide extra formatting since we are calling python from command line
                for setting,value in set_of_prec_mod_settings.items(): #strings will need extra set of quotes to maintain continuity
                    if type(value)==type(''): #since we need to match the input arguement formatting of LOCUST_run.py - which requires multiple quotes to format correctly
                        set_of_prec_mod_settings[setting]='"{value}"'.format(value=value)
                    else:
                        set_of_prec_mod_settings[setting]="'{value}'".format(value=value)

                LOCUST_run_flags="""\
--dir_LOCUST {dir_LOCUST} 
--dir_in {dir_input} 
--dir_out {dir_output} 
--dir_cache {dir_cache} 
--sys_name {system_name} 
--repo_URL {repo_URL} 
{commit} {commit_hash} 
{settings_prec_mod} {set_of_prec_mod_settings} 
{flags} {set_of_flags} 
                """.format(
                    dir_LOCUST=dir_LOCUST,
                    dir_input=dir_input,
                    dir_output=dir_output,
                    dir_cache=dir_cache,
                    system_name=system_name,
                    repo_URL=repo_URL,
                    commit='--commit' if commit_hash else '', #do not print if empty/not supplied
                    commit_hash=commit_hash if commit_hash else '',
                    settings_prec_mod='--settings_prec_mod ' if set_of_prec_mod_settings else '', #do not print if empty/not supplied
                    set_of_prec_mod_settings=' '.join(['{} {}'.format(setting,value) for setting,value in set_of_prec_mod_settings.items()]),
                    flags='--flags ' if set_of_flags else '', #do not print if empty/not supplied
                    set_of_flags=' '.join(['{}={}'.format(flag,value) if value else '{}'.format(flag) for flag,value in set_of_flags.items()])
                    ).replace('\n','')
                        
                #create a dummy batch script file
                batch_file="""\
{shebang}
{PBS_flags}

python3 {LOCUST_run} {LOCUST_run_flags}
                """.format(
                shebang=shebang,
                PBS_flags=PBS_flags,
                LOCUST_run=support.dir_run_scripts / 'LOCUST_run.py',
                LOCUST_run_flags=LOCUST_run_flags
                ).lstrip(' ')
 
                batch_file_path=support.dir_run_scripts / '.LOCUST_batch_submission_script_{}'.format(counter)
                with open(batch_file_path,'w') as file:
                    file.write(batch_file.lstrip(' '))

                #launch batch scripts here   
                try:
                    subprocess.run(['qsub',str(batch_file_path)])
                except subprocess.CalledProcessError as err:
                    raise(err)


            elif system_batch_name is None:
                pass
            else:
                print("ERROR: LOCUST_batch.launch given incorrect system_batch_name - options are {}\nreturning\n".format(self.system_batch_names_avail))
                return

#################################

##################################################################

###################################################################################################