#LOCUST_run.py

'''
Samuel Ward
28/09/2019
----
tools to run the LOCUST code
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

if __name__=='__main__':

    try:
        import context
    except:
        raise ImportError("ERROR: context.py could not be imported!\nreturning\n")
        sys.exit(1)

try:
    import run_scripts.LOCUST_environment
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/LOCUST_environment.py could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import run_scripts.LOCUST_build
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/LOCUST_build.py could not be imported!\nreturning\n")
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

class LOCUST_run:
    """
    class to run LOCUST

    notes:
        contains a LOCUST_build and a LOCUST_environment which store relevant settings for this run
        both this object and it's stored objects each have separate environments - in case one wants to run with a different environment that they built with
        when editing prec_mod.f90 source code using settings_mars_read, strings should be passed literally (see below)     
        python LOCUST_run.py --filepath_in 'some string formatted like this' --flags TOKAMAK=8 --settings_prec_mod a_string \'"should be formatted like this"\' a_variable 'like_this'
    usage:
        python LOCUST_run.py --flags TOKAMAK=8 BRELAX UNBOR=100 LEIID=8 OPENMESH OPENTRACK PFCMOD TOKHEAD JXB2 PROV GEQDSKFIX1 PITCHCUR EGSET=100000.0D0 EBASE UHST LNLBT B3D B3D_EX SPLIT NOINF SMALLEQ BP TIMAX=0.023D0 STDOUT
        some_run=LOCUST_run(system_name='TITAN',flags=flags) #by default this will clone LOCUST to the LOCUST_IO/LOCUST folder, run it on LOCUST_IO/data/input_files and dump to LOCUST_IO/data/output_files
        some_run.run() #this will do all stages of a LOCUST run including cloning, building, running and cleaning up afterwards
    """ 

    def __init__(self,dir_LOCUST=support.dir_locust,dir_input=support.dir_input_files,dir_output=support.dir_output_files,dir_cache=support.dir_cache_files,system_name='TITAN',repo_URL=settings.repo_URL_LOCUST,commit_hash=None,settings_prec_mod={},flags={}):
        """
        notes:
            most information stored in LOCUST_run.environment and LOCUST_run.build, most init args are to init these instances
        args:
            dir_LOCUST - directory to temporarily store source code
            dir_input - directory to read input data from 
            dir_output - directory to write results to 
            dir_cache - directory to read cache data from 
            system_name - string identifier to choose from selection of environments stored as class attributes 
            repo_URL - SSH address of remote target git repo 
            commit_hash - commit hash identifying this build - defaults to latest commit on settings.branch_default_LOCUST branch
            settings_prec_mod - dict denoting variable names and values to set them to within prec_mod.f90
            flags - dict denoting compile flags (no '-D' please e.g. STDOUT TOKAMAK=3)
        """
        
        #convert strings to paths just in case supplied at command line
        dir_LOCUST=pathlib.Path(dir_LOCUST)
        dir_input=pathlib.Path(dir_input)
        dir_output=pathlib.Path(dir_output)
        dir_cache=pathlib.Path(dir_cache)

        files_in_dir_LOCUST=dir_LOCUST.glob('*') #check for files in dir_LOCUST
        if [file_in_dir_LOCUST for file_in_dir_LOCUST in files_in_dir_LOCUST]:
            print("ERROR: new LOCUST_run found files in target dir_LOCUST - please specify empty dir!\nreturning\n")
            return

        self.dir_LOCUST=dir_LOCUST
        self.dir_input=dir_input
        self.dir_output=dir_output
        self.dir_cache=dir_cache
        self.environment=run_scripts.LOCUST_environment.LOCUST_environment(system_name=system_name) #create an environment for running
        self.build=run_scripts.LOCUST_build.LOCUST_build(system_name=system_name,repo_URL=repo_URL,commit_hash=commit_hash)
        self.build.flags_add(**flags) 
        self.build.settings_prec_mod_add(**settings_prec_mod) 

    def run(self):
        """
        notes:
        args:
        """

        #create output directory and directory to hold LOCUST if do not already exist
        if not self.dir_LOCUST.exists(): self.dir_LOCUST.mkdir()
        if not self.dir_output.exists(): self.dir_output.mkdir()

        #find where LOCUST will store its files - InputFiles, OutputFiles and CacheFiles
        if 'TOKHEAD' in self.build.flags: #user specified tokhead, which changes directory structure and must be accounted for
            tokhead='locust.'+self.build.tokamak
        else:
            tokhead='locust'

        #find root 
        if 'root' in self.build.settings_prec_mod: #variable is set in LOCUST source - prec_mod.f90
            root=pathlib.Path(self.build.settings_prec_mod['root'])
        else:
            root=pathlib.Path(settings.LOCUST_dir_root_default) #default root set within LOCUST prec_mod.f90

        #retrieve code
        self.build.clone(dir_LOCUST=self.dir_LOCUST)
        dir_LOCUST=self.dir_LOCUST / 'locust' #from now on all code is within second folder due to clone

        #compile (source files are edited within this step)
        self.build.make(dir_LOCUST=dir_LOCUST,clean=True)
        self.build.make(dir_LOCUST=dir_LOCUST,clean=False)

        #create root and child directories if do not exist
        for directory in [root,(root / settings.username),
                         (root / settings.username / tokhead),
                         (root / settings.username / tokhead / settings.LOCUST_dir_inputfiles_default),
                         (root / settings.username / tokhead / settings.LOCUST_dir_cachefiles_default),
                         (root / settings.username / tokhead / settings.LOCUST_dir_outputfiles_default)]:
            if not directory.exists():
                directory.mkdir()

        #copy input and cache files to correct location
        for file in self.dir_input.glob('*'): #move all input files to correct location
            subprocess.run(shlex.split('cp {file} {root}'.format(file=str(file),root=str(root / settings.username / tokhead / settings.LOCUST_dir_inputfiles_default))),shell=False)
        for file in self.dir_cache.glob('*'): #move all cache files to correct location
            subprocess.run(shlex.split('cp {file} {root}'.format(file=str(file),root=str(root / settings.username / tokhead / settings.LOCUST_dir_cachefiles_default))),shell=False)

        #run LOCUST
        command=' '.join([self.environment.create_command(),'; ./locust'])
        try:
            pass
            with subprocess.Popen(command,shell=True,cwd=str(dir_LOCUST)) as proc: #stdin=PIPE, stdout=PIPE, stderr=STDOUT
                pass
        except subprocess.CalledProcessError as err:
            print("ERROR: LOCUST_run failed to run LOCUST!\nreturning\n")
            #raise(err)

        #retrieve output
        subprocess.run(shlex.split('mv {prec_mod} {dir_output}'.format(prec_mod=str(dir_LOCUST/'prec_mod.f90'),dir_output=str(self.dir_output))),shell=False) #retain copy of prec_mod.f90            
        for file in (root / settings.username / tokhead / settings.LOCUST_dir_outputfiles_default).glob('*'):
            subprocess.run(shlex.split('mv {file} {dir_output}'.format(file=str(file),dir_output=str(self.dir_output))),shell=False)

        #remove source code, input file and cache file folders
        subprocess.run(shlex.split('rm -r {}'.format(str(root / settings.username))),shell=False) #delete folder up since added new child folder when cloning
        subprocess.run(shlex.split('rm -rf {}'.format(str(dir_LOCUST.parents))),shell=False) #delete folder holding LOCUST

if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser(description='run LOCUST from command line in Python!')

    parser.add_argument('--dir_LOCUST',type=str,action='store',default=support.dir_locust,dest='dir_LOCUST',help="directory to temporarily store source code",required=False)
    parser.add_argument('--dir_in',type=str,action='store',default=support.dir_input_files,dest='dir_input',help="directory to read input data from",required=False)
    parser.add_argument('--dir_out',type=str,action='store',default=support.dir_output_files,dest='dir_output',help="directory to write results to",required=False)
    parser.add_argument('--dir_cache',type=str,action='store',default=support.dir_cache_files,dest='dir_cache',help="directory to read cache data from",required=False)
    parser.add_argument('--sys_name',type=str,action='store',default='TITAN',dest='system_name',help="computer system currently running on e.g. \'TITAN\'",required=False)
    parser.add_argument('--repo_URL',type=str,action='store',default=settings.repo_URL_LOCUST,dest='repo_URL',help="URL of LOCUST git repository",required=False)
    parser.add_argument('--commit',type=str,action='store',default=None,dest='commit_hash',help="optional commit hash of specific build",required=False)
    parser.add_argument('--settings_prec_mod',nargs='+',type=str,action='store',default={},dest='settings_prec_mod',help="variable names and values to set them to within prec_mod.f90",required=False)
    parser.add_argument('--flags',nargs='+',type=str,action='store',default={},dest='flags',help="compile flags",required=False)
    
    args=parser.parse_args()

    #provide some extra parsing steps to flags, prec_mod settings and any similar dict-like input arguements
    parsed_flags=[flag.split('=') for flag in args.flags]    
    flags={}
    for flag in parsed_flags:
        if len(flag)>1:
            flags[flag[0]]=flag[1]
        else:
            flags[flag[0]]=None

    parsed_settings_prec_mod=[setting_prec_mod.split('=') for setting_prec_mod in args.settings_prec_mod]    
    settings_prec_mod={}
    for setting in parsed_settings_prec_mod:
        if len(setting)>1:
            settings_prec_mod[setting[0]]=setting[1]
        else:
            settings_prec_mod[setting[0]]=None


    this_run=LOCUST_run(system_name=args.system_name,repo_URL=args.repo_URL,commit_hash=args.commit_hash,dir_LOCUST=args.dir_LOCUST,dir_input=args.dir_input,dir_output=args.dir_output,dir_cache=args.dir_cache,settings_prec_mod=settings_prec_mod,flags=flags)
    this_run.run()

#################################

##################################################################

###################################################################################################