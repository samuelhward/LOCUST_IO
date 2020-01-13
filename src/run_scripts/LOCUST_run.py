#LOCUST_run.py

'''
Samuel Ward
28/09/2019
----
tools to run the LOCUST code
---
notes: 
    contains the tools to perform a single simple run of LOCUST
    custom steps can be added via add_command
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
    import run_scripts.LOCUST_build
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/LOCUST_build.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import run_scripts.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/utils.py could not be imported!\nreturning\n")
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

class LOCUST_run(run_scripts.workflow.Workflow):
    """
    defines a workflow which performs a simple LOCUST simulation

    notes:
        adds workflow components to workflow class for execution of a single standard LOCUST run
    usage:
        python LOCUST_run.py --flags TOKAMAK=8 BRELAX UNBOR=100 LEIID=8 OPENMESH OPENTRACK PFCMOD TOKHEAD JXB2 PROV GEQDSKFIX1 PITCHCUR EGSET=100000.0D0 EBASE UHST LNLBT B3D B3D_EX SPLIT NOINF SMALLEQ BP TIMAX=0.023D0 STDOUT
        some_run=LOCUST_run(system_name='TITAN',flags=flags) #by default this will clone LOCUST to the LOCUST_IO/LOCUST folder, run it on LOCUST_IO/data/input_files and dump to LOCUST_IO/data/output_files
        some_run.run() #this will execute all stages of a LOCUST run including cloning, building, running and cleaning up afterwards
    """ 

    def __init__(self,dir_LOCUST=support.dir_locust,dir_input=support.dir_input_files,dir_output=support.dir_output_files,dir_cache=support.dir_cache_files,system_name='TITAN',repo_URL=settings.repo_URL_LOCUST,commit_hash=None,settings_prec_mod={},flags={},*args,**kwargs):
        """
        notes:
            most information stored in LOCUST_run.environment and LOCUST_run.build, most init args are to init these instances
        args:
            dir_LOCUST - directory to temporarily store source code (must not already exist)
            dir_input - directory to read input data from 
            dir_output - directory to write results to 
            dir_cache - directory to read cache data from 
            system_name - string identifier to choose from selection of environments stored as class attributes 
            repo_URL - SSH address of remote target git repo 
            commit_hash - commit hash identifying this build - defaults to latest commit on settings.branch_default_LOCUST branch
            settings_prec_mod - dict denoting variable names and values to set them to within prec_mod.f90
            flags - dict denoting compile flags (no '-D' please e.g. STDOUT TOKAMAK=3)
        """

        #execute base class constructor to inherit required structures
        super().__init__()

        ################################# first generate class data that will be needed in workflow

        dir_LOCUST=pathlib.Path(dir_LOCUST) #convert strings to paths just in case supplied at command line
        dir_input=pathlib.Path(dir_input)
        dir_output=pathlib.Path(dir_output)
        dir_cache=pathlib.Path(dir_cache)
        if dir_LOCUST.is_dir(): #check if dir_LOCUST already exists
            print("ERROR: dir_LOCUST already exists for this run - please specify new dir!\nreturning\n")
            sys.exit(1)
        self.dir_LOCUST=dir_LOCUST if dir_LOCUST.is_absolute() else print("ERROR: LOCUST_run requires absolute dir_LOCUST path!")
        self.dir_input=dir_input if dir_input.is_absolute() else print("ERROR: LOCUST_run requires absolute dir_input path!")
        self.dir_output=dir_output if dir_output.is_absolute() else print("ERROR: LOCUST_run requires absolute dir_output path!")
        self.dir_cache=dir_cache if dir_cache.is_absolute() else print("ERROR: LOCUST_run requires absolute dir_cache path!")
        
        self.environment=run_scripts.environment.Environment(system_name=system_name) #create a runtime environment for this workflow
        
        self.build=run_scripts.LOCUST_build.LOCUST_build(system_name=system_name,repo_URL=repo_URL,commit_hash=commit_hash) #define a build for this workflow
        self.build.flags_add(**flags) 
        self.build.source_code_mods_add(source_code_filename='prec_mod.f90',**settings_prec_mod)

        #find where LOCUST will store its files - InputFiles, OutputFiles and CacheFiles
        if 'TOKHEAD' in self.build.flags: #user specified tokhead, which changes directory structure and must be accounted for
            self.tokhead='locust.'+self.build.tokamak
        else:
            self.tokhead='locust'
        #find root 
        if 'root' in self.build.source_code_mods['prec_mod.f90']: #variable to set in LOCUST source file prec_mod.f90
            self.root=pathlib.Path(self.build.source_code_mods['prec_mod.f90']['root'].strip('\'')) #user should have added extra quotes to properly edit prec_mod.f90 file
        else:
            self.root=pathlib.Path(settings.LOCUST_dir_root_default) #the default root set within LOCUST prec_mod.f90

        ################################# now make commands (defined below in this class) available to this workflow (and state position in execution order)
        
        self.add_command(command_name='mkdir',command_function=self.setup_LOCUST_dirs,position=1) #add new workflow stages
        self.add_command(command_name='get_code',command_function=self.get_code,position=2)
        self.add_command(command_name='make',command_function=self.make_code,position=3)
        self.add_command(command_name='get_input',command_function=self.get_inputs,position=4)
        self.add_command(command_name='run_code',command_function=self.run_code,position=5)
        self.add_command(command_name='get_output',command_function=self.get_outputs,position=6)
        self.add_command(command_name='get_cache',command_function=self.get_cache,position=7)
        self.add_command(command_name='cleanup',command_function=self.cleanup,position=8)

    def setup_LOCUST_dirs(self,*args,**kwargs):
        """
        LOCUST_run stage for initialising simulation directory structure

        notes:
        """

        #create output directory if does not already exist
        if not self.dir_output.is_dir(): self.dir_output.mkdir()
        #create LOCUST directory - must not already exis
        if self.dir_LOCUST.is_dir(): 
            print("ERROR: {workflow_name}.setup_LOCUST_dirs() found previous dir_LOCUST -  must not already exist!\nreturning\n".format(workflow_name=self.workflow_name))
            return
        else:
            self.dir_LOCUST.mkdir()

        #create self.root and child directories if do not exist
        for directory in [
                         (self.root / settings.username / self.tokhead / settings.LOCUST_dir_inputfiles_default),
                         (self.root / settings.username / self.tokhead / settings.LOCUST_dir_cachefiles_default),
                         (self.root / settings.username / self.tokhead / settings.LOCUST_dir_outputfiles_default)]:
            if not directory.exists():
                directory.mkdir(parents=True)

    def get_code(self,*args,**kwargs):
        """
        LOCUST_run stage for retrieving code

        notes:
        """

        #retrieve code
        try:
            self.build.clone(directory=self.dir_LOCUST)
        except:
            print("ERROR: {workflow_name}.get_code() could not clone to {directory}!\nreturning\n".format(directory=self.dir_LOCUST,workflow_name=self.workflow_name))
            return 

    def get_inputs(self,*args,**kwargs):
        """
        LOCUST_run stage for copying inputs to correct location

        notes:
        """

        #copy input and cache files to correct location
        for file in self.dir_input.glob('*'): #move all input files to correct location
            subprocess.run(shlex.split('cp {file} {inputdir}'.format(file=str(file),inputdir=str(self.root / settings.username / self.tokhead / settings.LOCUST_dir_inputfiles_default))),shell=False)
        for file in self.dir_cache.glob('*'): #move all cache files to correct location
            subprocess.run(shlex.split('cp {file} {cachedir}'.format(file=str(file),cachedir=str(self.root / settings.username / self.tokhead / settings.LOCUST_dir_cachefiles_default))),shell=False)

    def make_code(self,*args,**kwargs):
        """
        LOCUST_run stage for executing code build

        notes:
        """

        self.build.make(directory=self.dir_LOCUST / 'locust',clean=True) #compile (source files are edited within this step)
        self.build.make(directory=self.dir_LOCUST / 'locust',clean=False) #append /locust as this introduced in cloning stage

    def run_code(self,*args,**kwargs):
        """
        LOCUST_run stage for executing LOCUST

        notes:
        """

        #run LOCUST
        command=' '.join([self.environment.create_command_string(),'; ./locust'])
        try:
            pass
            subprocess.call(command,shell=True,cwd=str(self.dir_LOCUST / 'locust'))# as proc: #stdin=PIPE, stdout=PIPE, stderr=STDOUT

        except subprocess.CalledProcessError as err:
            print("ERROR: {workflow_name}.run_code() failed to run LOCUST!\nreturning\n".format(workflow_name=self.workflow_name))
            #raise(err)

    def get_outputs(self,*args,**kwargs):
        """
        LOCUST_run stage for retrieving run outputs

        notes:
        """

        #retrieve output
        subprocess.run(shlex.split('mv {prec_mod} {dir_output}'.format(prec_mod=str(self.dir_LOCUST / 'locust' / 'prec_mod.f90'),dir_output=str(self.dir_output))),shell=False) #retain copy of prec_mod.f90            
        for file in (self.root / settings.username / self.tokhead / settings.LOCUST_dir_outputfiles_default).glob('*'):
            subprocess.run(shlex.split('mv {file} {dir_output}'.format(file=str(file),dir_output=str(self.dir_output))),shell=False)

    def get_cache(self,*args,**kwargs):
        """
        LOCUST_run stage for retrieving run cache

        notes:
            only retrieves mesh data
        """

        #retrieve cache
        for file in (self.root / settings.username / self.tokhead / settings.LOCUST_dir_cachefiles_default).glob('MESH*'):
            subprocess.run(shlex.split('mv {file} {dir_cache}'.format(file=str(file),dir_cache=str(self.dir_cache))),shell=False)
        for file in (self.root / settings.username / self.tokhead / settings.LOCUST_dir_cachefiles_default).glob('LINK*'):
            subprocess.run(shlex.split('mv {file} {dir_cache}'.format(file=str(file),dir_cache=str(self.dir_cache))),shell=False)

    def cleanup(self,*args,**kwargs):
        """
        LOCUST_run stage for cleaning up temporary directory structure and code

        notes:
        """

        #remove source code and temporary input/cache file folders
        #currently works by deleting specific known target folders then deleting parent directories IFF now empty
        for directory in [ 
                         (self.root / settings.username / self.tokhead / settings.LOCUST_dir_inputfiles_default),
                         (self.root / settings.username / self.tokhead / settings.LOCUST_dir_cachefiles_default),
                         (self.root / settings.username / self.tokhead / settings.LOCUST_dir_outputfiles_default)]:
            try:
                subprocess.run(shlex.split('rm -r {}'.format(str(directory))),shell=False) #delete directories storing temporary InputFiles/OutputFiles/CacheFiles
            except:
                pass
        try:
            (self.root / settings.username / self.tokhead).rmdir() #delete parent folders IFF empty
        except:
            pass
        try:
            (self.root / settings.username).rmdir() 
        except:
            pass

        try: #now remove folder containing LOCUST repo
            subprocess.run(shlex.split('rm -rf {}'.format(str(self.dir_LOCUST / 'locust' / '.git'))),shell=False) #delete folder holding LOCUST git repo
        except:
            print("WARNING: .git repo not found in {directory} during {workflow_name}.cleanup!".format(directory=str(self.dir_LOCUST / 'locust'),workflow_name=self.workflow_name)) 

        subprocess.run(shlex.split('rm -r {}'.format(str(self.dir_LOCUST / 'locust'))),shell=False) #delete specific folder holding rest of LOCUST source   
        self.dir_LOCUST.rmdir() #delete original containing folder IFF empty

##################################################################

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

    #provide some extra parsing steps to flags, prec_mod settings and any similar dict-like input arguments
    args.settings_prec_mod=run_scripts.utils.command_line_arg_parse_dict(args.settings_prec_mod)
    args.flags=run_scripts.utils.command_line_arg_parse_dict(args.flags)

    this_run=LOCUST_run(system_name=args.system_name,repo_URL=args.repo_URL,commit_hash=args.commit_hash,dir_LOCUST=args.dir_LOCUST,dir_input=args.dir_input,dir_output=args.dir_output,dir_cache=args.dir_cache,settings_prec_mod=args.settings_prec_mod,flags=args.flags)
    this_run.run()

#################################

##################################################################

###################################################################################################