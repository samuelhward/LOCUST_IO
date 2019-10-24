#LOCUST_build.py

'''
Samuel Ward
27/09/2019
----
tools to retrieve and build the LOCUST code
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
    import run_scripts.build
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/build.py could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import run_scripts.environment
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/environment.py could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import run_scripts.LOCUST_edit_var
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/LOCUST_edit_var.py could not be imported!\nreturning\n")
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

class LOCUST_build(run_scripts.build.Build):
    """
    class to describe a LOCUST build

    notes:
        essentially a Build object but with git clone() functionality added
        LOCUST_build has its own instance of an environment
    usage:
        my_build=LOCUST_build('TITAN') #initialise a build using the TITAN environment
        my_build.environment.display() #environment is instance of Environment class 
        print(my_build.repo_URL) #if not specified at instantiation, default repo_URL and commit_hash are set - check LOCUST_IO/src/settings.py 
        print(my_build.commit_hash) #commit_hash is assumed to be latest hash of settings.branch_default_LOCUST branch
        my_build.flags_add(STDOUT=True,TOKAMAK=8,BRELAX=True,UNBOR=100,LEIID=8) #adds flags to this LOCUST_build
        my_build.source_code_mods_add(source_code_filename='prec_mod.f90',ThreadsPerBlock=56,some_string="'format like this'") #in prec_mod.f90 set ThreadsPerBlock to be declared=56 and string some_string to be initialised to 'format like this'
        my_build.clone(some_dir)
        my_build.make(directory=some_dir,clean=True) #execute 'make clean'
        my_build.make(directory=some_dir) #make with flags stored in LOCUST_build.flags
    """ 

    def __init__(self,system_name=None,repo_URL=settings.repo_URL_LOCUST,commit_hash=None):
        """
        args:
            system_name - optional specify system to choose environment e.g. TITAN
            repo_URL - SSH address of remote target git repo 
            commit_hash - commit hash identifying this build - defaults to latest commit on settings.branch_default_LOCUST branch
        notes:
        """

        #execute base class constructor to inherit required structures
        super().__init__(system_name=system_name)

        self.repo_URL=repo_URL
        if not commit_hash: 
            self.commit_hash=settings.commit_hash_default_LOCUST
            if not self.commit_hash: print("WARNING: commit_hash not set for LOCUST_build")

    def clone(self,directory=support.dir_locust):
        """
        notes:
            default behaviour is shallow clone from settings.repo_URL_LOCUST of settings.branch_default_LOCUST branch
            shallow clone is enabled if commit_hash has not been set or matches the latest commit hash of settings.branch_default_LOCUST branch
        args:
            directory - pathlib directory for where to clone LOCUST (new or empty dir), default to support.dir_locust
        """

        directory=pathlib.Path(directory) 
        if directory.exists():
            files_in_target_dir=directory.glob('*') #check for files in target directory
            files_in_target_dir=[file_in_target_dir for file_in_target_dir in files_in_target_dir] 
            if files_in_target_dir:
                print("ERROR: LOCUST_build.clone found files in target directory - please clone to empty dir!\nreturning\n")
                return
        else:
            directory.mkdir()

        if self.commit_hash==settings.commit_hash_default_LOCUST: shallow = True #shallow toggles whether or not to only clone latest commit

        command=[]
        command.append('git')
        command.append('clone')
        command.append(str(self.repo_URL))
        if shallow:
            command.append('--branch')
            command.append(settings.branch_default_LOCUST)
            command.append('--depth')
            command.append('1')

        subprocess.run(command,cwd=str(directory)) #code is now cloned
        directory_locust=directory / 'locust' #cloning adds additional folder
        
        if self.commit_hash and not shallow:
            command=['git','checkout','{commit_hash}'.format(commit_hash=self.commit_hash)]
            try:
                subprocess.run(command,shell=False,cwd=str(directory_locust))
            except subprocess.CalledProcessError as err:
                raise(err)

#################################

##################################################################

###################################################################################################