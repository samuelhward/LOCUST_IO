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

class LOCUST_build:
    """
    class to hold a LOCUST build

    notes:
        LOCUST_build has its own instance of a LOCUST_environment
    usage:
        my_build=LOCUST_build('TITAN') #initialise a build using the TITAN LOCUST_environment
        my_build.environment.display() #environment is instance of LOCUST_environment class 
        print(my_build.repo_URL) #if not specified at instantiation, default repo_URL and commit_hash are set - check LOCUST_IO/src/settings.py 
        print(my_build.commit_hash) #commit_hash is assumed to be latest hash of settings.branch_default_LOCUST branch
        my_build.add_flags(STDOUT=None,TOKAMAK=8,BRELAX=None,UNBOR=100,LEIID=8) #adds flags to this LOCUST_build
        my_build.make(clean=True) #execute 'make clean'
        my_build.make() #make with flags stored in LOCUST_build.flags
    """ 

    def __init__(self,system_name=None,repo_URL=settings.repo_URL_LOCUST,commit_hash=None):
        """
        args:
            system_name - optional specify system to choose environment e.g. TITAN
            repo_URL - SSH address of remote target git repo 
            commit_hash - commit hash identifying this build - defaults to latest commit on settings.branch_default_LOCUST branch
        notes:
        """

        self.environment=run_scripts.LOCUST_environment.LOCUST_environment(system_name=system_name) #create an environment for building
        self.repo_URL=repo_URL
        if not commit_hash: 
            self.commit_hash=settings.commit_hash_default_LOCUST
            if not self.commit_hash: print("WARNING: commit_hash not set for LOCUST_build")    

        if system_name: self.environment=run_scripts.LOCUST_environment.LOCUST_environment(system_name)

    def add_flags(self,**flags):
        """
        append LOCUST compile flags to flags currently stored in this LOCUST_build

        notes:    
            '-D' added by LOCUST_IO, so do not worry about this - just use the flag name e.g. STDOUT, TOKAMAK=1
        args:
            flags - set of kwargs denoting compile flags
        usage:
            LOCUST_build.add_flags('TOKAMAK'=1,'STDOUT'=None) #if flag does not need numerical value, set to None
        """
        if not hasattr(self,'flags'):
            self.flags={}
        for flag,value in flags.items():
            flag=''.join(['-D',flag])
            self.flags[flag]=value

    def clear_flags(self):
        """
        remove flags from this build 
        notes:
        """

        del(self.flags)

    def clone(self,directory=support.dir_locust):
        """
        notes:
            default behaviour is shallow clone from settings.repo_URL_LOCUST of settings.branch_default_LOCUST branch
            shallow clone is enabled if commit_hash has not been set or matches the latest commit hash of settings.branch_default_LOCUST branch
        args:
            directory - pathlib directory for where to clone LOCUST (make sure dir is empty), default to support.dir_locust
        """

        files_in_target_dir=pathlib.Path(directory).glob('*') #check for files in target directory
        files_in_target_dir=[file_in_target_dir for file_in_target_dir in files_in_target_dir] 
        if files_in_target_dir:
            print("ERROR: LOCUST_build.clone found files in target directory - please clone to empty dir!\nreturning\n")
            return

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
        directory=directory/'locust' #cloning adds additional folder
        
        if self.commit_hash and not shallow:
            command=['git','checkout','{commit_hash}'.format(commit_hash=self.commit_hash)]
            try:
                subprocess.run(command,shell=False,cwd=str(directory))
            except subprocess.CalledProcessError as err:
                raise(err)

    def make(self,directory=support.dir_locust,clean=False):
        """
        args:
            directory - path to directory holding code to make
            clean - toggle whether to execute make clean
        notes:
            set compile flags using LOCUST_build.add_flags
        """

        if clean:
            subprocess.run(shlex.split('make clean'),shell=False,cwd=str(directory))    

        else: #not making clean, so parse flags - two ways depending whether they hold numerical value e.g. -DSTDOUT vs -DTOKAMAK=1
            
            if not hasattr(self,'environment'):
                print("WARNING: LOCUST_build.make does not contain environment - using settings.environment_default!")                
                self.environment=run_scripts.LOCUST_environment.LOCUST_environment(settings.system_default)

            flags_=' '.join(['{}={}'.format(flag,value) if value else '{}'.format(flag) for flag,value in self.flags.items()])
            command=' '.join([self.environment.create_command(),'; make','FLAGS={}'.format(shlex.quote(flags_))])

            try:
                with subprocess.Popen(command,shell=True,cwd=str(directory)) as proc: #stdin=PIPE, stdout=PIPE, stderr=STDOUT
                    pass
            except subprocess.CalledProcessError as err:
                raise(err)

#################################

##################################################################

###################################################################################################