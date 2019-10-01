#MARS_builder_build.py

'''
Samuel Ward
2909/2019
----
tools to build the mars_builder code
---
notes: 
    have separated MARS_builder_build, MARS_builder_build etc. just because we need to be somewhat specific in some areas
        ideally we would have a generic code 'build' which can take arbitrary permutations of compile flags
---
'''


##################################################################
#Preamble

import context

try:
    import sys
    import subprocess
    import pathlib
    import shlex
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import run_scripts.LOCUST_environment
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/LOCUST_environment.py could not be imported!\nreturning\n")
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

class MARS_builder_build:
    """
    class to hold a MARS_builder build

    notes:
        MARS_builder_build has its own instance of a LOCUST_environment
        ignores the UPHASE, MPHASE and LPHASE settings - these can be set in LOCUST prec_mod.f90
    usage:
        my_build=MARS_builder_build('TITAN') #initialise a build using the TITAN LOCUST_environment
        my_build.environment.display() #environment is instance of LOCUST_environment class 
        print(my_build.repo_URL) #if not specified at instantiation, default repo_URL and commit_hash are set - check LOCUST_IO/src/settings.py 
        my_build.flags_add(PRODN=None) #adds flags to this MARS_builder_build
        my_build.settings_mars_read_add(file="'/home/my_directory'")
        my_build.make(clean=True) #execute 'make clean'
        my_build.make() #make with flags stored in MARS_builder_build.flags
    """ 

    def __init__(self,system_name=None):
        """
        args:
            system_name - optional specify system to choose environment e.g. TITAN
        notes:
        """

        if system_name: self.environment=run_scripts.LOCUST_environment.LOCUST_environment(system_name) #create an environment for building
        self.settings_mars_read={}

    def flags_add(self,**flags):
        """
        append compile flags to flags currently stored in this build

        notes:    
            '-D' added by LOCUST_IO automatically, so do not worry about this - just use the flag name e.g. STDOUT, TOKAMAK=1
        args:
            flags - set of kwargs denoting compile flags
        usage:
            MARS_builder_build.flags_add('TOKAMAK'=1,'STDOUT'=None) #if flag does not need numerical value, set to None
        """
        if not hasattr(self,'flags'):
            self.flags={}
        for flag,value in flags.items():
            flag=''.join(['-D',flag])
            self.flags[flag]=value

        if 'TOKAMAK' in flags:
            if int(flags['TOKAMAK'])==1:
                self.tokamak='ITER'
            if int(flags['TOKAMAK'])==2:
                self.tokamak='AUG'
            if int(flags['TOKAMAK'])==3:
                self.tokamak='MAST-U'
            if int(flags['TOKAMAK'])==4:
                self.tokamak='MAST'
            if int(flags['TOKAMAK'])==5:
                self.tokamak='JET'
            if int(flags['TOKAMAK'])==6:
                self.tokamak='W7-X'
            if int(flags['TOKAMAK'])==7:
                self.tokamak='CTF-YORK'
            if int(flags['TOKAMAK'])==8:
                self.tokamak='DIII-D'
            if int(flags['TOKAMAK'])==9:
                self.tokamak='ZOW_ANALYTIC'
            if int(flags['TOKAMAK'])==10:
                self.tokamak='STEP'

    def flags_clear(self):
        """
        remove flags from this build 
        notes:
        """

        del(self.flags)

    def settings_mars_read_add(self,**settings_mars_read):
        """
        append mars_read.f90 settings to those currently stored in this build

        notes:    
        args:
            settings_mars_read - set of kwargs denoting changes to make to mars_read.f90
        usage:
        """
        if not hasattr(self,'settings_mars_read'):
            self.settings_mars_read={}
        for setting_mars_read,value in settings_mars_read.items():
            self.settings_mars_read[setting_mars_read]=value

    def settings_mars_read_clear(self):
        """
        remove settings_mars_read from this build 
        notes:
        """

        del(self.settings_mars_read)

    def make(self,directory=support.dir_run_scripts / 'mars_builder',clean=False):
        """
        args:
            directory - path to directory holding code to make
            clean - toggle whether to execute make clean
        notes:
            set compile flags using MARS_builder_build.flags_add
        """

        if clean:
            subprocess.run(shlex.split('make clean'),shell=False,cwd=str(directory))    

        else: #not making clean, so parse flags - two ways depending whether they hold numerical value e.g. -DSTDOUT vs -DTOKAMAK=1
            
            run_scripts.LOCUST_edit_var.LOCUST_edit_var(filename_in=directory / 'mars_read.f90',filename_out=directory / 'mars_read.f90',**self.settings_mars_read) #perform source code edits

            if not hasattr(self,'environment'):
                print("WARNING: MARS_builder_build.make does not contain environment - using settings.environment_default!")                
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