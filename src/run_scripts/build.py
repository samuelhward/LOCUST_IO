#Build.py

'''
Samuel Ward
27/09/2019
----
tools to retrieve and build generic code
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

class Build:
    """
    class to hold a generic build of code

    notes:
        Build has its own instance of an environment
    usage:
        my_build=Build('TITAN') #initialise a build using the TITAN environment
        my_build.environment.display() #environment is instance of Environment class 
        my_build.flags_add(STDOUT=True,TOKAMAK=8,BRELAX=True,UNBOR=100,LEIID=8) #adds flags to this Build
        my_build.source_code_mods_add(source_code_filename='prec_mod.f90',ThreadsPerBlock=56,some_string="'format like this'") #in prec_mod.f90 set ThreadsPerBlock to be declared=56 and string some_string to be initialised to 'format like this'
        my_build.make(directory=some_dir,clean=True) #execute 'make clean'
        my_build.make(directory=some_dir) #make with flags stored in Build.flags
    """ 

    def __init__(self,system_name=None):
        """
        args:
            system_name - optional specify system to choose environment e.g. TITAN
        notes:
        """

        if system_name: self.environment=run_scripts.environment.Environment(system_name) #create an environment for building
        self.source_code_mods={}  

    def flags_add(self,**flags):
        """
        append compile flags to flags currently stored in this Build

        notes:    
            '-D' added by LOCUST_IO automatically, so do not worry about this - just use the flag name e.g. STDOUT, TOKAMAK=1
        args:
            flags - set of kwargs denoting compile flags
        usage:
            Build.flags_add('TOKAMAK'=1,'STDOUT'=True) #if flag does not need numerical value, set to True
        """

        if not hasattr(self,'flags'):
            self.flags={}
        for flag,value in flags.items():
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

    def source_code_mods_add(self,source_code_filename,**source_code_mods):
        """
        append source code modification settings to those currently stored in this Build

        notes:    
            must be careful that source code file does not contain a variable called source_code_filename!
        args:
            source_code_filename - target filename to apply source code modifications
            source_code_mods - set of kwargs denoting changes to make to individual files
        usage:
        """

        if source_code_filename not in self.source_code_mods.keys():
            self.source_code_mods[source_code_filename]={}
        for source_code_mod,setting in source_code_mods.items():
            self.source_code_mods[source_code_filename][source_code_mod]=setting

    def source_code_mods_clear(self):
        """
        remove source code modification settings from this build 
        notes:
        """

        del(self.source_code_mods)
        self.source_code_mods={}

    def make(self,directory,clean=False):
        """
        applies source code modifications before executing make

        args:
            directory - path to directory holding code to make
            clean - toggle whether to execute make clean
        notes:
            set compile flags using Build.flags_add
        """

        if clean:
            subprocess.run(shlex.split('make clean'),shell=False,cwd=str(directory))    

        else: #not making clean, so parse flags - two ways depending whether they hold numerical value e.g. -DSTDOUT vs -DTOKAMAK=1

            #perform source code edits if any requested
            for filename,mods in self.source_code_mods.items():
                run_scripts.LOCUST_edit_var.LOCUST_edit_var(filepath_in=directory / filename,filepath_out=directory / filename,**mods) #perform source code edits

            if not hasattr(self,'environment'):
                print("WARNING: Build.make does not contain environment - using settings.environment_default!")                
                self.environment=run_scripts.environment.Environment(settings.system_default)

            flags_=' '.join(['-D{}={}'.format(flag,value) if value is not True else '-D{}'.format(flag) for flag,value in self.flags.items()])
            command=' '.join([self.environment.create_command_string(),'; make','FLAGS={}'.format(shlex.quote(flags_))])

            try:
                subprocess.call(command,shell=True,cwd=str(directory))
            except subprocess.CalledProcessError as err:
                raise(err)

#################################

##################################################################

###################################################################################################