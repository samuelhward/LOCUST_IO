#MARS_builder_run.py

'''
Samuel Ward
30/09/2019
----
tools to run the MARS_builder code
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
    import os
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
    import run_scripts.MARS_builder_build
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/MARS_builder_build.py could not be imported!\nreturning\n")
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

class MARS_builder_run:
    """
    class to run MARS_builder

    notes:
        contains a MARS_builder_build and a LOCUST_environment which store relevant settings for this run
        both this object and it's stored objects each have separate environments - in case one wants to run with a different environment that they built with
        when editing mars_read.f90 source code using settings_mars_read, strings should be passed literally (see below) 
    usage:
        python MARS_builder_run.py --filepath_in 'some string formatted like this' --flags TOKAMAK=8 --settings_mars_read a_string \'"should be formatted like this"\' a_variable 'like_this'
        some_run=MARS_builder_run(filepath_input='here',system_name='TITAN',flags=flags) #by default this will run compile and run MARS_builder on filepath_input and dump to LOCUST_IO/data/input_files (since output of this function is still a LOCUST input)
        some_run.run() #this will do all stages of a mars_build run
    """ 

    def __init__(self,filepath_input,dir_output=support.dir_input_files,dir_MARS_builder=support.dir_run_scripts / 'mars_builder',system_name='TITAN',settings_mars_read={},flags={}):
        """
        notes:
            most information stored in MARS_builder_run.environment and MARS_builder_run.build, most init args are to init these instances
        args:
            filepath_input - default path to input file (can be overridden in MARS_builder_run.build.settings_mars_read)
            dir_output - directory to write results to (defaults to dir_input since output of this is a LOCUST input)
            dir_MARS_builder - directory to temporarily store source code
            system_name - string identifier to choose from selection of environments stored as class attributes 
            settings_mars_read - dict denoting variable names and values to set them to within mars_read.f90
            flags - dict denoting compile flags (no '-D' please e.g. STDOUT TOKAMAK=3)
        """
        
        #convert strings to paths just in case supplied at command line
        filepath_input=pathlib.Path(filepath_input)
        dir_output=pathlib.Path(dir_output)

        self.dir_MARS_builder=dir_MARS_builder
        self.filepath_input=filepath_input
        self.dir_output=dir_output
        self.environment=run_scripts.LOCUST_environment.LOCUST_environment(system_name=system_name) #create an environment for running
        self.build=run_scripts.MARS_builder_build.MARS_builder_build(system_name=system_name)
        self.build.flags_add(**flags) 
        self.build.settings_mars_read_add(**settings_mars_read) 

    def run(self):
        """
        notes:
        args:
        """

        #create output directory if does not exist
        if not self.dir_output.is_dir(): self.dir_output.mkdir()

        #perform default source code edits - set output directory
        dir_output=str(self.dir_output) if (str(self.dir_output)[-1]==str(os.sep) or str(self.dir_output)[-2:]==str(os.sep)) else ''.join([str(self.dir_output),str(os.sep)]) #add final separator character if stripped by pathlib
        default_edits={'file':"'{}'".format(self.filepath_input),'root':"'{}'".format(dir_output)}
        run_scripts.LOCUST_edit_var.LOCUST_edit_var(filename_in=self.dir_MARS_builder / 'mars_read.f90',filename_out=self.dir_MARS_builder / 'mars_read.f90',**default_edits)

        #compile (source files are edited within this step)
        self.build.make(dir_MARS_builder=self.dir_MARS_builder,clean=True)
        self.build.make(dir_MARS_builder=self.dir_MARS_builder,clean=False)

        #run MARS_builder
        command=' '.join([self.environment.create_command(),'; ./mars_build'])
        try:
            pass 
            with subprocess.Popen(command,shell=True,cwd=str(self.dir_MARS_builder)) as proc: #stdin=PIPE, stdout=PIPE, stderr=STDOUT
                pass
        except subprocess.CalledProcessError as err:
            print("ERROR: MARS_builder_run failed to run MARS_builder!\nreturning\n")
            #raise(err)

if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser(description='run MARS_builder from command line in Python!')

    parser.add_argument('--filepath_in',type=str,action='store',default=support.dir_input_files,dest='filepath_input',help="filepath to input data",required=True)
    parser.add_argument('--dir_out',type=str,action='store',default=support.dir_input_files,dest='dir_output',help="directory to write results to",required=False)
    parser.add_argument('--sys_name',type=str,action='store',default='TITAN',dest='system_name',help="computer system currently running on e.g. \'TITAN\'",required=False)
    parser.add_argument('--dir_MARS_builder',type=str,action='store',default=support.dir_run_scripts / 'mars_builder',dest='dir_MARS_builder',help="path to mars_builder directory source code",required=False)
    parser.add_argument('--settings_mars_read',nargs='+',type=str,action='store',default={},dest='settings_mars_read',help="variable names and values to set them to within mars_read.f90",required=False)
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

    parsed_settings_mars_read=[setting_prec_mod.split('=') for setting_prec_mod in args.settings_mars_read]    
    settings_mars_read={}
    for setting in parsed_settings_mars_read:
        if len(setting)>1:
            settings_mars_read[setting[0]]=setting[1]
        else:
            settings_mars_read[setting[0]]=None

    this_run=MARS_builder_run(filepath_input=args.filepath_input,dir_output=args.dir_output,system_name=args.system_name,dir_MARS_builder=args.dir_MARS_builder,settings_mars_read=settings_mars_read,flags=flags)
    this_run.run()

#################################

##################################################################

###################################################################################################