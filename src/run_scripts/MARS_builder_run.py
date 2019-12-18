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
    import run_scripts.build
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/build.py could not be imported!\nreturning\n")
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

class MARS_builder_run(run_scripts.workflow.Workflow):
    """
    defines a workflow which performs a simple LOCUST simulation

    notes:
        contains a Build and an Environment which store relevant settings for this run
        both this object and it's stored objects each have separate environments - in case one wants to run with a different environment that they built with
        when editing mars_read.f90 source code using settings_mars_read, strings should be passed literally (see below) 
    usage:
        python MARS_builder_run.py --filepath_in 'some string formatted like this' --flags TOKAMAK=8 --settings_mars_read a_string="'should be formatted like this'" a_number="like_this"
        some_run=MARS_builder_run(filepath_input='here',system_name='TITAN',flags=flags) #by default this will run compile and run MARS_builder on filepath_input and dump to LOCUST_IO/data/input_files (since output of this function is still a LOCUST input)
        some_run.run() #this will execute all stages of a mars_build run
    """ 

    def __init__(self,filepath_input,dir_output=support.dir_input_files,dir_MARS_builder=support.dir_run_scripts / 'mars_builder' / 'mars_builder_temp',system_name='TITAN',settings_mars_read={},flags={}):
        """
        notes:
            most information stored in MARS_builder_run.environment and MARS_builder_run.build, most init args are to init these instances
        args:
            filepath_input - default path to input file 
            dir_output - directory to write results to (defaults to dir_input since output of this is a LOCUST input)
            dir_MARS_builder - directory to temporarily store source code
            system_name - string identifier to choose from selection of environments stored as class attributes 
            settings_mars_read - dict denoting variable names and values to set them to within mars_read.f90
            flags - dict denoting compile flags (no '-D' please e.g. STDOUT TOKAMAK=3)
        """
        
        #execute base class constructor to inherit required structures
        super().__init__()

        ################################# first generate class data that will be needed in workflow

        self.dir_MARS_builder=pathlib.Path(dir_MARS_builder)
        self.filepath_input=pathlib.Path(filepath_input)
        self.dir_output=pathlib.Path(dir_output)

        self.environment=run_scripts.environment.Environment(system_name=system_name) #create an environment for running
        
        if not settings_mars_read: settings_mars_read={} #turn __init__ args such as filepath_input into the corresponding source code edits
        settings_mars_read['file']="'{}'".format(str(self.filepath_input))
        dir_output=str(self.dir_output) if (str(self.dir_output)[-1]==str(os.sep) or str(self.dir_output)[-2:]==str(os.sep)) else ''.join([str(self.dir_output),str(os.sep)]) #add final separator character if stripped by pathlib
        settings_mars_read['root']="'{}'".format(str(dir_output))

        if not flags: #set default compile flags 
            flags={}
            flags['TOKAMAK']=1

        self.build=run_scripts.build.Build(system_name=system_name)
        self.build.flags_add(**flags)         
        self.build.source_code_mods_add(source_code_filename='mars_read.f90',**settings_mars_read)

        ################################# now make commands (defined below in this class) available to this workflow (and state position in execution order)

        self.add_command(command_name='mkdir',command_function=self.setup_MARS_builder_dirs,position=1) #add new workflow stages
        self.add_command(command_name='get_code',command_function=self.get_code,position=2)
        self.add_command(command_name='make',command_function=self.build_code,position=3)
        self.add_command(command_name='run_code',command_function=self.run_code,position=4)
        self.add_command(command_name='cleanup',command_function=self.clean_up_code,position=5)

    def setup_MARS_builder_dirs(self,*args,**kwargs):
        """
        notes:
        """

        #create output and MARS_builder directories if do not exist
        if not self.dir_output.is_dir(): self.dir_output.mkdir()
        if not self.dir_MARS_builder.is_dir(): self.dir_MARS_builder.mkdir()

    def get_code(self,*args,**kwargs):
        """
        notes:
        """

        #copy over source code
        for file in (support.dir_run_scripts / 'mars_builder').glob('*'):
            subprocess.run(shlex.split('cp {file} {dir_MARS_builder}'.format(file=str(file),dir_MARS_builder=self.dir_MARS_builder)),shell=False)


    def build_code(self,*args,**kwargs):
        """
        notes:
        """

        #execute make
        self.build.make(directory=self.dir_MARS_builder,clean=True)
        self.build.make(directory=self.dir_MARS_builder,clean=False)

    def run_code(self,*args,**kwargs):
        """
        notes:
        args:
        """

        #run MARS_builder
        command=' '.join([self.environment.create_command_string(),'; ./mars_build'])
        try:
            with subprocess.Popen(command,shell=True,cwd=str(self.dir_MARS_builder)) as proc: #stdin=PIPE, stdout=PIPE, stderr=STDOUT
                pass
        except subprocess.CalledProcessError as err:
            print("ERROR: MARS_builder_run failed to run MARS_builder!\nreturning\n")

    def clean_up_code(self,*args,**kwargs):
        """
        notes:
        """

        for directory in [ 
                         (self.dir_MARS_builder)]:
            subprocess.run(shlex.split('rm -r {}'.format(str(directory))),shell=False) #delete temporary directories


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

    #provide some extra parsing steps to flags, prec_mod settings and any similar dict-like input arguments
    args.settings_mars_read=run_scripts.utils.command_line_arg_parse_dict(args.settings_mars_read)
    args.flags=run_scripts.utils.command_line_arg_parse_dict(args.flags)

    this_run=MARS_builder_run(filepath_input=args.filepath_input,dir_output=args.dir_output,system_name=args.system_name,dir_MARS_builder=args.dir_MARS_builder,settings_mars_read=args.settings_mars_read,flags=args.flags)
    this_run.run()

#################################

##################################################################

###################################################################################################