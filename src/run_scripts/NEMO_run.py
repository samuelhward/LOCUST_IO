#NEMO_run.py

'''
Samuel Ward
01/01/2019
----
tools to run the NEMO code in python actor form
---
notes: 
    contains the tools to perform a single run of NEMO
    custom steps can be added via add_command

    to get nemo all ready:
    git clone ssh://git@git.iter.org/heat/nemo.git
    cd nemo
    git checkout develop
    (source correct environment)
    make 
    kepler_load <some_kepler_installation_name> - in my case this is new_lowlevel
    make actor
    DONE!
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
    import numpy as np
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
    import run_scripts.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/utils.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import run_scripts.environment
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/environment.py could not be imported!\nreturning\n")
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

class NEMO_run(run_scripts.workflow.Workflow):
    """
    defines a workflow which performs a simple NEMO simulation

    notes:
        adds workflow components to workflow class for execution of a single standard NEMO run
        since calling from within Python, compatible environment must already be loaded before calling this workflow (call from command line using argparse if necessary) 
    usage:
        python NEMO_run.py --run_in 1 --shot_in 1 --run_out 2 --shot_out 2 --xml_settings nmarker 1000 
        some_run=NEMO_run(run_in=1,shot_in=1,run_out=2,shot_out=2,xml_settings=some_dict) 
        some_run.run() #this will execute all stages of a NEMO run
    """ 

    def __init__(self,
                dir_NEMO,
                shot_in,
                shot_out,
                run_in,
                run_out,
                username=settings.username,
                imasdb=settings.imasdb,
                imas_version=settings.imas_version,
                xml_settings=None,
                commands=['write_xml','run_code'],
                *args,**kwargs):
        """
        notes:
        args:
            dir_NEMO - directory where NEMO source is stored
            shot_in - input IDS shot number
            shot_out - output IDS shot number
            run_in - input IDS run number
            run_out - output IDS run number
            username - IMAS username
            imasdb - local IMAS database name set by 'imasdb' command, sometimes called 'tokamak name'
            imas_version - string denoting IMAS major version e.g. '3'
            xml_settings - dict containing settings to hard code in input XML
            commands - optional list of strings specifying order of subcommands to execute by workflow
        """

        #execute base class constructor to inherit required structures
        super().__init__()

        ################################# first generate class data that will be needed in workflow

        self.dir_NEMO=pathlib.Path(dir_NEMO)
        self.shot_in=shot_in
        self.shot_out=shot_out
        self.run_in=run_in
        self.run_out=run_out
        self.username=username
        self.imasdb=imasdb
        self.imas_version=imas_version
        self.xml_settings=xml_settings

        ################################# now make commands (defined below in this class) available to this workflow (and state position in execution order)

        self.commands_available={}
        self.commands_available['write_xml']=self.write_xml_settings
        self.commands_available['run_code']=self.call_NEMO_actor
        for command in commands:
            self.add_command(command_name=command,command_function=self.commands_available[command]) #add all workflow stages

    def write_xml_settings(self,*args,**kwargs):
        """
                

        notes:
        """

        if self.xml_settings:
            try:
                import xml.etree.ElementTree
            except:
                raise ImportError("ERROR: NEMO_run.write_xml_settings could not import xml module!\nreturning\n")
                return 

            nemo_settings_xml_filepath = self.dir_NEMO / 'input' / 'input_nemo_sa_imas.xml'
            tree=xml.etree.ElementTree.parse(nemo_settings_xml_filepath) #open/grab xml data
            root=tree.getroot()
            for key,value in self.xml_settings.items():
                next(root.iter(key)).text=str(value) #do a search for each setting and apply value
            tree.write(nemo_settings_xml_filepath,method='xml',xml_declaration=True,encoding='utf-8')

    def call_NEMO_actor(self,*args,**kwargs):
        """
        NEMO_run stage for executing NEMO python actor

        notes:
            assumes you have already cloned and compiled the nemo source and actor into dir_NEMO 
            assumes target IDS already created
        """

        try:
            import imas 
        except:
            raise ImportError("ERROR: NEMO_run.call_NEMO_actor could not import IMAS module!\nreturning\n")
            return 


        # IMPORT MODULE(S) FOR SPECIFIC PHYSICS CODE(S)
        actor_path=pathlib.Path('/') / 'home' / 'ITER' / f'{self.username}' / 'public' / 'imas_actors' / 'nemo'
        sys.path.insert(1,str(actor_path))
        from nemo.wrapper import nemo_actor as nemo
        time=0.0

        # OPEN INPUT DATAFILE TO GET DATA FROM IMAS SCENARIO DATABASE
        print('=> Open input datafile')
        input = imas.ids(self.shot_in,self.run_in,0,0)
        input.open_env(self.username,self.imasdb,self.imas_version)

        # READ INPUT IDS'S AND CLOSE INPUT FILE
        print('=> Read input IDSs')
        print('   ---> equilibrium')
        input.equilibrium.getSlice(time,1)
        print('   ---> core_profiles')
        input.core_profiles.getSlice(time,1)
        print('   ---> nbi')
        input.nbi.getSlice(time,1)
        print('   ---> wall')
        input.wall.getSlice(time,1)
        print('   ---> distribution_sources')
        input.distribution_sources.getSlice(time,1)
        print('   ---> distributions')
        input.distributions.getSlice(time,1)
        input.close()

        # IF DISTRIBUTIONS HAS NEVER BEEN FILLED
        # (IF THIS IS THE FIRST TIME STEP OF THE SIMULATION)
        # THEN HOMOGENEOUS_TIME HAS NEVER BEEN FILLED
        # => SYSTEMATICALLY FILLED HERE THEN
        input.distributions.ids_properties.homogeneous_time = 1

        # OPEN OUTPUT OBJECT, IN VIEW OF SAVING RESULTS TO LOCAL DB
        print('=> Create output datafile')
        output=imas.ids(self.shot_out,self.run_out,0)
        output.open_env(self.username,self.imasdb,self.imas_version)
        idx_out=output.distribution_sources.getPulseCtx()

        # EXECUTE PHYSICS CODE
        print('=> Execute NEMO')
        output.distribution_sources=nemo(input.equilibrium,
                                            input.core_profiles,
                                            input.nbi,
                                            input.wall,
                                            input.distribution_sources,
                                            pathlib.Path(self.dir_NEMO) / 'input' / 'input_nemo_sa_imas.xml')

        # CREATE OUTPUT DATAFILE AND EXPORT RESULTS TO LOCAL DATABASE
        print('=> Export output IDSs to local database')
        output.distribution_sources.setPulseCtx(idx_out)
        output.distribution_sources.ids_properties.homogeneous_time=1
        output.distribution_sources.time=np.array([0.0])
        output.distribution_sources.put()
        output.close()
        print('Done exporting.')

    def call_NEMO_actor_command_line(self):
        """
        helper function for calling NEMO python wrapper from command line - calls this file

        notes:
            this stage is NOT intended for use in workflows via add_command - otherwise recursion will occur
            motivation for this function is to be able to define runtime environment whilst calling NEMO python actor from command line
        """

        NEMO_run_args={}
        NEMO_run_args['dir_NEMO']=[self.dir_NEMO]
        NEMO_run_args['shot_in']=[self.shot_in]
        NEMO_run_args['shot_out']=[self.shot_out]
        NEMO_run_args['run_in']=[self.run_in]
        NEMO_run_args['run_out']=[self.run_out]
        NEMO_run_args['username']=[self.username]
        NEMO_run_args['imasdb']=[self.imasdb]
        NEMO_run_args['imas_version']=[self.imas_version]
        NEMO_run_args['xml_settings']=[self.xml_settings]

        NEMO_run_environment=run_scripts.environment.Environment('TITAN_NEMO')
        command=' '.join([NEMO_run_environment.create_command_string(),
                                '; python {}'.format(str(support.dir_run_scripts / 'NEMO_run.py')),
                                run_scripts.utils.command_line_arg_parse_generate_string(**NEMO_run_args)])

        try:
            subprocess.call(command,shell=True)# as proc: #stdin=PIPE, stdout=PIPE, stderr=STDOUT

        except subprocess.CalledProcessError as err:
            print("ERROR: {workflow_name}.call_NEMO_actor_command_line() failed to run NEMO!\nreturning\n".format(workflow_name=self.workflow_name))
            #raise(err)

##################################################################

if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser(description='run NEMO from command line in Python!')

    parser.add_argument('--dir_NEMO',type=str,action='store',default=support.dir_nemo,dest='dir_NEMO',help="directory where NEMO source is stored",required=False)
    parser.add_argument('--shot_in',type=int,action='store',dest='shot_in',help="input IDS shot number",required=True)
    parser.add_argument('--shot_out',type=int,action='store',dest='shot_out',help="output IDS shot number",required=True)
    parser.add_argument('--run_in',type=int,action='store',dest='run_in',help="input IDS run number",required=True)
    parser.add_argument('--run_out',type=int,action='store',dest='run_out',help="output IDS run number",required=True)
    parser.add_argument('--username',type=str,action='store',default=settings.username,dest='username',help="IMAS username",required=False)
    parser.add_argument('--imasdb',type=str,action='store',default=settings.imasdb,dest='imasdb',help="local IMAS database name set by imasdb command, sometimes called 'tokamak name'",required=False)
    parser.add_argument('--imas_version',type=str,action='store',default=settings.imas_version,dest='imas_version',help="string denoting IMAS major version e.g. '3'",required=False)
    parser.add_argument('--xml_settings',nargs='+',type=str,action='store',default={},dest='xml_settings',help="settings contained in input XML file e.g. nmarker=64",required=False)
    parser.add_argument('--commands',type=str,action='store',dest='commands',help="", required=False)

    args=parser.parse_args()
    args.xml_settings=run_scripts.utils.command_line_arg_parse_dict(args.xml_settings)
    args.commands=run_scripts.utils.literal_eval(args.commands)

    this_run=NEMO_run(**{key:arg for key,arg in args._get_kwargs()})
    this_run.run()

#################################

##################################################################

###################################################################################################