#NEMO_run.py

'''
Samuel Ward
28/09/2019
----
tools to run the NEMO code in python actor form
---
notes: 
    contains the tools to perform a single run of NEMO
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
    import os

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
        python NEMO_run.py --run_in 1 --shot_in 1 --run_out 2 --shot_out 2 --nmarker 1000 
        some_run=NEMO_run(run_in=1,shot_in=1,run_out=2,shot_out=2,nmarker=1000) 
        some_run.run() #this will execute all stages of a NEMO run
    """ 

    workflow_name='NEMO_run'

    def __init__(self,
                dir_NEMO,
                shot_in,
                shot_out,
                run_in,
                run_out,
                username=settings.username,
                imasdb=settings.imasdb,
                imas_version=settings.imas_version,
                nmarker=int(1.e6),
                fokker_flag=0,
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
            nmarker - number of Monte Carlo markers to generate
            fokker_flag - 
        """

        #execute base class constructor to inherit required structures
        super().__init__()

        ################################# first generate class data that will be needed in workflow

        self.dir_NEMO=dir_NEMO
        self.shot_in=shot_in
        self.shot_out=shot_out
        self.run_in=run_in
        self.run_out=run_out
        self.username=username
        self.imasdb=imasdb
        self.imas_version=imas_version
        self.nmarker=nmarker
        self.fokker_flag=fokker_flag

        ################################# now make commands (defined below in this class) available to this workflow (and state position in execution order)
        
        self.add_command(command_name='run_code',command_function=self.call_NEMO_actor,position=1) #add new workflow 

    def call_NEMO_actor(self,*args,**kwargs):
        """
        NEMO_run stage for executing NEMO python actor

        notes:
            assumes you have already cloned and compiled the nemo source and actor into dir_NEMO 
            adapted from run_NEMO.py within NEMO source code
        """

        try:
            import imas 
        except:
            raise ImportError("ERROR: NEMO_run.call_NEMO_actor could not import IMAS module!\nreturning\n")
            return 

        # IMPORT MODULE(S) FOR SPECIFIC PHYSICS CODE(S)
        actor_path = os.path.join(os.getenv('KEPLER'), 'imas/src/org/iter/imas/python')
        list_of_actors = ['nemo']
        for name in list_of_actors:
            sys.path.insert(1,str(os.path.join(actor_path,name)))
            globals()[name] = getattr(__import__(name), name)

        # OPEN INPUT DATAFILE TO GET DATA FROM IMAS SCENARIO DATABASE
        print('=> Open input datafile')
        input = imas.ids(self.shot_in,self.run_in,0,0)
        input.open_env(self.username,self.imasdb,self.imas_version)

        # READ INPUT IDS'S AND CLOSE INPUT FILE
        print('=> Read input IDSs')
        print('   ---> equilibrium')
        input.equilibrium.get()
        print('   ---> core_profiles')
        input.core_profiles.get()
        print('   ---> nbi')
        input.nbi.get()
        print('   ---> distribution_sources')
        input.distribution_sources.get()
        print('   ---> distributions')
        input.distributions.get()
        input.close()

        # IF DISTRIBUTIONS HAS NEVER BEEN FILLED
        # (IF THIS IS THE FIRST TIME STEP OF THE SIMULATION)
        # THEN HOMOGENEOUS_TIME HAS NEVER BEEN FILLED
        # => SYSTEMATICALLY FILLED HERE THEN
        input.distributions.ids_properties.homogeneous_time = 1

        # OPEN OUTPUT OBJECT, IN VIEW OF SAVING RESULTS TO LOCAL DB
        print('=> Create output datafile')
        output = imas.ids(self.shot_out,self.run_out,0,0)

        # EXECUTE PHYSICS CODE
        print('=> Execute NEMO')
        output.distribution_sources = nemo(input.equilibrium,input.core_profiles,input.nbi,
                                           input.distribution_sources,self.fokker_flag,self.nmarker,pathlib.Path(self.dir_NEMO) / 'input' / 'input_nemo_reg_imas.xml')

        # CREATE OUTPUT DATAFILE AND EXPORT RESULTS TO LOCAL DATABASE
        print('=> Export output IDSs to local database')
        output.create_env(self.user,self.tokamakname,self.version)
        output.distribution_sources.put()
        output.close()
        print('Done exporting.')

##################################################################

if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser(description='run NEMO from command line in Python!')

    parser.add_argument('--dir_NEMO',type=str,action='store',dest=dir_NEMO,help="directory where NEMO source is stored",required=True)
    parser.add_argument('--shot_in',type=int,action='store',dest=shot_in,help="input IDS shot number",required=True)
    parser.add_argument('--shot_out',type=int,action='store',dest=shot_out,help="output IDS shot number",required=True)
    parser.add_argument('--run_in',type=int,action='store',dest=run_in,help="input IDS run number",required=True)
    parser.add_argument('--run_out',type=int,action='store',dest=run_out,help="output IDS run number",required=True)
    parser.add_argument('--username',type=str,action='store',default=settings.username,dest=username,help="IMAS username",required=False)
    parser.add_argument('--imasdb',type=str,action='store',default=settings.imasdb,dest=imasdb,help="local IMAS database name set by imasdb command, sometimes called 'tokamak name'",required=False)
    parser.add_argument('--imas_version',type=str,action='store',default=settings.imas_version,dest=imas_version,help="string denoting IMAS major version e.g. '3'",required=False)
    parser.add_argument('--nmarker',type=int,action='store',default=int(1.e6),dest=nmarker,help="number of Monte Carlo markers to generate",required=False)
    parser.add_argument('--fokker_flag',type=int,action='store',default=0,dest=fokker_flag,help="",required=False)

    args=parser.parse_args()

    this_run=NEMO_run(dir_NEMO=args.dir_NEMO,
                    shot_in=args.shot_in,
                    shot_out=args.shot_out,
                    run_in=args.run_in,
                    run_out=args.run_out,
                    username=args.username,
                    imasdb=args.imasdb,
                    imas_version=args.imas_version,
                    nmarker=args.nmarker,
                    fokker_flag=args.fokker_flag):
    this_run.run()

#################################

##################################################################

###################################################################################################