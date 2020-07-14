#BBNBI_run.py

'''
Samuel Ward
13/07/2020
----
tools to run the BBNBI code in python actor form
---
notes: 
    contains the tools to perform a single run of BBNBI
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

class BBNBI_run(run_scripts.workflow.Workflow):
    """
    defines a workflow which performs a simple BBNBI simulation

    notes:
        adds workflow components to workflow class for execution of a single standard BBNBI run
        since calling from within Python, compatible environment must already be loaded before calling this workflow (call from command line using argparse if necessary) 
    usage:
        python BBNBI_run.py --run_in 1 --shot_in 1 --run_out 2 --shot_out 2 --xml_settings nmarker 1000 
        some_run=BBNBI_run(run_in=1,shot_in=1,run_out=2,shot_out=2,xml_settings=some_dict) 
        some_run.run() #this will execute all stages of a BBNBI run
    """ 

    def __init__(self,
                dir_BBNBI,
                shot_in,
                shot_out,
                run_in,
                run_out,
                username=settings.username,
                imasdb=settings.imasdb,
                imas_version=settings.imas_version,
                xml_settings=None,
                number_particles=10000,
                number_processors=1,
                *args,**kwargs):
        """
        notes:
        args:
            dir_BBNBI - directory where BBNBI source is stored
            shot_in - input IDS shot number
            shot_out - output IDS shot number
            run_in - input IDS run number
            run_out - output IDS run number
            username - IMAS username
            imasdb - local IMAS database name set by 'imasdb' command, sometimes called 'tokamak name'
            imas_version - string denoting IMAS major version e.g. '3'
            xml_settings - dict containing settings to hard code in input XML
            number_particles - set number of markers to generate
            number_processors - set number of CPUs
        """

        #execute base class constructor to inherit required structures
        super().__init__()

        ################################# first generate class data that will be needed in workflow

        self.dir_BBNBI=pathlib.Path(dir_BBNBI)
        self.shot_in=shot_in
        self.shot_out=shot_out
        self.run_in=run_in
        self.run_out=run_out
        self.username=username
        self.imasdb=imasdb
        self.imas_version=imas_version
        self.xml_settings=xml_settings
        self.number_particles=number_particles
        self.number_processors=number_processors

        ################################# now make commands (defined below in this class) available to this workflow (and state position in execution order)
        
        self.add_command(command_name='write_xml',command_function=self.write_xml_settings,position=1) #add new workflow 
        self.add_command(command_name='run_code',command_function=self.call_BBNBI_actor,position=2) 

    def write_xml_settings(self,*args,**kwargs):
        """
                

        notes:
        """

        if self.xml_settings:
            try:
                import xml.etree.ElementTree
            except:
                raise ImportError("ERROR: BBNBI_run.write_xml_settings could not import xml module!\nreturning\n")
                return 

            BBNBI_settings_xml_filepath = self.dir_BBNBI / 'actors' / 'bbnbi_codeparam.xml'
            tree=xml.etree.ElementTree.parse(BBNBI_settings_xml_filepath) #open/grab xml data
            root=tree.getroot()
            for key,value in self.xml_settings.items():
                next(root.iter(key)).text=str(value) #do a search for each setting and apply value
            tree.write(BBNBI_settings_xml_filepath,method='xml',xml_declaration=True,encoding='utf-8')

    def call_BBNBI_actor(self,*args,**kwargs):
        """
        BBNBI_run stage for executing BBNBI python actor

        notes:
            assumes you have already cloned and compiled the BBNBI source and actor into dir_BBNBI (or at least store XML in equivalent structure)
            assumes target IDS already created
        """

        try:
            import imas
            from bbnbi.wrapper          import bbnbi_actor          as bbnbi
            from ascot4serial.wrapper   import ascot4serial_actor   as ascot4serial
            from ascot4parallel.wrapper import ascot4parallel_actor as ascot4parallel
        except:
            raise ImportError("ERROR: BBNBI_run.call_BBNBI_actor could not import IMAS modules!\nreturning\n")
            return 

        dt_required = 1.
        time_slice=-1. #simulate first time slice
        ntimes=1 #simulate one time step
        run_bbnbi      = True
        run_ascot      = False
        ascot_parallel = False

        def find_nearest(a, a0):
            "Element in nd array `a` closest to the scalar value `a0`"
            idx = np.abs(a - a0).argmin()
            return a.flat[idx],idx

        input = imas.ids(self.shot_in,self.run_in,0,0)
        input.open_env(self.username,self.imasdb,self.imas_version)
        time_array = input.equilibrium.partialGet('time')

        if time_slice == 0:
            ntimes = len(time_array)
            it = 0
        else:
          # FIND INDEX OF NEAREST TIME SLICE IN TIME ARRAY
            [tc,it] = find_nearest(time_array,time_slice)
        idx_in = input.distributions.getPulseCtx()

        # CREATE OUTPUT DATAFILE
        print('=> Create output datafile')
        output = imas.ids(self.shot_out,self.run_out,0,0)
        output.open_env(self.username,self.imasdb,self.imas_version)
        idx_out = output.distributions.getPulseCtx()

        # DEFINED A WORK IMAS STRUCTURE
        work = imas.ids(0,0)

        first_time = True

        for itime in range(it,it+ntimes):

            # READ INPUT IDSS FOR THIS SPECIFIC TIME SLICE
            print('=> Read input IDSs for time slice =',str(time_array[itime]),'s')
            print('   ---> equilibrium')
            input.distributions.setPulseCtx(idx_in)
            input.equilibrium.getSlice(time_array[itime],1)
            print('   ---> core_profiles')
            input.core_profiles.getSlice(time_array[itime],1)
            print('   ---> wall')
            input.wall.getSlice(time_array[itime],1)
            #input.wall.copyValues(wall)
            print('   ---> nbi')
            input.nbi.getSlice(time_array[itime],1)
            print('   ---> distribution_sources')
            input.distribution_sources.getSlice(time_array[itime],1)

            # IF THE DISTRIBUTIONS IDS EXISTS IN THE INPUT, READ IT, OTHERWISE LEAVE IT EMPTY 
            # (STILL WITH HOMOGENEOUS_TIME=1)
            print('   ---> distributions')
            if input.distributions.ids_properties.homogeneous_time>=0 and first_time:
                input.distributions.getSlice(time_array[itime],1)
            input.distributions.ids_properties.homogeneous_time = 1

            # COPY THE DISTRIBUTIONS IDS TO THE LAST_DISTRIBUTIONS OBJECT
            if first_time:
                work.distributions.copyValues(input.distributions)
                work.distributions.ids_properties.homogeneous_time=1
                work.distributions.time=input.core_profiles.time

            # RUN THE BBNBI ACTOR
            if run_bbnbi:
                print('=> Execute BBNBI')
                output.distribution_sources = bbnbi(self.number_particles,input.nbi,input.wall,input.core_profiles,input.equilibrium,self.dir_BBNBI / 'actors' / 'bbnbi_codeparam.xml')
            else:
                output.distribution_sources.copyValues(input.distribution_sources)

            # RUN THE ASCOT ACTOR
            if run_ascot:
                print('=> Execute ASCOT')
                if ascot_parallel:
                    output.distributions = ascot4parallel(input.core_profiles,input.equilibrium,input.wall,output.distribution_sources,work.distributions,self.dir_BBNBI / 'actors' / 'ascot4parallel_codeparam.xml','mpi_local',mpi_processes=self.number_processors)
                else:
                    output.distributions = ascot4serial(input.core_profiles,input.equilibrium,input.wall,output.distribution_sources,work.distributions,self.dir_BBNBI / 'actors' / 'ascot4serial_codeparam.xml')

            output.distributions.setPulseCtx(idx_out)
            output.distribution_sources.setPulseCtx(idx_out)
            output.distributions.time=np.array([time_array[itime]])
            output.distributions.ids_properties.homogeneous_time=1

            # WRITE THE OUTPUT TIME SLICE
            if first_time:
                output.distribution_sources.put() # WRITE STATIC DATA
                output.distributions.put()
            else:
                output.distribution_sources.putSlice()
                output.distributions.putSlice()

            # COPY THE OUTPUT DISTRIBUTIONS TO USE IS BACK TO THE INPUT
            work.distributions.copyValues(output.distributions)

            print('-------------------------------------')
            print('Output time = ',output.distributions.time[0])
            print('-------------------------------------')

            first_time = False

        input.close()
        output.close()
        print('Done exporting.')

    def call_BBNBI_actor_command_line(self):
        """
        helper function for calling BBNBI python wrapper from command line - calls this file

        notes:
            this stage is NOT intended for use in workflows via add_command - otherwise recursion will occur
            motivation for this function is to be able to define runtime environment whilst calling BBNBI python actor from command line
        """

        BBNBI_run_args={}
        BBNBI_run_args['dir_BBNBI']=[self.dir_BBNBI]
        BBNBI_run_args['shot_in']=[self.shot_in]
        BBNBI_run_args['shot_out']=[self.shot_out]
        BBNBI_run_args['run_in']=[self.run_in]
        BBNBI_run_args['run_out']=[self.run_out]
        BBNBI_run_args['username']=[self.username]
        BBNBI_run_args['imasdb']=[self.imasdb]
        BBNBI_run_args['imas_version']=[self.imas_version]
        BBNBI_run_args['xml_settings']=[self.xml_settings]

        BBNBI_run_environment=run_scripts.environment.Environment('TITAN_BBNBI')
        command=' '.join([BBNBI_run_environment.create_command_string(),
                                '; python {}'.format(str(support.dir_run_scripts / 'BBNBI_run.py')),
                                run_scripts.utils.command_line_arg_parse_generate_string(**BBNBI_run_args)])

        try:
            pass
            subprocess.call(command,shell=True,cwd=str(self.dir_BBNBI))# as proc: #stdin=PIPE, stdout=PIPE, stderr=STDOUT

        except subprocess.CalledProcessError as err:
            print("ERROR: {workflow_name}.call_BBNBI_actor_command_line() failed to run BBNBI!\nreturning\n".format(workflow_name=self.workflow_name))
            #raise(err)

##################################################################

if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser(description='run BBNBI from command line in Python!')

    parser.add_argument('--dir_BBNBI',type=str,action='store',default=support.dir_bbnbi,dest='dir_BBNBI',help="directory where BBNBI source is stored",required=False)
    parser.add_argument('--shot_in',type=int,action='store',dest='shot_in',help="input IDS shot number",required=True)
    parser.add_argument('--shot_out',type=int,action='store',dest='shot_out',help="output IDS shot number",required=True)
    parser.add_argument('--run_in',type=int,action='store',dest='run_in',help="input IDS run number",required=True)
    parser.add_argument('--run_out',type=int,action='store',dest='run_out',help="output IDS run number",required=True)
    parser.add_argument('--username',type=str,action='store',default=settings.username,dest='username',help="IMAS username",required=False)
    parser.add_argument('--imasdb',type=str,action='store',default=settings.imasdb,dest='imasdb',help="local IMAS database name set by imasdb command, sometimes called 'tokamak name'",required=False)
    parser.add_argument('--imas_version',type=str,action='store',default=settings.imas_version,dest='imas_version',help="string denoting IMAS major version e.g. '3'",required=False)
    parser.add_argument('--xml_settings',nargs='+',type=str,action='store',default={},dest='xml_settings',help="settings contained in input XML file e.g. nmarker=64",required=False)
    parser.add_argument('--number_particles',type=int,action='store',default=10000,dest='number_particles',help="set number of markers to generate",required=False)
    parser.add_argument('--number_processors',type=int,action='store',default=1,dest='number_processors',help="set number of CPUs",required=False)

    args=parser.parse_args()
    args.xml_settings=run_scripts.utils.command_line_arg_parse_dict(args.xml_settings)

    this_run=BBNBI_run(dir_BBNBI=args.dir_BBNBI,
                    shot_in=args.shot_in,
                    shot_out=args.shot_out,
                    run_in=args.run_in,
                    run_out=args.run_out,
                    username=args.username,
                    imasdb=args.imasdb,
                    imas_version=args.imas_version,
                    xml_settings=args.xml_settings,
                    number_particles=args.number_particles,
                    number_processors=args.number_processors)
    this_run.run()

#################################

##################################################################

###################################################################################################