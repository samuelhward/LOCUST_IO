#LOCUST_environment.py

'''
Samuel Ward
27/09/2019
----
LOCUST environments and related tools for various machines running LOCUST
---
notes: 
---
'''


##################################################################
#Preamble

try:
    import sys
    import shlex
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

##################################################################
#Main

class LOCUST_environment:
    """
    class to hold a LOCUST runtime environment

    notes:
    usage:
        my_env=LOCUST_environment('TITAN')
        my_env.display(command_types=['export']) #display all necessary export commands
        print(my_env.command_types) #show the different commands used
        print(my_env.required_modules) #show the required modules
    """ 

    ##################################################################
    #common environments to be made available globally
    
    ################################# GA
    LOCUST_environment_GA={}
    LOCUST_environment_GA['export']={}
    LOCUST_environment_GA['export']['NO_AT_BRIDGE']=1
    LOCUST_environment_GA['export']['OMP_NUM_THREADS']=1
    LOCUST_environment_GA['export']['OMP_STACKSIZE']=102400
    LOCUST_environment_GA['export']['FFLAGS']="-I$HDF5_DIR/include"
    LOCUST_environment_GA['export']['LDFLAGS']="-L$HDF5_DIR/lib"
    LOCUST_environment_GA['module load']=[
                                'pgf/18.3',
                                'hdf5/1.8.19-mpich3.2-pgf18.3',
                                'cuda/9.0',
                                'gcc7/default'  ]
    LOCUST_environment_GA['module unload']=[
                                'gcc7/default'  ]
    LOCUST_environment_GA['module switch']=[]
    LOCUST_environment_GA['misc']={}
    LOCUST_environment_GA['misc']['ulimit']="-s 2000000"

    ################################# TITAN
    LOCUST_environment_TITAN={}
    LOCUST_environment_TITAN['export']={}
    LOCUST_environment_TITAN['export']['FFLAGS']='"-I$HDF5_DIR/include"'
    LOCUST_environment_TITAN['export']['LDFLAGS']='"-L$HDF5_DIR/lib"'
    LOCUST_environment_TITAN['export']['NO_AT_BRIDGE']=1
    LOCUST_environment_TITAN['export']['OMP_NUM_THREADS']=8
    LOCUST_environment_TITAN['export']['OMP_STACKSIZE']=102400
    LOCUST_environment_TITAN['export']['CUDA_CACHE_DISABLE']=1
    LOCUST_environment_TITAN['module load']=[
                                'IMAS/3.23.3-4.1.3',
                                'CUDA/10.1.105'  ]
    LOCUST_environment_TITAN['module unload']=[]
    LOCUST_environment_TITAN['module switch']=[
                                'PGI/19.4-GCC-6.4.0-2.28',
                                'HDF5/1.10.5-PGI-19.4-GCC-6.4.0-2.28'  ]
    LOCUST_environment_TITAN['misc']={}
    LOCUST_environment_TITAN['misc']['ulimit']="-s 2000000"


    ################################# NEMO TITAN
    NEMO_environment_TITAN={}
    NEMO_environment_TITAN['module load']=['IMAS/3.21.0-3.8.11',
                                           'Kepler/2.5p4-2.1.5',
                                           'FC2K/4.4.0',
                                           'PyUAL/1.0.0-foss-2018a-Python-3.6.4',
                                           'sh/1.12.14-foss-2018a-Python-3.6.4',
                                           'lxml/4.2.0-foss-2018a-Python-3.6.4',
                                           'FRUIT/3.4.3-intel-2018a-Ruby-2.5.1',
                                           'FRUIT_processor/3.4.3-intel-2018a-Ruby-2.5.1',
                                           'interpos/8.2.1-ifort',
                                           'XMLlib/3.1.0-intel-2018a',
                                           'PSPLINE/20181008-intel-2018a']
    NEMO_environment_TITAN['module unload']=['Anaconda3']
    NEMO_environment_TITAN['module switch']=['Python/3.6.4-foss-2018a',
                                             'matplotlib/2.1.2-foss-2018a-Python-3.6.4',
                                             'PyYAML/3.12-foss-2018a-Python-3.6.4',
                                             'UDA/2.2.5-foss-2018a',
                                             'IDStools/1.0.9-Python-3.6.4',
                                             'PostgreSQL/10.3-foss-2018a-Python-3.6.4',
                                             'SWIG/3.0.12-foss-2018a-Python-3.6.4',
                                             'HDF5/1.10.1-foss-2018a',
                                             'MDSplus/7.46.1-foss-2018a',
                                             'MDSplus-Python/7.46.1-foss-2018a-Python-3.6.4',
                                             'Boost/1.66.0-foss-2018a # ?',
                                             'Tkinter/3.6.4-foss-2018a-Python-3.6.4']
    NEMO_environment_TITAN['export']={}
    NEMO_environment_TITAN['export']['MD_ACCESS']='no'
    NEMO_environment_TITAN['export']['DIAG_INFO']='-DNO_DIAG_INFO'
    NEMO_environment_TITAN['misc']={}
    NEMO_environment_TITAN['misc']['ulimit']='-s unlimited'

    ################################# MARFE
    LOCUST_environment_MARFE={}
    LOCUST_environment_MARFE['export']={}
    LOCUST_environment_MARFE['export']['PATH']='$PATH:/usr/local/cuda-9.1/bin'
    LOCUST_environment_MARFE['export']['PATH']='$PATH:/usr/local/hdf5/bin'
    LOCUST_environment_MARFE['export']['HDF5_DIR']='/usr/local/hdf5'
    LOCUST_environment_MARFE['export']['LD_LIBRARY_PATH']='/usr/local/hdf5/lib'
    LOCUST_environment_MARFE['export']['FFLAGS']="-I$\"{HDF5_DIR}\"/include"
    LOCUST_environment_MARFE['export']['LDFLAGS']="-L$\"{HDF5_DIR}\"/lib"
    LOCUST_environment_MARFE['export']['OMP_NUM_THREADS']=1
    LOCUST_environment_MARFE['export']['OMP_STACKSIZE']=102400
    LOCUST_environment_MARFE['export']['NO_AT_BRIDGE']=1
    LOCUST_environment_MARFE['module load']=[
                                '/usr/local/pgi/modulefiles/pgi/18.4']
    LOCUST_environment_MARFE['module unload']=[]
    LOCUST_environment_MARFE['module switch']=[]
    LOCUST_environment_MARFE['misc']={}
    LOCUST_environment_MARFE['misc']['ulimit']="-s 2000000"

    ################################# GPU
    LOCUST_environment_GPU={}
    ################################# CUMULUS
    LOCUST_environment_CUMULUS={}

    system_names=['GA','TITAN','TITAN_NEMO','MARFE','GPU','CUMULUS']

    def __init__(self,system_name=None):
        """
        initialise an environment

        args:
            system_name - string identifier to choose from selection of environments stored as class attributes 
        notes:
            options:
                GA 
                TITAN
                TITAN_NEMO
                MARFE
                GPU
                CUMULUS
        """ 

        if system_name=='GA':
            self.environment=LOCUST_environment.LOCUST_environment_GA
        elif system_name=='TITAN':
            self.environment=LOCUST_environment.LOCUST_environment_TITAN
        elif system_name=='TITAN_NEMO':
            self.environment=LOCUST_environment.NEMO_environment_TITAN
        elif system_name=='MARFE':
            self.environment=LOCUST_environment.LOCUST_environment_MARFE
        elif system_name=='GPU':
            self.environment=LOCUST_environment.LOCUST_environment_GPU
        elif system_name=='CUMULUS':
            self.environment=LOCUST_environment.LOCUST_environment_CUMULUS
        else:
            print("ERROR: LOCUST_environment initialised with incorrect system_name - options are {}".format(LOCUST_environment.system_names))
            self.environment=None

        self.system_name=system_name
        self.command_types=[command_type for command_type in self.environment]
        self.required_modules=[module_name for module_name in self.environment['module load']]

    def create_command(self):
        """
        create single command string to load/source a chosen environment

        notes:
            also returns command string 
        """ 

        commands=[]

        for command in ['module load','module unload','module switch','misc','export']: #user can control order of execution of commands here
            things_to_command=self.environment[command]

            if command=='export':
                commands.extend([' '.join([command,''.join([thing_to_command,'=',str(value)])]) for thing_to_command,value in things_to_command.items()])
            elif command=='module load':
                commands.extend(['module load {}'.format(module) for module in things_to_command])
            elif command=='module switch':
                commands.extend(['module switch {}'.format(module) for module in things_to_command])
            elif command=='module unload':
                commands.extend(['module unload {}'.format(module) for module in things_to_command])
            elif command=='misc':
                commands.extend([' '.join([thing_to_command,str(self.environment[command][thing_to_command])]) for thing_to_command in things_to_command])

        commands=' ; '.join(command for command in commands)
        self.commands=commands
        return commands

    def display(self,string=False,command_types=None):
        """
        print chosen environment

        notes:
        args:
            string - print in single-line string form
            command_types - choose which commands to print
        """ 

        commands=[]
        if not command_types: command_types=self.command_types

        for command in command_types: #user can control order of execution of commands here
            things_to_command=self.environment[command]

            if command=='export':
                commands.extend([' '.join([command,thing_to_command+'=',str(value)]) for thing_to_command,value in things_to_command.items()])
            elif command=='module load':
                commands.extend(['module load {}'.format(module) for module in things_to_command])
            elif command=='module switch':
                commands.extend(['module switch {}'.format(module) for module in things_to_command])
            elif command=='module unload':
                commands.extend(['module unload {}'.format(module) for module in things_to_command])
            elif command=='misc':
                commands.extend([' '.join([thing_to_command+'=',str(self.environment[command][thing_to_command])]) for thing_to_command in things_to_command])

        if not string:
            for command in commands:
                print(command)
        else:
            commands=' ; '.join(command for command in commands)
            print(commands)

#################################

##################################################################

###################################################################################################