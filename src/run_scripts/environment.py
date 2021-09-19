#environment.py

'''
Samuel Ward
27/09/2019
----
LOCUST environments and related tools for various machines running LOCUST and related tools
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

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/src/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import constants
except:
    raise ImportError("ERROR: LOCUST_IO/src/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import settings
except:
    raise ImportError("ERROR: LOCUST_IO/src/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)

##################################################################
#Main

class Environment:
    """
    class to hold a LOCUST runtime environment

    notes:
    usage:
        my_env=Environment('TITAN')
        my_env.display(command_types=['export']) #display all necessary export commands
        print(my_env.command_types) #show the different commands used
        print(my_env.required_modules) #show the required modules
    """ 

    ##################################################################
    #common environments to be made globally available by default 
    
    environments={}

    ################################# GA
    environments['GA']={}
    environments['GA']['export']={}
    environments['GA']['export']['NO_AT_BRIDGE']=1
    environments['GA']['export']['OMP_NUM_THREADS']=1
    environments['GA']['export']['OMP_STACKSIZE']=102400
    environments['GA']['export']['FFLAGS']="-I$HDF5_DIR/include"
    environments['GA']['export']['LDFLAGS']="-L$HDF5_DIR/lib"
    environments['GA']['module load']=[
                                'pgf/18.3',
                                'hdf5/1.8.19-mpich3.2-pgf18.3',
                                'cuda/9.0',
                                'gcc7/default' ]
    environments['GA']['module unload']=[
                                'gcc7/default' ]
    environments['GA']['module switch']=[]
    environments['GA']['misc']={}
    environments['GA']['misc']['ulimit']="-s 2000000"
    ################################# TITAN
    environments['TITAN']={}
    environments['TITAN']['export']={}
    environments['TITAN']['export']['FFLAGS']='"-I$HDF5_DIR/include"'
    environments['TITAN']['export']['LDFLAGS']='"-L$HDF5_DIR/lib"'
    environments['TITAN']['export']['NO_AT_BRIDGE']=1
    environments['TITAN']['export']['OMP_NUM_THREADS']=8
    environments['TITAN']['export']['OMP_STACKSIZE']=102400
    environments['TITAN']['export']['CUDA_CACHE_DISABLE']=1
    environments['TITAN']['module load']=[
                                'IMAS/3.32.1-4.9.1',
                                'CUDA/10.1.105']
    environments['TITAN']['module unload']=[
                                'matplotlib',
                                'GCCcore/10.2.0',
                                            ]
    environments['TITAN']['module switch']=[
                                'PGI/19.4-GCC-6.4.0-2.28',
                                'HDF5/1.10.5-PGI-19.4-GCC-6.4.0-2.28',
                                'Python/3.8.6-GCCcore-10.2.0'
                                ]
    environments['TITAN']['misc']={}
    environments['TITAN']['misc']['ulimit']="-s 2000000"
    ################################# NEMO TITAN
    environments['TITAN_NEMO']={}
    environments['TITAN_NEMO']['module load']=['IMAS/3.28.1-4.7.2',
                                                'FC2K',
                                                'FRUIT/3.4.3-intel-2018a-Ruby-2.5.1',
                                                'FRUIT_processor/3.4.3-intel-2018a-Ruby-2.5.1',
                                                'INTERPOS/9.0.2-intel-2018a',
                                                'XMLlib/3.2.0-intel-2018a',
                                                'PyAL/1.3.1-intel-2018a-Python-3.6.4']
    environments['TITAN_NEMO']['module unload']=[
                                                'matplotlib']
    environments['TITAN_NEMO']['module switch']=[]
    environments['TITAN_NEMO']['export']={}
    environments['TITAN_NEMO']['misc']={}
    environments['TITAN_NEMO']['misc']['ulimit']='-Ss unlimited'
    environments['TITAN_NEMO']['misc']['imasdb']=f'{settings.imasdb}'
    ################################# BBNBI TITAN
    environments['TITAN_BBNBI']={}
    environments['TITAN_BBNBI']['module load']=['ASCOT']
    environments['TITAN_BBNBI']['module unload']=[
                                                'matplotlib']
    environments['TITAN_BBNBI']['module switch']=[]
    environments['TITAN_BBNBI']['export']={}
    environments['TITAN_BBNBI']['misc']={}
    environments['TITAN_BBNBI']['misc']['ulimit']="-s 2000000"
    ################################# MARFE
    environments['MARFE']={}
    environments['MARFE']['export']={}
    environments['MARFE']['export']['PATH']='$PATH:/usr/local/cuda-9.1/bin'
    environments['MARFE']['export']['PATH']='$PATH:/usr/local/hdf5/bin'
    environments['MARFE']['export']['HDF5_DIR']='/usr/local/hdf5'
    environments['MARFE']['export']['LD_LIBRARY_PATH']='/usr/local/hdf5/lib'
    environments['MARFE']['export']['FFLAGS']="-I$\"{HDF5_DIR}\"/include"
    environments['MARFE']['export']['LDFLAGS']="-L$\"{HDF5_DIR}\"/lib"
    environments['MARFE']['export']['OMP_NUM_THREADS']=1
    environments['MARFE']['export']['OMP_STACKSIZE']=102400
    environments['MARFE']['export']['NO_AT_BRIDGE']=1
    environments['MARFE']['module load']=[
                                '/usr/local/pgi/modulefiles/pgi/18.4']
    environments['MARFE']['module unload']=[]
    environments['MARFE']['module switch']=[]
    environments['MARFE']['misc']={}
    environments['MARFE']['misc']['ulimit']="-s 2000000"
    ################################# GPU5
    environments['GPU5']={}
    environments['GPU5']['export']={}
    environments['GPU5']['export']['OMP_NUM_THREADS']=1
    environments['GPU5']['export']['OMP_STACKSIZE']=102400
    environments['GPU5']['export']['NO_AT_BRIDGE']=1
    environments['GPU5']['export']['CUDA_CACHE_DISABLE']=1
    environments['GPU5']['module load']=[
                            'python/3.7',
                            'cuda/9.1',
                            'pgi/20.1',
                            'hdf5/1.8.20',
                            'hdf5-devel/1.8.20']
    environments['GPU5']['module unload']=[]
    environments['GPU5']['module switch']=[]
    environments['GPU5']['misc']={}
    environments['GPU5']['misc']['ulimit']="-s 2000000"
    ################################# GPU6
    environments['GPU6']={}
    environments['GPU6']['export']={}
    environments['GPU6']['export']['OMP_NUM_THREADS']=8
    environments['GPU6']['export']['OMP_STACKSIZE']=102400
    environments['GPU6']['export']['NO_AT_BRIDGE']=1
    environments['GPU6']['export']['CUDA_CACHE_DISABLE']=1
    environments['GPU6']['module load']=[
                            'python/3.7',
                            'cuda/9.1',
                            'pgi/18.1',
                            'hdf5/1.8.20',
                            'hdf5-devel/1.8.20']
    environments['GPU6']['module unload']=[]
    environments['GPU6']['module switch']=[]
    environments['GPU6']['misc']={}
    environments['GPU6']['misc']['ulimit']="-s 2000000"
    ################################# GPU7
    environments['GPU7']={}
    environments['GPU7']['export']={}
    environments['GPU7']['export']['OMP_NUM_THREADS']=8
    environments['GPU7']['export']['OMP_STACKSIZE']=102400
    environments['GPU7']['export']['NO_AT_BRIDGE']=1
    environments['GPU7']['export']['CUDA_CACHE_DISABLE']=1
    environments['GPU7']['module load']=[
                            'python/3.7',
                            'cuda/9.1',
                            'pgi/18.1',
                            'hdf5/1.8.20',
                            'hdf5-devel/1.8.20']
    environments['GPU7']['module unload']=[]
    environments['GPU7']['module switch']=[]
    environments['GPU7']['misc']={}
    environments['GPU7']['misc']['ulimit']="-s 2000000"
    ################################# GPU8
    environments['GPU8']={}
    environments['GPU8']['export']={}
    environments['GPU8']['export']['OMP_NUM_THREADS']=16
    environments['GPU8']['export']['OMP_STACKSIZE']=102400
    environments['GPU8']['export']['NO_AT_BRIDGE']=1
    environments['GPU8']['export']['CUDA_CACHE_DISABLE']=1
    environments['GPU8']['module load']=[
                            'python/3.7',
                            'cuda/9.1',
                            'pgi/18.1',
                            'hdf5/1.8.20',
                            'hdf5-devel/1.8.20']
    environments['GPU8']['module unload']=[]
    environments['GPU8']['module switch']=[]
    environments['GPU8']['misc']={}
    environments['GPU8']['misc']['ulimit']="-s 2000000"
    ################################# GPU9
    environments['GPU9']={}
    environments['GPU9']['export']={}
    environments['GPU9']['export']['OMP_NUM_THREADS']=16
    environments['GPU9']['export']['OMP_STACKSIZE']=102400
    environments['GPU9']['export']['NO_AT_BRIDGE']=1
    environments['GPU9']['export']['CUDA_CACHE_DISABLE']=1
    environments['GPU9']['module load']=[
                            'python/3.7',
                            'cuda/9.1',
                            'pgi/18.1',
                            'hdf5/1.8.20',
                            'hdf5-devel/1.8.20']
    environments['GPU9']['module unload']=[]
    environments['GPU9']['module switch']=[]
    environments['GPU9']['misc']={}
    environments['GPU9']['misc']['ulimit']="-s 2000000"
    ################################# GPU10
    environments['GPU10']={}
    environments['GPU10']['export']={}
    environments['GPU10']['export']['OMP_NUM_THREADS']=16
    environments['GPU10']['export']['OMP_STACKSIZE']=102400
    environments['GPU10']['export']['NO_AT_BRIDGE']=1
    environments['GPU10']['export']['CUDA_CACHE_DISABLE']=1
    environments['GPU10']['module load']=[
                            'python/3.7',
                            'cuda/9.1',
                            'pgi/18.1',
                            'hdf5/1.8.20',
                            'hdf5-devel/1.8.20']
    environments['GPU10']['module unload']=[]
    environments['GPU10']['module switch']=[]
    environments['GPU10']['misc']={}
    environments['GPU10']['misc']['ulimit']="-s 2000000"
    ################################# CUMULUS
    environments['CUMULUS']={}
    environments['CUMULUS']['export']={}
    environments['CUMULUS']['export']['OMP_NUM_THREADS']=8
    environments['CUMULUS']['export']['OMP_STACKSIZE']=102400
    environments['CUMULUS']['export']['CUDA_CACHE_DISABLE']=1
    environments['CUMULUS']['export']['NO_AT_BRIDGE']=1
    environments['CUMULUS']['module load']=[
                            'cuda/10.1',
                            'pgi/20.1',
                            'hdf5/1.12.0',
                            'hdf5-devel/1.12.0',
                            'python/3.6.4'
                            ]
    environments['CUMULUS']['module unload']=[]
    environments['CUMULUS']['module switch']=[]
    environments['CUMULUS']['misc']={}
    environments['CUMULUS']['misc']['ulimit']="-s 2000000"
    ################################# VIKING
    environments['VIKING']={}
    environments['VIKING']['export']={}
    environments['VIKING']['export']['OMP_NUM_THREADS']=1
    environments['VIKING']['export']['OMP_STACKSIZE']=102400
    environments['VIKING']['export']['NO_AT_BRIDGE']=1
    environments['VIKING']['export']['FFLAGS']='"-I$HDF5_DIR/include"'
    environments['VIKING']['export']['LDFLAGS']='"-L$HDF5_DIR/lib"'
    environments['VIKING']['module load']=[
                            'compiler/PGI/19.10-GCC-7.3.0-2.30',
                            'system/CUDA/9.2.88-GCC-7.3.0-2.30',
                            'lang/Python/3.7.0-foss-2018b',
                            'data/HDF5/1.10.2-PGI-19.10-GCC-7.3.0-2.30']
    environments['VIKING']['module unload']=[]
    environments['VIKING']['module switch']=[]
    environments['VIKING']['misc']={}
    environments['VIKING']['misc']['ulimit']="-s 2000000"

    environments_avail=environments.keys()
    for environment in environments: environments[environment]['module purge']=[] #add empty module purge to all environments

    def __init__(self,environment_name=None):
        """
        initialise an environment

        args:
            environment_name - string identifier to choose from selection of environments stored as class attributes 
        notes:
            available options are held in environment.environments_avail
        """ 

        try:
            self.environment=Environment.environments[environment_name]
        except:
            print("WARNING: environment failed to initialise - environment_name options are {}".format(Environment.environments_avail))
            self.environment=None

        self.environment_name=environment_name
        self.command_types=[command_type for command_type in self.environment]
        self.required_modules=[module_name for module_name in self.environment['module load']]

    def create_command_string(self,listed=False):
        """
        create single command string to load/source a chosen environment
        args:
            listed - return string as list of individual commands
        notes:
            also returns command string 
        """ 

        commands=[]

        for command in ['module purge','module load','module switch','module unload','misc','export']: #user can control order of execution of commands here
            things_to_command=self.environment[command]

            if command=='module load':
                commands.extend(['module load {}'.format(module) for module in things_to_command])
            elif command=='module switch':
                commands.extend(['module switch {}'.format(module) for module in things_to_command])
            elif command=='module unload':
                commands.extend(['module unload {}'.format(module) for module in things_to_command])
            elif command=='module purge':
                commands.extend(['module purge'])
            elif command=='export':
                commands.extend([' '.join([command,''.join([thing_to_command,'=',str(value)])]) for thing_to_command,value in things_to_command.items()])
            elif command=='misc':
                commands.extend([' '.join([thing_to_command,str(self.environment[command][thing_to_command])]) for thing_to_command in things_to_command])

        if listed:
            commands=list([command for command in commands])
        else:
            commands=' ; '.join(command for command in commands)
        return commands

    def display(self,string=False,command_types=None):
        """
        print chosen environment

        notes:
        args:
            string - print in single-line string form
            command_types - choose which commands to print, default to printing all commands
        """ 

        commands=[]
        if not command_types: command_types=self.command_types

        for command in command_types: #control order of execution of commands here
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