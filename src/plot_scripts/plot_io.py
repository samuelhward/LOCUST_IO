#plot_io.py

'''
Samuel Ward
06/11/2019
----
plot a LOCUST_IO object from command line
---
notes: 
    essentially wraps argparse to feed command line args straight to a quick instantiation/plotting of a LOCUST_IO object
---
'''


##################################################################
#Preamble

import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway

if __name__=='__main__':

    try:
        import context
    except:
        raise ImportError("ERROR: context.py could not be imported!\nreturning\n")
        sys.exit(1)

try:
    import ast
    import numpy as np
    import matplotlib.pyplot as plt
    import pathlib
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import classes.input_classes.beam_deposition
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/beam_deposition.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import classes.input_classes.equilibrium
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/equilibrium.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import classes.input_classes.number_density
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/number_density.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import classes.input_classes.perturbation
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/perturbation.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import classes.input_classes.rotation
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/rotation.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import classes.input_classes.temperature
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/temperature.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import classes.input_classes.wall
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/wall.py could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import classes.output_classes.distribution_function
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/output_classes/distribution_function.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import classes.output_classes.moments
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/output_classes/moments.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import classes.output_classes.orbits
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/output_classes/orbits.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import classes.output_classes.particle_list
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/output_classes/particle_list.py could not be imported!\nreturning\n")
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
    import constants
except:
    raise ImportError("ERROR: LOCUST_IO/src/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from settings import *
except:
    raise ImportError("ERROR: LOCUST_IO/src/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)


##################################################################

if __name__=='__main__': 

    dispatch_table_objects={} #create dispatch table to easily assign user plot choice to corresponding LOCUST_IO class type
    dispatch_table_objects['beam_deposition']=classes.input_classes.beam_deposition.Beam_Deposition
    dispatch_table_objects['equilibrium']=classes.input_classes.equilibrium.Equilibrium
    dispatch_table_objects['number_density']=classes.input_classes.number_density.Number_Density
    dispatch_table_objects['perturbation']=classes.input_classes.perturbation.Perturbation
    dispatch_table_objects['rotation']=classes.input_classes.rotation.Rotation
    dispatch_table_objects['temperature']=classes.input_classes.temperature.Temperature
    dispatch_table_objects['wall']=classes.input_classes.wall.Wall
    dispatch_table_objects['distribution_function']=classes.output_classes.distribution_function.Distribution_Function
    dispatch_table_objects['moments']=classes.output_classes.moments.Moments
    dispatch_table_objects['orbits']=classes.output_classes.orbits.Orbits
    dispatch_table_objects['particle_list']=classes.output_classes.particle_list.Final_Particle_List
    data_types_available=[key for key in dispatch_table_objects.keys()]

    import argparse
    parser=argparse.ArgumentParser(description='')

    parser.add_argument('--data_type',type=str,action='store',dest='data_type',help="target data type - available options: {} ".format(data_types_available),default=None,required=True) #in addition to the usual arguements given within Python, 'data_type' is needed to select the correct LOCUST_IO class
    parser.add_argument('--read_settings',nargs='+',type=str,action='store',dest='read_settings',help="settings for LOCUST_IO object instantiation e.g. data_format='GEQDSK'",default={})
    parser.add_argument('--plot_settings',nargs='+',type=str,action='store',dest='plot_settings',help="settings for LOCUST_IO object plot method e.g. fill=False",default={})
    
    args=parser.parse_args()

    read_settings=run_scripts.utils.command_line_arg_parse_dict(args.read_settings) #these are dict-style objects so need an extra layer of parsing
    plot_settings=run_scripts.utils.command_line_arg_parse_dict(args.plot_settings)

    if args.data_type not in data_types_available:
        print("ERROR: {data_type_supplied} is unavailable --data_type - available options: {data_types_available} ".format(data_type_supplied=args.data_type,data_types_available=data_types_available))
    else:    

        if 'filename' in read_settings and 'input_files' in str(read_settings['filename']): #in case user supplies full filepath to target
            read_settings['filename']=pathlib.Path(read_settings['filename']).relative_to(support.dir_input_files)
        if 'axes' in plot_settings: #apply some extra formatting if defining plot axes at command line such as axes="['R','Z']"
            plot_settings['axes']=ast.literal_eval(plot_settings['axes'])
        read_settings['ID']='' #include dummy ID to assign the object to be plotted, since this is stipulated in LOCUST_IO

        try:
            obj=dispatch_table_objects[args.data_type](**read_settings)
            obj.plot(**plot_settings)
        except:
            print("ERROR: plot_io could not plot")

#################################

##################################################################

###################################################################################################