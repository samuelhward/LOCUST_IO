#run_scripts.ASCOT_2_LOCUST.py
 
"""
Samuel Ward
26/03/2019
----
read ASCOT inputs and generate LOCUST inputs 
---
usage:
    see README.md for usage
 
notes:         
---
"""

###################################################################################################
#Preamble

import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway

if __name__=='__main__':

    try:
        import context
    except:
        raise ImportError("ERROR: context.py could not be imported!\nreturning\n")
        sys.exit(1)

try:
    import numpy as np
    import argparse
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import run_scripts.utils
except:
    raise ImportError("ERROR: LOCUST_IO/run_scripts/utils.py could not be imported!\nreturning\n")
    sys.exit(1)  

try:
    import classes.input_classes.equilibrium
except:
    raise ImportError("ERROR: LOCUST_IO/classes/input_classes/equilibrium.py could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from constants import *
except:
    raise ImportError("ERROR: LOCUST_IO/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from settings import *
except:
    raise ImportError("ERROR: LOCUST_IO/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)

###################################################################################################
#Main

#read in all relevant ASCOT inputs and outputs
def ASCOT_2_LOCUST(path_LOCUST='',path_ASCOT='',beam_depo_GC=True,filename_ASCOT_equilibrium='ASCOT_GEQDSK',GEQDSKFIX=0,species_numbers=[1],wall_type='2D',tag=''):
    """
    read all LOCUST inputs from ASCOT run data  

    notes:
        remember to make sure the directory structure already exists!
        XXX assumes equilibrium stored as GEQDSK in path_ASCOT since no way of assimilating ASCOT equilibrium for now
        takes wall from ASCOT input, not output
    args:
        path_LOCUST - path to LOCUST files in input_files dir (input_files/path_LOCUST...)
        path_ASCOT - path to ASCOT files in input_files dir (input_files/path_ASCOT...)
        beam_depo_GC - toggle dumping birth list at guiding-centre or particle position
        filename_ASCOT_equilibrium - filename of equilibrium stored AS GEQDSK in path_ASCOT
        GEQDSKFIX - LOCUST-equivalent flag to optionally flip fields in GEQDSK
        species_numbers - species number labels from input particle list to read in
        wall_type - set wall type to '2D' or '3D'
        tag - optional identifier tag for each set of run files produced
    """

    try:
        temperature_array,density_array,temperature_e,density_e,beam_deposition,wall=run_scripts.utils.read_inputs_ASCOT(input_path=path_ASCOT,beam_depo_GC=beam_depo_GC,species_numbers=species_numbers,wall_type=wall_type)
        equilibrium=classes.input_classes.equilibrium.Equilibrium(ID='made using ASCOT_2_LOCUST',data_format='GEQDSK',filename=path_ASCOT+filename_ASCOT_equilibrium,GEQDSKFIX=GEQDSKFIX)
    except:
        print("ERROR: ASCOT_2_LOCUST could not read_inputs_ASCOT from LOCUST_IO/input_files/{}\n".format(path_ASCOT))
        return 

    #try: #run dump_inputs_LOCUST without ion temperature - dump these separately due to possible multiple species
    run_scripts.utils.dump_inputs_LOCUST(temperature_e=temperature_e,density_e=density_e,equilibrium=equilibrium,beam_deposition=beam_deposition,wall=wall,beam_depo_GC=beam_depo_GC,beam_depo_weighted=True,BCHECK=False,wall_type=wall_type,input_path=path_LOCUST,tag=tag)
        
    for temperature in temperature_array:
        temperature.dump_data(data_format='LOCUST',filename='profile_Ti{}.dat'.format(temperature.properties['species_number']))
    
    for density in density_array:
        density.dump_data(data_format='LOCUST',filename='profile_ni{}.dat'.format(density.properties['species_number']))

    #except:
    #    print("ERROR: ASCOT_2_LOCUST could not dump_inputs_LOCUST to LOCUST_IO/input_files/{}\n".format(path_LOCUST))
    #    return         

if __name__=='__main__':

    import argparse #take user command line input here
    parser=argparse.ArgumentParser(description='convert ASCOT files to LOCUST input files')
    parser.add_argument('--path_LOCUST',type=str,action='store',default='',dest='path_LOCUST',help="source path to ASCOT inputs within input_files",required=False)
    parser.add_argument('--path_ASCOT',type=str,action='store',default='',dest='path_ASCOT',help="target path to LOCUST inputs within input_files",required=False)
    parser.add_argument('--GC',action='store_true',default=False,dest='beam_depo_GC',help="toggle guiding centre particle list (default False)",required=False)
    parser.add_argument('--ASCOT_eq',type=str,action='store',default='ASCOT_GEQDSK',dest='filename_ASCOT_equilibrium',help="filename of equilibrium stored as GEQDSK in path_ASCOT",required=False)
    parser.add_argument('--GEQDSKFIX',action='store',default=0,dest='GEQDSKFIX',help="LOCUST-equivalent flag to optionally flip fields in GEQDSK",required=False)    
    parser.add_argument('--species_numbers',nargs='+',type=int,action='store',default=1,dest='species_numbers',help="species number labels from input particle list to read in e.g. --species_numbers 1 3 5",required=False)
    parser.add_argument('--wall_type',type=str,action='store',default='2D',dest='wall_type',help="set wall type to '2D' or '3D'",required=False)
    parser.add_argument('--tag',type=str,action='store',default='',dest='tag',help="optional identifier tag for each set of run files produced",required=False)    
    
    args=parser.parse_args()

    #convert 
    ASCOT_2_LOCUST(run_ID=args.run_ID,shot_number=args.shot_number,path_LOCUST=args.path_LOCUST,path_ASCOT=args.path_ASCOT,beam_depo_GC=args.beam_depo_GC,GEQDSKFIX=args.GEQDSKFIX,beam_depo_weighted=args.beam_depo_weighted,tag=args.tag)

#################################
 
##################################################################
 
###################################################################################################