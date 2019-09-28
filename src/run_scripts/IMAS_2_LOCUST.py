#run_scripts.IMAS_2_LOCUST.py
 
"""
Samuel Ward
06/08/2019
----
read IMAS inputs from IMAS IDSs and generates corresponding LOCUST inputs 
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
    import pathlib
    import argparse
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
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

###################################################################################################
#Main

#read in all relevant IMAS inputs and outputs
def IMAS_2_LOCUST(shot,run,path_LOCUST=pathlib.Path(''),beam_depo_GC=True,GEQDSKFIX=0,tag=''):
    """
    read all LOCUST inputs from IMAS run data  

    notes:
        remember to make sure the directory structure already exists!
        assumes first time slice
        dumps kinetic profiles for all available species        
    args:
        shot - IDS shot number 
        run - IDS run number
        path_LOCUST - path to LOCUST files in input_files dir (input_files/path_LOCUST...)
        beam_depo_GC - toggle dumping birth list at guiding-centre or particle position
        GEQDSKFIX - LOCUST-equivalent flag to optionally flip fields in GEQDSK
        tag - optional identifier tag for each set of run files produced
    """

    path_LOCUST=pathlib.Path(path_LOCUST)

    try:
        temperature_array,density_array,perturbation_array,temperature_e,density_e,beam_deposition,wall,equilibrium=run_scripts.utils.read_inputs_IMAS(shot=shot,run=run,GEQDSKFIX=GEQDSKFIX)
    except:
        print("ERROR: IMAS_2_LOCUST could not read_inputs_IMAS from IDS (shot - {shot}, run - {run})".format(shot=shot,run=run))
        return 

    try: #run dump_inputs_LOCUST without ion temperature - dump these separately due to possible multiple species

        run_scripts.utils.dump_inputs_LOCUST(temperature_e=temperature_e,density_e=density_e,equilibrium=equilibrium,beam_deposition=beam_deposition,wall=wall,beam_depo_GC=beam_depo_GC,beam_depo_weighted=True,BCHECK=False,wall_type='2D',input_path=path_LOCUST,tag=tag)
            
        for temperature in temperature_array:
            temperature.dump_data(data_format='LOCUST',filename=path_LOCUST / 'profile_Ti_Z={}.dat'.format(temperature.properties['Z']))
        
        for density in density_array:
            density.dump_data(data_format='LOCUST',filename=path_LOCUST / 'profile_ni_Z={}.dat'.format(density.properties['Z']))

        for perturbation in perturbation_array:
            perturbation.dump_data(data_format='LOCUST',filename=path_LOCUST / 'BPLASMA_n{}.dat'.format(perturbation.mode_number))

    except:
        print("ERROR: IMAS_2_LOCUST could not dump_inputs_LOCUST to LOCUST_IO/input_files/{}\n".format(path_LOCUST))
        return         

if __name__=='__main__':

    import argparse #take user command line input here
    parser=argparse.ArgumentParser(description='convert IMAS IDSs files to LOCUST input files')

    parser.add_argument('--shot',type=int,action='store',dest='shot',help="IDS shot number",required=True)
    parser.add_argument('--run',type=int,action='store',dest='run',help="IDS run number",required=True)
    parser.add_argument('--path_LOCUST',type=str,action='store',default='',dest='path_LOCUST',help="source path to IMAS inputs within input_files",required=False)
    parser.add_argument('--GC',action='store_true',default=False,dest='beam_depo_GC',help="toggle guiding centre particle list (default False)",required=False)
    parser.add_argument('--GEQDSKFIX',action='store',default=0,dest='GEQDSKFIX',help="LOCUST-equivalent flag to optionally flip fields in GEQDSK",required=False)    
    parser.add_argument('--tag',type=str,action='store',default='',dest='tag',help="optional identifier tag for each set of run files produced",required=False)    
    
    args=parser.parse_args()

    #convert 
    IMAS_2_LOCUST(shot=args.shot,run=args.run,path_LOCUST=args.path_LOCUST,beam_depo_GC=args.beam_depo_GC,GEQDSKFIX=args.GEQDSKFIX,tag=args.tag)

#################################
 
##################################################################
 
###################################################################################################