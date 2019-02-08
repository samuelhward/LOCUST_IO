#run_scripts.TRANSP_2_LOCUST.py
 
"""
Samuel Ward
31/01/2019
----
read TRANSP I/O and generate LOCUST inputs 
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

#read in all relevant TRANSP inputs and outputs

def TRANSP_2_LOCUST(run_ID,shot_number,input_path='',output_path='',beam_depo_GC=True,beam_depo_weighted=True,tag=''):
    """
    read all LOCUST inputs from TRANSP run data  

    notes:
    args:
        run_ID - TRANSP run_ID e.g. W01
        shot_number - TRANSP shot number e.g. 29034
        input_path - path to target in input_files dir (input_files/path/)
        output_path - path to target in output_files dir (output_files/path/)
        beam_depo_GC - toggle dumping birth list at guiding-centre or particle position
        beam_depo_weighted - toggle dumping weighted birth list 
        tag - optional identifier tag for each set of run files produced
    """

    try:
        temperature_i,temperature_e,density_e,equilibrium,beam_deposition,wall=run_scripts.utils.read_inputs_TRANSP(run_ID=run_ID,shot_number=shot_number,input_path=input_path,output_path=output_path,beam_depo_GC=beam_depo_GC)
    except:    

    try:
        run_scripts.utils.dump_inputs_LOCUST(temperature_i=temperature_i,temperature_e=temperature_e,density_e=density_e,equilibrium=equilibrium,beam_deposition=beam_deposition,wall=wall,beam_depo_GC=beam_depo_GC,beam_depo_weighted=beam_depo_weighted,BCHECK=False,wall_type='2D',tag=tag)
    except:

if __name__=='__main__':

    import argparse #take user command line input here
    parser=argparse.ArgumentParser(description='convert TRANSP files to LOCUST input files')
    parser.add_argument('--run_ID',type=str,action='store',dest='run_ID',help='TRANSP run ID',required=True)
    parser.add_argument('--shot_number',type=str,action='store',dest='shot_number',help='TRANSP shot number',required=True)
    parser.add_argument('--path_in',type=str,action='store',dest='input_path',help='path to TRANSP inputs within input_files',required=False)
    parser.add_argument('--path_out',type=str,action='store',dest='output_path',help='path to TRANSP outputs within output_files',required=False)
    parser.add_argument('--GC',action='store_true',default=False,dest='beam_depo_GC',help='toggle guiding centre particle list',required=False)
    parser.add_argument('--weight',action='store_true',default=False,dest='beam_depo_weighted',help='toggle weighted particle list',required=False)    
    parser.add_argument('--tag',type=str,action='store',dest='tag',help='optional identifier tag for each set of run files produced',required=False)    
    
    args=parser.parse_args()

    #convert 
    TRANSP_2_LOCUST(run_ID=args.run_ID,shot_number=args.shot_number,input_path=args.input_path,output_path=args.output_path,beam_depo_GC=args.beam_depo_GC,beam_depo_weighted=args.beam_depo_weighted,tag=args.tag)

#################################
 
##################################################################
 
###################################################################################################