#get_fbm.py
 
"""
Samuel Ward
20/02/2019
----
execute get_fbm command line program 
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


if __name__=='__main__':

    import argparse #take user command line input here
    parser=argparse.ArgumentParser(description='convert TRANSP files to ASCOT input files')
    parser.add_argument('--run_ID',type=str,action='store',dest='run_ID',help='TRANSP run ID',required=True)
    parser.add_argument('--shot_number',type=str,action='store',dest='shot_number',help='TRANSP shot number',required=True)
    parser.add_argument('--number_files',type=int,action='store',default=1,dest='number_files',help='total number of .DATA# files to extract CDF from',required=True)       
    parser.add_argument('--GC',action='store_true',default=False,dest='GC',help='toggle whether to generate set of CDFs at particle position or guiding centre',required=False)
    parser.add_argument('--path_TRANSP',type=str,action='store',default='',dest='path_TRANSP',help='path to TRANSP files in output_files dir (output_files/path_TRANSP...)',required=False)
    parser.add_argument('--device',type=str,action='store',default='d3d',dest='device',help='device code for machine under study',required=False)    

    args=parser.parse_args()

    if args.GC:
        particle_position=False
        guiding_centre=True
    else:
        particle_position=True
        guiding_centre=False

    run_scripts.utils.TRANSP_get_fbm_FI_CDF(run_ID=args.run_ID,shot_number=args.shot_number,number_files=args.number_files,particle_position=particle_position,guiding_centre=guiding_centre,path_TRANSP=args.path_TRANSP,device=args.device)

#################################
 
##################################################################
 
###################################################################################################