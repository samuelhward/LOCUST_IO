#run_scripts.TRANSP_2_LOCUST.py
 
'''
Samuel Ward
31/01/2019
----
read TRANSP I/O and generate LOCUST inputs 
---
usage:
    see README.md for usage
 
notes:         
---
'''

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
    import settings
except:
    raise ImportError("ERROR: LOCUST_IO/src/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)

###################################################################################################
#Main

#read in all relevant TRANSP inputs and outputs

def TRANSP_2_LOCUST(run_ID,shot_number,path_LOCUST=pathlib.Path(''),path_TRANSP=pathlib.Path(''),beam_depo_GC=True,beam_depo_number=None,GEQDSKFIX1=0,GEQDSKFIX2=0,beam_depo_weighted=True,tag=''):
    """
    read all LOCUST inputs from TRANSP run data  

    notes:
        remember to make sure the directory structure already exists!
    args:
        run_ID - TRANSP run_ID e.g. W01
        shot_number - TRANSP shot number e.g. 29034
        path_LOCUST - path to LOCUST files in input_files dir (input_files/path_LOCUST...)
        path_TRANSP - path to TRANSP files in input_files dir (input_files/path_TRANSP...)
        beam_depo_GC - toggle dumping birth list at guiding-centre or particle position
        beam_depo_number - integer number of beam depo file to read elif None then read all available beam depositions and combine into single object 
        GEQDSKFIX1 - LOCUST-equivalent flag to optionally flip fields in GEQDSK
        GEQDSKFIX2 - LOCUST-equivalent flag to optionally flip fields in GEQDSK
        beam_depo_weighted - toggle dumping weighted birth list 
        tag - optional identifier tag for each set of run files produced
    """

    try:
        temperature_i,temperature_e,density_e,equilibrium,beam_deposition,wall=run_scripts.utils.read_inputs_TRANSP(run_ID=run_ID,shot_number=shot_number,input_path=path_TRANSP,beam_depo_GC=beam_depo_GC,beam_depo_number=beam_depo_number,GEQDSKFIX1=GEQDSKFIX1,GEQDSKFIX2=GEQDSKFIX2)
    except:
        print("ERROR: TRANSP_2_LOCUST could not read_inputs_TRANSP from LOCUST_IO/data/input_files/{}\n".format(path_TRANSP))
        return 

    try:
        run_scripts.utils.dump_inputs_LOCUST(temperature_i=temperature_i,temperature_e=temperature_e,density_e=density_e,equilibrium=equilibrium,beam_deposition=beam_deposition,wall=wall,beam_depo_GC=beam_depo_GC,beam_depo_weighted=beam_depo_weighted,BCHECK=False,wall_type='2D',input_path=path_LOCUST,tag=tag)
    except:
        print("ERROR: TRANSP_2_LOCUST could not dump_inputs_LOCUST to LOCUST_IO/data/input_files/{}\n".format(path_LOCUST))
        return         

if __name__=='__main__':

    import argparse #take user command line input here
    parser=argparse.ArgumentParser(description='convert TRANSP files to LOCUST input files')
    parser.add_argument('--run_ID',type=str,action='store',dest='run_ID',help="TRANSP run ID",required=True)
    parser.add_argument('--shot_number',type=str,action='store',dest='shot_number',help="TRANSP shot number",required=True)
    parser.add_argument('--path_LOCUST',type=str,action='store',default='',dest='path_LOCUST',help="source path to TRANSP inputs within input_files",required=False)
    parser.add_argument('--path_TRANSP',type=str,action='store',default='',dest='path_TRANSP',help="target path to LOCUST inputs within input_files",required=False)
    parser.add_argument('--GC',action='store_true',default=False,dest='beam_depo_GC',help="toggle guiding centre particle list (default False)",required=False)
    parser.add_argument('--beam_depo_number',action='store',default=None,dest='beam_depo_number',help="integer number of beam depo file to read elif None then read all available beam depositions and combine into single object",required=False)    
    parser.add_argument('--GEQDSKFIX1',action='store',default=0,dest='GEQDSKFIX1',help="LOCUST-equivalent flag to optionally flip fields in GEQDSK",required=False)    
    parser.add_argument('--GEQDSKFIX2',action='store',default=0,dest='GEQDSKFIX2',help="LOCUST-equivalent flag to optionally flip fields in GEQDSK",required=False)    
    parser.add_argument('--weight',action='store_true',default=False,dest='beam_depo_weighted',help="toggle weighted particle list (default False)",required=False)    
    parser.add_argument('--tag',type=str,action='store',default='',dest='tag',help="optional identifier tag for each set of run files produced",required=False)    
    
    args=parser.parse_args()

    path_LOCUST=pathlib.Path(args.path_LOCUST)
    path_TRANSP=pathlib.Path(args.path_TRANSP)

    #convert 
    TRANSP_2_LOCUST(run_ID=args.run_ID,shot_number=args.shot_number,path_LOCUST=path_LOCUST,path_TRANSP=path_TRANSP,beam_depo_GC=args.beam_depo_GC,beam_depo_number=args.beam_depo_number,GEQDSKFIX1=args.GEQDSKFIX1,beam_depo_weighted=args.beam_depo_weighted,tag=args.tag)

#################################
 
##################################################################
 
###################################################################################################