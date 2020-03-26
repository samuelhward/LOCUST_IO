#get_filepaths_single.py
 
"""
Samuel Ward
18/04/2019
----
script to grab filepaths for MAST V and V - single set of runs
---
usage:
 
notes:         

to do:
    could have it so that ASCOT and LOCUST time slices are all stored in a single folder, however their dimension and time are inferred from looking at the flags in the run_data files
        THIS WOULD WARRANT A LOCUST run class! this could read a run_data file and extract the flags and main parameters for example
    assumes TRANSP always outputs <11 time slices - ptherwise need to change how TRANSP files are sorted to a more robust way
---
"""

import context
import glob
import support
import processing.utils
import numpy as np
import pathlib

def get_filepaths_single(shot='29034',run='W03',dimension='2D',extra_info='NUBEAM_birth_GC',code='LOCUST'):
    """
    grabs all the filenames for the specified run, given the same directory structure
    """

    #now each code directory structure begins to diverge
    #find all sub directories for this run, for example directories containing time slices in the case for LOCUST

    time_slices=[] #time slices that we have written DFNs out at - for ascot and locust we should infer these from the directory names e.g. 2D_30
    files=[]
    
    if code == 'TRANSP':
        dirs=support.dir_output_files / code / (shot+run)
        files=list(dirs.glob('*_fi_*')) #for TRANSP, all time slice files are contained in one directory
        files.sort() #sort alphanumerically i.e. numerically since only changing parameter is index of the time slice
        files.append(files[0]) #10 is sorted to beginning
        del(files[0])
        files=[str(file).split('output_files')[1][1:] for file in files]
    else:
        dirs=list((support.dir_output_files / code / (shot+run) / extra_info).glob('{}*'.format(dimension)))

        for dir_ in dirs: #cycling through time slices                
            time_slice=int(dir_.parts[-1][len(dimension)+1:]) #get time from name of directory holding simulation e.g. 2D_10
            time_slices=np.append(time_slices,time_slice)

            #get the filename of the distribution function in this directory
            if code == 'LOCUST':
                files.append(list(dir_.glob('*.dfn'))[0])
            elif code == 'ASCOT':
                files.append(list(dir_.glob('*.h5'))[0])

        time_slices,files=processing.utils.sort_arrays(np.asarray(time_slices),np.asarray(files))
        files=[str(file).split('output_files')[1][1:] for file in files] #remove initial directories before output_files, [1:] to get rid of '/' at start of path

    return files

#################################
 
##################################################################
 
###################################################################################################
