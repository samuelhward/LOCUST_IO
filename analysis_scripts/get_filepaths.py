#get_filepaths.py
 
"""
Samuel Ward
18/04/2019
----
script to grab filepaths for MAST V and V
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

def get_filepaths(shot='29034',run='W03',dimension_LOCUST='2D',dimension_ASCOT='2D',extra_info_LOCUST='NUBEAM_birth_GC',extra_info_ASCOT='NUBEAM_birth_GC'):
    """
    grabs all the LOCUST, ASCOT and TRANSP filenames for the specified run, given the same directory structure
    """

    #now each code directory structure begins to diverge
    #find all sub directories for this run, for example directories containing time slices in the case for LOCUST

    dirs_locust=list((support.dir_output_files / 'LOCUST' / (shot+run) / extra_info_LOCUST).glob('{}*'.format(dimension_LOCUST)))
    dir_transp=support.dir_output_files / 'TRANSP' / (shot+run)
    dirs_ascot=list((support.dir_output_files / 'ASCOT'/ (shot+run) / extra_info_ASCOT).glob('{}*'.format(dimension_ASCOT)))
    time_slices_locust=[] #time slices that we have written DFNs out at
    time_slices_ascot=[] #for ascot and locust we should infer these from the directory names e.g. 2D_30

    files_locust=[]
    files_ascot=[]
    files_transp=list(dir_transp.glob('*_fi_*')) #for TRANSP, all time slice files are contained in one directory
    files_transp.sort() #sort alphanumerically i.e. numerically since only changing parameter is index of the time slice
    files_transp.append(files_transp[0]) #10 is sorted to beginning
    del(files_transp[0])

    for dir_locust,dir_ascot in zip(dirs_locust,dirs_ascot): #cycling through time slices
            
        time_slice_locust=int(dir_locust.parts[-1][len(dimension_LOCUST)+1:]) #get time from name of directory holding LOCUST simulation e.g. 2D_10
        time_slice_ascot=int(dir_ascot.parts[-1][len(dimension_ASCOT)+1:])

        time_slices_locust=np.append(time_slices_locust,time_slice_locust)
        time_slices_ascot=np.append(time_slices_ascot,time_slice_ascot)

        #get the filename of the distribution function in this directory
        files_locust.append(list(dir_locust.glob('*.dfn'))[0])
        files_ascot.append(list(dir_ascot.glob('*.h5'))[0])

    time_slices_locust,files_locust=processing.utils.sort_arrays(np.asarray(time_slices_locust),np.asarray(files_locust))
    time_slices_ascot,files_ascot=processing.utils.sort_arrays(np.asarray(time_slices_ascot),np.asarray(files_ascot))

    files_locust=[str(file_locust).split('output_files')[1][1:] for file_locust in files_locust] #remove initial directories before output_files, [1:] to get rid of '/' at start of path
    files_ascot=[str(file_ascot).split('output_files')[1][1:] for file_ascot in files_ascot]
    files_transp=[str(file_transp).split('output_files')[1][1:] for file_transp in files_transp]

    return files_locust,files_ascot,files_transp

#################################
 
##################################################################
 
###################################################################################################
