#support.py

'''
Samuel Ward
13/12/2017
----
File which holds supporting functions and variables for the LOCUST-IO package
---
notes: 
---
'''


##################################################################
#Preamble
import sys
import os
import numpy as np
import pathlib


##################################################################
#Project directory paths

#software version number 
LOCUST_IO_version=str('LOCUST_IO version 1.2.0')
numpy_version=np.version.version.split('.') 
python_version=sys.version_info

#look for input_files and output_files directories
thisdirectory=pathlib.Path.cwd()
parts=list(thisdirectory.parts)
parts.reverse()
for counter,level in enumerate(parts): #find LOCUST_IO directory by looking for licence, then if user deletes licence code will not work
    if 'LOCUST_IO_LICENCE.md' in [str(list(path.parts)[-1]) for path in thisdirectory.parents[counter].glob('*')]:
        dir_locust_io=thisdirectory.parents[counter]
        break
dir_input_files=dir_locust_io / 'data' / 'input_files'
dir_output_files=dir_locust_io / 'data' / 'output_files'
dir_cache_files=dir_locust_io / 'data' / 'cache_files'
dir_classes=dir_locust_io / 'classes'

#data which must exist to be able to run LOCUST
required_equilibrium=['rdim','zdim','rcentr','rleft','zmid','rmaxis','zmaxis','simag','sibry','bcentr',
'current','simag','rmaxis','zmaxis','sibry','fpol','pres','ffprime','pprime','psirz','qpsi','lcfs_n','limitr','lcfs_r','lcfs_z','rlim','zlim']
required_beam_deposition=['R','phi','Z','V_R','V_tor','V_Z','weight','absorption_fraction','absorption_scaling']
required_temperature=['flux_pol_norm','T']
required_number_density=['flux_pol_norm','n']
required_perturbation=['R_2D','Z_2D','dB_field_R_real','dB_field_R_imag','dB_field_Z_real','dB_field_Z_imag','dB_field_tor_real','dB_field_tor_imag']
required_rotation=['rotation']

#################################

##################################################################

###################################################################################################