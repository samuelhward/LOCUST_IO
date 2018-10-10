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


##################################################################
#Project directory paths





#software version number 
LOCUST_IO_version=str('LOCUST_IO version 1.2.0')
numpy_version=np.version.version.split('.')
python_version=sys.version_info



#get the paths to folders within the package for importing modules, opening files etc
pwd=os.path.dirname(os.path.abspath(__file__)) #get the directory this script is in
dir_locust_io=os.path.dirname(pwd) #go one level up from that
dir_input_files=os.path.join(dir_locust_io,'input_files/') #add the / at the end so the user only needs to append filenames
dir_output_files=os.path.join(dir_locust_io,'output_files/')
dir_classes=os.path.join(dir_locust_io,'classes/')



#data which must exist to be able to run LOCUST
required_equilibrium=['rdim','zdim','rcentr','rleft','zmid','rmaxis','zmaxis','simag','sibry','bcentr',
'current','simag','rmaxis','zmaxis','sibry','fpol','pres','ffprime','pprime','psirz','qpsi','lcfs_n','limitr','lcfs_r','lcfs_z','rlim','zlim']
required_beam_deposition=['R','phi','Z','V_R','V_tor','V_Z','absorption_fraction','absorption_scaling']
required_temperature=['flux_pol_norm','T']
required_number_density=['flux_pol_norm','n']
required_perturbation=['R_2D','Z_2D','B_field_R_real','B_field_R_imag','B_field_Z_real','B_field_Z_imag','B_field_tor_real','B_field_tor_imag']

#################################

##################################################################

###################################################################################################