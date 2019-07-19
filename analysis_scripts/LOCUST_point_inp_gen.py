#LOCUST_point_inp_gen.py
"""
Samuel Ward
28/03/2019
----
check supplied 3D field components
---
usage:
    see README.md for usage
 
notes:         
	uses magnetic field components given at location of the particles to check 3D field and orientation
---
"""

import context
from classes.input_classes.beam_deposition import Beam_Deposition as bd 
from classes.input_classes.perturbation import Perturbation as pert
import numpy as np

testbd=bd(ID='ASCOT ITPA beam deposition to test field components for -DBCHECK mode',filename='ASCOT/input_nbi8.particles',data_format='ASCOT_FO')
perty=pert(ID='perturbation object to hold particle marker positions to then dump to -DBCHECK file') #generate a blank object for now to fill with data in a sec

perty.set(R_point_data=testbd['R'],Z_point_data=testbd['Z'],phi_point_data=testbd['phi'],time_point_data=np.zeros(len(testbd['phi'])))
perty.dump_data(data_format='point_data',filename='point_data.inp',BCHECK=1)


