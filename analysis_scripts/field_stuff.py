#LOCUST_data_prep.py
"""
Samuel Ward
26/03/2019
----
convert from supplied MARS-F fields and ASCOT input data to LOCUST input data using LOCUST_IO
---
usage:
    see README.md for usage
 
notes:         
LOCUST_IO version - cee6151f87ad6450f631f7f713d063c0752312a8
---
"""


import context
from classes.input_classes.perturbation import Perturbation as pert
from run_scripts.ASCOT_2_LOCUST import ASCOT_2_LOCUST as A2L
import numpy as np

path_ASCOT='ASCOT/'
path_MARSF='MARSF/'
path_LOCUST='LOCUST/'

n2_vacuum_L=pert(ID='',data_format='LOCUST',filename=path_LOCUST+'BPLASMA_n2_vac_0.dat')
n2_vacuum_M=pert(ID='',data_format='MARSF',filename=path_MARSF+'B_field_n2_vacuum_33143_2730_0.txt')


n2_vacuum_L.set(B_R=np.sqrt(n2_vacuum_L['B_field_R_real']**2+n2_vacuum_L['B_field_R_imag']**2))
n2_vacuum_L.set(B_Z=np.sqrt(n2_vacuum_L['B_field_Z_real']**2+n2_vacuum_L['B_field_Z_imag']**2))
n2_vacuum_L.set(B_tor=np.sqrt(n2_vacuum_L['B_field_tor_real']**2+n2_vacuum_L['B_field_tor_imag']**2))
n2_vacuum_L.set(B=np.sqrt(n2_vacuum_L['B_tor']**2+n2_vacuum_L['B_R']**2+n2_vacuum_L['B_Z']**2))

n2_vacuum_M.set(B_R=np.sqrt(n2_vacuum_M['B_field_R_real']**2+n2_vacuum_M['B_field_R_imag']**2))
n2_vacuum_M.set(B_Z=np.sqrt(n2_vacuum_M['B_field_Z_real']**2+n2_vacuum_M['B_field_Z_imag']**2))
n2_vacuum_M.set(B_tor=np.sqrt(n2_vacuum_M['B_field_tor_real']**2+n2_vacuum_M['B_field_tor_imag']**2))
n2_vacuum_M.set(B=np.sqrt(n2_vacuum_M['B_tor']**2+n2_vacuum_M['B_R']**2+n2_vacuum_M['B_Z']**2))

n2_vacuum_L.plot(key='B_R')
n2_vacuum_M.plot(key='B_R')

n2_vacuum_L.look()
n2_vacuum_M.look()

