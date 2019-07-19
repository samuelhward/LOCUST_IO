#LOCUST_sample_data_prep.py
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
from classes.input_classes.beam_deposition import Beam_Deposition as bd
from run_scripts.ASCOT_2_LOCUST import ASCOT_2_LOCUST as A2L
import numpy as np
import pathlib

path_ASCOT=pathlib.Path('ASCOT')
path_MARSF=pathlib.Path('MARSF')
path_LOCUST=pathlib.Path('LOCUST')

n2_total_0=pert(ID='',data_format='MARSF',filename=path_MARSF+'B_field_n2_total_33143_2730_0.txt')
n2_total_340=pert(ID='',data_format='MARSF',filename=path_MARSF+'B_field_n2_total_33143_2730_340.txt')
n2_vacuum_0=pert(ID='',data_format='MARSF',filename=path_MARSF+'B_field_n2_vacuum_33143_2730_0.txt')
n2_vacuum_340=pert(ID='',data_format='MARSF',filename=path_MARSF+'B_field_n2_vacuum_33143_2730_340.txt')
n6_total=pert(ID='',data_format='MARSF',filename=path_MARSF+'B_field_n6_total_33143_2730_340.txt')
n6_vacuum=pert(ID='',data_format='MARSF',filename=path_MARSF+'B_field_n6_vacuum_33143_2730_340.txt')

'''
n2_vacuum_0.set(B_R=np.sqrt(n2_vacuum_0['B_field_R_real']**2+n2_vacuum_0['B_field_R_imag']**2))
n2_vacuum_0.set(B_Z=np.sqrt(n2_vacuum_0['B_field_Z_real']**2+n2_vacuum_0['B_field_Z_imag']**2))
n2_vacuum_0.set(B_tor=np.sqrt(n2_vacuum_0['B_field_tor_real']**2+n2_vacuum_0['B_field_tor_imag']**2))
n2_vacuum_0.set(B=np.sqrt(n2_vacuum_0['B_tor']**2+n2_vacuum_0['B_R']**2+n2_vacuum_0['B_Z']**2))
n2_vacuum_0.plot(key='B_field_R_real')
'''

n2_total_0.dump_data(data_format='LOCUST',filename=path_LOCUST+'BPLASMA_n2_resp_0.dat')
n2_total_340.dump_data(data_format='LOCUST',filename=path_LOCUST+'BPLASMA_n2_resp_340.dat')
n2_vacuum_0.dump_data(data_format='LOCUST',filename=path_LOCUST+'BPLASMA_n2_vac_0.dat')
n2_vacuum_340.dump_data(data_format='LOCUST',filename=path_LOCUST+'BPLASMA_n2_vac_340.dat')
n6_total.dump_data(data_format='LOCUST',filename=path_LOCUST+'BPLASMA_n6_resp.dat')
n6_vacuum.dump_data(data_format='LOCUST',filename=path_LOCUST+'BPLASMA_n6_vac.dat')

A2L(path_LOCUST='LOCUST/',path_ASCOT='ASCOT/',beam_depo_GC=False,filename_ASCOT_equilibrium='g033143.02730',GEQDSKFIX=0,species_numbers=[1,2],wall_type='2D',tag='')

mybd=bd('ITPA beam deposition',data_format='LOCUST_FO_weighted',filename=path_LOCUST+'ptcles.dat')
mybd.dump_data(data_format='LOCUST_FO',filename=path_LOCUST+'ptcles_unweighted.dat') #need to also generate a non-weighted version of the particle list due to -DSPLIT mode being incompatible with weighted lists (need -DSPLIT to generate a final particle list and prompt loss map)
