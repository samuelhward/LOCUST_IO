#LOCUST_field_data_gen.py
"""
Samuel Ward
07/04/2019
----
generate multiple LOCUST fields with upper/lower coils rotated
---
usage:
    see README.md for usage
 
notes:         
	assumes that n2 and n6 perturbations will rotate together
---
"""


import context
from classes.input_classes.perturbation import Perturbation as pert  
import numpy as np
import pathlib


phases=np.linspace(0,180,10)

n2=True
n6=True

filepaths=[]
harmonics=[]


if n2:
    filepaths.append(pathlib.Path('ASCOT')/'dryan_data'/'33143_2730_mars_data')
    harmonics.append('n2')
if n6:
    filepaths.append(pathlib.Path('ASCOT')/'dryan_data'/'33143_2730_n_6/')
    harmonics.append('n6')

for filepath,harmonic in zip(filepaths,harmonics): #cycle through harmonics
	for phase in phases: #for each harmonic vary the phase between the upper and lower coils

		perturb=pert(ID='harmonic - {harmonic}, phase - {phase}'.format(harmonic=harmonic,phase=phase),data_format='MARSF_bplas',filename=filepath,response=True,ideal=False,phase_shift=phase)
		perturb.dump_data(data_format='LOCUST',filename='LOCUST/coil_scans/BPLASMA_{harmonic}_{phase}_resist_resp.dat'.format(harmonic=harmonic,phase=phase))
		perturb.plot()

