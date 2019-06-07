#GENERATE_DUMMY_IDS.py

'''
Samuel Ward
07/06/2019
----
Run this file to generate an IDS filled with example LOCUST inputs for testing
---
notes:
    some parts are not fully generated and require external file
---
'''



##################################################################
#Preamble

import context
import numpy as np

from classes.input_classes.temperature import Temperature as T
from classes.input_classes.number_density import Number_Density as ND
from classes.input_classes.beam_deposition import Beam_Deposition as BD
from classes.input_classes.equilibrium import Equilibrium as EQ
from classes.input_classes.perturbation import Perturbation as P
from classes.input_classes.wall import Wall as W

##################################################################
#Main

use_temperature=True
use_number_density=True
use_beam_deposition=True
use_equilibrium=False
use_perturbation=False
use_wall=False

shot=2
run=2

if use_temperature:
    DUMMY_profile=np.linspace(1,100,100)
    DUMMY_flux_pol=np.linspace(0,1,100)
    temp_i=T(ID='',species='ions')
    temp_e=T(ID='',species='electrons')
    temp_i.set(flux_pol=DUMMY_flux_pol)
    temp_i.set(T=DUMMY_profile)
    temp_e.set(flux_pol=DUMMY_flux_pol)
    temp_e.set(T=DUMMY_profile)
    temp_i.dump_data(data_format='IDS',shot=shot,run=run,species='ions')
    temp_e.dump_data(data_format='IDS',shot=shot,run=run,species='electrons')

if use_number_density:
    DUMMY_profile=np.linspace(1,100,100)
    DUMMY_flux_pol=np.linspace(0,1,100)
    numd_i=ND(ID='',species='ions')
    numd_e=ND(ID='',species='electrons')
    numd_i.set(flux_pol=DUMMY_flux_pol)
    numd_i.set(T=DUMMY_profile)
    numd_e.set(flux_pol=DUMMY_flux_pol)
    numd_e.set(T=DUMMY_profile)
    numd_i.dump_data(data_format='IDS',shot=shot,run=run,species='ions')
    numd_e.dump_data(data_format='IDS',shot=shot,run=run,species='electrons')

if use_beam_deposition:
    beam=BD(ID='')
    beam.set(R=np.random.uniform(size=10))
    beam.set(phi=np.random.uniform(size=10))
    beam.set(Z=np.random.uniform(size=10))
    beam.set(V_R=np.random.uniform(size=10))
    beam.set(V_tor=np.random.uniform(size=10))
    beam.set(V_Z=np.random.uniform(size=10))
    beam.set(weight=np.random.uniform(size=10))
    beam.dump_data(data_format='IDS',shot=shot,run=run)

if use_equilibrium:
    filename_equilibrium='g157418'
    equil=EQ(ID='',filename=filename_equilibrium,data_format='GEQDSK')
    equil.dump_data(data_format='IDS',shot=shot,run=run)

if use_perturbation:
    filename_perturbation='BPLASMA_n3'
    pert=P(ID='')

if use_wall:
    filename_wall=''
    wall=W(ID='')

#################################

##################################################################

###################################################################################################