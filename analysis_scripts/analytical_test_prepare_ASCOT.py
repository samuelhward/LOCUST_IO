#analytical test

#generate particle list for ASCOT from LOCUST inputs
#test is for the -DTEST mode built into LOCUST


import context
from classes.input_classes.beam_deposition import Beam_Deposition
from classes.input_classes.equilibrium import Equilibrium
from classes.input_classes.temperature import Temperature
from classes.input_classes.number_density import Number_Density
from classes.input_classes.wall import Wall
import processing.utils
import run_scripts.utils
import constants
import numpy as np


#INPUTS



equi=Equilibrium('',data_format='GEQDSK',filename='locust/g157418.03000',GEQDSKFIX=1)

beam_deposition=Beam_Deposition('locust particle list',data_format='LOCUST_FO',filename='locust/ptcles_15-01-2019_21-57-56.430.dat')
beam_deposition['X'],beam_deposition['Y']=processing.utils.RphiZ_to_XYZ(beam_deposition['R'],beam_deposition['phi'])
beam_deposition.set(weight=np.full(len(beam_deposition['X']),1.0))

temperature_i=Temperature('locust ion temperature',data_format='LOCUST',filename='locust/profile_Ti.dat',species='ions')
temperature_e=Temperature('locust electron temperature',data_format='LOCUST',filename='locust/profile_Te.dat',species='electrons')

density_e=Number_Density('locust electron density',data_format='LOCUST',filename='locust/profile_ne.dat',species='electrons')

rotation_toroidal=np.zeros(len(temperature_e['T'])) #set rotation = 0

wall=Wall('R=1.5 limiter wall','LOCUST_2D','locust/LCFS_DIII-D.dat')


#OUTPUTS

run_scripts.utils.ASCOT_run_gen(temperature_e=temperature_e,temperature_i=temperature_i,density_e=density_e,density_i=density_e,rotation_toroidal=rotation_toroidal,equilibrium=equi,wall=wall,beam_deposition=beam_deposition,guiding_centre=False)