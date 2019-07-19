#generate the equivalent ASCOT run for 157418U46/47 TRANSP run 

import numpy as np
import context
import run_scripts.utils
from classes.input_classes.beam_deposition import Beam_Deposition as bd
from classes.input_classes.equilibrium import Equilibrium as eq
from classes.input_classes.temperature import Temperature
from classes.input_classes.number_density import Number_Density



#XXX STUCK BECAUSE MY ASCOT_DUMP_GC NEEDS FULL 3D VELOCITY VECTOR STILL
#NEED TO MAKE A GET_FBM PTCLES LIST THAT I CAN RUN ASCOT AND LOCUST OFF

#XXX EDIT, ASCOT INPUT FOR GC DOES NOT REQUIRE 3D VELOCITY VECTOR - THIS IS POSSIBLE NOW FOR GC RUNS

filename_eq='g157418.03000'
equi=eq('','GEQDSK',filename_eq)

filename_beam_depo='U46U47_ptcles.dat'
my_beam_deposition=bd('U46/U47 beam deposition at guiding centre','LOCUST',filename_beam_depo)

filename_Ti='profile_Ti.dat'
my_temperature_i=Temperature('U46/U47 ion temperature','LOCUST',filename_Ti)
filename_Te='profile_Te.dat'
my_temperature_e=Temperature('U46/U47 electron temperature','LOCUST',filename_Te)

filename_ni='profile_ni.dat'
my_density_i=Number_Density('U46/U47 ion density','LOCUST',filename_ni)
filename_ne='profile_ne.dat'
my_density_e=Number_Density('U46/U47 electron density','LOCUST',filename_ne)

my_rotation_toroidal=np.zeros(len(my_temperature_e['T'])) #set rotation = 0

filename_ASCOT_profile='input.plasma_1d'
dump_profiles_ASCOT(filename_ASCOT_profile,my_temperature_i,my_temperature_e,my_density_i,my_density_e,my_rotation_toroidal):
filename_ASCOT_wall='input.wall_2d'
dump_wall_ASCOT(filename_ASCOT_wall,equi)
filename_ASCOT_beam_deposition='input.particles'
my_beam_deposition.dump_data('ASCOT_gc',filename_ASCOT_beam_deposition,equilibrium=equi) #need to supply equilibrium to write out interpolated B field at particle location


ASCOT_run_gen(run_file='ascot4.cmd',initialdir=None,output_file='ascot.out',max_proc=50,min_proc=25,error_file='ascot.err',executable='test_ascot'):