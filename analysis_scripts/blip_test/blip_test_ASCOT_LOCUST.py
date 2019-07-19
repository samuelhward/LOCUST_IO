#blip test

#generate particle list for ASCOT and LOCUST that are injected along the magnetic axis for testing dfn normalisation and core slowing-down


import context
from classes.input_classes.beam_deposition import Beam_Deposition
from classes.input_classes.equilibrium import Equilibrium
import processing.utils
import run_scripts.utils
import constants
import numpy as np
import pathlib

filename_eq='g157418.03000'
equi=Equilibrium(ID=filename_eq,data_format='GEQDSK',filename=pathlib.Path('blip_test') / filename_eq,GEQDSKFIX=1)

#need the R Z phi V_par V and weight of particles in our beam deposition
number_of_particles=100
energy_of_particles_ev=80000 #energy of injected particles in keV
velocity_of_particles=np.sqrt(constants.species_charge*energy_of_particles_ev*2./constants.species_mass)
my_beam_depo=Beam_Deposition(ID='blip test') #make a blank Beam_Deposition object
my_beam_depo.set(R=np.full(number_of_particles,equi['rmaxis']+0.1),
                 Z=np.full(number_of_particles,equi['zmaxis']), #shift particles slightly away from axis to avoid problems
    phi=np.linspace(0.,2.*constants.pi,number_of_particles),
    V=np.full(number_of_particles,velocity_of_particles),
    weight=np.full(number_of_particles,1.0),
    V_R=np.full(number_of_particles,0.0),
    V_tor=np.full(number_of_particles,-1.0*velocity_of_particles),
    V_Z=np.full(number_of_particles,0.0),
    E=np.full(number_of_particles,energy_of_particles_ev)
    )

my_beam_depo.set(V_pitch=processing.utils.pitch_calc_2D(my_beam_depo,equi))
my_beam_depo['X'],my_beam_depo['Y']=processing.utils.RphiZ_to_XYZ(my_beam_depo['R'],my_beam_depo['phi'])
my_beam_depo.dump_data(data_format='ASCOT_GC',filename=pathlib.Path('blip_test') / 'input.particles',equilibrium=equi)
my_beam_depo.dump_data(data_format='LOCUST_FO',filename=pathlib.Path('blip_test') / 'ptcles.dat',equilibrium=equi)
