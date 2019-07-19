import context
from classes.input_classes.beam_deposition import Beam_Deposition as bd
from classes.input_classes.equilibrium import Equilibrium as equi
import processing.utils
from constants import *
import pathlib


shot_time_path=pathlib.Path('29034.360')

filename_BBNBI=shot_time_path / 'mast_29034_360ms_injSS_BBNBI_ward.particles'
filename_eq=shot_time_path / 'g29034'

eq=equi('',data_format='GEQDSK',filename=filename_eq,GEQDSKFIX=1)
mybd=bd('',data_format='ASCOT_FO',filename=filename_BBNBI)
mybd.set(V_pitch=processing.utils.pitch_calc_2D(mybd,eq))
mybd.set(E=0.5*species_mass*(mybd['V_R']**2+mybd['V_Z']**2+mybd['V_tor']**2)/e_charge)


mybd.dump_data(data_format='LOCUST_FO_weighted',filename=shot_time_path / 'ptcles.dat_BBNBI_weighted_FO')
mybd.dump_data(data_format='LOCUST_FO',filename=shot_time_path / 'ptcles.dat_BBNBI_FO')
mybd.dump_data(data_format='LOCUST_GC_weighted',filename=shot_time_path / 'ptcles.dat_BBNBI_weighted_GC',equilibrium=eq)

mybd.set(R=mybd['R'][0:100000])
mybd.dump_data(data_format='ASCOT_FO',filename=shot_time_path / 'input.particles_BBNBI_FO',equilibrium=eq)
mybd.dump_data(data_format='ASCOT_GC',filename=shot_time_path / 'input.particles_BBNBI_GC',equilibrium=eq)

#mybd.plot(real_scale=True,number_bins=400)
