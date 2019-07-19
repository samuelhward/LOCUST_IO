#script to generate beam deposition input for LOCUST

import context
from classes.input_classes.beam_deposition import Beam_Deposition as bd
from classes.input_classes.equilibrium import Equilibrium as eq
from classes.input_classes.wall import Wall as W
import constants
import numpy as np


runs=['U69','U70','U71','U72','U73','U74']
walls=['LCFS_DIII-D.dat_1.05','LCFS_DIII-D.dat_1.10','LCFS_DIII-D.dat_1.20','LCFS_DIII-D.dat_1.30','LCFS_DIII-D.dat_1.40','LCFS_DIII-D.dat_1.50']


filename_eq='g157418.03000'
equi=eq('','GEQDSK',filename_eq)

for run,wall in zip(runs,walls): #add the files

    filename='157418{}_birth.cdf1'.format(run)
    beam_depo=bd(run,'TRANSP_birth_gc',filename)
    limiter=W(wall,'LOCUST_2D',wall)
    beam_depo.dump_data(data_format='LOCUST_weighted',filename='{}_ptcles.dat'.format(run),equilibrium=equi)


    #scale up the weights of ASCOT simulation such that they add up to Pdep in LOCUST (1W in these runs)

    Pdep_desired=1.
    Pdep_current=np.sum(beam_depo['weight']*beam_depo['E']*constants.e_charge)
	beam_depo['weight']*=Pdep_desired/Pdep_current

   	beam_depo.dump_data(filename='input_particles_{}'.format(run),data_format='ASCOT_gc',equilibrium=equi)
    beam_depo.plot(some_equilibrium=limiter,limiters=True,axes=['X','Y'],number_bins=300)