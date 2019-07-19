#!/usr/bin/env python

"""
Samuel Ward
09/06/2018
----
run script for converting NEMO beam depositions to LOCUST input files

---
notes:

"""

##################################################################
#Preamble

import sys
import numpy as np
import os

import context
from classes.input_classes.equilibrium import Equilibrium
from classes.input_classes.beam_deposition import Beam_Deposition 
from classes.input_classes.temperature import Temperature
from classes.input_classes.number_density import Number_Density
from processing import process_input as prin


##################################################################
#Main Code

pi=np.pi #define pi
au=1.660539e-27 #define atomic mass unit
species_charge=1.60217662e-19 #electron charge
shot=44
runs=[221,331,551]
runs=[221]

#######################################################################################################################################
#read in all the data

for run in runs:

    #read in all the required data
    equilibrium=Equilibrium(ID='ITER_NEMO Benchmark shot=44 run={} - EQ'.format(run),data_format='IDS',shot=shot,run=run)
    particles=Beam_Deposition(ID='ITER_NEMO Benchmark shot=44 run={} - ptcles'.format(run),data_format='IDS',shot=shot,run=run)
    temperature_i=Temperature(ID='ITER_NEMO Benchmark shot=44 run={} - Ti'.format(run),data_format='IDS',shot=shot,run=run,species='ions')
    temperature_e=Temperature(ID='ITER_NEMO Benchmark shot=44 run={} - Te'.format(run),data_format='IDS',shot=shot,run=run,species='electrons')
    density_i=Number_Density(ID='ITER_NEMO Benchmark shot=44 run={} - ni'.format(run),data_format='IDS',shot=shot,run=run,species='ions')
    density_e=Number_Density(ID='ITER_NEMO Benchmark shot=44 run={} - ne'.format(run),data_format='IDS',shot=shot,run=run,species='electrons')

    #######################################################################################################################################
    #get stuff ready to input to LOCUST by calculating poloidal flux profiles and dummy/missing data

    for locust_input in temperature_i,temperature_e,density_i,density_e:
        #calculate toroidal flux using the B_t on-axis from the equilibrium IDS, since we're only given toroidal flux coordinate 
        locust_input.set(flux_tor=equilibrium['bcentr']*(locust_input['flux_tor_coord']**2)/2.) #/2. to convert to Wb/rad (pi cancels)
        #calculate missing poloidal fluxes for LOCUST format
        locust_input['flux_pol']=prin.QTP_calc(T=locust_input['flux_tor'],Q=locust_input['q']) #gives flux_pol in Wb/rad
        locust_input['flux_pol_norm']=(locust_input['flux_pol']-locust_input['flux_pol'][0])/(equilibrium['sibry']-locust_input['flux_pol'][0]) #sibry already in Wb/rad

    ##################################################################

    #calculate missing equilibrium data
    equilibrium.set(xdum=np.array(0))
    equilibrium.set(limitr=np.array(0))
    rmag,zmag=prin.mag_axis_calc(equilibrium)
    equilibrium.set(zcentr=zmag)

    ##################################################################

    #rename the fields in the beam deposition keys in the particles_ objects to LOCUST_IO-friendly names
    particles['V_R']=particles.data.pop('V_R')
    particles['V_tor']=particles.data.pop('V_PHI')
    particles['V_Z']=particles.data.pop('V_Z')
    particles['R']=particles.data.pop('R')
    particles['phi']=particles.data.pop('phi')
    particles['Z']=particles.data.pop('Z')
    particles['E']=particles.data.pop('Energy')
    particles['V_pitch']=particles.data.pop('Pitch angle')
    particles['X']=particles['R']*np.cos(particles['phi'])
    particles['Y']=particles['R']*np.sin(particles['phi'])

    #######################################################################################################################################
    #dump all the data to the LOCUST formats



    #currently normalises poloidal flux coordinates for the 1D profiles
    
    equilibrium.dump_data(data_format='GEQDSK',filename='LOCUST_NEMO_{}.eqdsk'.format(run))
    particles.dump_data(data_format='LOCUST_FO_weighted',filename='LOCUST_NEMO_ptcles_{}.dat'.format(run))
    temperature_i.dump_data(data_format='LOCUST',filename='LOCUST_NEMO_profile_Ti_{}.dat'.format(run))
    temperature_e.dump_data(data_format='LOCUST',filename='LOCUST_NEMO_profile_Te_{}.dat'.format(run))
    density_i.dump_data(data_format='LOCUST',filename='LOCUST_NEMO_profile_ni_{}.dat'.format(run))
    density_e.dump_data(data_format='LOCUST',filename='LOCUST_NEMO_profile_ne_{}.dat'.format(run))

    ###################################################################################################
    #plot some of the data

    #particles.plot(some_equilibrium=equilibrium,style='histogram',number_bins=50,axes=['R'])
    #particles.plot(some_equilibrium=equilibrium,style='histogram',number_bins=50,axes=['E','V_pitch'],LCFS=True,real_scale=True)
    #equilibrium.plot()
    #temperature_i.plot()
    #density_i.plot()

#################################

##################################################################

###################################################################################################
