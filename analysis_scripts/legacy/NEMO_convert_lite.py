#!/usr/bin/env python

"""
Samuel Ward
09/06/2018
----
run script for converting NEMO beam depositions to LOCUST input files

using run 44 shot 221

shots 225,335 and 555 were all runs I performed with NEMO after Mireille added new feature to dump to LOCUST coordinates
shot 221 is same as 225 (I think) but after I hacked NEMO to write out real toroidal angle of marker
---
notes:
    ---
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
from processing import plot_input as plin
from processing import process_input as prin


##################################################################
#Main Code

pi=np.pi #define pi
au=1.660539e-27 #define atomic mass unit
e_charge=1.60217662e-19 #electron charge




#######################################################################################################################################
#read in all the data

#read in all the required data
equilibrium_221=Equilibrium(ID='ITER_NEMO Benchmark 221 - EQ',data_format='IDS',shot=44,run=221)
particles_221=Beam_Deposition(ID='ITER_NEMO Benchmark 221 - ptcles',data_format='IDS',shot=44,run=221)
temperature_221_i=Temperature(ID='ITER_NEMO Benchmark 221 - Ti',data_format='IDS',shot=44,run=221,species='ions')
temperature_221_e=Temperature(ID='ITER_NEMO Benchmark 221 - Te',data_format='IDS',shot=44,run=221,species='electrons')
density_221_i=Number_Density(ID='ITER_NEMO Benchmark 221 - ni',data_format='IDS',shot=44,run=221,species='ions')
density_221_e=Number_Density(ID='ITER_NEMO Benchmark 221 - ne',data_format='IDS',shot=44,run=221,species='electrons')


#######################################################################################################################################
#get stuff ready to input to LOCUST by calculating poloidal flux profiles and dummy/missing data

#calculate toroidal flux using the B_t on-axis from the equilibrium IDS, since we're only given toroidal flux coordinate
temperature_221_i['flux_tor']=equilibrium_221['bcentr']*(temperature_221_i['flux_tor_coord']**2)/2. #/2. to convert to Wb/rad
temperature_221_e['flux_tor']=equilibrium_221['bcentr']*(temperature_221_e['flux_tor_coord']**2)/2.
density_221_i['flux_tor']=equilibrium_221['bcentr']*(density_221_i['flux_tor_coord']**2)/2.
density_221_e['flux_tor']=equilibrium_221['bcentr']*(density_221_e['flux_tor_coord']**2)/2.


##################################################################

#calculate missing poloidal fluxes for LOCUST format

temperature_221_i['flux_pol']=prin.QTP_calc(T=temperature_221_i['flux_tor'],Q=temperature_221_i['q']) #gives flux_pol in Wb/rad
temperature_221_e['flux_pol']=prin.QTP_calc(T=temperature_221_e['flux_tor'],Q=temperature_221_e['q'])
density_221_i['flux_pol']=prin.QTP_calc(T=density_221_i['flux_tor'],Q=density_221_i['q'])
density_221_e['flux_pol']=prin.QTP_calc(T=density_221_e['flux_tor'],Q=density_221_e['q'])

temperature_221_i['flux_pol_norm']=temperature_221_i['flux_pol']/equilibrium_221['sibry'] #sibry already in Wb/rad
temperature_221_e['flux_pol_norm']=temperature_221_e['flux_pol']/equilibrium_221['sibry']
density_221_i['flux_pol_norm']=density_221_i['flux_pol']/equilibrium_221['sibry']
density_221_e['flux_pol_norm']=density_221_e['flux_pol']/equilibrium_221['sibry']

##################################################################

#calculate missing equilibrium data
equilibrium_221.set(xdum=np.array(0))
equilibrium_221.set(limitr=np.array(0))
#rmag,zmag=prin.mag_axis_calc(equilibrium_221)
#equilibrium_221.set(zcentr=zmag)
equilibrium_221.set(zcentr=np.array(0.6)) #XXX could probably interpolate this

##################################################################

#rename the fields in the beam deposition keys in the particles_ objects to LOCUST_IO-friendly names
#particles_221['V_R']=particles_221.data.pop('V_R')
#particles_221['V_tor']=particles_221.data.pop('V_PHI')
#particles_221['V_Z']=particles_221.data.pop('V_Z')
#particles_221['R']=particles_221.data.pop('R')
#particles_221['phi']=particles_221.data.pop('phi')
#particles_221['Z']=particles_221.data.pop('Z')


particles_221['X']=particles_221['R']*np.cos(particles_221['phi'])
particles_221['Y']=particles_221['R']*np.sin(particles_221['phi'])



#######################################################################################################################################
#dump all the data to the LOCUST formats


'''
#currently normalises poloidal flux coordinates for the 1D profiles
equilibrium_221.dump_data(data_format='GEQDSK',filename='LOCUST_NEMO_221.eqdsk')
#particles_221.dump_data(data_format='LOCUST',filename='LOCUST_NEMO_221_ptcles.dat')
temperature_221_i.dump_data(data_format='LOCUST',filename='LOCUST_NEMO_221_profile_Ti.dat')
temperature_221_e.dump_data(data_format='LOCUST',filename='LOCUST_NEMO_221_profile_Te.dat')
density_221_i.dump_data(data_format='LOCUST',filename='LOCUST_NEMO_221_profile_ni.dat')
density_221_e.dump_data(data_format='LOCUST',filename='LOCUST_NEMO_221_profile_ne.dat')
'''




###################################################################################################
#plot some of the data

plin.plot_beam_deposition(particles_221,some_equilibrium=equilibrium_221,type='histogram',number_bins=50,axes=['X','Y'],LCFS=True,real_scale=True)


#################################

##################################################################

###################################################################################################