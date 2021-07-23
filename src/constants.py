#constants.py

'''
Samuel Ward
18/10/2018
----
File which holds physical constants and other data
---
notes: 
    - should replace with scipy.constants
---
'''


##################################################################
#Preamble
import numpy


##################################################################


################################# fundamental numbers

pi=numpy.pi

################################# physical constants

epsilon_0=8.854187817e-12 #electrical permitivity of free space
c=2.99792458e+08 #speed of light in vacuum
h=6.62606896e-34 #planck constant
k_B=1.3806504e-23 #boltzmann constant
amu=1.66053904e-27 

################################# specific constants

mass_electron_amu=0.00054858
mass_e=9.10938215e-31
charge_e=1.60217662e-19

mass_neutron_amu=1.0013783
mass_neutron=mass_neutron_amu*amu
charge_neutron=0.

mass_deuteron_amu=2.013553212
mass_deuteron=mass_deuteron_amu*amu
charge_deuterium=charge_e

mass_triton_amu=2.9937181
mass_triton=mass_triton_amu*amu
charge_triton=charge_e

mass_helium3_amu=2.9931529
mass_helium3=mass_helium3_amu*amu
charge_helium3=2.*charge_e

################################# defaults

species_mass=mass_deuteron #set defaults
species_mass_amu=mass_deuteron_amu
species_charge=charge_e

#################################

##################################################################

###################################################################################################
