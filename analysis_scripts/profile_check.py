#profile_check.py

'''
Samuel Ward
16/05/2019
----
Plot kinetic profiles from shot 157418
---
notes: 
---
'''


##################################################################
#Preamble

import context
from classes.input_classes.number_density import Number_Density as ND 
from classes.input_classes.temperature import Temperature as T

import matplotlib.pyplot as plt 

##################################################################

ND_e=ND(ID='',data_format='LOCUST',filename='profile_ne.dat',species='electrons')
T_e=T(ID='',data_format='LOCUST',filename='profile_Te.dat',species='electrons')
T_i=T(ID='',data_format='LOCUST',filename='profile_Ti.dat',species='ions')


plt.plot(T_e['flux_pol_norm'],T_i['T']/T_e['T'],'r')

plt.plot(T_e['flux_pol_norm'],T_i['T']/ND_e['n'],'b')

plt.plot(T_e['flux_pol_norm'],ND_e['n']/T_e['T']**(2/3),'b')

plt.legend(['Ti/Te','Ti/Ni','Ne/Te^2/3 - electron heating'])

plt.show()


#################################

##################################################################

###################################################################################################