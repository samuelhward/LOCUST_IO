#!/usr/bin/env python

#input_classes_testing.py

'''
Samuel Ward
02/11/2017
----
Unit tests for the input_classes.py methods
---

---
'''



##################################################################
#Preamble

import sys

try:
	from LOCUST_IO.classes import input_classes 
except:
	raise ImportError("ERROR: input_classes.py could not be imported!\nreturning\n")
	sys.exit(1)

try:
	import matplotlib.pyplot as plt
	import numpy
except:
	raise ImportError("ERROR: supporting modules not found!\nreturning\n")
	sys.exit(1)



##################################################################
#Main Code


#specify test GEQDSK data
ID='eq_ID_here'
input_filename="test.eqdsk"
data_format='GEQDSK'
test=input_classes.Equilibrium(ID,input_filename,data_format)

print(test.input_filename)
print(test.data_format)
print(test.input_type)
print(test.ID)

#plot the GEQDSK data

#0D data
'''
plt.plot(test.data['nh'])		#number of points in R (x or width)
plt.plot(test.data['nw'])		#number of points in Z (y or height)
plt.plot(test.data['idum']) 	#number of spatial dimensions?
plt.plot(test.data['rdim'])		#size of the R dimension in m
plt.plot(test.data['zdim'])		#size of the Z dimension in m
plt.plot(test.data['rcentr'])	#reference value of R
plt.plot(test.data['bcentr'])	#vacuum toroidal magnetic field at rcentr
plt.plot(test.data['rleft'])	#R at left (inner) boundary
plt.plot(test.data['zmid'])		#Z at middle of domain
plt.plot(test.data['rmaxis'])	#R at magnetic axis (O-point)
plt.plot(test.data['zmaxis'])	#Z at magnetic axis (O-point)
plt.plot(test.data['simag'])	#poloidal flux psi at magnetic axis (Weber / rad)
plt.plot(test.data['sibry'])	#poloidal flux at plasma boundary (Weber / rad)
plt.plot(test.data['current'])	#plasma current [Amps]   
plt.plot(test.data['xdum'])		#dummy variable - just contains zero
plt.plot(test.data['nbbbs'])	#plasma boundary
plt.plot(test.data['limitr'])	#wall boundary

#1D data
plt.plot(test.data['fpol'])		#poloidal current function on uniform flux grid (1D array of f(psi)=R*Bt  [meter-Tesla])
plt.plot(test.data['pres'])		#plasma pressure in nt/m^2 on uniform flux grid (1D array of p(psi) [Pascals])
plt.plot(test.data['ffprime'])	#workk1
plt.plot(test.data['pprime'])	#workk1
plt.plot(test.data['qpsi'])		#q values on uniform flux grid
plt.plot(test.data['rlim']) 	#r wall boundary
plt.plot(test.data['zlim'])		#z wall boundary
plt.plot(test.data['rbbbs'])	#r plasma boundary
plt.plot(test.data['zbbbs'])	#z plasma boundary

#2D data
plt.plot(test.data['psirz']) #2D array (nx,ny) of poloidal flux
'''

#plt.plot(test.data['psirz']) #2D array (nx,ny) of poloidal flux
#plt.show()




#################################

##################################################################

###################################################################################################