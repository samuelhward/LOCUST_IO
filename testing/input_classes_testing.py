#!/usr/bin/env python

#input_classes_testing.py

'''
Samuel Ward
02/11/2017
----
Unit tests for the input_classes.py methods
---
notes:
---
'''



##################################################################
#Preamble

import sys

try:
	import context
except:
	raise ImportError("ERROR: local context.py not found - please create one and add it to this directory\nreturning\n")
	sys.exit(1)

try: 
	import input_classes
except:
	raise ImportError("ERROR: LOCUST_IO/classes/input_classes.py could not be imported!\nreturning\n")
	sys.exit(1)





##################################################################
#Main Code

#first test is to open an Equilibrium GEQDSK and test that reading in upon initialisation is OK and that all member data is set correctly
test1_ID='test1_ID'
test1_input_filename='test.eqdsk'
test1_data_format='GEQDSK'
test_1=input_classes.Equilibrium(test1_ID,test1_input_filename,test1_data_format)
print(test_1.input_filename)
print(test_1.data_format)
print(test_1.LOCUST_input_type)
print(test_1.ID)


#next test to see if we can initialise a blank Equilibrium before populating with GEQDSK formatted data using read_data
test2_ID='test2_ID'
test2_input_filename='test.eqdsk'
test2_data_format='GEQDSK'
test_2=input_classes.Equilibrium(test2_ID)
test_2.read_data(test2_input_filename,test2_data_format)
print(test_2.input_filename)
print(test_2.data_format)
print(test_2.LOCUST_input_type)
print(test_2.ID)


#next test to see if we can dump an Equilibrium to GEQDSK
test3_ID='test3_ID'
test3_input_filename='test.eqdsk'
test3_data_format='GEQDSK'
test3=input_classes.Equilibrium(test3_ID,test3_input_filename,test3_data_format)
test3.dump_data('test3_output.eqdsk','GEQDSK')




















#examine the data - uncomment relevant sections

'''
#0D data

print(test.data['nh'])		#number of points in R (x or width)
print(test.data['nw'])		#number of points in Z (y or height)
print(test.data['idum']) 	#number of spatial dimensions?
print(test.data['rdim'])	#size of the R dimension in m
print(test.data['zdim'])	#size of the Z dimension in m
print(test.data['rcentr'])	#reference value of R
print(test.data['bcentr'])	#vacuum toroidal magnetic field at rcentr
print(test.data['rleft'])	#R at left (inner) boundary
print(test.data['zmid'])	#Z at middle of domain
print(test.data['rmaxis'])	#R at magnetic axis (O-point)
print(test.data['zmaxis'])	#Z at magnetic axis (O-point)
print(test.data['simag'])	#poloidal flux psi at magnetic axis (Weber / rad)
print(test.data['sibry'])	#poloidal flux at plasma boundary (Weber / rad)
print(test.data['current'])	#plasma current [Amps]   
print(test.data['xdum'])	#dummy variable - just contains zero
print(test.data['nbbbs'])	#plasma boundary
print(test.data['limitr'])	#wall boundary

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
import matplotlib.pyplot as plt
from matplotlib import cm #get colourmaps
import numpy
from mpl_toolkits import mplot3d #import 3D plotting axes
fig=plt.figure()

X=numpy.arange(test.data['nh']) #make a mesh
Y=numpy.arange(test.data['nw'])
X,Y=numpy.meshgrid(X,Y)
Z=test.data['psirz'] #2D array (nx,ny) of poloidal flux
ax=plt.axes(projection='3d')
#ax.view_init(elev=90, azim=None) #rotate the camera
ax.plot_surface(X,Y,Z,rstride=1,cstride=1,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*numpy.amin(Z),vmax=1.01*numpy.amax(Z))
plt.show()
'''



#################################

##################################################################

###################################################################################################