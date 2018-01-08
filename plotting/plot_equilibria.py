#plot_equilibria.py

'''
Samuel Ward
15/12/2017
----
Plotting routines for LOCUST equilibria data
---
notes:
	NOTE here can define all the plotting member functions for the input/output classes	
	https://stackoverflow.com/questions/9455111/python-define-method-outside-of-class-definition
---
'''

'''



#GEQDSK variable names:


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
