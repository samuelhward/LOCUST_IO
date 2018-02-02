#plot_input.py

'''
Samuel Ward
15/12/2017
----
plotting routines for LOCUST input data
---
notes:
---
'''

##################################################################
#Preamble
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d #import 3D plotting axes
from matplotlib import cm #get colourmaps



##################################################################
#Main Code





def plot_number_density(some_number_density):
    """
    simple number density plot
     
    notes:
    """

    plt.plot(some_number_density['psi'],some_number_density['n'])
    plt.show()







def plot_equilibrium(some_equilibrium,key=None,boundary=None):
    """
    simple equilibrium plot 
     
    notes:
    """

    #default plot mode (plot psirz grid)
    if key==None:

        fig=plt.figure()

        X=np.arange(some_equilibrium['nw']) #make a mesh
        Y=np.arange(some_equilibrium['nh']) 
        X,Y=np.meshgrid(X,Y) #NOTE change this to R_mesh Z_mesh after issue 23 is fixed
        Z=some_equilibrium['psirz'].T 
        #plt.contour(X,Y,Z,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(Z),vmax=1.01*np.amax(Z))
        plt.pcolormesh(X,Y,Z,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(Z),vmax=1.01*np.amax(Z))

        plt.colorbar()

    else:
        
        #0D data
        if some_equilibrium[key].ndim==0:
            print(some_equilibrium[key])
        
        #1D data
        elif some_equilibrium[key].ndim==1:
            plt.plot(some_equilibrium[key])

        #2D data
        elif some_equilibrium[key].ndim==2:

            fig=plt.figure()

            X=np.arange(some_equilibrium['nw']) #make a mesh
            Y=np.arange(some_equilibrium['nh']) 
            X,Y=np.meshgrid(X,Y) #NOTE change this to R_mesh Z_mesh after issue 23 is fixed
            Z=some_equilibrium['psirz'].T #2D array (nw,nh) of poloidal flux
            
            #2D plot
            #plt.contour(X,Y,Z,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(Z),vmax=1.01*np.amax(Z))
            plt.pcolormesh(X,Y,Z,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(Z),vmax=1.01*np.amax(Z))

            #3D plot
            #ax=plt.axes(projection='3d')
            #ax.view_init(elev=90, azim=None) #rotate the camera
            #ax.plot_surface(X,Y,Z,rstride=1,cstride=1,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(Z),vmax=1.01*np.amax(Z))
            
            plt.colorbar()


    #add boundaries if requested
    if boundary=='limiters':
        pass
        #plt.plot(some_equilibrium['rlim'],some_equilibrium['zlim']) #NOTE waiting for issue 23 fix

    elif boundary=='plasma':
        pass
        #plt.plot(some_equilibrium['rbbbs'],some_equilibrium['zbbbs']) #NOTE waiting for issue 23 fix


    plt.show()




#################################

##################################################################

###################################################################################################