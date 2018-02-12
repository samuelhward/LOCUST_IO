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
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d #import 3D plotting axes
from matplotlib import cm #get colourmaps



##################################################################
#Main Code



def plot_number_density(some_number_density,axis=None):
    """
    plots number density
     
    notes:
    """

    if axis is None:
        plt.plot(some_number_density['flux_pol'],some_number_density['n'])
    else:
        plt.plot(some_number_density[axis],some_number_density['n'])

    plt.show()



def plot_temperature(some_temperature,axis=None):
    """
    plots number density

    notes:

    """

    if axis is None:
        plt.plot(some_temperature['flux_pol'],some_temperature['T'])
    else:
        plt.plot(some_temperature[axis],some_temperature['T'])

    plt.show()



def plot_beam_deposition(some_beam_depo,ndim=1,number_bins=50,axes=None):
    """
    plots beam deposition

    notes:
        ndim - set number of dimensions to plot
        number_bins - set bin for histogram
        axes - list of strings specifying which axes should be plotted
    """
    
    if ndim==1:
        if axes is None:
            some_beam_depo_binned,some_beam_depo_binned_edges=np.histogram(some_beam_depo['r'],bins=number_bins)
            some_beam_depo_binned_centres=(some_beam_depo_binned_edges[:-1]+some_beam_depo_binned_edges[1:])*0.5
            plt.plot(some_beam_depo_binned_centres,some_beam_depo_binned)
        else:
            some_beam_depo_binned,some_beam_depo_binned_edges=np.histogram(some_beam_depo[axes],bins=number_bins)
            some_beam_depo_binned_centres=(some_beam_depo_binned_edges[:-1]+some_beam_depo_binned_edges[1:])*0.5
            plt.plot(some_beam_depo_binned_centres,some_beam_depo_binned)

    elif ndim==2:

        if axes is None: #some_beam_depo_binned_x and some_beam_depo_binned_x are first edges then converted to centres
            some_beam_depo_binned,some_beam_depo_binned_x,some_beam_depo_binned_y=np.histogram2d(some_beam_depo['r'],some_beam_depo['z'],bins=number_bins) 
            some_beam_depo_binned_x=(some_beam_depo_binned_x[:-1]+some_beam_depo_binned_x[1:])*0.5
            some_beam_depo_binned_y=(some_beam_depo_binned_y[:-1]+some_beam_depo_binned_y[1:])*0.5
            some_beam_depo_binned_x,some_beam_depo_binned_y=np.meshgrid(some_beam_depo_binned_x,some_beam_depo_binned_y)
            plt.pcolormesh(some_beam_depo_binned_x,some_beam_depo_binned_y,some_beam_depo_binned,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(some_beam_depo_binned),vmax=1.01*np.amax(some_beam_depo_binned))
            #plt.contourf(some_beam_depo_binned_x,some_beam_depo_binned_y,some_beam_depo_binned,levels=np.linspace(0.99*np.amin(some_beam_depo_binned),1.01*np.amax(some_beam_depo_binned),num=20),cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(some_beam_depo_binned),vmax=1.01*np.amax(some_beam_depo_binned))

        else:
            some_beam_depo_binned,some_beam_depo_binned_x,some_beam_depo_binned_y=np.histogram2d(some_beam_depo[axes[0]],some_beam_depo[axes[1]],bins=number_bins)
            some_beam_depo_binned_x=(some_beam_depo_binned_x[:-1]+some_beam_depo_binned_x[1:])*0.5
            some_beam_depo_binned_y=(some_beam_depo_binned_y[:-1]+some_beam_depo_binned_y[1:])*0.5
            some_beam_depo_binned_x,some_beam_depo_binned_y=np.meshgrid(some_beam_depo_binned_x,some_beam_depo_binned_y)
            plt.pcolormesh(some_beam_depo_binned_x,some_beam_depo_binned_y,some_beam_depo_binned,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(some_beam_depo_binned),vmax=1.01*np.amax(some_beam_depo_binned))
            #plt.contourf(some_beam_depo_binned_x,some_beam_depo_binned_y,some_beam_depo_binned,levels=np.linspace(0.99*np.amin(some_beam_depo_binned),1.01*np.amax(some_beam_depo_binned),num=20),cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(some_beam_depo_binned),vmax=1.01*np.amax(some_beam_depo_binned))
    
    plt.show()




def plot_equilibrium(some_equilibrium,key=None,boundary=None):
    """
    plots equilibrium
     
    notes:
    """

    #default plot mode (plot psirz grid)
    if key==None:

        fig=plt.figure()

        X=some_equilibrium['R_1D'] #make a mesh
        Y=some_equilibrium['Z_1D'] 
        X,Y=np.meshgrid(X,Y)
        Z=some_equilibrium['psirz'].T 
        plt.contourf(X,Y,Z,levels=np.linspace(0.99*np.amin(Z),1.01*np.amax(Z),num=20),cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(Z),vmax=1.01*np.amax(Z))
        #plt.pcolormesh(X,Y,Z,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(Z),vmax=1.01*np.amax(Z))

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

            X=some_equilibrium['R_1D'] #make a mesh
            Y=some_equilibrium['Z_1D'] 
            X,Y=np.meshgrid(X,Y)
            Z=some_equilibrium['psirz'].T #2D array (nw,nh) of poloidal flux
            
            #2D plot
            plt.contourf(X,Y,Z,levels=np.linspace(0.99*np.amin(Z),1.01*np.amax(Z),num=20),cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(Z),vmax=1.01*np.amax(Z))
            #plt.pcolormesh(X,Y,Z,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(Z),vmax=1.01*np.amax(Z))

            #3D plot
            #ax=plt.axes(projection='3d')
            #ax.view_init(elev=90, azim=None) #rotate the camera
            #ax.plot_surface(X,Y,Z,rstride=1,cstride=1,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(Z),vmax=1.01*np.amax(Z))
            
            plt.colorbar()

    if boundary is not None:
        #add boundaries if requested
        if 'limiters'in boundary:
            plt.plot(some_equilibrium['rlim'],some_equilibrium['zlim'],'k-') 

        if 'plasma' in boundary:
            plt.plot(some_equilibrium['rbbbs'],some_equilibrium['zbbbs'],'m-') 


    plt.show()




#################################

##################################################################

###################################################################################################