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
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm #get colourmaps

pi=np.pi #define pi

##################################################################
#Main Code



def plot_number_density(some_number_density,axis=None):
    """
    plots number density
     
    notes:
        axis - selects x axis of plot
    """

    if axis is None:
        plot_number_density(some_number_density=some_number_density,axis='flux_pol') #default plot
    else:
        plt.plot(some_number_density[axis],some_number_density['n'])

    plt.show()



def plot_temperature(some_temperature,axis=None):
    """
    plots number density

    notes:
        axis - selects x axis of plot
    """

    if axis is None:
        plot_temperature(some_temperature=some_temperature,axis='flux_pol') #default plot
    else:
        plt.plot(some_temperature[axis],some_temperature['T'])

    plt.show()



def plot_beam_deposition(some_beam_depo,ndim=1,number_bins=50,axes=None,plasma=False,plasma_r=None,plasma_z=None,real_scale=False,R_1D=None,Z_1D=None):
    """
    plots beam deposition

    notes:
        ndim - set number of dimensions to plot
        number_bins - set bin for histogram
        axes - list of strings specifying which axes should be plotted
        plasma - toggles whether plasma boundary is included
        real_scale - sets r,z scale to real tokamak cross section
    """
    
    if ndim==1:

        if axes is None:

            plot_beam_deposition(some_beam_depo=some_beam_depo,ndim=ndim,number_bins=number_bins,axes='r') #default plot

        else:
            some_beam_depo_binned,some_beam_depo_binned_edges=np.histogram(some_beam_depo[axes],bins=number_bins)
            some_beam_depo_binned_centres=(some_beam_depo_binned_edges[:-1]+some_beam_depo_binned_edges[1:])*0.5
            plt.plot(some_beam_depo_binned_centres,some_beam_depo_binned)

    elif ndim==2:

        if axes is None: 

            plot_beam_deposition(some_beam_depo=some_beam_depo,ndim=ndim,number_bins=number_bins,axes=['r','z']) #default plot

        else: #some_beam_depo_binned_x and some_beam_depo_binned_x are first edges then converted to centres
            some_beam_depo_binned,some_beam_depo_binned_x,some_beam_depo_binned_y=np.histogram2d(some_beam_depo[axes[0]],some_beam_depo[axes[1]],bins=number_bins)
            some_beam_depo_binned_x=(some_beam_depo_binned_x[:-1]+some_beam_depo_binned_x[1:])*0.5
            some_beam_depo_binned_y=(some_beam_depo_binned_y[:-1]+some_beam_depo_binned_y[1:])*0.5
            some_beam_depo_binned_x,some_beam_depo_binned_y=np.meshgrid(some_beam_depo_binned_x,some_beam_depo_binned_y)
            plt.pcolormesh(some_beam_depo_binned_x,some_beam_depo_binned_y,some_beam_depo_binned,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(some_beam_depo_binned),vmax=1.01*np.amax(some_beam_depo_binned))
            #plt.contourf(some_beam_depo_binned_x,some_beam_depo_binned_y,some_beam_depo_binned,levels=np.linspace(0.99*np.amin(some_beam_depo_binned),1.01*np.amax(some_beam_depo_binned),num=20),cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(some_beam_depo_binned),vmax=1.01*np.amax(some_beam_depo_binned))
        
        if real_scale is True: #set x and y plot limits to real scales 
            plt.axis('scaled')
            plt.xlim(np.min(R_1D),np.max(R_1D))
            plt.ylim(np.min(Z_1D),np.max(Z_1D))
            
        if plasma is True: #plot plasma boundary
            plt.plot(plasma_r,plasma_z,'m-') 

    plt.show()




def plot_equilibrium(some_equilibrium,key=None,boundary=None,number_contours=20):
    """
    plots equilibrium
     
    notes:
        key - selects which data in equilibrium to plot
        boundary - toggles plasma boundary on/off in 2D plots
        number_contours - set fidelity of 2D contour plot
    """

    #default plot mode (plot psirz grid)
    if key==None:

        plot_equilibrium(some_equilibrium,key='psirz',boundary=boundary,number_contours=number_contours)

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
            Z=some_equilibrium[key].T #2D array (nw,nh) of poloidal flux
            
            #2D plot
            contour=plt.contourf(X,Y,Z,levels=np.linspace(0.99*np.amin(Z),1.01*np.amax(Z),num=number_contours),cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(Z),vmax=1.01*np.amax(Z))
            #plt.pcolormesh(X,Y,Z,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(Z),vmax=1.01*np.amax(Z))
            for c in contour.collections:
                c.set_edgecolor("face")

            #3D plot
            #ax=plt.axes(projection='3d')
            #ax.view_init(elev=90, azim=None) #rotate the camera
            #ax.plot_surface(X,Y,Z,rstride=1,cstride=1,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=0.99*np.amin(Z),vmax=1.01*np.amax(Z))
            
            plt.colorbar()
            plt.axis('scaled')
            plt.xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
            plt.ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))

            if boundary is not None:
                #add boundaries if requested
                if 'limiters' in boundary:
                    plt.plot(some_equilibrium['rlim'],some_equilibrium['zlim'],'k-') 

                if 'plasma' in boundary:
                    plt.plot(some_equilibrium['rbbbs'],some_equilibrium['zbbbs'],'m-') 

    plt.show()


def plot_field_line(B_field,R_1D,Z_1D,mag_axis_R,mag_axis_Z,plasma_boundary_r,plasma_boundary_z,number_field_lines=1,number_points=100):
    """
    plots random field lines for pi radians around the tokamak

    notes:
        essentially uses the Euler method of integration
        number_field_lines - the number of field lines to plot
        number_points - the number of points to plot the field line over
    """

    dl=mag_axis_R*pi/number_points #set integration path step size

    R_points=np.array([])
    Z_points=np.array([])
    tor_points=np.array([])

    for line in range(number_field_lines): 

        R_point=np.random.uniform(mag_axis_R,0.5*max(plasma_boundary_r))    #pick a random starting point
        Z_point=np.random.uniform(mag_axis_Z,0.5*max(plasma_boundary_z))    #NOTE : more rigorous check to see if we are in plasma boundary is needed
        tor_point=np.random.uniform(0.0,2.0*pi*R_point)

        R_points=np.append(R_points,R_point) #add this position to our array of points along trajectory
        Z_points=np.append(Z_points,Z_point)
        tor_points=np.append(tor_points,tor_point)     

        while np.abs(tor_point)<pi: #keep going until we rotate half way around the tokamak
        
            B_field_R=interpolate_2D(R_point,Z_point,R_1D,Z_1D,B_field[:,:,0])  #calculate the B field at this point
            B_field_Z=interpolate_2D(R_point,Z_point,R_1D,Z_1D,B_field[:,:,2])
            B_field_tor=interpolate_2D(R_point,Z_point,R_1D,Z_1D,B_field[:,:,1])

            B_field_mag=np.sqrt(B_field_R**2+B_field_Z**2+B_field_tor**2)   #normalise the vector magnetic field
            B_field_R/=B_field_mag
            B_field_Z/=B_field_mag
            B_field_tor/=B_field_mag

            R_point+=B_field_R*dl
            Z_point+=B_field_Z*dl
            tor_point+=B_field_tor*dl

            R_points=np.append(R_points,R_point)
            Z_points=np.append(Z_points,Z_point)
            tor_points=np.append(tor_points,tor_point)           

    X_points=R_points*np.cos(tor_points)
    Y_points=R_points*np.sin(tor_points)

    X_points=X_points.reshape(number_field_lines,len(X_points)/number_field_lines)
    Y_points=Y_points.reshape(number_field_lines,len(Y_points)/number_field_lines)
    Z_points=Z_points.reshape(number_field_lines,len(Z_points)/number_field_lines)

    fig = plt.figure()
    ax = fig.gca(projection='3d') #plot the result
    for line in range(number_field_lines):
        ax.plot(X_points[line,:],Y_points[line,:],Z_points[line,:],color=cm.viridis(np.random.uniform()))
    plt.show()


def plot_B_field(B_field,R_1D,Z_1D):
    """
    stream plot of magnetic field in R,Z plane

    notes:
    """

    R_2D,Z_2D=np.meshgrid(R_1D,Z_1D) #set up domain
    B_mag=np.sqrt(B_field[:,:,0]**2+B_field[:,:,2]**2) #calculate poloidal field magnitude



#################################

##################################################################

###################################################################################################