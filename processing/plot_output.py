#plot_output.py

'''
Samuel Ward
25/1/2018
----
plotting routines for LOCUST output data
---
notes:
---
'''

##################################################################
#Preamble
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

cmap_viridis=matplotlib.cm.get_cmap('viridis') #set colourmap
pi=np.pi #define pi

##################################################################
#Main Code


def plot_orbits(some_orbits,some_equilibrium=None,particles=[0],axes=['r','z'],plasma=False,real_scale=False,start_mark=False,ax=False):
    """
    simple orbits plot in the R,Z plane
     
    notes:
        particles - iterable list of particle numbers
        axes - define plot axes
        plasma - show plasma boundary outline
        real_scale - plot to Tokamak scale (requires equilibrium arguement)
        start_mark - include marker to show birth point
        ax - take input axes (can be used to stack plots)
    """
    
    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True
    
    ndim=len(axes)
    if ax_flag is False: #if user has not externally supplied axes, generate them #initialise plot
        fig = plt.figure() #initialise plot
        ax = fig.add_subplot(111)
    ax.set_aspect('equal')

    if ndim==2: #2D plotting
        
        if axes==['r','z']: #if we're plotting along poloidal projection, then give options to include full cross-section and plasma boundary
           
            if real_scale is True:
                ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))
            
            if plasma is True: #plot plasma boundary
                ax.plot(some_equilibrium['rbbbs'],some_equilibrium['zbbbs'],'m-') 

            for particle in particles: #plot all the particle tracks one by one
                ax.plot(some_orbits['orbits'][1::2,0,particle],some_orbits['orbits'][1::2,2,particle],color=cmap_viridis(np.random.uniform()),linewidth=0.5) #plot every other position along trajectory
                if start_mark: #show birth point
                    ax.plot(some_orbits['orbits'][0,0,particle],some_orbits['orbits'][0,2,particle],color='red',marker='o',markersize=1)

        elif axes==['x','y']: #plotting top-down
            
            if real_scale is True:
                ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))

            if plasma: #plot concentric rings to show inboard/outboard plasma boundaries
                plasma_max_R=np.max(some_equilibrium['rbbbs'])
                plasma_min_R=np.min(some_equilibrium['rbbbs'])
                ax.plot(plasma_max_R*np.cos(np.linspace(0,2.0*pi,100)),plasma_max_R*np.sin(np.linspace(0.0,2.0*pi,100)),'m-')
                ax.plot(plasma_min_R*np.cos(np.linspace(0,2.0*pi,100)),plasma_min_R*np.sin(np.linspace(0.0,2.0*pi,100)),'m-')

            for particle in particles: #plot all the particle tracks one by one
                x_points=some_orbits['orbits'][1::2,0,particle]*np.cos(some_orbits['orbits'][1::2,1,particle]) #calculate using every other position along trajectory
                y_points=some_orbits['orbits'][1::2,0,particle]*np.sin(some_orbits['orbits'][1::2,1,particle])   
                ax.plot(x_points,y_points,color=cmap_viridis(np.random.uniform()),linewidth=1) 
                
                if start_mark: #show birth point
                    ax.plot(x_points[0],y_points[0],color='red',marker='o',markersize=1)

        ax.set_xlabel(axes[0])
        ax.set_ylabel(axes[1])


    elif ndim==3: #3D plotting

        if ax_flag is False:
            ax = fig.gca(projection='3d')

        if real_scale:
            ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
            ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
            ax.set_zlim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))

        if plasma: #plot periodic poloidal cross-sections in 3D 
            for angle in np.linspace(0.0,2.0*pi,4,endpoint=False):
                x_points=some_equilibrium['rbbbs']*np.cos(angle)
                y_points=some_equilibrium['rbbbs']*np.sin(angle)
                z_points=some_equilibrium['zbbbs']
                ax.plot(x_points,y_points,zs=z_points,color='m')

        for particle in particles: #plot all the particle tracks one by one
            x_points=some_orbits['orbits'][1::2,0,particle]*np.cos(some_orbits['orbits'][1::2,1,particle]) #calculate using every other position along trajectory
            y_points=some_orbits['orbits'][1::2,0,particle]*np.sin(some_orbits['orbits'][1::2,1,particle])   
            z_points=some_orbits['orbits'][1::2,2,particle]
            ax.plot(x_points,y_points,zs=z_points,color=cmap_viridis(np.random.uniform())) 
            
            if start_mark: #show birth point
                ax.scatter(x_points[0],y_points[0],z_points[0],color='red',marker='o',s=10)

        ax.set_xlabel(axes[0])
        ax.set_ylabel(axes[1])
        ax.set_zlabel(axes[2])

    if ax_flag is False:
        plt.show()


#################################

##################################################################

###################################################################################################