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

import sys

try:
    import scipy
    import numpy as np
    import matplotlib
    from matplotlib import cm
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d #import 3D plotting axes
    from mpl_toolkits.mplot3d import Axes3D
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from constants import *
except:
    raise ImportError("ERROR: LOCUST_IO/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from settings import *
except:
    raise ImportError("ERROR: LOCUST_IO/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)



##################################################################
#Main Code

def plot_number_density(some_number_density,axis='flux_pol_norm',ax=False,fig=False):
    """
    plots number density
     
    notes:
        axis - selects x axis of plot
        ax - take input axes (can be used to stack plots)
        fig - take input fig (can be used to add colourbars etc)
    """

    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if fig is False:
        fig_flag=False
    else:
        fig_flag=True

    if fig_flag is False:
        fig = plt.figure() #if user has not externally supplied figure, generate
    
    if ax_flag is False: #if user has not externally supplied axes, generate them
        ax = fig.add_subplot(111)

    ax.plot(some_number_density[axis],some_number_density['n'])
    ax.set_xlabel(axis)
    ax.set_ylabel('number density [m^-3]')
    ax.set_title(some_number_density.ID)
    
    if ax_flag is False and fig_flag is False:
        plt.show()

def plot_temperature(some_temperature,axis='flux_pol_norm',ax=False,fig=False):
    """
    plots number density

    notes:
        axis - selects x axis of plot
        ax - take input axes (can be used to stack plots)
        fig - take input fig (can be used to add colourbars etc)
    """

    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if fig is False:
        fig_flag=False
    else:
        fig_flag=True

    if fig_flag is False:
        fig = plt.figure() #if user has not externally supplied figure, generate
    
    if ax_flag is False: #if user has not externally supplied axes, generate them
        ax = fig.add_subplot(111)
   
    ax.plot(some_temperature[axis],some_temperature['T'])
    ax.set_xlabel(axis)
    ax.set_ylabel('temperature [eV]')
    ax.set_title(some_temperature.ID)

    if ax_flag is False and fig_flag is False:
        plt.show()

def plot_beam_deposition(some_beam_depo,some_equilibrium=False,grid=False,style='histogram',weight=True,number_bins=50,axes=['R','Z'],LCFS=False,limiters=False,real_scale=False,colmap=cmap_default,ax=False,fig=False):
    """
    plots beam deposition

    notes:
        some_equilibrium - corresponding equilibrium for plotting plasma boundary, scaled axes etc.
        grid - grid-like object containing same 'axes' to bin against e.g. distribution_function object with ['R'] and ['Z'] data
        style - choose from scatter or histogram
        weight - toggle whether to include marker weights in histograms
        number_bins - set bin for histogram
        axes - list of strings specifying which axes should be plotted
        LCFS - toggles whether plasma boundary is included (requires equilibrium arguement)
        limiters - toggles limiters on/off in 2D plots
        real_scale - sets r,z scale to real tokamak cross section
        colmap - set the colour map (use get_cmap names)
        ax - take input axes (can be used to stack plots)
        fig - take input fig (can be used to add colourbars etc)
    """

    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if fig is False:
        fig_flag=False
    else:
        fig_flag=True

    if fig_flag is False:
        fig = plt.figure() #if user has not externally supplied figure, generate
    
    if ax_flag is False: #if user has not externally supplied axes, generate them
        ax = fig.add_subplot(111)

    ndim=len(axes) #infer how many dimensions user wants to plot
    if ndim==1: #plot 1D histograms
        if weight:
            some_beam_depo_binned,some_beam_depo_binned_edges=np.histogram(some_beam_depo[axes[0]],bins=number_bins,weights=some_beam_depo['weight'])
        else:
            some_beam_depo_binned,some_beam_depo_binned_edges=np.histogram(some_beam_depo[axes[0]],bins=number_bins)
        some_beam_depo_binned_centres=(some_beam_depo_binned_edges[:-1]+some_beam_depo_binned_edges[1:])*0.5
        ax.plot(some_beam_depo_binned_centres,some_beam_depo_binned)
        ax.set_xlabel(axes[0])
        ax.set_title(some_beam_depo.ID)

    elif ndim==2: #plot 2D histograms

        if axes==['R','Z']: #check for commonly-used axes
            if real_scale is True: #set x and y plot limits to real scales
                if some_equilibrium:
                    ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                    ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))
                ax.set_aspect('equal')
            else:
                ax.set_aspect('auto')

        elif axes==['X','Y']:          
            if real_scale is True: 
                if some_equilibrium:
                    ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                    ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_aspect('equal')
            else:
                ax.set_aspect('auto')

        if style=='histogram':
            if grid is not False: #bin according to pre-defined grid
                if weight:
                    some_beam_depo_binned,some_beam_depo_binned_x,some_beam_depo_binned_y=np.histogram2d(some_beam_depo[axes[0]],some_beam_depo[axes[1]],bins=[grid[axes[0]],grid[axes[1]]],weights=some_beam_depo['weight'])
                else:
                    some_beam_depo_binned,some_beam_depo_binned_x,some_beam_depo_binned_y=np.histogram2d(some_beam_depo[axes[0]],some_beam_depo[axes[1]],bins=[grid[axes[0]],grid[axes[1]]])
            else:
                if weight:
                    some_beam_depo_binned,some_beam_depo_binned_x,some_beam_depo_binned_y=np.histogram2d(some_beam_depo[axes[0]],some_beam_depo[axes[1]],bins=number_bins,weights=some_beam_depo['weight'])
                else:
                    some_beam_depo_binned,some_beam_depo_binned_x,some_beam_depo_binned_y=np.histogram2d(some_beam_depo[axes[0]],some_beam_depo[axes[1]],bins=number_bins)
            #some_beam_depo_binned_x and some_beam_depo_binned_x are first edges then converted to centres
            some_beam_depo_binned_x=(some_beam_depo_binned_x[:-1]+some_beam_depo_binned_x[1:])*0.5
            some_beam_depo_binned_y=(some_beam_depo_binned_y[:-1]+some_beam_depo_binned_y[1:])*0.5
            some_beam_depo_binned_y,some_beam_depo_binned_x=np.meshgrid(some_beam_depo_binned_y,some_beam_depo_binned_x)
            
            ax.set_facecolor(colmap(np.amin(some_beam_depo_binned)))
            mesh=ax.pcolormesh(some_beam_depo_binned_x,some_beam_depo_binned_y,some_beam_depo_binned,cmap=colmap,vmin=np.amin(some_beam_depo_binned),vmax=np.amax(some_beam_depo_binned))
            #mesh=ax.contourf(some_beam_depo_binned_x,some_beam_depo_binned_y,some_beam_depo_binned,levels=np.linspace(np.amin(some_beam_depo_binned),np.amax(some_beam_depo_binned),num=20),cmap=colmap,edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(some_beam_depo_binned),vmax=np.amax(some_beam_depo_binned))
            if fig_flag is False:    
                fig.colorbar(mesh,ax=ax,orientation='horizontal')

        elif style=='scatter':
            ax.scatter(some_beam_depo[axes[0]],some_beam_depo[axes[1]],color='red',marker='x',s=1)

        if axes==['R','Z']:
            if real_scale is True: #set x and y plot limits to real scales
                if some_equilibrium:
                    ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                    ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))
                ax.set_aspect('equal')
            else:
                ax.set_aspect('auto')
            if LCFS is True: #plot plasma boundary
                ax.plot(some_equilibrium['lcfs_r'],some_equilibrium['lcfs_z'],plot_style_LCFS) 
            if limiters is True: #add boundaries if desired
                ax.plot(some_equilibrium['rlim'],some_equilibrium['zlim'],plot_style_limiters)

        elif axes==['X','Y']:
            if real_scale is True: #set x and y plot limits to real scales
                if some_equilibrium:
                    ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                    ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_aspect('equal')
            else:
                ax.set_aspect('auto')
            if LCFS is True: #plot plasma boundary
                plasma_max_R=np.max(some_equilibrium['lcfs_r'])
                plasma_min_R=np.min(some_equilibrium['lcfs_r'])
                ax.plot(plasma_max_R*np.cos(np.linspace(0,2.0*pi,100)),plasma_max_R*np.sin(np.linspace(0.0,2.0*pi,100)),plot_style_LCFS)
                ax.plot(plasma_min_R*np.cos(np.linspace(0,2.0*pi,100)),plasma_min_R*np.sin(np.linspace(0.0,2.0*pi,100)),plot_style_LCFS)          
            if limiters is True: #add boundaries if desired
                ax.plot(some_equilibrium['rlim'],some_equilibrium['zlim'],plot_style_limiters)    
        
        if ax_flag is True or fig_flag is True: #return the plot object
            if 'mesh' in locals():
                return mesh

        ax.set_xlabel(axes[0])
        ax.set_ylabel(axes[1])
        ax.set_title(some_beam_depo.ID)
       
    elif ndim==3: #plot 3D scatter - assume X,Y,Z

        if style!='scatter':
            print("ERROR: plot_beam_deposition() can only plot scatter style in 3D!")
            return

        if ax_flag is False and len(axes)==3:
            ax = fig.gca(projection='3d')
        
        if LCFS: #plot periodic poloidal cross-sections in 3D
            for angle in np.linspace(0.0,2.0*pi,4,endpoint=False):
                x_points=some_equilibrium['lcfs_r']*np.cos(angle)
                y_points=some_equilibrium['lcfs_r']*np.sin(angle)
                z_points=some_equilibrium['lcfs_z']
                ax.plot(x_points,y_points,zs=z_points,color=plot_style_LCFS)

        if limiters: #plot periodic poloidal cross-sections in 3D
            for angle in np.linspace(0.0,2.0*pi,4,endpoint=False):
                x_points=some_equilibrium['rlim']*np.cos(angle)
                y_points=some_equilibrium['rlim']*np.sin(angle)
                z_points=some_equilibrium['zlim']
                ax.plot(x_points,y_points,zs=z_points,color=plot_style_limiters)

        if real_scale is True:
            ax.set_aspect('equal')
            if some_equilibrium:
                ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_zlim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D'])) 

        ax.scatter(some_beam_depo[axes[0]],some_beam_depo[axes[1]],some_beam_depo[axes[2]],color=colmap(np.random.uniform()),s=0.1)
    
    if ax_flag is False and fig_flag is False:
        plt.show()

def plot_equilibrium(some_equilibrium,key='psirz',LCFS=False,limiters=False,number_contours=20,contour_fill=True,colmap=cmap_default,ax=False,fig=False):
    """
    plots equilibrium
    
    notes:
        
    args:
        key - selects which data in equilibrium to plot
        LCFS - toggles plasma boundary on/off in 2D plots (requires equilibrium arguement)
        limiters - toggles limiters on/off in 2D plots
        number_contours - set fidelity of 2D contour plot
        contour_fill - toggle contour fill on 2D plots
        colmap - set the colour map (use get_cmap names)
        ax - take input axes (can be used to stack plots)
        fig - take input fig (can be used to add colourbars etc)
    """

    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if fig is False:
        fig_flag=False
    else:
        fig_flag=True

    #0D data
    if some_equilibrium[key].ndim==0:
        print(some_equilibrium[key])
        return
    
    #>0D data is plottable
    if fig_flag is False:
        fig = plt.figure() #if user has not externally supplied figure, generate
    
    if ax_flag is False: #if user has not externally supplied axes, generate them
        ax = fig.add_subplot(111)

    #1D data
    if some_equilibrium[key].ndim==1:
        ax.plot(some_equilibrium[key])
        ax.set_ylabel(key)

    #2D data
    elif some_equilibrium[key].ndim==2:

        X=some_equilibrium['R_1D'] #make a mesh
        Y=some_equilibrium['Z_1D'] 
        Y,X=np.meshgrid(Y,X) #swap since things are defined r,z 
        Z=some_equilibrium[key] #2D array (nR_1D,nZ_1D) of poloidal flux
        
        #2D plot
        if contour_fill is True:
            mesh=ax.contourf(X,Y,Z,levels=np.linspace(np.amin(Z),np.amax(Z),num=number_contours),cmap=colmap,edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(Z),vmax=np.amax(Z))
            for c in mesh.collections: #for use in contourf
                c.set_edgecolor("face")
        else:
            mesh=ax.contour(X,Y,Z,levels=np.linspace(np.amin(Z),np.amax(Z),num=number_contours),cmap=colmap,edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(Z),vmax=np.amax(Z))
        #mesh=ax.pcolormesh(X,Y,Z,cmap=colmap,edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(Z),vmax=np.amax(Z))

        #3D plot
        #ax=ax.axes(projection='3d')
        #ax.view_init(elev=90, azim=None) #rotate the camera
        #ax.plot_surface(X,Y,Z,rstride=1,cstride=1,cmap=colmap,edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(Z),vmax=np.amax(Z))
        
        if fig_flag is False:    
            fig.colorbar(mesh,ax=ax,orientation='horizontal')
        ax.set_aspect('equal')
        ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
        ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))
        ax.set_xlabel('R [m]')
        ax.set_ylabel('Z [m]')
        ax.set_title(some_equilibrium.ID)

        if LCFS is True:
            ax.plot(some_equilibrium['lcfs_r'],some_equilibrium['lcfs_z'],plot_style_LCFS) 
        if limiters is True: #add boundaries if desired
            ax.plot(some_equilibrium['rlim'],some_equilibrium['zlim'],plot_style_limiters) 

        if ax_flag is True or fig_flag is True: #return the plot object
            return mesh

    if ax_flag is False and fig_flag is False:
        plt.show()


def plot_B_field_line(some_equilibrium,axes=['X','Y','Z'],LCFS=True,limiters=False,number_field_lines=1,angle=2.0*pi,plot_full=False,start_mark=False,colmap=cmap_default,ax=False,fig=False):
    """
    plots random field lines for 'angle' radians around the tokamak

    notes:
        essentially uses the Euler method of integration
    args:
        axes - list of strings specifying which axes should be plotted
        LCFS - show plasma boundary outline (requires equilibrium arguement)
        limiters - toggles limiters on/off in 2D plots
        number_field_lines - the number of field lines to plot
        angle - plot field line for this many radians around central column
        plot_full - choose whether each field line will trace the whole plasma topology (see below also)
        start_mark - add marker for start of the field line
        colmap - set the colour map (use get_cmap names)
        ax - take input axes (can be used to stack plots)
        fig - take input fig (can be used to add colourbars etc)
    """
    
    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if fig is False:
        fig_flag=False
    else:
        fig_flag=True

    if 'B_field' not in some_equilibrium.data: #check we have a B field first
        print("ERROR: 'B_field' missing in equilibrium object (calculate first with B_calc)")
        return

    if fig_flag is False:
        fig = plt.figure() #if user has not externally supplied figure, generate
    
    if ax_flag is False: #if user has not externally supplied axes, generate them
        ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    
    dr=np.abs(some_equilibrium['R_1D'][1]-some_equilibrium['R_1D'][0])
    dz=np.abs(some_equilibrium['Z_1D'][1]-some_equilibrium['Z_1D'][0])
    
    if plot_full is True: #set integration path step size
        #if this option, then dl chosen to cause slight numerical drift that such that field line strays outward onto
        #different flux surfaces - thus tracing the whole field topology
        dl=3.0*np.sqrt(dr**2+dz**2) 
        angle=pi*200.0
        number_field_lines=1
        ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
        ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D'])) 
    else:
        #if this option chosen, then dl is reduced to stop numerical drift - thus a single flux surface is plotted
        dl=0.005*np.sqrt(dr**2+dz**2)
        if ax_flag is False and len(axes)==3:
            ax = fig.gca(projection='3d')
        
    print('plot_B_field_line - generating B field interpolators')
    B_field_R_interpolator=processing.utils.interpolate_2D(some_equilibrium['R_1D'],some_equilibrium['Z_1D'],some_equilibrium['B_field'][:,:,0])
    B_field_Z_interpolator=processing.utils.interpolate_2D(some_equilibrium['R_1D'],some_equilibrium['Z_1D'],some_equilibrium['B_field'][:,:,2])
    B_field_tor_interpolator=processing.utils.interpolate_2D(some_equilibrium['R_1D'],some_equilibrium['Z_1D'],some_equilibrium['B_field'][:,:,1])
    print('plot_B_field_line - finished generating B field interpolators')


    for line in range(number_field_lines): 

        R_points=np.array([]) #reset the arrays that will hold marker coordinates along a single field line
        Z_points=np.array([])
        tor_points=np.array([])

        R_point=float(some_equilibrium['rmaxis']) #pick some starting points                                                       
        tor_point=np.random.uniform(0.0,2.0*pi*R_point) 
        if plot_full is True: 
            Z_point=float(1.05*some_equilibrium['zmaxis'])
        else:
            Z_point=np.random.uniform(1.05*some_equilibrium['zmaxis'],0.9*np.max(some_equilibrium['lcfs_z']))      
            
        R_points=np.append(R_points,R_point) #add this position to our array of points along trajectory
        Z_points=np.append(Z_points,Z_point)
        tor_points=np.append(tor_points,tor_point)     

        tor_point_start=tor_point #remember where we started toroidally

        while np.abs(tor_point-tor_point_start)<some_equilibrium['rmaxis']*angle: #keep going until we rotate 'angle' radians around the tokamak
            
            B_field_R=B_field_R_interpolator(R_point,Z_point)
            B_field_tor=B_field_tor_interpolator(R_point,Z_point)
            B_field_Z=B_field_Z_interpolator(R_point,Z_point)

            #could save computing power by just dividing by largest B component*some_constant - since we only need to divide the B vector by something we know is going to be larger than than the largest component
            #and remember (in 3D) the magnitude cannot ever be more than sqrt(3) times larger than the largest component of the field, so if we divide by at least sqrt(3)*np.max(field) then we always know our normalised values will be < 1.0
            B_field_mag=np.sqrt(B_field_R**2+B_field_Z**2+B_field_tor**2)
            
            B_field_R/=B_field_mag #normalise the vector magnetic field
            B_field_Z/=B_field_mag
            B_field_tor/=B_field_mag

            Z_point+=B_field_Z*dl
            tor_point+=B_field_tor*dl
            R_point+=B_field_R*dl 

            R_points=np.append(R_points,R_point)
            Z_points=np.append(Z_points,Z_point)
            tor_points=np.append(tor_points,tor_point)    

        R_points=R_points[1::2] #take every other value to help the visuals
        Z_points=Z_points[1::2]
        tor_points=tor_points[1::2]
        X_points=R_points*np.cos(tor_points/R_points) #transform to cartesian
        Y_points=R_points*np.sin(tor_points/R_points) 

        if plot_full is True: #if wanting to trace the flux surfaces, then plot in r,z plane
            ax.plot(R_points,Z_points,color=colmap(np.random.uniform()))
            if start_mark: 
                ax.scatter(R_points[0],Z_points[0],color='red',s=10)
        else:
            
            if axes==['R','Z']: #poloidal plot
                ax.plot(R_points,Z_points,color=colmap(np.random.uniform()))
                if start_mark: 
                    ax.scatter(R_points[0],Z_points[0],color='red',s=10)
                if LCFS is True: #plot plasma boundary
                    ax.plot(some_equilibrium['lcfs_r'],some_equilibrium['lcfs_z'],plot_style_LCFS) 
                if limiters is True: #add boundaries if desired
                    ax.plot(some_equilibrium['rlim'],some_equilibrium['zlim'],plot_style_limiters)       
        
                ax.set_xlabel('R [m]')
                ax.set_ylabel('Z [m]')
                ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))
                ax.set_title(some_equilibrium.ID)

            elif axes==['X','Y']: #top-down plot
                ax.plot(X_points,Y_points,color=colmap(np.random.uniform()))
                if start_mark: 
                    ax.scatter(X_points[0],Y_points[0],color='red',s=10)
                if LCFS is True: #plot plasma boundary
                    plasma_max_R=np.max(some_equilibrium['lcfs_r'])
                    plasma_min_R=np.min(some_equilibrium['lcfs_r'])
                    ax.plot(plasma_max_R*np.cos(np.linspace(0,2.0*pi,100)),plasma_max_R*np.sin(np.linspace(0.0,2.0*pi,100)),plot_style_LCFS)
                    ax.plot(plasma_min_R*np.cos(np.linspace(0,2.0*pi,100)),plasma_min_R*np.sin(np.linspace(0.0,2.0*pi,100)),plot_style_LCFS) 
                if limiters is True: #add boundaries if desired
                    ax.plot(some_equilibrium['rlim'],some_equilibrium['zlim'],plot_style_limiters)
                
                ax.set_xlabel('X [m]')
                ax.set_ylabel('Y [m]')
                ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_title(some_equilibrium.ID)
            
            else: #3D plot
                if start_mark: 
                    ax.scatter(X_points[0],Y_points[0],Z_points[0],color='red',s=10)
                if LCFS: #plot periodic poloidal cross-sections in 3D
                    for angle in np.linspace(0.0,2.0*pi,4,endpoint=False):
                        x_points=some_equilibrium['lcfs_r']*np.cos(angle)
                        y_points=some_equilibrium['lcfs_r']*np.sin(angle)
                        z_points=some_equilibrium['lcfs_z']
                        ax.plot(x_points,y_points,zs=z_points,color=plot_style_LCFS)
                if limiters: #plot periodic poloidal cross-sections in 3D
                    for angle in np.linspace(0.0,2.0*pi,4,endpoint=False):
                        x_points=some_equilibrium['rlim']*np.cos(angle)
                        y_points=some_equilibrium['rlim']*np.sin(angle)
                        z_points=some_equilibrium['zlim']
                        ax.plot(x_points,y_points,zs=z_points,color=plot_style_limiters)

                ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_zlim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D'])) 
                ax.plot(X_points,Y_points,zs=Z_points,color=colmap(np.random.uniform()))

    if ax_flag is False and fig_flag is False:
        plt.show()


def plot_B_field_stream(some_equilibrium,colmap=cmap_default,ax=False,fig=False):
    """
    stream plot of magnetic field in R,Z plane

    notes:
        take transpose due to streamplot index convention
        colmap - set the colour map (use get_cmap names)
        ax - take input axes (can be used to stack plots)
        fig - take input fig (can be used to add colourbars etc)
    """

    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if fig is False:
        fig_flag=False
    else:
        fig_flag=True

    if fig_flag is False:
        fig = plt.figure() #if user has not externally supplied figure, generate
    
    if ax_flag is False: #if user has not externally supplied axes, generate them
        ax = fig.add_subplot(111)
    ax.set_aspect('equal')
        
    if 'B_field' not in some_equilibrium.data: #check we have a B field first
        print("ERROR: B_field missing in equilibrium object (calculate first with B_calc)")
        return

    B_mag=np.sqrt(some_equilibrium['B_field'][:,:,0]**2+some_equilibrium['B_field'][:,:,2]**2) #calculate poloidal field magnitude
    strm = ax.streamplot(some_equilibrium['R_1D'],some_equilibrium['Z_1D'],some_equilibrium['B_field'][:,:,0].T,some_equilibrium['B_field'][:,:,2].T, color=B_mag.T, linewidth=1, cmap=colmap)

    if fig_flag is False:    
        fig.colorbar(strm.lines,ax=ax,orientation='horizontal')
    ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
    ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))

    ax.set_xlabel('R [m]')
    ax.set_ylabel('Z [m]')
    ax.set_title(some_equilibrium.ID)

    if ax_flag is False and fig_flag is False:
        plt.show()
    

#################################

##################################################################

###################################################################################################