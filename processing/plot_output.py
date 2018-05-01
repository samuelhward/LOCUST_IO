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

try:
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    import copy
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)
try:
    from processing import process_output
except:
    raise ImportError("ERROR: LOCUST_IO/processing/process_output.py could not be imported!\nreturning\n")
    sys.exit(1)

cmap_viridis=matplotlib.cm.get_cmap('viridis') #set colourmap
pi=np.pi #define pi

##################################################################
#Main Code


def plot_orbits(some_orbits,some_equilibrium=None,particles=[0],axes=['R','Z'],LCFS=False,real_scale=False,start_mark=False,ax=False):
    """
    simple orbits plot in the R,Z/X,Y planes
     
    notes:
        some_equilibrium - corresponding equilibrium for plotting plasma boundary, scaled axes etc.
        particles - iterable list of particle numbers
        axes - define plot axes
        LCFS - show plasma boundary outline
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
        
        if axes==['R','Z']: #if we're plotting along poloidal projection, then give options to include full cross-section and plasma boundary
           
            if real_scale is True:
                ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))
            
            if LCFS is True: #plot plasma boundary
                ax.plot(some_equilibrium['lcfs_r'],some_equilibrium['lcfs_z'],'m-') 

            for particle in particles: #plot all the particle tracks one by one
                ax.plot(some_orbits['orbits'][1::2,0,particle],some_orbits['orbits'][1::2,2,particle],color=cmap_viridis(np.random.uniform()),linewidth=0.5) #plot every other position along trajectory
                if start_mark: #show birth point
                    ax.plot(some_orbits['orbits'][0,0,particle],some_orbits['orbits'][0,2,particle],color='red',marker='o',markersize=1)

        elif axes==['X','Y']: #plotting top-down
            
            if real_scale is True:
                ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))

            if LCFS is True: #plot concentric rings to show inboard/outboard plasma boundaries
                plasma_max_R=np.max(some_equilibrium['lcfs_r'])
                plasma_min_R=np.min(some_equilibrium['lcfs_r'])
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
        ax.set_title(some_orbits.ID)


    elif ndim==3: #3D plotting

        if ax_flag is False:
            ax = fig.gca(projection='3d')

        if real_scale:
            ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
            ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
            ax.set_zlim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))

        if LCFS is True: #plot periodic poloidal cross-sections in 3D 
            for angle in np.linspace(0.0,2.0*pi,4,endpoint=False):
                x_points=some_equilibrium['lcfs_r']*np.cos(angle)
                y_points=some_equilibrium['lcfs_r']*np.sin(angle)
                z_points=some_equilibrium['lcfs_z']
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
        ax.set_title(some_orbits.ID)

    if ax_flag is False:
        plt.show()


def plot_final_particle_list(some_final_particle_list,some_equilibrium=None,some_dfn=None,type='histogram',number_bins=50,axes=['R','Z'],LCFS=False,real_scale=False,status_flags=['PFC_intercept'],weight=1.0,ax=False):
    """
    plot the final particle list
     
    notes:
        some_equilibrium - corresponding equilibrium for plotting plasma boundary, scaled axes etc.
        some_dfn - corresponding distribution_function for matching binning axes etc. 
        type - choose from scatter or histogram
        number_bins - set bin for histogram
        axes - define plot axes
        LCFS - show plasma boundary outline
        real_scale - plot to Tokamak scale (requires equilibrium arguement)
        status_flags - plot particles with these statuses
        weight - weight per marker (must be constant)
        ax - take input axes (can be used to stack plots)
    """

    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if ax_flag is False: #if user has not externally supplied axes, generate them #initialise plot
        fig = plt.figure() #initialise plot
        ax = fig.add_subplot(111)

    ndim=len(axes)
    if ndim==1: #plot 1D histograms

        for status in status_flags:
            p=np.where(some_final_particle_list['status_flag']==some_final_particle_list['status_flags'][status]) #find the particle indices which have the desired status_flag
            some_final_particle_list_binned,some_final_particle_list_binned_edges=np.histogram(some_final_particle_list[axes[0]][p],bins=number_bins)
            some_final_particle_list_binned_centres=(some_final_particle_list_binned_edges[:-1]+some_final_particle_list_binned_edges[1:])*0.5
            ax.plot(some_final_particle_list_binned_centres,weight*some_final_particle_list_binned)

            ax.set_xlabel(axes[0])
            ax.set_title(some_final_particle_list.ID)

    elif ndim==2: #plot 2D histograms

        if axes==['R','Z']: #check for common axes
            if real_scale is True: #set x and y plot limits to real scales
                ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))
                ax.set_aspect('equal')
            else:
                ax.set_aspect('auto')

        elif axes==['X','Y']:          
            if real_scale is True: 
                ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_aspect('equal')
            else:
                ax.set_aspect('auto')

        for status in status_flags: #XXX THIS MIGHT BE CAUSING THE BUG FOR PLOTTING MULTIPLE STATUS FLAGS, AS AXES COULD BE RESET BETWEEN EACH PLOT
            p=np.where(some_final_particle_list['status_flag']==some_final_particle_list['status_flags'][status]) #find the particle indices which have the desired status_flag
            
            if type=='histogram':
                if some_dfn is not None:
                    some_final_particle_list_binned,some_final_particle_list_binned_x,some_final_particle_list_binned_y=np.histogram2d(some_final_particle_list[axes[0]][p],some_final_particle_list[axes[1]][p],bins=[some_dfn['R'],some_dfn['Z']])
                else:  
                    some_final_particle_list_binned,some_final_particle_list_binned_x,some_final_particle_list_binned_y=np.histogram2d(some_final_particle_list[axes[0]][p],some_final_particle_list[axes[1]][p],bins=number_bins)
                #some_final_particle_list_binned_x and some_final_particle_list_binned_x are first edges then converted to centres
                some_final_particle_list_binned_x=(some_final_particle_list_binned_x[:-1]+some_final_particle_list_binned_x[1:])*0.5
                some_final_particle_list_binned_y=(some_final_particle_list_binned_y[:-1]+some_final_particle_list_binned_y[1:])*0.5
                some_final_particle_list_binned_y,some_final_particle_list_binned_x=np.meshgrid(some_final_particle_list_binned_y,some_final_particle_list_binned_x)
                ax=plt.axes(facecolor=cmap_viridis(np.amin(some_final_particle_list_binned)))
                mesh=ax.pcolormesh(some_final_particle_list_binned_x,some_final_particle_list_binned_y,weight*some_final_particle_list_binned,cmap='viridis',edgecolor='face',linewidth=0,antialiased=True,vmin=np.amin(weight*some_final_particle_list_binned),vmax=np.amax(weight*some_final_particle_list_binned))
                #ax.contourf(some_final_particle_list_binned_x,some_final_particle_list_binned_y,some_final_particle_list_binned,levels=np.linspace(np.amin(some_final_particle_list_binned),np.amax(some_final_particle_list_binned),num=20),cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(some_final_particle_list_binned),vmax=np.amax(some_final_particle_list_binned))
                plt.colorbar(mesh)

            elif type=='scatter':
                ax.scatter(some_final_particle_list[axes[0]][p],some_final_particle_list[axes[1]][p],color=cmap_viridis(np.random.uniform()),marker='x',s=1)

        if axes==['R','Z']:
            if real_scale is True: #set x and y plot limits to real scales
                ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))
                ax.set_aspect('equal')
            else:
                ax.set_aspect('auto')
            if LCFS is True: #plot plasma boundary
                ax.plot(some_equilibrium['lcfs_r'],some_equilibrium['lcfs_z'],'m-') 

        elif axes==['X','Y']:
            if real_scale is True: #set x and y plot limits to real scales
                ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_aspect('equal')
            else:
                ax.set_aspect('auto')
            if LCFS is True: #plot plasma boundary
                plasma_max_R=np.max(some_equilibrium['lcfs_r'])
                plasma_min_R=np.min(some_equilibrium['lcfs_r'])
                ax.plot(plasma_max_R*np.cos(np.linspace(0,2.0*pi,100)),plasma_max_R*np.sin(np.linspace(0.0,2.0*pi,100)),'m-')
                ax.plot(plasma_min_R*np.cos(np.linspace(0,2.0*pi,100)),plasma_min_R*np.sin(np.linspace(0.0,2.0*pi,100)),'m-')          
        
        ax.set_xlabel(axes[0])
        ax.set_ylabel(axes[1])
        ax.set_title(some_final_particle_list.ID)
            
    if ax_flag is False:
        plt.show()


def plot_distribution_function(some_distribution_function,some_equilibrium=None,key='dfn',axes=['R','Z'],LCFS=False,real_scale=False,integrate=['velocity'],ax=False):
    """
    plot the distribution function

    notes:
        some_equilibrium - corresponding equilibrium for plotting plasma boundary, scaled axes etc.
        key - select which data to plot
        axes - define plot axes
        LCFS - show plasma boundary outline
        real_scale - plot to Tokamak scale (requires equilibrium arguement)
        integrate - when plotting the actual distribution function, toggle to integrate over space/velocity to get to particles/cell instead of SI units
        ax - take input axes (can be used to stack plots)
        
        assumes that the dfn data has not been collapsed yet
    """


    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    #0D data
    if some_distribution_function[key].ndim==0:
        print(some_distribution_function[key])
        return
    
    #>0D data is plottable
    if ax_flag is False: #if user has not externally supplied axes, generate them
        fig = plt.figure() #initialise plot
        ax = fig.add_subplot(111)

    #1D data
    if some_distribution_function[key].ndim==1:
        ax.plot(some_distribution_function[key])
        ax.set_xlabel(key)
        ax.set_title(some_distribution_function.ID)

    #plot distribution function
    elif key=='dfn':

        if len(axes)>2:
            print("ERROR: plot_distribution_function given too many plot axes")
            return
        
        dfn_copy=copy.deepcopy(some_distribution_function.data) #deep copy the object data dictionary (remember dicitonaries act like classes with a __getitem__ method)
        if integrate:
            dfn_copy[key]=process_output.dfn_integrate(dfn_copy,integrate) #integrate dfn over velocity to get number of particles per m^3
        dfn_copy[key]=process_output.dfn_collapse(dfn_copy,coordinates=axes) #collapse dfn down to 2D

        if 'V' in axes: #convert velocity to keV if velocity if a chosen axis
            dfn_copy['V']=(0.5*2.0*1.67e-27*dfn_copy['V']**2)/(1000.*1.6e-19)

        if real_scale is True: #set x and y plot limits to real scales
            ax.set_xlim(np.min(dfn_copy[axes[0]]),np.max(dfn_copy[axes[0]]))
            ax.set_ylim(np.min(dfn_copy[axes[1]]),np.max(dfn_copy[axes[1]]))
            ax.set_aspect('equal')
        else:
            ax.set_aspect('auto')

        X=dfn_copy[axes[0]] #make a mesh
        Y=dfn_copy[axes[1]]
        Y,X=np.meshgrid(Y,X) #dfn is r,z so need to swap order here
        ax=plt.axes(facecolor=cmap_viridis(np.amin(dfn_copy[key])))
        mesh=ax.pcolormesh(X,Y,dfn_copy[key],cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(dfn_copy[key]),vmax=np.amax(dfn_copy[key]))
        #mesh=ax.contourf(X,Y,dfn_copy[key],levels=np.linspace(np.amin(dfn_copy[key]),np.amax(dfn_copy[key]),num=number_contours),cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(dfn_copy[key]),vmax=np.amax(dfn_copy[key]))
        '''for c in mesh.collections: #for use in contourf
            c.set_edgecolor("face")'''        
        plt.colorbar(mesh)
        ax.set_xlabel(axes[0])
        ax.set_ylabel(axes[1])
        ax.set_title(some_distribution_function.ID)
        
        if real_scale is True: #set x and y plot limits to real scales
            ax.set_xlim(np.min(dfn_copy[axes[0]]),np.max(dfn_copy[axes[0]]))
            ax.set_ylim(np.min(dfn_copy[axes[1]]),np.max(dfn_copy[axes[1]]))
            ax.set_aspect('equal')
        else:
            ax.set_aspect('auto')
        if LCFS is True: #plot plasma boundary
            ax.plot(some_equilibrium['lcfs_r'],some_equilibrium['lcfs_z'],'m-') 
        
    if ax_flag is False:
        plt.show()

#################################

##################################################################

###################################################################################################