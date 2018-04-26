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

try:
    import scipy
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib import cm #get colourmaps
    from mpl_toolkits import mplot3d #import 3D plotting axes
    from mpl_toolkits.mplot3d import Axes3D
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)
try:
    from processing import process_input
except:
    raise ImportError("ERROR: LOCUST_IO/processing/process_input.py could not be imported!\nreturning\n")
    sys.exit(1)

cmap_viridis=matplotlib.cm.get_cmap('viridis') #set colourmap
pi=np.pi #define pi

##################################################################
#Main Code



def plot_number_density(some_number_density,axis='flux_pol_norm',ax=False):
    """
    plots number density
     
    notes:
        axis - selects x axis of plot
        ax - take input axes (can be used to stack plots)
    """

    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if ax_flag is False: #if user has not externally supplied axes, generate them
        fig = plt.figure() #initialise plot
        ax = fig.add_subplot(111)
    ax.plot(some_number_density[axis],some_number_density['n'])
    
    if ax_flag is False:
        plt.show()


def plot_temperature(some_temperature,axis='flux_pol_norm',ax=False):
    """
    plots number density

    notes:
        axis - selects x axis of plot
        ax - take input axes (can be used to stack plots)
    """

    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if ax_flag is False: #if user has not externally supplied axes, generate them
        fig = plt.figure() #initialise plot
        ax = fig.add_subplot(111)
    ax.plot(some_temperature[axis],some_temperature['T'])

    if ax_flag is False:
        plt.show()



def plot_beam_deposition(some_beam_depo,some_equilibrium=None,type='histogram',number_bins=50,axes=['R','Z'],LCFS=False,real_scale=False,ax=False):
    """
    plots beam deposition

    notes:
        some_equilibrium - corresponding equilibrium for plotting plasma boundary, scaled axes etc.
        type - choose from scatter or histogram
        number_bins - set bin for histogram
        axes - list of strings specifying which axes should be plotted
        LCFS - toggles whether plasma boundary is included
        real_scale - sets r,z scale to real tokamak cross section
        ax - take input axes (can be used to stack plots)
    """

    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if ax_flag is False: #if user has not externally supplied axes, generate them
        fig = plt.figure() #initialise plot
        ax = fig.add_subplot(111)

    ndim=len(axes) #infer how many dimensions user wants to plot
    if ndim==1: #plot 1D histograms

        some_beam_depo_binned,some_beam_depo_binned_edges=np.histogram(some_beam_depo[axes[0]],bins=number_bins)
        some_beam_depo_binned_centres=(some_beam_depo_binned_edges[:-1]+some_beam_depo_binned_edges[1:])*0.5
        ax.plot(some_beam_depo_binned_centres,some_beam_depo_binned)
        ax.set_xlabel(axes[0])

    elif ndim==2: #plot 2D histograms

        if axes==['R','Z']: #check for common axes
            if real_scale is True: #set x and y plot limits to real scales
                ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))
                ax.set_aspect('equal')

        elif axes==['X','Y']:          
            if real_scale is True: 
                ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_aspect('equal')

        if type=='histogram':
            #some_beam_depo_binned_x and some_beam_depo_binned_x are first edges then converted to centres
            some_beam_depo_binned,some_beam_depo_binned_x,some_beam_depo_binned_y=np.histogram2d(some_beam_depo[axes[0]],some_beam_depo[axes[1]],bins=number_bins)
            some_beam_depo_binned_x=(some_beam_depo_binned_x[:-1]+some_beam_depo_binned_x[1:])*0.5
            some_beam_depo_binned_y=(some_beam_depo_binned_y[:-1]+some_beam_depo_binned_y[1:])*0.5
            some_beam_depo_binned_y,some_beam_depo_binned_x=np.meshgrid(some_beam_depo_binned_y,some_beam_depo_binned_x)
            axes=plt.axes(facecolor=cmap_viridis(np.amin(some_beam_depo_binned)))
            mesh=ax.pcolormesh(some_beam_depo_binned_x,some_beam_depo_binned_y,some_beam_depo_binned,cmap='viridis',edgecolor='face',linewidth=0,antialiased=True,vmin=np.amin(some_beam_depo_binned),vmax=np.amax(some_beam_depo_binned))
            #mesh=ax.contourf(some_beam_depo_binned_x,some_beam_depo_binned_y,some_beam_depo_binned,levels=np.linspace(np.amin(some_beam_depo_binned),np.amax(some_beam_depo_binned),num=20),cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(some_beam_depo_binned),vmax=np.amax(some_beam_depo_binned))
            plt.colorbar(mesh)

        elif type=='scatter':
            ax.scatter(some_beam_depo[axes[0]],some_beam_depo[axes[1]],color='red',marker='x',s=1)

        if axes==['R','Z']:
            if real_scale is True: #set x and y plot limits to real scales
                ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))
                ax.set_aspect('equal')
            if LCFS is True: #plot plasma boundary
                ax.plot(some_equilibrium['lcfs_r'],some_equilibrium['lcfs_z'],'m-') 

        elif axes==['X','Y']:
            if real_scale is True: #set x and y plot limits to real scales
                ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
                ax.set_aspect('equal')
            if LCFS is True: #plot plasma boundary
                plasma_max_R=np.max(some_equilibrium['lcfs_r'])
                plasma_min_R=np.min(some_equilibrium['lcfs_r'])
                ax.plot(plasma_max_R*np.cos(np.linspace(0,2.0*pi,100)),plasma_max_R*np.sin(np.linspace(0.0,2.0*pi,100)),'m-')
                ax.plot(plasma_min_R*np.cos(np.linspace(0,2.0*pi,100)),plasma_min_R*np.sin(np.linspace(0.0,2.0*pi,100)),'m-')          
        
        ax.set_xlabel(axes[0])
        ax.set_ylabel(axes[1])
        
    if ax_flag is False:
        plt.show()




def plot_equilibrium(some_equilibrium,key='psirz',LCFS=None,limiters=None,number_contours=20,ax=False):
    """
    plots equilibrium
     
    notes:
        key - selects which data in equilibrium to plot
        LCFS - toggles plasma boundary on/off in 2D plots
        limiters - toggles limiters on/off in 2D plots
        number_contours - set fidelity of 2D contour plot
        ax - take input axes (can be used to stack plots)
    """

    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    #0D data
    if some_equilibrium[key].ndim==0:
        print(some_equilibrium[key])
        return
    
    #>0D data is plottable
    if ax_flag is False: #if user has not externally supplied axes, generate them
        fig = plt.figure() #initialise plot
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
        mesh=ax.contourf(X,Y,Z,levels=np.linspace(np.amin(Z),np.amax(Z),num=number_contours),cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(Z),vmax=np.amax(Z))
        #mesh=ax.pcolormesh(X,Y,Z,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(Z),vmax=np.amax(Z))
        for c in mesh.collections: #for use in contourf
            c.set_edgecolor("face")

        #3D plot
        #ax=ax.axes(projection='3d')
        #ax.view_init(elev=90, azim=None) #rotate the camera
        #ax.plot_surface(X,Y,Z,rstride=1,cstride=1,cmap='viridis',edgecolor='none',linewidth=0,antialiased=True,vmin=np.amin(Z),vmax=np.amax(Z))
        
        plt.colorbar(mesh)
        ax.set_aspect('equal')
        ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
        ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))
        ax.set_xlabel('R [m]')
        ax.set_ylabel('Z [m]')
        ax.set_title(key)

        if limiters is True: #add boundaries if desired
            ax.plot(some_equilibrium['rlim'],some_equilibrium['zlim'],'k-') 
        if LCFS is True:
            ax.plot(some_equilibrium['lcfs_r'],some_equilibrium['lcfs_z'],'m-') 

    if ax_flag is False:
        plt.show()


def plot_B_field_line(some_equilibrium,axes=['X','Y','Z'],number_field_lines=1,angle=2.0*pi,plot_full=False,start_mark=False,boundary=False,ax=False):
    """
    plots random field lines for 'angle' radians around the tokamak

    notes:
        essentially uses the Euler method of integration
        number_field_lines - the number of field lines to plot
        angle - plot field line for this many radians around central column
        plot_full - choose whether each field line will trace the whole plasma topology (see below also)
        start_mark - add marker for start of the field line
        boundary - show plasma boundary outline
        ax - take input axes (can be used to stack plots)
    """
    
    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if 'B_field' not in some_equilibrium.data: #check we have a B field first
        print("ERROR: 'B_field' missing in equilibrium object (calculate first with B_calc)")
        return

    if ax_flag is False: #if user has not externally supplied axes, generate them
        fig = plt.figure() #initialise plot
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
        if ax_flag is False:
            ax = fig.gca(projection='3d')
        ax.set_xlim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
        ax.set_ylim(-1.0*np.max(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
        ax.set_zlim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D'])) 


    print('plot_B_field_line - generating B field interpolators')
    B_field_R_interpolator=process_input.interpolate_2D(some_equilibrium['R_1D'],some_equilibrium['Z_1D'],some_equilibrium['B_field'][:,:,0])
    B_field_Z_interpolator=process_input.interpolate_2D(some_equilibrium['R_1D'],some_equilibrium['Z_1D'],some_equilibrium['B_field'][:,:,2])
    B_field_tor_interpolator=process_input.interpolate_2D(some_equilibrium['R_1D'],some_equilibrium['Z_1D'],some_equilibrium['B_field'][:,:,1])
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
            B_field_mag=np.sqrt(B_field_R**2+B_field_Z**2+B_field_tor**2)   #calculate the magnitude
            
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
            ax.plot(R_points,Z_points,color=cmap_viridis(np.random.uniform()))
            if start_mark: 
                ax.scatter(R_points[0],Z_points[0],color='red',s=10)
        else:
            if start_mark: 
                ax.scatter(X_points[0],Y_points[0],Z_points[0],color='red',s=10)
            if axes==['R','Z']:
                ax.plot(R_points,Z_points,color=cmap_viridis(np.random.uniform()))
            elif axes==['X','Y']:
                ax.plot(X_points,Y_points,color=cmap_viridis(np.random.uniform()))
            else:
                ax.plot(X_points,Y_points,zs=Z_points,color=cmap_viridis(np.random.uniform()))

    if boundary: #plot periodic poloidal cross-sections in 3D
        print('plot_B_field_line - plotting plasma boundary')
        for angle in np.linspace(0.0,2.0*pi,4,endpoint=False):
            x_points=some_equilibrium['lcfs_r']*np.cos(angle)
            y_points=some_equilibrium['lcfs_r']*np.sin(angle)
            z_points=some_equilibrium['lcfs_z']
            ax.plot(x_points,y_points,zs=z_points,color='m')

    if ax_flag is False:
        plt.show()


def plot_B_field_stream(some_equilibrium,ax=False):
    """
    stream plot of magnetic field in R,Z plane

    notes:
        take transpose due to streamplot index convention
        ax - take input axes (can be used to stack plots)
    """

    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if ax_flag is False: #if user has not externally supplied axes, generate them #initialise plot
        fig = plt.figure() #initialise plot
        ax = fig.add_subplot(111)
    ax.set_aspect('equal')
        
    if 'B_field' not in some_equilibrium.data: #check we have a B field first
        print("ERROR: B_field missing in equilibrium object (calculate first with B_calc)")
        return

    B_mag=np.sqrt(some_equilibrium['B_field'][:,:,0]**2+some_equilibrium['B_field'][:,:,2]**2) #calculate poloidal field magnitude
    strm = ax.streamplot(some_equilibrium['R_1D'],some_equilibrium['Z_1D'],some_equilibrium['B_field'][:,:,0].T,some_equilibrium['B_field'][:,:,2].T, color=B_mag.T, linewidth=1, cmap='viridis')
    plt.colorbar(strm.lines)
    ax.set_xlim(np.min(some_equilibrium['R_1D']),np.max(some_equilibrium['R_1D']))
    ax.set_ylim(np.min(some_equilibrium['Z_1D']),np.max(some_equilibrium['Z_1D']))

    if ax_flag is False:
        plt.show()
    

#################################

##################################################################

###################################################################################################
