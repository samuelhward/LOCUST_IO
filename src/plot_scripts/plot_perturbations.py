#plot_perturbations.py

'''
Samuel Ward
20/08/2020
----
plots multiple perturbations
---
notes:
    essentially wrapping individual perturbation.plot() method 
---
'''


##################################################################
#Preamble

import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway

if __name__=='__main__':

    try:
        import context
    except:
        raise ImportError("ERROR: context.py could not be imported!\nreturning\n")
        sys.exit(1)

try:
    import numpy as np
    import matplotlib.pyplot as plt
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import run_scripts.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/utils.py could not be imported!\nreturning\n")
    sys.exit(1)
    
try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/src/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import constants
except:
    raise ImportError("ERROR: LOCUST_IO/src/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import settings
except:
    raise ImportError("ERROR: LOCUST_IO/src/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)


##################################################################

def plot_perturbations(perturbations,key='dB_field_R_real',axes=['R','Z'],LCFS=False,limiters=False,number_bins=20,fill=True,vminmax=None,i3dr=-1,phase=0.,colmap=settings.cmap_default,colmap_val=np.random.uniform(),line_style=settings.plot_line_style,gridlines=False,label='',ax=False,fig=False):
    """
    plots a perturbation
    
    notes:
        
    args:
        key - selects which data in perturbation to plot
        axes - list of strings specifying which axes should be plotted (current options only XY or RZ)
        LCFS - toggles plasma boundary on/off in 2D plots (requires equilibrium argument)
        limiters - object which contains limiter data rlim and zlim
        number_bins - set number of bins or levels
        fill - toggle contour fill on 2D plots
        vminmax - set mesh Vmin/Vmax values
        i3dr - flip definition of phi (+1 anti-clockwise, -1 clockwise)
        phase - global field phase shift (of field origin with respect to locust origin) (radians, anti-clockwise)
        colmap - set the colour map (use get_cmap names)
        colmap_val - optional numerical value for defining single colour plots 
        line_style - set 1D line style
        gridlines - toggle gridlines on plot
        label - plot label for legends
        ax - take input axes (can be used to stack plots)
        fig - take input fig (can be used to add colourbars etc)
    """

    #do some preliminary parsing of variables in case supplied as strings from command line etc.
    axes,number_bins,vminmax,i3dr,phase,colmap_val=run_scripts.utils.literal_eval(axes,number_bins,vminmax,i3dr,phase,colmap_val)

    import scipy
    import matplotlib
    from matplotlib import cm
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d #import 3D plotting axes
    from mpl_toolkits.mplot3d import Axes3D

    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if fig is False:
        fig_flag=False
    else:
        fig_flag=True
    
    #>0D data is plottable
    if fig_flag is False:
        fig = plt.figure() #if user has not externally supplied figure, generate
    
    if ax_flag is False: #if user has not externally supplied axes, generate them
        polar=True if axes==['phi','R'] else False
        ax = fig.add_subplot(111,polar=polar)
    
    ax.grid(False) if gridlines is False else ax.grid(True)
    if gridlines is False: ax.set_yticklabels([])

    #1D data

    if key in perturbations[0].data.keys():
        if perturbations[0][key].ndim==1:
            ax.plot(np.sum(perturbation[key] for perturbation in perturbations),color=colmap(colmap_val),label=label)
            ax.set_ylabel(key)

    #2D data
    if len(axes)==2:
        if axes==['R','Z']:

            R=perturbations[0]['R_1D'] #make a mesh
            Z=perturbations[0]['Z_1D'] 
            dr,dz=R[1]-R[0],Z[1]-Z[0]
            ax.set_xticks(R) #set axes ticks
            ax.set_yticks(Z)

            Z,R=np.meshgrid(Z-dz/2.,R-dr/2.) #offset ticks onto bin centres

            if key not in perturbations[0].data.keys():

                R_poloidal,Z_poloidal=np.meshgrid(perturbations[0]['R_1D'],perturbations[0]['Z_1D']) 
                R_flat,Z_flat=R_poloidal.flatten(),Z_poloidal.flatten()
                phi_flat=np.full(len(R_flat),0.) #XXX this zero should be phi parameter for tomographic slices

                dB_R=np.zeros(len(R_flat))
                dB_tor=np.zeros(len(R_flat))
                dB_Z=np.zeros(len(R_flat))

                for perturbation in perturbations:
                    dB=perturbation.evaluate(R=R_flat,phi=phi_flat,Z=Z_flat,phase=phase,i3dr=i3dr,mode_number=perturbation.mode_number) #evaluate poloidally
                    dB_R+=dB[0]
                    dB_tor+=dB[1]
                    dB_Z+=dB[2]       

                if key=='dB_field_R':
                    values=dB_R
                elif key=='dB_field_tor':
                    values=dB_tor
                elif key=='dB_field_Z':
                    values=dB_Z
                elif key=='dB_field_mag':
                    values=np.sqrt(dB_R**2+dB_tor**2+dB_Z**2)
                else:
                    print(f"WARNING: perturbation.plot() could not plot key={key} in {perturbation.ID} - plotting dB_field_mag!\n")
                    key='dB_field_mag'
                    values=np.sqrt(dB_R**2+dB_tor**2+dB_Z**2)
            
                values=values.reshape(len(perturbations[0]['Z_1D']),len(perturbations[0]['R_1D'])).T

            else:
                values=np.sum([perturbation[key] for perturbation in perturbations])
 
            if vminmax:
                vmin=vminmax[0]
                vmax=vminmax[1]
            else:
                vmin=np.min(values)
                vmax=np.max(values)

            #2D plot
            if fill is True:
                mesh=ax.contourf(R,Z,values,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',antialiased=True,vmin=vmin,vmax=vmax)
                for c in mesh.collections: #for use in contourf
                    c.set_edgecolor("face")
            else:
                mesh=ax.contour(R,Z,values,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',antialiased=True,vmin=vmin,vmax=vmax)
                if settings.plot_contour_labels:
                    ax.clabel(mesh,inline=1,fontsize=10)
                
            #mesh=ax.pcolormesh(R,Z,values,colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',antialiased=True,vmin=np.amin(values),vmax=np.amax(values))

            #3D plot
            #ax=ax.axes(projection='3d')
            #ax.view_init(elev=90, azim=None) #rotate the camera
            #ax.plot_surface(R,Z,values,rstride=1,cstride=1,colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',antialiased=True,vmin=np.amin(values),vmax=np.amax(values))
            
            ax.set_aspect('equal')
            ax.set_xlim(np.min(perturbations[0]['R_1D']),np.max(perturbations[0]['R_1D']))
            ax.set_ylim(np.min(perturbations[0]['Z_1D']),np.max(perturbations[0]['Z_1D']))
            ax.set_xlabel('R [m]')
            ax.set_ylabel('Z [m]')

            if LCFS:
                ax.plot(LCFS['lcfs_r'],LCFS['lcfs_z'],color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS') 
            if limiters: #add boundaries if desired
                ax.plot(limiters['rlim'],limiters['zlim'],color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall') 

        elif axes==['phi','R']:

            R=perturbations[0]['R_1D'] #make a mesh
            phi=np.linspace(0.,2.*np.pi,100) 
            nR,nphi=len(R),len(phi)
            phi,R=np.meshgrid(phi,R)
            R_flat,phi_flat=R.flatten(),phi.flatten()
            Z_flat=np.full(len(phi_flat),0.) #XXX this zero should be Z parameter for tomographic slices

            dB_R=np.zeros(len(R_flat))
            dB_tor=np.zeros(len(R_flat))
            dB_Z=np.zeros(len(R_flat))

            if key not in perturbations[0].data.keys():

                for perturbation in perturbations:
                    dB=perturbation.evaluate(R_flat,phi_flat,Z_flat,mode_number=perturbation.mode_number,i3dr=i3dr,phase=phase)
                    dB_R+=dB[0]
                    dB_tor+=dB[1]
                    dB_Z+=dB[2]     

                if key=='dB_field_R':
                    values=dB_R
                elif key=='dB_field_tor':
                    values=dB_tor
                elif key=='dB_field_Z':
                    values=dB_Z
                elif key=='dB_field_mag':
                    values=np.sqrt(dB_R**2+dB_tor**2+dB_Z**2)
                else:
                    print(f"WARNING: perturbation.plot() could not plot key={key} in {perturbation.ID} - plotting dB_field_mag!\n")
                    key='dB_field_mag'
                    values=np.sqrt(dB_R**2+dB_tor**2+dB_Z**2)
                
                values=values.reshape(nR,nphi)

            else:
                values=np.sum([perturbation[key] for perturbation in perturbations])

            if vminmax:
                vmin=vminmax[0]
                vmax=vminmax[1]
            else:
                vmin=np.min(values)
                vmax=np.max(values)

            #2D plot
            if fill is True:
                mesh=ax.contourf(phi,R,values,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',antialiased=True,vmin=vmin,vmax=vmax)
                for c in mesh.collections: #for use in contourf
                    c.set_edgecolor("face")
            else:
                mesh=ax.contour(phi,R,values,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',antialiased=True,vmin=vmin,vmax=vmax)
                if settings.plot_contour_labels:
                    ax.clabel(mesh,inline=1,fontsize=10)

            if LCFS: #plot plasma boundary
                plasma_max_R=np.max(LCFS['lcfs_r'])
                plasma_min_R=np.min(LCFS['lcfs_r'])
                ax.plot(np.linspace(0,2.0*constants.pi,100),np.full(plasma_max_R,100),color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS')
                ax.plot(np.linspace(0,2.0*constants.pi,100),np.full(plasma_min_R,100),color=settings.plot_colour_LCFS,linestyle=settings.plot_line_style_LCFS,label='LCFS')
                ax.set_rmax(1.1*plasma_max_R)

            if limiters: #add boundaries if desired
                limiters_max_R=np.max(limiters['rlim'])
                limiters_min_R=np.min(limiters['rlim'])
                ax.plot(np.linspace(0,2.0*constants.pi,100),np.full(limiters_max_R,100),color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall')
                ax.plot(np.linspace(0,2.0*constants.pi,100),np.full(limiters_min_R,100),color=settings.plot_colour_limiters,linestyle=settings.plot_line_style_limiters,label='wall')
                ax.set_rmax(1.1*limiters_max_R)
                
            ax.set_rmin(0.0)

        elif axes==['X','Y']:

            R=perturbation['R_1D'] #make a mesh
            phi=np.linspace(0.,2.*np.pi,100) 
            nR,nphi=len(R),len(phi)
            phi,R=np.meshgrid(phi,R)
            X,Y=processing.utils.RphiZ_to_XYZ(R,phi,RH=True)
            R_flat,phi_flat=R.flatten(),phi.flatten()
            Z_flat=np.zeros(len(phi_flat))

            if key not in perturbations[0].data.keys():

                dB_R=np.zeros(len(R_flat))
                dB_tor=np.zeros(len(R_flat))
                dB_Z=np.zeros(len(R_flat))

                for perturbation in perturbations:
                    dB=perturbation.evaluate(R_flat,phi_flat,Z_flat,mode_number=perturbation.mode_number,i3dr=i3dr,phase=phase)
                    dB_R+=dB[0]
                    dB_tor+=dB[1]
                    dB_Z+=dB[2]  

                if key=='dB_field_R':
                    values=dB_R
                elif key=='dB_field_tor':
                    values=dB_tor
                elif key=='dB_field_Z':
                    values=dB_Z
                elif key=='dB_field_mag':
                    values=np.sqrt(dB_R**2+dB_tor**2+dB_Z**2)
                else:
                    print(f"WARNING: perturbation.plot() could not plot key={key} in {perturbation.ID} - plotting dB_field_mag!\n")
                    key='dB_field_mag'
                    values=np.sqrt(dB_R**2+dB_tor**2+dB_Z**2)
                
                values=values.reshape(nR,nphi)

            else:
                values=np.sum([perturbation[key] for perturbation in perturbations])

            if vminmax:
                vmin=vminmax[0]
                vmax=vminmax[1]
            else:
                vmin=np.min(values)
                vmax=np.max(values)

            #2D plot
            if fill is True:
                mesh=ax.contourf(X,Y,values,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',antialiased=True,vmin=vmin,vmax=vmax)
                for c in mesh.collections: #for use in contourf
                    c.set_edgecolor("face")
            else:
                mesh=ax.contour(X,Y,values,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',antialiased=True,vmin=vmin,vmax=vmax)
                if settings.plot_contour_labels:
                    ax.clabel(mesh,inline=1,fontsize=10)

        if fig_flag is False:    
            fig.colorbar(mesh,ax=ax,orientation='horizontal')

        if ax_flag is True or fig_flag is True: #return the plot object
            return mesh

    if ax_flag is False and fig_flag is False:
        plt.show()

#################################

##################################################################

###################################################################################################