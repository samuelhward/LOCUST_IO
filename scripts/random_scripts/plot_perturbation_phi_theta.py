
import context
import matplotlib.pyplot as plt  
import numpy as np 
from classes.input_classes.equilibrium import Equilibrium as eq 
from classes.input_classes.perturbation import Perturbation as pert 
from classes.output_classes.particle_list import Final_Particle_List as fpl 
import processing.utils
import settings

def plot_perturbation_phi_theta(perturbations,equilibrium,psi,vminmax=None,number_bins=100,colmap=settings.cmap_default,ax=None,fig=None):

    # extract contour coordinates at psi in equilibrium
    X=equilibrium['R_1D'] #make a mesh
    Y=equilibrium['Z_1D']
    dx,dy=X[1]-X[0],Y[1]-Y[0]
    Y,X=np.meshgrid(Y-dy/2.,X-dx/2.) #offset ticks onto bin centres
    Z=(equilibrium['psirz']-equilibrium['simag'])/(equilibrium['sibry']-equilibrium['simag']) 
    levels=[
            psi
            ]
    fig_,ax_=plt.subplots(1)
    #ax.contourf(X,Y,Z,cmap='Greys')
    mesh=ax_.contour(X,Y,Z,levels=levels)

    if ax is None:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if fig is None:
        fig_flag=False
    else:
        fig_flag=True

    if fig_flag is False:
        fig = plt.figure() #if user has not externally supplied figure, generate

    if ax_flag is False: #if user has not externally supplied axes, generate them
        ax = fig.add_subplot(111)

    levels_coords=[]
    #total number of points to interpolate onto
    #this will be spread out over the contour's sub-paths (if any) in proportion to their length
    number_coords_new=400

    for contour in mesh.collections: 
        contour_x=[]
        contour_y=[]
        contour_coords_x=[]
        contour_coords_y=[]
        path_lengths=[]
        #scan each separate path to find longest
        #assume longest is closed flux surface and get rid of rest
        for path in contour.get_paths(): 
            v = path.vertices
            path_x=v[:,0]
            path_y=v[:,1]        
            coords=np.array([path_x,path_y])
            distance_total=np.cumsum(np.sqrt(np.sum(np.diff(coords,axis=1)**2,axis=0)))[-1]
            path_lengths.append(distance_total)
        for path_number,path in enumerate(contour.get_paths()): 
            if path_lengths[path_number]==np.max(np.array(path_lengths)):
                v = path.vertices
                contour_x.extend(v[:,0])
                contour_y.extend(v[:,1])
                path_x=v[:,0]
                path_y=v[:,1]        
                #interpolate along trajectory
                coords=np.array([path_x,path_y])
                distance=np.cumsum(np.sqrt(np.sum(np.diff(coords,axis=1)**2,axis=0)))
                distance=np.insert(distance,0,0.) #first distance is 0
                distance=(distance-distance[0])/(distance[-1]-distance[0]) #normalise
                interpolator=processing.utils.interpolate_1D(distance,coords,function='linear',type='interp1d',smooth=0)
                distance_coords_new=np.linspace(0,1,number_coords_new+1)[:-1] #this loops back onto itself, so add  extra point and remove the end
                coords_new=interpolator(distance_coords_new)
                contour_coords_x.extend(coords_new[0,:])
                contour_coords_y.extend(coords_new[1,:])
                #ax.scatter(*coords_new,color='y')
        levels_coords.append([contour_coords_x,contour_coords_y])
        #ax.scatter(contour_x,contour_y,color='r')
        # so now levels_coords[n,0,p] contains the x coordinate for point p along contour n
        # so now levels_coords[n,1,p] contains the y coordinate for point p along contour n
    #theta=np.linspace(0,2.*np.pi)
    #phi=np.linspace(0,2.*np.pi)
    #theta_2D,phi_2D=np.meshgrid(theta,phi)

    phi=np.linspace(0,2.*np.pi,450)

    levels_coords=np.array(levels_coords)
    R_contour=levels_coords[0,0,:]
    Z_contour=levels_coords[0,1,:]
    R_flat=np.array(list(R_contour)*len(phi)).flatten()
    Z_flat=np.array(list(Z_contour)*len(phi)).flatten()
    phi_flat=np.array([np.full(len(R_contour),angle) for angle in phi]).flatten()

    dB_R=np.zeros(len(R_flat))
    dB_tor=np.zeros(len(R_flat))
    dB_Z=np.zeros(len(R_flat))
    for perturbation in perturbations:
        dB=perturbation.evaluate(R=R_flat,phi=phi_flat,Z=Z_flat,phase=0,i3dr=-1,mode_number=perturbation.mode_number) #evaluate poloidally
        dB_R+=dB[0]
        dB_tor+=dB[1]
        dB_Z+=dB[2]

    dB_R=dB_R.reshape(len(R_contour),len(phi),order='F')
    dB_tor=dB_tor.reshape(len(R_contour),len(phi),order='F')
    dB_Z=dB_Z.reshape(len(R_contour),len(phi),order='F')
    values=np.sqrt(dB_R**2+dB_tor**2+dB_Z**2)

    theta=processing.utils.angle_pol(R=R_contour,Z=Z_contour,R_major=0.639277243e1,Z_major=0.597384943)
    
    theta[theta>np.pi]-=2.*np.pi

    if vminmax is None: vminmax=[np.min(values),np.max(values)]
    vmin,vmax=(_ for _ in vminmax)

    mesh=ax.contourf(phi,theta,values,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',antialiased=True,vmin=vmin,vmax=vmax)

    if fig_flag is False:    
        fig.colorbar(mesh,ax=ax,orientation='horizontal')

    if ax_flag is True or fig_flag is True: #return the plot object
        return mesh

    if ax_flag is False and fig_flag is False:
        plt.show()


if __name__=='__main__':

    import matplotlib
    from matplotlib import cm
    e=eq('','GEQDSK','geqdsk_case2')
    p=pert('','MARSF','BPLASMA_n3',mode_number=-3)
    f=fpl('','LOCUST','ptcl_cache.dat')
    f['theta']=processing.utils.angle_pol(R=f['R'],Z=f['Z'],R_major=0.639277243e1,Z_major=0.597384943)
    f['theta'][f['theta']>np.pi]-=2.*np.pi
    f['phi'][f['phi']<0]+=2.*np.pi

    fig,ax=plt.subplots(1)
    plot_perturbation_phi_theta(perturbations=[p],equilibrium=e,psi=1.03,colmap=matplotlib.cm.get_cmap('Greys'),ax=ax,fig=fig)
    f.plot(axes=['phi','theta'],colfield='V_phi',ax=ax,fig=fig,style='scatter')
    plt.show()
