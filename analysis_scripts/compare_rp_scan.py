import context
from classes.output_classes.particle_list import Final_Particle_List 
from classes.input_classes.equilibrium import Equilibrium 
import processing.utils
import settings
import support
from scipy.io import readsav
from scipy.spatial import ConvexHull
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

#define some colourmaps
cmap_r=settings.colour_custom([194,24,91,1])
cmap_g=settings.colour_custom([76,175,80,1])
cmap_g_=settings.colour_custom([205,220,57,1])

filename_sav='compare_rp.sav'
filepath_sav=str(support.dir_output_files/filename_sav)
sav=readsav(filepath_sav)
status_flags={}
status_flags['PFC_intercept_2D']=0.

fpl_locust_3D=Final_Particle_List('')
fpl_locust_3D.set(
    R_start=sav['b'].rpzi[0][sav['wloc'],0],
    V_pitch_start=sav['b'].pitches[0][sav['wloc'],0],
    status_flag=np.full(len(sav['b'].rpzi[0][sav['wloc'],0]),0.0),
    status_flags=status_flags
    )

fpl_orb=Final_Particle_List('')
fpl_orb.set(
    R_start=sav['b'].rpzi[0][sav['wlo'],0],
    V_pitch_start=sav['b'].pitches[0][sav['wlo'],0],
    status_flag=np.full(len(sav['b'].rpzi[0][sav['wlo'],0]),0.0),
    status_flags=status_flags
    )

phase_space_R,phase_space_pitch=np.meshgrid(sav['b'].rpzi[0][:,0],sav['b'].pitches[0][:,0])
phase_space_R,phase_space_pitch=sav['b'].rpzi[0][:,0],sav['b'].pitches[0][:,0]

eq=Equilibrium(ID='',filename='g157418.03000',data_format='GEQDSK',GEQDSKFIX1=True)
fpl_locust_2D=Final_Particle_List(ID='',filename='ptcl_cache.dat',data_format='LOCUST') #this is different since ran these simulations later
fpl_locust_2D.set(R_start=phase_space_R)
fpl_locust_2D.set(V_pitch_start=phase_space_pitch)
fpl_locust_2D.set(V_pitch=processing.utils.pitch_calc_2D(particle_list=fpl_locust_2D,equilibrium=eq))

theta=[]
for R,Z in zip(fpl_locust_2D['R'],fpl_locust_2D['Z']): theta.append(processing.utils.angle_pol(R_major=eq['rmaxis'],Z_major=eq['zmaxis'],R=R,Z=Z))
fpl_locust_2D.set(theta=np.array(theta)-2.*np.pi)
PFC_intercept_2D=np.where(fpl_locust_2D['status_flag']==fpl_locust_2D['status_flags']['PFC_intercept_2D'])[0] #find indices of markers which intercept wall

fig=plt.figure()
axgrid=gridspec.GridSpec(2,2)
ax1=plt.subplot(axgrid[0,:])
scatter=True
if scatter:
    ax1.scatter(sav['b'].rpzi[0][:,0],sav['b'].pitches[0][:,0],color='k',marker='.',s=5,label='not lost') #start
    ax1.scatter(sav['b'].rpzi[0][sav['wloc'],0],sav['b'].pitches[0][sav['wloc'],0],color=cmap_g(0.),marker='o',s=60,label='LOCUST 3D') #locust 3d
    ax1.scatter(sav['b'].rpzi[0][sav['wlo'],0],sav['b'].pitches[0][sav['wlo'],0],color=cmap_r(0.),marker='x',s=40,label='ORB 2D') #orb
    ax1.scatter(fpl_locust_2D['R_start'][PFC_intercept_2D],fpl_locust_2D['V_pitch_start'][PFC_intercept_2D],color=cmap_g_(0.),marker='+',s=20,label='LOCUST 2D') #locust 2d
    ax1.set_facecolor('w')
else: #XXX this code is sketchy since perimeter hull is concave
    ax1.scatter(sav['b'].rpzi[0][:,0],sav['b'].pitches[0][:,0],color='k',marker='.',s=5) #start
    #fpl_locust_2D.plot(axes=['R_start','V_pitch_start'],status_flags=['PFC_intercept_2D'],colmap=cmap_g_,fill=True,number_bins=4)
    for particle_list,color in zip([fpl_orb,fpl_locust_2D,fpl_locust_3D],[settings.cmap_r,settings.cmap_g,settings.cmap_m]):
        hull=ConvexHull([[R,pitch] for R,pitch in zip(particle_list['R_start'],particle_list['V_pitch_start'])])
        ax1.plot(particle_list['R_start'][hull.vertices],particle_list['V_pitch_start'][hull.vertices],color=color(0.5))
    fpl_locust_2D.plot(axes=['R_start','V_pitch_start'],status_flags=['PFC_intercept_2D'],colmap=settings.cmap_g,ax=ax1,fig=fig,fill=False,number_bins=3)
    fpl_locust_3D.plot(axes=['R_start','V_pitch_start'],status_flags=['PFC_intercept_2D'],colmap=settings.cmap_r,ax=ax1,fig=fig,fill=False,number_bins=2)
    fpl_orb.plot(axes=['R_start','V_pitch_start'],status_flags=['PFC_intercept_2D'],colmap=settings.cmap_m,ax=ax1,fig=fig,fill=False,number_bins=2)

ax2=plt.subplot(axgrid[1,:])
ax2.set_xlim([-np.pi,np.pi])
ax2.set_ylim([-np.pi,np.pi])
ax2.scatter(sav['phiw'][sav['wlo']],sav['thetw'][sav['wlo']],color=cmap_r(0.),s=5)
ax2.scatter(sav['phiw_loc'],sav['thetw_loc'],color=cmap_g(0.),s=5)
ax2.scatter(fpl_locust_2D['phi'][PFC_intercept_2D],fpl_locust_2D['theta'][PFC_intercept_2D],color=cmap_g_(0.),s=5)



ax1.set_xlabel('R [m]')
ax1.set_ylabel('$\lambda$')
ax1.legend(markerscale=2)
ax1.set_title("a) Prompt losses in initial position phase space",pad=20)
ax2.set_title("b) Prompt losses in final position phase space",pad=20)
ax2.set_xlabel('Toroidal angle [rad]')
ax2.set_ylabel('Poloidal angle [rad]')
plt.show()


print('LOCUST Loss 3D: '+str(sav['lossfrac_loc'])[0:5])
print('LOCUST Loss 2D: '+str(len(np.where(fpl_locust_2D['status_flag']==fpl_locust_2D['status_flags']['PFC_intercept_2D'])[0])/len(fpl_locust_2D['status_flag'])))
print('DIII-D Loss: '+str(sav['lossfrac_orb'])[0:5])