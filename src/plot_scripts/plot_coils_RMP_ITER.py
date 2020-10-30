#plot_coils_RMP_ITER.py

'''
Samuel Ward
14/05/2020
----
plots the ITER RMP coils
---
notes: 
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
    import matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d #import 3D plotting axes
    from mpl_toolkits.mplot3d import Axes3D
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import run_scripts.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/processing/utils.py could not be imported!\nreturning\n")
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

def plot_coils_RMP_ITER(axes=['R','Z'],shot=1180,run=17,username='public',imasdb='ITER_MD',imas_version='3',imas_entry=0,plot_centres=True,colmap=settings.cmap_default,colmap_val=np.random.uniform(),line_style=settings.plot_line_style,ax=False,fig=False):
    """
    plot the ITER RMP coils

    notes:
        this code only works where a relatively high number of points lie on the part of the coil we care about i.e. the centre of the points must lie sufficiently close to the coil
    args:
        shot - IDS shot number
        run - IDS run number 
        username - IMAS username
        imasdb - local IMAS database name set by 'imasdb' command, sometimes called 'tokamak name'
        imas_version - string denoting IMAS major version e.g. '3'
        imas_entry - IDS entry number
        plot_centres - toggle plotting coil centres
        colmap - plotted line colour
        colmap_val - optional numerical value for defining single colour plots 
        line_style - set 1D line style
        ax - ax object to add plot to
        fig - take input fig (can be used to add colourbars etc)
    """

    axes,shot,run,imas_entry,plot_centres,colmap_val=run_scripts.utils.literal_eval(axes,shot,run,imas_entry,plot_centres,colmap_val)

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
        polar=True if axes==['phi','R'] else False
        ax = fig.add_subplot(111,polar=polar)

        if len(axes)==3: ax = fig.add_subplot(111,projection='3d')

    try:
        import imas 
    except:
        raise ImportError("ERROR: plot_coils_RMP_ITER() could not import IMAS module!\nreturning\n")
        return

    print("plot_coils_RMP_ITER()")

    IDS=imas.ids(int(shot),int(run)) #initialise new blank IDS
    IDS.open_env(username,imasdb,imas_version)
    IDS.coils_non_axisymmetric.get(imas_entry) #open the file and get all the data from it

    R_major=6.412033474019571 #ITER dimension
    Z_major=0.5699986946776257

    
    for counter,coil in enumerate(IDS.coils_non_axisymmetric.coil): 
        
        coil_data={}
        coil_data['R']=coil.conductor[0].elements.start_points.r
        coil_data['Z']=coil.conductor[0].elements.start_points.z
        coil_data['phi']=np.remainder(coil.conductor[0].elements.start_points.phi,2.*np.pi)
        coil_data['X']=coil_data['R']*np.cos(coil_data['phi'])
        coil_data['Y']=coil_data['R']*np.sin(coil_data['phi'])
        coil_data['r_min']=np.array([processing.utils.minor_radius(R_major,R,Z,Z_major) for R,Z in zip(coil_data['R'],coil_data['Z'])])
        coil_data['theta']=np.remainder(np.array([processing.utils.angle_pol(R_major,R,Z,Z_major) for R,Z in zip(coil_data['R'],coil_data['Z'])]),2.*np.pi) 
        
        for counter,(phi_,theta_) in enumerate(zip(coil_data['phi'],coil_data['theta'])): #move coils split about 0/2pi 
            if phi_>+2.*np.pi-0.5: coil_data['phi'][counter]-=2.*np.pi
            if theta_>+np.pi: coil_data['theta'][counter]-=2.*np.pi
        
        mask_old=np.where(coil_data['r_min']<1.4*np.min(coil_data['r_min']))[0] #preprocess slightly with some coarse filtering out of points at high coil_data['r_min']
        tolerance_theta=2.35 #hand-tuned factor denoting how far from initial coil centre all coil points lie
        tolerance_phi=2.2
        tolerance_r_min=2.3
        
        converged=False
        while not converged:
            mask_new=np.where( #initially mask out parts of coil further than average distance to centre
            (np.abs(coil_data['theta']-np.mean(coil_data['theta'][mask_old]))<tolerance_theta*np.mean(np.abs(coil_data['theta'][mask_old]-np.mean(coil_data['theta'][mask_old])))) 
            & (np.abs(coil_data['phi']-np.mean(coil_data['phi'][mask_old]))<tolerance_phi*np.mean(np.abs(coil_data['phi'][mask_old]-np.mean(coil_data['phi'][mask_old]))))
            & (np.abs(coil_data['r_min']-np.mean(coil_data['r_min'][mask_old]))<tolerance_r_min*np.mean(np.abs(coil_data['r_min'][mask_old]-np.mean(coil_data['r_min'][mask_old]))))
            )[0]         
            converged=True if np.all(mask_new==mask_old) else False
            mask_old=mask_new

        items_to_delete=list(set(mask_new)^set(np.arange(len(coil_data['theta']))))
        for key,value in coil_data.items(): coil_data[key]=np.delete(coil_data[key],items_to_delete)

        if plot_centres: ax.scatter(*[np.mean(coil_data[variable]) for variable in axes],color=settings.colour_start_mark,marker=settings.marker_start_mark,s=settings.markersize_start_mark,label='coil centres')
        ax.plot(*[coil_data[variable] for variable in axes],color=colmap(colmap_val),label='ITER RMP coils',linestyle=line_style)

    if ax_flag is False and fig_flag is False:
        plt.show() 


    print("finished plot_coils_RMP_ITER()")


if __name__=='__main__': #plot collision operator and expected steady state distribution function (diffusion disabled)

    import argparse
    parser=argparse.ArgumentParser(description='plot coulomb logarithm')

    parser.add_argument('--axes',type=str,action='store',default="['phi','theta']",dest='axes',help="IDS shot number",required=False)
    parser.add_argument('--shot',type=int,action='store',default=1180,dest='shot',help="IDS shot number",required=False)
    parser.add_argument('--run',type=int,action='store',default=17,dest='run',help="IDS run number",required=False)
    parser.add_argument('--username',type=str,action='store',default='public',dest='username',help="IMAS username",required=False)
    parser.add_argument('--imasdb',type=str,action='store',default='ITER_MD',dest='imasdb',help="local IMAS database name set by 'imasdb' command, sometimes called 'tokamak name'",required=False)
    parser.add_argument('--imas_version',type=str,action='store',default='3',dest='imas_version',help="string denoting IMAS major version e.g. '3'",required=False)
    parser.add_argument('--imas_entry',type=int,action='store',default=0,dest='imas_entry',help="IDS imas_entry number",required=False)
    parser.add_argument('--plot_centres',type=bool,action='store',default=True,dest='plot_centres',help="toggle plotting coil centres",required=False)
    parser.add_argument('--colmap_val',type=float,action='store',default=np.random.uniform(),dest='colmap_val',help="optional numerical value for defining single colour plots",required=False)

    args=parser.parse_args()

    plot_coils_RMP_ITER(**{key:arg for key,arg in args._get_kwargs() if arg is not None},colmap=settings.cmap_default)
    plt.show()

#################################

##################################################################

###################################################################################################