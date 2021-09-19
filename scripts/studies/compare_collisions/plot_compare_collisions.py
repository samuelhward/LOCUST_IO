#plot_compare_collisions.py
 
"""
Samuel Ward
24/08/20
----
plot collision type comparison
---
 
notes:         
---
"""

###################################################################################################
#Preamble
 
import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway

try:
    import numpy as np
    import pathlib
    import copy
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    import os
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

if __name__=='__main__':

    try:
        import context
    except:
        raise ImportError("ERROR: context.py could not be imported!\nreturning\n")
        sys.exit(1)

try:
    from classes.input_classes.equilibrium import Equilibrium
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/equilibrium.py could not be imported!\nreturning\n") 
    sys.exit(1)

try:
    import processing.utils
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

try:
    cwd=pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    sys.path.append(str(cwd.parents[1]))
    import templates.plot_mod
except:
    raise ImportError("ERROR: templates/template_mod.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## 
#Main 

import compare_collisions_launch as batch_data

#cycle through trajectories and poincare maps
poincares=templates.plot_mod.get_output_files(batch_data,'poinc')
particles_to_read=[7]
particles_to_plot=[0]
orbits2D=templates.plot_mod.get_output_files(batch_data,'orbit2D',particles=particles_to_read)
orbits3D=templates.plot_mod.get_output_files(batch_data,'orbit3D',particles=particles_to_read)

for run_number,(orbit2D,orbit3D,poincare) in enumerate(zip(orbits2D,orbits3D,poincares)):

    axes=['R','Z']
    #axes=['phi','theta']
    #axes=['X','Y']
    #axes=['X','Y','Z']
    if len(axes)==3:
        fig=plt.figure()
        from mpl_toolkits import mplot3d #import 3D plotting axes
        ax=fig.gca(projection='3d')
    else:
        fig,ax=plt.subplots(1)

    equilibrium=Equilibrium(ID='',data_format='GEQDSK',filename=batch_data.args_batch['RMP_study__filepath_equilibrium'][run_number],GEQDSKFIX1=True,GEQDSKFIX2=True)
    if orbit3D:
        for quant in ['R','phi','Z','X','Y']:   
            orbit3D[quant]=orbit3D[quant][int(0.9*len(orbit3D[quant])):] #separate out last 10% of trajectory
        orbit3D['theta']=np.mod(processing.utils.angle_pol(R_major=equilibrium['rmaxis'],R=orbit3D['R'],Z=orbit3D['Z'],Z_major=equilibrium['zmaxis'])+np.pi,2.*np.pi)
        if axes==['R','Z']: 
            poincare.plot(colmap=settings.cmap_inferno_r,style='histogram',ax=ax,fig=fig)
            equilibrium.plot(fill=False,LCFS=True,limiters=True,ax=ax,fig=fig)
        orbit3D.plot(particles=particles_to_plot,axes=axes,start_mark=True,colmap=settings.cmap_b,real_scale=True,label='3D',ax=ax,fig=fig)
        if orbit2D:
            orbit2D['theta']=np.mod(processing.utils.angle_pol(R_major=equilibrium['rmaxis'],R=orbit2D['R'],Z=orbit2D['Z'],Z_major=equilibrium['zmaxis'])+np.pi,2.*np.pi)
            orbit2D.plot(particles=particles_to_plot,axes=axes,start_mark=True,colmap=settings.cmap_g,real_scale=True,label='2D',ax=ax,fig=fig)
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_title("")
    ax.legend()
    plt.show()

#################################
    
##################################################################
 
###################################################################################################