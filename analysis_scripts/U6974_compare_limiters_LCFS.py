import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm #get colourmaps

import context
from classes.input_classes.equilibrium import Equilibrium
from classes.input_classes.wall import Wall 

radii=['1.05','1.10','1.20','1.30','1.40','1.50'] #radii for limiter profiles
wall_file='LCFS_DIII-D.dat_'

filename_eq='g157418.03000'
equi=Equilibrium(filename_eq,'GEQDSK',filename_eq)

fig,(ax1)=plt.subplots(1)

for radius in radii:
    wall=Wall(radius,'LOCUST_2D',wall_file+radius)
    ax1.plot(wall['rlim'],wall['zlim'],'k-')

ax1.plot(equi['lcfs_r'],equi['lcfs_z'],'m-')
ax1.set_aspect('equal')

plt.show()