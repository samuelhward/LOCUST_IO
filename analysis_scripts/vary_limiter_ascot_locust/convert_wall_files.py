#convert wall files

import context
from classes.input_classes.wall import Wall as w
import matplotlib.pyplot as plt

head='OMF157418.LIM_'
wall_files=['0.95','1.00','1.01']

fig,ax=plt.subplots(1)
for wall_file in wall_files:
    file=head+wall_file
    wall=w('','UFILE',file)       
    wall.plot(ax=ax,fig=fig,colmap='r')
    wall.dump_data('LOCUST_2D','LCFS_DIII-D.dat_'+wall_file)
    wall.dump_data('ASCOT_2D_input','input.wall_2d_'+wall_file)

plt.show()