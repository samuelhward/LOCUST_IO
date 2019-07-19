import context
from classes.input_classes.wall import Wall as w
from classes.input_classes.equilibrium import Equilibrium




ASCOT_input=w('ASCOT input','ASCOT_2D_input','input.wall_2d_1.05')
ASCOT_output=w('ASCOT output','ASCOT_2D_output','ascot_freia_8849968.h5')
LOCUST_input=w('LOCUST input','LOCUST_2D','LCFS_DIII-D.dat_1.05')
TRANSP_input=w('TRANSP input','ufile','OMF157418.LIM_1.05')

ASCOT_input.dump_data('LOCUST_2D','test_wall_locust')
LOCUST_test_wall=w('test-wall','LOCUST_2D','test_wall_locust')

filename_eq='g157418.03000'
equi=Equilibrium(filename_eq,'GEQDSK',filename_eq)


import matplotlib.pyplot as plt

fig,(ax1,ax2)=plt.subplots(1,2)

TRANSP_input.plot(colmap='k',ax=ax1,fig=fig,real_scale=True,LCFS=equi)
ASCOT_output.plot(colmap='b',ax=ax1,fig=fig,real_scale=True)
ASCOT_input.plot(colmap='r',ax=ax1,fig=fig,real_scale=True)
LOCUST_input.plot(colmap='g',ax=ax1,fig=fig,real_scale=True)

ASCOT_output.plot(colmap='b',ax=ax2,fig=fig,real_scale=True)
LOCUST_test_wall.plot(colmap='g',ax=ax2,fig=fig,real_scale=True)

plt.show()