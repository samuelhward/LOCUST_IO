import context
import input_classes as ic
import process_input as prin
import output_classes as oc 
import plot_input as plin
import plot_output as plou
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d #import 3D plotting axes
from mpl_toolkits.mplot3d import Axes3D

#these are 1us no flags
orbit_filenames_1=['ORBIT_05-04-2018_22-52-55.712.dat','ORBIT_05-04-2018_23-00-17.585.dat','ORBIT_05-04-2018_23-02-22.484.dat','ORBIT_05-04-2018_23-03-13.312.dat']

#these are 10us BFLIP and GEQDSKFIX1
orbit_filenames_2=['ORBIT_06-04-2018_20-14-19.910.dat','ORBIT_06-04-2018_20-18-51.734.dat','ORBIT_06-04-2018_20-20-18.702.dat','ORBIT_06-04-2018_20-21-23.441.dat']

#these are 10us with GEQDSKFIX1
orbit_filenames_2=['ORBIT_08-04-2018_06-34-09.140.dat','ORBIT_08-04-2018_06-35-52.322.dat','ORBIT_08-04-2018_06-38-02.390.dat','ORBIT_08-04-2018_06-39-11.531.dat']

#these are 3us with GEQDSKFIX1
orbit_filenames_2=['ORBIT_09-04-2018_06-18-54.661.dat','ORBIT_09-04-2018_06-20-12.100.dat','ORBIT_09-04-2018_06-21-03.962.dat','ORBIT_09-04-2018_06-21-54.182.dat']

#orbit_filenames_1=[]
#orbit_filenames_1=[]
#orbit_filenames_1=[]
#orbit_filenames_1=[]



#read in the two sets of orbits
ORBIT_1a=oc.Orbits('ORBIT_1a','ASCII',orbit_filenames_1[0])
ORBIT_2a=oc.Orbits('ORBIT_2a','ASCII',orbit_filenames_1[1])
ORBIT_3a=oc.Orbits('ORBIT_3a','ASCII',orbit_filenames_1[2])
ORBIT_4a=oc.Orbits('ORBIT_4a','ASCII',orbit_filenames_1[3])
ORBIT_1b=oc.Orbits('ORBIT_1b','ASCII',orbit_filenames_2[0])
ORBIT_2b=oc.Orbits('ORBIT_2b','ASCII',orbit_filenames_2[1])
ORBIT_3b=oc.Orbits('ORBIT_3b','ASCII',orbit_filenames_2[2])
ORBIT_4b=oc.Orbits('ORBIT_4b','ASCII',orbit_filenames_2[3])






#import GEQDSK and calculate magnetic field
eq_filename='g157418.03000'
eq=ic.Equilibrium('ID','GEQDSK',eq_filename)
eq.set(fpolrz=prin.fpolrz_calc(eq))
eq.set(B_field=prin.B_calc(eq))







#2D - try plotting with plot_equilibrium?

#ORBIT 1 poloidal
fig, some_axes = plt.subplots(1, 1)
plou.plot_orbits(ORBIT_1a,eq,real_scale=True,start_mark=True,plasma=True,axes=['r','z'],ax=some_axes)
plou.plot_orbits(ORBIT_1b,eq,real_scale=True,start_mark=True,plasma=True,axes=['r','z'],ax=some_axes)
plt.show()

#ORBIT 1 top down
fig, some_axes = plt.subplots(1, 1)
plou.plot_orbits(ORBIT_1a,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y'],ax=some_axes)
plou.plot_orbits(ORBIT_1b,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y'],ax=some_axes)
plt.show()




#ORBIT 2 poloidal
fig, some_axes = plt.subplots(1, 1)
plou.plot_orbits(ORBIT_2a,eq,real_scale=True,start_mark=True,plasma=True,axes=['r','z'],ax=some_axes)
plou.plot_orbits(ORBIT_2b,eq,real_scale=True,start_mark=True,plasma=True,axes=['r','z'],ax=some_axes)
plt.show()

#ORBIT 2 top down
fig, some_axes = plt.subplots(1, 1)
plou.plot_orbits(ORBIT_2a,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y'],ax=some_axes)
plou.plot_orbits(ORBIT_2b,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y'],ax=some_axes)
plt.show()





#ORBIT 3 poloidal
fig, some_axes = plt.subplots(1, 1)
plou.plot_orbits(ORBIT_3a,eq,real_scale=True,start_mark=True,plasma=True,axes=['r','z'],ax=some_axes)
plou.plot_orbits(ORBIT_3b,eq,real_scale=True,start_mark=True,plasma=True,axes=['r','z'],ax=some_axes)
plt.show()

#ORBIT 3 top down
fig, some_axes = plt.subplots(1, 1)
plou.plot_orbits(ORBIT_3a,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y'],ax=some_axes)
plou.plot_orbits(ORBIT_3b,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y'],ax=some_axes)
plt.show()




#ORBIT 4 poloidal
fig, some_axes = plt.subplots(1, 1)
plou.plot_orbits(ORBIT_4a,eq,real_scale=True,start_mark=True,plasma=True,axes=['r','z'],ax=some_axes)
plou.plot_orbits(ORBIT_4b,eq,real_scale=True,start_mark=True,plasma=True,axes=['r','z'],ax=some_axes)
plt.show()

#ORBIT 4 top down
fig, some_axes = plt.subplots(1, 1)
plou.plot_orbits(ORBIT_4a,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y'],ax=some_axes)
plou.plot_orbits(ORBIT_4b,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y'],ax=some_axes)
plt.show()









'''
#3D

#ORBIT 1
fig = plt.figure() #initialise plot
some_axes = fig.gca(projection='3d')
plin.plot_B_field_line(eq,start_mark=True,plot_full=False,angle=2.0,boundary=True,ax=some_axes)
plou.plot_orbits(ORBIT_1a,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y','z'],ax=some_axes)
plou.plot_orbits(ORBIT_1b,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y','z'],ax=some_axes)
plt.show()


#ORBIT 2
fig = plt.figure() #initialise plot
some_axes = fig.gca(projection='3d')
plin.plot_B_field_line(eq,start_mark=True,plot_full=False,angle=2.0,boundary=True,ax=some_axes)
plou.plot_orbits(ORBIT_2a,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y','z'],ax=some_axes)
plou.plot_orbits(ORBIT_2b,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y','z'],ax=some_axes)
plt.show()



#ORBIT 3
fig = plt.figure() #initialise plot
some_axes = fig.gca(projection='3d')
plin.plot_B_field_line(eq,start_mark=True,plot_full=False,angle=2.0,boundary=True,ax=some_axes)
plou.plot_orbits(ORBIT_3a,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y','z'],ax=some_axes)
plou.plot_orbits(ORBIT_3b,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y','z'],ax=some_axes)
plt.show()



#ORBIT 4
fig = plt.figure() #initialise plot
some_axes = fig.gca(projection='3d')
plin.plot_B_field_line(eq,start_mark=True,plot_full=False,angle=2.0,boundary=True,ax=some_axes)
plou.plot_orbits(ORBIT_4a,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y','z'],ax=some_axes)
plou.plot_orbits(ORBIT_4b,eq,real_scale=True,start_mark=True,plasma=True,axes=['x','y','z'],ax=some_axes)
plt.show()
'''



