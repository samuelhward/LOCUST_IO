import context
from plot_scripts.plot_collision_operator import plot_collision_operator
import matplotlib.pyplot as plt
import numpy as np
from constants import *



fh = open('ASCOT_coeffs.dat','r')
fh.readline()
cf = np.zeros([120,12])
for i in range(0,120):
    cf[i,:] = [float(x) for x in fh.readline().split()]
fh.close()
e = 1.602176565e-19


fig,ax=plt.subplots(1)
E,drag=plot_collision_operator(At=mass_deuteron_amu,Ai=[mass_deuteron_amu,mass_electron_amu],Zt=1,Zi=[1,1],Ti=[9419.,4119.],ni=[5.96e19,5.96e19],Einj=80000,Pdep=1,Bmod=1.5,type='LOCUST',fig=fig,ax=ax)
ax.plot(cf[:,0]/e/1e3,cf[:,10]/1.e-13,color='C0',linestyle='--',label='ASCOT, electron')
ax.plot(cf[:,0]/e/1e3,cf[:,11]/1.e-13,color='C1',linestyle='--',label='ASCOT, ion')
ax.legend()
ax.set_title('Energy drift (D ion on D background)')
ax.grid()
plt.show()


fig,ax=plt.subplots(1)
ax.plot(E/e/1e3,(drag[0]+drag[1])/1.e-13,color='C0',linestyle='--',label='LOCUST, total')
ax.plot(cf[:,0]/e/1e3,(cf[:,10]+cf[:,11])/1.e-13,color='C1',linestyle='--',label='ASCOT total')
ax.legend()
plt.show()


fig,ax=plt.subplots(1)
ax.plot(E/e/1e3,1./-(drag[0]+drag[1]),color='C1',linestyle='--',label='LOCUST, total')
ax.plot(cf[:,0]/e/1e3,-1/(cf[:,10]),color='C0',linestyle='--',label='ASCOT total')
ax.legend()
plt.show()

