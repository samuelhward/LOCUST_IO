#quick script for checking magnetic field in g157418.03000 to make sure we map Rhomax to the same point as the artificial limiter in the new LOCUST runs 
import context
from classes.input_classes.equilibrium import Equilibrium as eq 
from processing.utils import interpolate_2D as i2 
import numpy as np

myeq=eq('','GEQDSK','g157418.03000',GEQDSKFIX=1)
interpolator=i2(myeq['R_1D'],myeq['Z_1D'],myeq['psirz'])
print('huasdkas')

r=[2.25535e+00,2.28561e+00,2.29166e+00] #0.95,1.00,1.01
z=[9.21947e-03,9.70470e-03,9.80175e-03]

#r=[2.2691] #makes a difference of 0.5% whether it is actually at Z=0 or not
#z=[0]

for r_,z_ in zip(r,z): #rho defined as sqrt of normalised poloidal flux
    rho=interpolator(r_,z_)[0][0]
    rho=(rho-myeq['simag'])/(myeq['sibry']-myeq['simag'])
    print(rho)
