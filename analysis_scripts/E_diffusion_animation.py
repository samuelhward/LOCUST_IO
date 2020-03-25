#script to plot 1D energy distribution of TRANSP runs

import context
from classes.input_classes.equilibrium import Equilibrium as eq
from run_scripts import utils

run_1='U46'
run_2='U47'

filename_eq='g157418.03000'
equi=eq('','GEQDSK',filename_eq)

number_files=13 #number files per run
FI_array=[]

for file in range(number_files): #add the files
    file+=1 #starts from file 1
    filename='157418{}_fi_{}_gc.cdf'.format(run_1,file)
    FI_array.append(utils.TRANSP_output_FI(filename))

for file in range(number_files): #add the next files
    file+=1 #starts from file 1
    filename='157418{}_fi_{}_gc.cdf'.format(run_2,file)
    FI_array.append(utils.TRANSP_output_FI(filename))

FI_array[0].plot(axes=['E','time'],TRANSP_output_FI_list=[FI for FI in FI_array])

#FI_1.plot(axes=['R','Z'],some_equilibrium=equi,LCFS=True,limiters=True,real_scale=True)