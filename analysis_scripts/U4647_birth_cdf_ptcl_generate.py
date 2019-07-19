#script to generate beam deposition input for LOCUST

import context
from classes.input_classes.beam_deposition import Beam_Deposition as bd
from classes.input_classes.equilibrium import Equilibrium as eq
from processing import plot_input as plin
import matplotlib.pyplot as plt
import numpy as np

run_1='U46'
run_2='U47'

filename_eq='g157418.03000'
equi=eq('','GEQDSK',filename_eq)

number_files=13 #number files per run
beam_depo_array=[]

for file in range(number_files): #add the files
    file+=1 #starts from file 1
    filename='157418{}_birth.cdf{}'.format(run_1,file)
    beam_depo_array.append(bd(str(file),'TRANSP_birth_gc',filename))

for file in range(number_files): #add the next files
    file+=1 #starts from file 1
    filename='157418{}_birth.cdf{}'.format(run_2,file)
    beam_depo_array.append(bd(str(file),'TRANSP_birth_gc',filename))


R_out=np.array([])
Phi_out=np.array([])
Z_out=np.array([])
E_out=np.array([])
V_pitch_out=np.array([])
weight_out=np.array([])
X_out=np.array([])
Y_out=np.array([])

for beam_depo in beam_depo_array:
    R_out=np.append(R_out,beam_depo['R'])
    Phi_out=np.append(Phi_out,beam_depo['phi'])
    Z_out=np.append(Z_out,beam_depo['Z'])
    E_out=np.append(E_out,beam_depo['E'])
    V_pitch_out=np.append(V_pitch_out,beam_depo['V_pitch'])
    weight_out=np.append(weight_out,beam_depo['weight'])
    X_out=np.append(X_out,beam_depo['X'])
    Y_out=np.append(Y_out,beam_depo['Y'])

blank_beam_depot=bd('COMBINED')
blank_beam_depot.set(R=R_out,Z=Z_out,phi=Phi_out,E=E_out,V_pitch=V_pitch_out,weight=weight_out,X=X_out,Y=Y_out)

blank_beam_depot.dump_data('LOCUST_weighted','{}{}_ptcles.dat'.format(run_1,run_2))
