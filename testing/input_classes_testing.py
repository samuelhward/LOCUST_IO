#!/usr/bin/env python

#input_classes_testing.py

"""
Samuel Ward
02/11/2017
----
Small examples to test the input_classes.py methods
---
notes:
---
"""



##################################################################
#Preamble

import sys
import numpy as np
import matplotlib.pyplot as plt

try:
	import context
except:
	raise ImportError("ERROR: local context.py not found - please create one and add it to this directory\nreturning\n")
	sys.exit(1)

try: 
	import input_classes
except:
	raise ImportError("ERROR: LOCUST_IO/classes/input_classes.py could not be imported!\nreturning\n")
	sys.exit(1)





##################################################################
#Main Code


#################################first test is to open an Equilibrium GEQDSK and test that reading in upon initialisation is OK and that all member data is set correctly
test1_ID='test1_ID'
test1_input_filename='test.eqdsk'
test1_data_format='GEQDSK'
test_1=input_classes.Equilibrium(test1_ID,test1_input_filename,test1_data_format)
print(test_1.input_filename)
print(test_1.data_format)
print(test_1.LOCUST_input_type)
print(test_1.ID)


#################################next test to see if we can initialise a blank Equilibrium before populating with GEQDSK formatted data using read_data
test2_ID='test2_ID'
test2_input_filename='test.eqdsk'
test2_data_format='GEQDSK'
test_2=input_classes.Equilibrium(test2_ID)
test_2.read_data(test2_input_filename,test2_data_format)
print(test_2.input_filename)
print(test_2.data_format)
print(test_2.LOCUST_input_type)
print(test_2.ID)


#################################next test to see if we can dump an Equilibrium to GEQDSK
test3_ID='test3_ID'
test3_input_filename='test.eqdsk'
test3_data_format='GEQDSK'
test3=input_classes.Equilibrium(test3_ID,test3_input_filename,test3_data_format)
test3.dump_data('test3_output.eqdsk','GEQDSK')


#################################now need to check the data we just outputted in test 3 to see if it's the same - subtract all the member data away from eachother to see if it's = 0
test4_ID='test4_ID'
test4_input_filename='test3_output.eqdsk'
test4_data_format='GEQDSK'
test4=input_classes.Equilibrium(test4_ID,test4_input_filename,test4_data_format)
print(np.subtract(test4.data['fpol'],test3.data['fpol']))







#QTP_calc() tests

#analytical
P=np.linspace(0,8,num=1000) 
T=np.sin(P)
Q=np.cos(P)

#QTP results
P_recover=prin.QTP_calc(T=T,Q=Q)
Q_recover=prin.QTP_calc(T=T,P=P)
T_recover=prin.QTP_calc(Q=Q,P=P)

plt.plot(P_recover-P,'r-')
plt.plot(T_recover-T,'b-')
plt.plot(Q_recover-Q,'g-')
plt.show()


#################################

##################################################################

###################################################################################################