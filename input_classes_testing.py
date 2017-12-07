#!/usr/bin/env python

#input_classes_testing.py

'''
Samuel Ward
02/11/2017
----
Unit tests for the input_classes.py methods
---

---
'''



##################################################################
#Preamble

import sys
#sys.path.insert(1,'/home/ITER/wards2/random_code_tests/')

try:
	import input_classes 
except:
	raise ImportError("ERROR: input_classes.py could not be imported!\nreturning\n")
	sys.exit(1)

try:
	import numpy
except:
	raise ImportError("ERROR: numpy not found!\nreturning\n")
	sys.exit(1)



##################################################################
#Main Code


ID='eq_ID_here'
input_filename="test.eqdsk"
data_format='GEQDSK'

test=input_classes.Equilibrium(ID,input_filename,data_format)

print(test.input_type)
print(test.data)
print(test.input_filename)
print(test.data_format)

#################################

##################################################################

###################################################################################################