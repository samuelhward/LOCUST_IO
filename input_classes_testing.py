#!/work/imas/opt/EasyBuild/software/Anaconda2/4.4.0/bin/python

#input_classes_testing.py

'''
Samuel Ward
02/11/2017
----
Tests for the input_classes.py methods
---

---
'''



##################################################################
#Preamble

import sys
sys.path.insert(0,'/home/ITER/wards2/random_code_tests/')

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


#initialise some object here and test the class methods 

ID='hm'

input_filename="input_file.txt"
data_format='GEQDSK'
test=input_classes.Equilibrium(ID,input_filename,data_format)

print(test.data) #this should be 10 but is None
print(test.input_filename)
print(test.data_format)

print("Do we get here?\n")

data_format='IDS_equilibrium'
test2=input_classes.Equilibrium(ID,input_filename,data_format)
print(test2.data) #this should be 20
print(test2.input_filename)
print(test2.data_format)

#################################

##################################################################

###################################################################################################