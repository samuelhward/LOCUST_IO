#support.py

'''
Samuel Ward
13/12/2017
----
File which holds supporting functions and variables for the LOCUST-IO package
---
notes: 
	TODO this could hold all non-class functions? 
	TODO could also have another file called read_write functions to hold all the file methods too
---
'''


##################################################################
#Preamble
import sys
import os


##################################################################
#Project directory paths

#get the paths to folders within the package for importing modules, opening files etc
pwd=os.path.dirname(os.path.abspath(__file__)) #get the directory this script is in

dir_locust_io=os.path.dirname(pwd) #go one level up from that

dir_input_files=os.path.join(dir_locust_io,'input_files/') #add the / at the end so the user can just append with filename strings
dir_output_files=os.path.join(dir_locust_io,'output_files/')
dir_classes=os.path.join(dir_locust_io,'classes/')





#sys.path.insert(1,dir_classes) #add that to python path (must specify the exact module directory, python does not search recursively)

#################################

##################################################################

###################################################################################################