#support.py

'''
Samuel Ward
13/12/2017
----
File which holds supporting functions and variables for the LOCUST-IO package
---
notes: 
	TODO this could hold all the read/write functions as well as their supporting functions? 
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

dir_input_files=os.path.join(dir_locust_io,'input_files/') #add the / at the end so the user only needs to append filenames
dir_output_files=os.path.join(dir_locust_io,'output_files/')
dir_classes=os.path.join(dir_locust_io,'classes/')






#software version number 
LOCUST_IO_version=str('LOCUST_IO version 0.99')






#################################

##################################################################

###################################################################################################