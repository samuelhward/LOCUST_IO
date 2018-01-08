#context.py

'''
Samuel Ward
13/12/2017
----
Add this file to any LOCUST_IO/sub_dirs so that context and module importing is automatically handled
---
notes:
	when using a project interactively, just import this module before importing any classes
---
'''



##################################################################
#Preamble

import sys
import os

##################################################################
#Main


#append Python sys path
thisdirectory=os.path.dirname(os.path.abspath(__file__)) #get the directory this script is in
thisdirectory_up1=os.path.dirname(thisdirectory) #go one level up from that
thisdirectory_up1_classes=os.path.join(thisdirectory_up1,'classes') #into classes/ directory
thisdirectory_up1_plotting=os.path.join(thisdirectory_up1,'plotting') #into plotting/ directory
sys.path.insert(1,thisdirectory_up1_classes) #add LOCUST_IO/classes/ directory to python path
sys.path.insert(1,thisdirectory_up1_plotting) #add LOCUST_IO/plotting/ directory to python path

#################################

##################################################################

###################################################################################################