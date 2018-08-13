#context.py

'''
Samuel Ward
13/12/2017
----
Add this file to any LOCUST_IO projects so visibility is automatically handled
---
notes:
    when using a project interactively, just import this file first
---
'''



##################################################################
#Preamble

import sys
import os

##################################################################
#Main

thisdirectory=os.path.dirname(os.path.abspath(__file__)) #get the directory this script is in
thisdirectory_up1=os.path.dirname(thisdirectory) #go one level up from that (should be LOCUST_IO directory)
sys.path.insert(1,thisdirectory_up1) #append Python sys path

#################################

##################################################################

###################################################################################################