#context.py

'''
Samuel Ward
13/12/2017
----
Testing scripts need to import this file in order to append the python path to include other folders 
---
notes:
---
'''



##################################################################
#Preamble

import sys
import os
thisdirectory=os.path.dirname(os.path.abspath(__file__)) #get the directory this script is in
thisdirectory_up1=os.path.dirname() #go one level up from that
#thisdirectory_up1_classes=os.path.join(thisdirectory_up1,'classes/') #add classes/ to the path
sys.path.insert(1,thisdirectory_up1) #add LOCUST_IO directory to python path

#################################

##################################################################

###################################################################################################