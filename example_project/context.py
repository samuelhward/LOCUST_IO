#context_updated.py

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

tail=''
thisdirectory=os.path.dirname(os.path.abspath(__file__)) #get the directory this script is in
LOCUST_IO_dir=thisdirectory

while tail!='LOCUST_IO':
	LOCUST_IO_dir=os.path.dirname(LOCUST_IO_dir) #keep moving up levels until LOCUST_IO dir is found
	tail=LOCUST_IO_dir[-9:]
sys.path.insert(1,LOCUST_IO_dir) #append Python sys path

#################################

##################################################################

###################################################################################################