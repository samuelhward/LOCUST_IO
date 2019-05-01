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
import pathlib

##################################################################
#Main

thisdirectory=pathlib.Path.cwd()
parts=list(thisdirectory.parts)
parts.reverse()
for counter,level in enumerate(parts): #find LOCUST_IO directory by looking for licence, then if user deletes licence code will not work
    if 'LOCUST_IO_LICENCE.md' in [str(list(path.parts)[-1]) for path in thisdirectory.parents[counter].glob('*')]:
        LOCUST_IO_dir=thisdirectory.parents[counter]
        sys.path.insert(1,str(LOCUST_IO_dir)) #append Python sys path
        break
#################################

##################################################################

###################################################################################################