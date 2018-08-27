#orbits.py

"""
Samuel Ward
15/01/2018
----
class to handle LOCUST orbit output data
---
usage:
    see README.md for usage

notes:         
---
"""


###################################################################################################
#Preamble

import sys #have global imports --> makes less modular (no "from output_classes import x") but best practice to import whole output_classes module anyway
try:
    import numpy as np
    import copy
    import re
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)
try:
    import imas 
except:
    print("WARNING: IMAS module could not be imported!\n")
try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/processing/ could not be imported!\nreturning\n")
    sys.exit(1)  
try:
    from classes import base_output 
except:
    raise ImportError("ERROR: base_output.py could not be imported!\nreturning\n")
    sys.exit(1) 
try:
    from classes import support #import support module from this directory
except:
    raise ImportError("ERROR: support.py could not be imported!\nreturning\n") 
    sys.exit(1)


np.set_printoptions(precision=5,threshold=3) #set printing style of numpy arrays

################################################################## Orbits functions

def read_orbits_LOCUST(filepath):
    """
    reads orbits stored in LOCUST format - r phi z

    notes:
        reads in a headerline for number of particles
        reads in a footerline for number of time steps
    """

    print("reading orbits from LOCUST")

    with open(filepath) as file:
        
        lines=file.readlines() #return lines as list
        if not lines: #check to see if the file opened
            raise IOError("ERROR: read_orbits_LOCUST() cannot read from "+filepath)
    
        number_particles=int(lines[0]) #extract number of particles
        number_timesteps=int(lines[-1])-1 #extract number of time steps of each trajectory
        number_coords=int(3)

        del(lines[0])
        del(lines[-1])

        input_data = {} #initialise the dictionary to hold the dat
        input_data['orbits']=np.array([[float(value) for value in line.split()] for line in lines]).reshape((number_timesteps+1,number_coords,number_particles)) #read all data and reshape accordingly
        input_data['number_particles']=np.asarray(number_particles)
        input_data['number_timesteps']=np.asarray(number_timesteps)
   
    print("finished reading orbits from LOCUST")

    return input_data


'''
def dump_orbits_LOCUST(output_data,filepath): 
    """
    writes orbits to LOCUST format - r phi z
    
    notes:
        writes out a headerline for number of particles
        writes out a footerline for number of time steps
    """

    print("writing orbits to LOCUST")

    with open(filepath,'w') as file: #open file

        file.write("{}\n".format(output_data['number_particles'].size)) #re-insert line containing number of particles


        #some stuff here


        file.write("{}".format(output_data['number_timesteps'].size)) #re-insert line containing number of time steps

    print("finished writing orbits to LOCUST")
'''
'''
def dump_orbits_vtk(output_data,filepath)
    """
    notes:
    
    here just need to print x,y,z \n locations at all times for particle 1 THEN all times for particle 2 etc
    """

    printF, 1, '# vtk DataFile Version 2.0'
    printF, 1, 'Unstructured Grid ORBIT data'
    printF, 1, 'ASCII'
    printF, 1, 'DATASET UNSTRUCTURED_GRID'
    printF, 1, 'POINTS '+strcompress(string(npoints),/re)+' float'
'''

################################################################## Orbits class

class Orbits(base_output.LOCUST_output):
    """
    class describing orbits output for LOCUST
    
    inheritedfrom LOCUST_output:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all output data in dictionary object
        self.LOCUST_output_type     string which holds this class' output type, this case = 'orbits'
    class data
        self.data_format            data format of original data e.g. LOCUST
        self.filename               name of file in output_files folder
        self.filepath               full path of file in output_files folder  
        self.properties             data to hold additional class-specific information e.g. ion species
        key, value                  key for data dictionary to specify data entry holding value
        target                      external object to copy from
        filename                    name of file to write to
        filepath                    full path to output file in output_files folder

    notes:
        data is stored such that coordinate i at time t for particle p is my_orbit['orbits'][t,i,p]
        in this way, a single particle's trajectory is my_orbit['orbits'][:,i,p] where i=0,1,2=r,phi,z
    """

    LOCUST_output_type='orbits'

    def read_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        read orbits from file 

        notes:
        """

        if processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot read_data() - data_format required\n",data_format): #must always have data_format if reading in data
            pass

        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot read_data() from LOCUST - filename required\n",filename): #must check we have all info required for reading

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files+filename
                self.properties={**properties}
                self.data=read_orbits_LOCUST(self.filepath) #read the file
        else:
            print("ERROR: cannot read_data() - please specify a compatible data_format (LOCUST)\n")            

    def dump_data(self,data_format=None,filename=None,shot=None,run=None):
        """
        write orbits to file

        notes: 
        """
        if processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot dump_data() - self.data and compatible data_format required\n",self.data,data_format):
            pass
        
        elif data_format=='LOCUST':
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: cannot dump_data() to LOCUST - filename required\n",filename):
                filepath=support.dir_output_files+filename
                dump_orbits_LOCUST(self.data,filepath)
        else:
            print("ERROR: cannot dump_data() - please specify a compatible data_format (LOCUST)\n")



#################################

##################################################################

###################################################################################################