#rundata.py

'''
Samuel Ward
04/02/2020
----
class to handle LOCUST output rundata message files
---
usage:
    see README.md for usage

notes:         
---
'''


###################################################################################################
#Preamble

import sys #have global imports --> makes less modular (no "from output_classes import x") but best practice to import whole output_classes module anyway

try:
    import numpy as np
    import pathlib
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)  

try:
    import classes.base_output 
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/base_output.py could not be imported!\nreturning\n")
    sys.exit(1) 

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/src/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import constants
except:
    raise ImportError("ERROR: LOCUST_IO/src/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    import settings
except:
    raise ImportError("ERROR: LOCUST_IO/src/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## Orbits functions

def read_rundata_LOCUST(filepath):
    """
    reads generic rundata file output from LOCUST

    notes:
    """

    print("reading rundata from LOCUST")

    def escape_ansi(line):
        import re
        ansi_escape =re.compile(r'(\x9B|\x1B\[)[0-?]*[ -\/]*[@-~]')
        return ansi_escape.sub('', line)

    file=open(filepath,'r')
    lines=file.readlines()

    rundata={}

    for line in lines: #pre-parse the output data according to calling function, message or variable and value if any
        try:
            split_line=escape_ansi(line).strip('\n').split(':')
            for column_number,column in enumerate(split_line):
                if all([char==' ' for char in column]): del(split_line[column_number])

            calling_function=split_line[0].strip()
            variable_or_message=split_line[1].strip()
            
            if calling_function not in rundata:
                rundata[calling_function]={}

            try:
                values=[value.strip() for value in split_line[2:]]
                if variable_or_message not in rundata[calling_function]:
                    rundata[calling_function][variable_or_message]=[values] 
                else: #if key already occupied then append value
                    rundata[calling_function][variable_or_message].append(values) 
            except:
                rundata[calling_function][variable_or_message]=None
        except:            
            pass

    for calling_function in rundata:
        for variable_or_message in rundata[calling_function]:
            if len(rundata[calling_function][variable_or_message])==1: 
                rundata[calling_function][variable_or_message]=rundata[calling_function][variable_or_message][0]

    input_data={} #initialise data dictionary

    #extract power flux data
    input_data['PFC_power']={}
    input_data['PFC_power']['component']={}

    for message in rundata['write_vtk']:

        if 'Integrated power' in message:
            if 'to component' in message:
                for component_power in rundata['write_vtk'][message]:
                    input_data['PFC_power']['component'][str(component_power[0])]=float(component_power[1][:-2])*1.e6
            elif 'to PFCs' in message:
                #print(message)
                #print(rundata['write_vtk'][message])
                input_data['PFC_power']['total']=float(rundata['write_vtk'][message][0][:-2])*1.e6

    #extract loss channel data
    input_data['loss_channels']={}

    for message in rundata['loss_chk']:
        try:
            input_data['loss_channels'][message]=rundata['loss_chk'][message][0].replace(' ','').replace('MW','')
            if '%' in input_data['loss_channels'][message]:
                input_data['loss_channels'][message]=input_data['loss_channels'][message].replace('%','')
            else:
                input_data['loss_channels'][message]=float(input_data['loss_channels'][message])*1.e6
        except:
            pass

    #look for kernel end time
    input_data['time_total']=None
    input_data['time_total']=float(rundata['ctrk']['End of Kernel wrapper'][0])

    print("finished reading rundata from LOCUST")

    return input_data

################################################################## Orbits class

class Rundata(classes.base_output.LOCUST_output):
    """
    class describing a generic rundata output for LOCUST
    
    inherited from LOCUST_output:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all output data in dictionary object
        self.LOCUST_output_type     string which holds this class' output type, this case = 'rundata'
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
    """

    LOCUST_output_type='rundata'

    def read_data(self,data_format=None,filename=None,**properties):
        """
        read rundata from file 

        notes:
        """

        if processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot read_data() - data_format required\n".format(self.ID),data_format): #must always have data_format if reading in data
            pass

        elif data_format=='LOCUST': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot read_data() from LOCUST - filename required\n".format(self.ID),filename): #must check we have all info required for reading

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files / filename
                self.properties={**properties}
                self.data=read_rundata_LOCUST(self.filepath,**properties)

        else:
            print("ERROR: {} cannot read_data() - please specify a compatible data_format (LOCUST)\n".format(self.ID))            

#################################

##################################################################

###################################################################################################
