#input_classes.py

'''
Samuel Ward
02/11/2017
----
python module for classes which hold LOCUST's input data
contains methods for reading/writing/converting/plotting/manipulating LOCUST input data
---
usage:
    see README.md for usage

notes: 
    TODO GEQDSK read in seems ok, may need to have a look online at certain functions and need testing. now need to check GEQDSK dump and add IDS functionality
    TODO need to decide how to standardise data in dicts
    TODO need to read up on readline() to see if there is some global counter keeps track of which file line the entire script is currently on
    i.e. whether two separate calls to readline() will read the same line or different line due to some global current line tracker. That will help explain
    the file_numbers function somewhat and whether all file lines are held by the thing that it returns when its called in the main code body
    TODO please check how get_next() works and whether it just returns one value at a time (this is what I think)
    TODO integrate with JET SAL API for instant access to JET data (will need error handling for use of this module on systems without JET SAL)

    NOTE use argparse to input command line arguments etc to the python code which calls from input_classes. e.g. pass file names etc via command line
    and argparse to a python script which implements this module
    NOTE use fabric/paramiko for remote host handling stuff https://dtucker.co.uk/hack/ssh-for-python-in-search-of-api-perfection.html (fab and spur both built on paramiko and simplifies it, although more options obvs then with paramiko)
---
'''


###################################################################################################
#Preamble

import sys

try:
    import support #import support module from this directory
    import numpy as np
    import copy
    import re
    import imas
except:
    raise ImportError("ERROR: initial imported modules not found!\nreturning\n") #TODO give each of these its own error
    sys.exit(1)






























###################################################################################################
#Main Code




################################################################## supporting functions
def none_check(ID,input_type,error_message,*args):
    '''
    generic function for checking if a None value appears in *args
    '''

    for arg in args: #loop through and return this error if any element in args is None, otherwise return False 
        if arg == None:
            print("Warning: none_check returned True (input_type={input_type}, ID={ID}): {message}".format(input_type=input_type,ID=ID,message=error_message))
            return True
        
    return False



def file_numbers(ingf):#instance of generators are objects hence can use .next() method on them
    '''
    generator to read numbers in a file, originally written by Ben Dudson
    '''
    toklist = []
    while True:#runs until break statement below (i.e. until what is read is not a line)
        line = ingf.readline() #read in line from INGF
        if not line: break #break if readline() doesnt return a line i.e. end of file
        line = line.replace("NaN","-0.00000e0") #replaces NaNs with 0s
        pattern = r'[+-]?\d*[\.]?\d+(?:[Ee][+-]?\d+)?' #regular expression to find numbers
        toklist = re.findall(pattern,line) #toklist now holds all the numbers in a file line (is essentially a single file line)
        for tok in toklist: #runs whilst True and so yields every line in a file
            yield tok #tok is each number in a file line



def get_next(obj):
    '''
    generic object iterator
    '''
    pyVer = sys.version_info[0]
    if pyVer == 2:
        return obj.next() #XX how exactly does .next() work? will it accept tok in toklist values from file_numbers() and allow file_numbers to then read in the next bunch (since yield will pause the While True loop?)
    else:
        return next(obj)



def read_GEQDSK(ID,input_filepath): 
    ''' 
    read a G-EQDSK formatted equilibrium file, originally written by Ben Dudson and edited by Nick Walkden

    0D data
        nh          #number of points in R (x or width)
        nw          #number of points in Z (y or height)
        idum        #number of spatial dimensions?
        rdim        #size of the R dimension in m
        zdim        #size of the Z dimension in m
        rcentr      #reference value of R
        bcentr      #vacuum toroidal magnetic field at rcentr
        rleft       #R at left (inner) boundary
        zmid        #Z at middle of domain
        rmaxis      #R at magnetic axis (O-point)
        zmaxis      #Z at magnetic axis (O-point)
        simag       #poloidal flux psi at magnetic axis (Weber / rad)
        sibry       #poloidal flux at plasma boundary (Weber / rad)
        current     #plasma current [Amps]   
        xdum        #dummy variable - just contains zero
        nbbbs       #plasma boundary
        limitr      #wall boundary
    1D data
        fpol        #poloidal current function on uniform flux grid (1D array of f(psi)=R*Bt  [meter-Tesla])
        pres        #plasma pressure in nt/m^2 on uniform flux grid (1D array of p(psi) [Pascals])
        ffprime     #workk1
        pprime      #workk1
        qpsi        #q values on uniform flux grid
        rlim        #r wall boundary
        zlim        #z wall boundary
        rbbbs       #r plasma boundary
        zbbbs       #z plasma boundary
    2D data
        psirz       #array (nx,ny) of poloidal flux (array of arrays)   

    '''

    input_data = {}
    flags = {'loaded' : False } #XXX might not need this variable now

    file = open(input_filepath) #open file, assumed input_filepath is object
    
    line = file.readline() #first line should be case, id number and dimensions
    if not line:
        raise IOError("ERROR: Cannot read from file"+input_filepath)
    
    #extract case, id number and dimensions  
    conts = line.split()    #split by white space (no arguement in .split())
    input_data['nh'] = int(conts[-1]) #same as nyefit or height dimension
    input_data['nw'] = int(conts[-2]) #same as nxefit or width dimension
    input_data['idum'] = int(conts[-3])
    flags['case'] = conts[0:-4]

    #now use generator to read numbers from remaining lines in file
    token = file_numbers(file) #token now holds all lines in the file, containing all values. each get_next() call will grab the next number in a line (since get_next() returns the next value in the yield loop? check this plz)
    
    float_keys = [
    'rdim','zdim','rcentr','rleft','zmid',
    'rmaxis','zmaxis','simag','sibry','bcentr',
    'current','simag','xdum','rmaxis','xdum',
    'zmaxis','xdum','sibry','xdum','xdum']
    
    #read in all floats
    for key in float_keys:                              
        input_data[key] = float(get_next(token)) #these are all single value floats
    
    #now read arrays  
    def read_1d(n): 
        input_data = np.zeros([n]) #initialise blank lists
        for i in np.arange(n): #instead of using linspace or something makes a temporary numpy array of dimension n to iterate through
            input_data[i] = float(get_next(token))
        return input_data

    def read_2d(nx,ny):
        input_data = np.zeros([nx,ny])
        for i in np.arange(nx): #same technique here, iterate through one dimension and call read_1d along the other dimension
            input_data[i,:] = read_1d(ny)
        return input_data
    
    input_data['fpol'] = read_1d(input_data['nw']) #remember input_data['nw'] holds width, input_data['nh'] holds height
    input_data['pres'] = read_1d(input_data['nw'])
    input_data['ffprime'] = read_1d(input_data['nw'])
    input_data['pprime'] = read_1d(input_data['nw'])
    input_data['psirz'] = read_2d(input_data['nw'],input_data['nh'])
    input_data['qpsi'] = read_1d(input_data['nw'])

    #now deal with boundaries
    input_data['nbbbs'] = int(get_next(token))
    input_data['limitr'] = int(get_next(token))

    def read_bndy(nb,nl): #number of boundaries and limiters
        if nb > 0: #read in boundaries
            rb = np.zeros(nb)
            zb = np.zeros(nb)
            for i in np.arange(nb):
                rb[i] = float(get_next(token)) #read in R,Z pairs
                zb[i] = float(get_next(token))
        else:
            rb = [0]
            zb = [0]
    
        if nl > 0: #read in limiters
            rl = np.zeros(nl)
            zl = np.zeros(nl)
            for i in np.arange(nl):
                rl[i] = float(get_next(token))
                zl[i] = float(get_next(token))
        else:
            rl = [0]
            zl = [0]

        return rb,zb,rl,zl

    input_data['rbbbs'],input_data['zbbbs'],input_data['rlim'],input_data['zlim'] = read_bndy(input_data['nbbbs'],input_data['limitr'])
    flags['loaded'] = True
    
    return input_data



def dump_GEQDSK(ID,data,output_filepath):
    '''
    generic function for writing a GEQDSK to file, originally written by Ben Dudson and edited by Nick Walkden

    0D data
        nh          #number of points in R (x or width)
        nw          #number of points in Z (y or height)
        idum        #number of spatial dimensions?
        rdim        #size of the R dimension in m
        zdim        #size of the Z dimension in m
        rcentr      #reference value of R
        bcentr      #vacuum toroidal magnetic field at rcentr
        rleft       #R at left (inner) boundary
        zmid        #Z at middle of domain
        rmaxis      #R at magnetic axis (O-point)
        zmaxis      #Z at magnetic axis (O-point)
        simag       #poloidal flux psi at magnetic axis (Weber / rad)
        sibry       #poloidal flux at plasma boundary (Weber / rad)
        current     #plasma current [Amps]   
        xdum        #dummy variable - just contains zero
        nbbbs       #plasma boundary
        limitr      #wall boundary
    1D data
        fpol        #poloidal current function on uniform flux grid (1D array of f(psi)=R*Bt  [meter-Tesla])
        pres        #plasma pressure in nt/m^2 on uniform flux grid (1D array of p(psi) [Pascals])
        ffprime     #workk1
        pprime      #workk1
        qpsi        #q values on uniform flux grid
        rlim        #r wall boundary
        zlim        #z wall boundary
        rbbbs       #r plasma boundary
        zbbbs       #z plasma boundary
    2D data
        psirz       #array (nx,ny) of poloidal flux (array of arrays)   
    '''

    import time #TODO check this
    from itertools import cycle
    cnt = cycle([0,1,2,3,4]) #counter

    def write_number(file,number,counter):
        if number < 0:
            seperator = "-"
            number = np.abs(number)
        else:
            seperator = " "
        if get_next(counter) == 4:
            last = "\n"
        else:
            last = ""
        
        string = '%.10E'%number
        #mant,exp = string.split('E')
        file.write(seperator+string+last)

    def write_1d(file,array,counter):
        for num in array:
            write_number(file,num,counter)

    def write_2d(file,array,counter):
        nx = array.shape[0]
        for i in np.arange(nx):
            write_1d(file,array[i],counter)
    
    def write_bndry(file,R,Z,counter):
        for i in np.arange(len(list(R))):
            write_number(file,R[i],counter)
            write_number(file,Z[i],counter)
        file.write("\n")
    
    with open(output_filepath,'w') as file:
        line = " pyEquilibrium "+time.strftime("%d/%m/%Y")+" # 0 0 "+str(output_data['nw'])+" "+str(output_data['nh'])+"\n"
        file.write(line)

        float_keys = [
        'rdim','zdim','rcentr','rleft','zmid',
        'rmaxis','zmaxis','simag','sibry','bcentr',
        'current','simag','xdum','rmaxis','xdum',
        'zmaxis','xdum','sibry','xdum','xdum']
        for key in float_keys:
            write_number(file,output_data[key],cnt)
        write_1d(file,output_data['fpol'],cnt)
        write_1d(file,output_data['pres'],cnt)
        write_1d(file,output_data['ffprime'],cnt)
        write_1d(file,output_data['pprime'],cnt)
        write_2d(file,output_data['psirz'],cnt)
        #No qpsi in eq object for now
        write_1d(file,np.zeros(output_data['nw']),cnt)    
        file.write("\n"+str(len(list(output_data['rbbbs'])))+"\t"+str(len(list(output_data['rlim'])))+"\n")
        write_bndry(file,output_data['rbbbs'],output_data['zbbbs'],cnt)
        write_bndry(file,output_data['rlim'],output_data['zlim'],cnt)



def read_IDS_equilibrium(ID,shot,run): 
    '''
    generic function for reading an IDS equilibrium file and returning it as a dictionary
    '''

    input_IDS=imas.ids(shot,run) #initialise new blank IDS
    input_IDS.create() 
 
    return data



def dump_IDS_equilibrium(ID,data,shot,run):
    '''
    generic function for reading a GEQDSK file and returning as a dictionary
    '''

    output_IDS=imas.ids(shot,run)
    output_IDS.equilibrium.put() #finally "put" the data
    output_IDS.close()





























################################################################## base class
class LOCUST_input:
    """
    base class for a LOCUST input object

    self.input_type - string which holds this class' input type, this case = 'equilibrium'
    self.ID - unique object identifier, good convention to fill these for error handling etc
    self.input_filename - name of file in input_files folder
    self.input_filepath - full path of file in input_files folder  
    self.data_format - data format of input_file e.g. G-EQDSK
    self.data - holds all input data in dictionary object

    notes:
    """

    input_type='base_input'

    def __init__(self,ID,input_filename=None,data_format=None,*args,**kwargs): #this is common to all children (not overloaded)

        self.ID=ID
        self.data_format=data_format
        self.input_filename=input_filename
        self.input_filepath=support.dir_input_files+input_filename #prepend the path to the file as calculated within the support module

        if none_check(self.ID,self.input_type,'blank LOCUST_input initialised - input_filename or data_format is missing\n',self.input_filename,self.data_format): #check to make sure input_filename and data_format are specified, if either are blank then do nothing, can still read later in program using read_data() directly
            pass #TODO check whether this generates blank object
        else: #read input data if sufficient arguements are given
            self.read_data(self.input_filename,self.data_format)

    def read_data(self,input_filename=None,data_format=None): #bad practice to change overridden method signatures, so retain all method arguements             
        self.data=None #read_data definition is to be overloaded in children classes 




























################################################################## main classes
class Equilibrium(LOCUST_input):
    """
    class describing the equilibrium input for LOCUST

    inherited from LOCUST_input:
        self.input_type - string which holds this class' input type, this case = 'equilibrium'
        self.ID - unique object identifier, good convention to fill these for error handling etc
        self.input_filename - name of file in input_files folder
        self.input_filepath - full path of file in input_files folder  
        self.data_format - data format of input_file e.g. G-EQDSK
        self.data - holds all input data in dictionary object
    Equilibrium data
        key - key in data dictionary to specify data entries
        target - external object to copy from

    notes:
    """

    input_type='equilibrium'

    def __enter__(self,*args,**kwargs):
        return self

    def __exit__(self,*args,**kwargs):
        pass

    def read_data(self,input_filename=None,data_format=None):
        """
        read equilibrium from file 

        notes:
        """

        if not none_check(self.ID,self.input_type,'cannot read_data, both input_filename and data_format required\n',input_filename,data_format): #check to make sure input_filename and data_format are specified

            self.data_format=data_format
            self.input_filename=input_filename
            self.input_filepath=support.dir_input_files+input_filename

            if self.data_format=='GEQDSK':
                self.data=read_GEQDSK(self.ID,self.input_filepath)
               
            elif self.data_format=='IDS_equilibrium':
                self.data=read_IDS_equilibrium(self.ID,self.input_filepath)

            
    def dump_data(self,output_filename="test.eqdsk",data_format='GEQDSK'):
        """
        write data to file

        output_filename - name of file

        notes: 
        """

        if none_check(self.ID,self.input_type,'cannot dump_data, data=None\n',self.data): #check to make sure input_filename and data_format are specified
            pass
        else:
            if data_format=='GEQDSK':
                dump_GEQDSK(self.ID,self.data,output_filename)

            elif data_format=='IDS':
                dump_IDS_equilibrium(self.ID,self.data,output_filename)


    def copy(self,target,key=None):
        """
        copy two equilibrium objects 

        notes:
            -if target data is empty this function does nothing
            -if no key supplied then copy everything
            -if key supplied then copy/append dictionary data accordingly
        """

        if target.data == None:    #return warning if target is empty
            raise RuntimeError("Warning (ID=%s): tried to copy from blank equilibrium\n" % self.ID) #TODO need to python3 this + print statement here?
        elif key == None: #if no key supplied then copy everything
            self.data_format=target.data_format
            self.input_filename=target.input_filename
            self.data=target.data
        else: #if key supplied then overwrite self.input_filename to show data has been edited
            self.input_filename='MIXED' #don't append here because many copy() calls will create huge input_filename
            self.data[key]=target.data[key]


    def set(self,key=None,value=None):
        """
        set equilibrium object data 

        notes:
        """
        
        if none_check(self.ID,self.input_type,'tried to call set() with empty key or value\n',key,value):
            pass
        else:
            self.data[key]=value






























#################################
'''
def Ptcles(LOCUST_input):
    """
    Neutral beam deposition profile

    NOTE this needs to store in the "distribution_sources" IDS which stores marker information
    """
    class_level_attribute=attribute_here 

    def __init__(self,*args,**kwargs):
        instance data, methods etc
    def __enter__(self,*args,**kwargs):
        some_things
    def __exit__(self,*args,**kwargs):
        some_things
    def class_methods(self,*args,**kwargs):
        some_things
'''


#################################
'''
def T(LOCUST_input):
    """
    Species temperature as a function of Psi
    """
    class_level_attribute=attribute_here 

    def __init__(self,*args,**kwargs):
        instance data, methods etc
    def __enter__(self,*args,**kwargs):
        some_things
    def __exit__(self,*args,**kwargs):
        some_things
    def class_methods(self,*args,**kwargs):
        some_things
'''

#################################
'''
def N(LOCUST_input):
    """
    Particle number density as a function of Psi
    """
    class_level_attribute=attribute_here 

    def __init__(self,*args,**kwargs):
        instance data, methods etc
    def __enter__(self,*args,**kwargs):
        some_things
    def __exit__(self,*args,**kwargs):
        some_things
    def class_methods(self,*args,**kwargs):
        some_things
'''






























#NOTE post-processing IDL script for LOCUST to be re-written in Python?

'''
pro locust_orplot

device, decompose=0

file    = '/tmp/sward/locust.AUG/OutputFiles/ORBIT.dat'

n       = 0L
npoints = 0

SPAWN, 'tail -1 '+file, STR, /SH

m       = fix(STR[0],type=3)

openr, 1, file

readf, 1, n

DAT = fltarr(n,3,m)
DUM = fltarr(n,3)

for i=0,m-1 do begin
    readf, 1, DUM

    DAT[*,*,i]=DUM
endfor

close, 1

i = where(DAT[*,0,*] gt 0.0)

!X.range=[min(DAT[*,0,*]),max(DAT[*,0,*])]
!Y.range=[min(DAT[*,2,*]),max(DAT[*,2,*])]

plot, [0.,1.], [0.,1.], /iso, /nodata, xtitle='R [m]', ytitle='Z [m]'

for i=0,n-1 do begin
    j = where( (DAT[i,0,*] eq SHIFT(DAT[i,0,*],-1)) and                       $
               (DAT[i,1,*] eq SHIFT(DAT[i,1,*],-1)) and                       $
               (DAT[i,2,*] eq SHIFT(DAT[i,2,*],-1)) ) & j=j[0]

;-  CHECK THIS FOR ORBITS THAT WORKED....

    if( j eq -1 )then j = m-1

    print, DAT[i,0,0], DAT[i,0,j], j

    npoints = npoints + j+1
    print, i, npoints
    oplot, DAT[i,0,0:j], DAT[i,2,0:j], color=7*16-1
endfor

COL  = fltarr(npoints)
COL2 = fltarr(npoints)

openW,  1, '/tmp/sward/ORBITS.vtk'
printF, 1, '# vtk DataFile Version 2.0'
printF, 1, 'Unstructured Grid ORBIT data'
printF, 1, 'ASCII'
printF, 1, 'DATASET UNSTRUCTURED_GRID'
printF, 1, 'POINTS '+strcompress(string(npoints),/re)+' float'

for i=0,n-1 do begin

    j = where( (DAT[i,0,*] eq SHIFT(DAT[i,0,*],-1)) and                       $
               (DAT[i,1,*] eq SHIFT(DAT[i,1,*],-1)) and                       $
               (DAT[i,2,*] eq SHIFT(DAT[i,2,*],-1)) ) & j=j[0]

;- CHECK

    if( j eq -1 )then j=m-1


    for k=0,j do begin
        printF, 1, DAT[i,0,k]*cos(DAT[i,1,k]), DAT[i,0,k]*sin(DAT[i,1,k]),    $
                   DAT[i,2,k]
    endfor
endfor

printF, 1, 'CELLS '+strcompress(string(n),/re)+' '+                           $
                    strcompress(string(npoints+n),/re)
ioff = 0

for i=0,n-1 do begin

    j = where( (DAT[i,0,*] eq SHIFT(DAT[i,0,*],-1)) and                       $
               (DAT[i,1,*] eq SHIFT(DAT[i,1,*],-1)) and                       $
               (DAT[i,2,*] eq SHIFT(DAT[i,2,*],-1)) ) & j=j[0]
;- CHECK

    if( j eq -1 )then j=m-1

    printF, 1, j+1, indgen(j+1)+ioff
    COL [indgen(j+1)+ioff] = findgen(j+1) + 1
    COL2[indgen(j+1)+ioff] = DAT[i,0,0]
    ioff = ioff + j+1
endfor

printF, 1, 'CELL_TYPES '+strcompress(string(n),/re)

for i=0,n-1 do begin
    printF, 1, 4
endfor

printF, 1, 'POINT_DATA '+strcompress(string(npoints),/re)
printF, 1, 'SCALARS INDX float 2'
printF, 1, 'LOOKUP_TABLE default'

for i=0,n_elements(COL)-1 do begin
    printF, 1, COL[i], COL2[i]
endfor

close, 1

stop
end

'''

#################################

##################################################################

###################################################################################################