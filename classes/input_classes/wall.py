#wall.py
 
"""
Samuel Ward
04/11/2018
----
class to handle LOCUST wall input data
---
usage:
    see README.md for usage
 
notes:         
---
"""


###################################################################################################
#Preamble
 
import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway

try:
    import numpy as np
    import re
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/processing/utils.py could not be imported!\nreturning\n")
    sys.exit(1)  

try:
    import classes.base_input 
except:
    raise ImportError("ERROR: LOCUST_IO/classes/base_input.py could not be imported!\nreturning\n")
    sys.exit(1) 

try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from constants import *
except:
    raise ImportError("ERROR: LOCUST_IO/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from settings import *
except:
    raise ImportError("ERROR: LOCUST_IO/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## Wall read functions

def read_wall_LOCUST_3D(filepath):
    """
    notes:
    """
    pass

def read_wall_LOCUST_2D(filepath):
    """
    notes:
    """
    pass

def read_wall_GEQDSK(filepath):
    """
    notes:
    """

    buffer_data={}
    input_data={}

    def file_numbers(ingf):#instance of generators are objects hence can use .next() method on them
        """
        generator to read numbers in a file
     
        notes:
            originally written by Ben Dudson
     
            calling get_next() on a generator produced with this function will permanently advance the iterator in this generator
            when generators reach the end, that's permanently it!
        """
        toklist = []
        while True: #runs forever until break statement (until what is read is not a line) below will cause a return - which permanently stops a generator
            line = ingf.readline() #read in line from INGF
            if not line: break #break if readline() doesnt return a line i.e. end of file
            line = line.replace("NaN","-0.00000e0") #replaces NaNs with 0s
            pattern = r'[+-]?\d*[\.]?\d+(?:[Ee][+-]?\d+)?' #regular expression to find numbers
            toklist = re.findall(pattern,line) #toklist now holds all the numbers in a file line (is essentially a single file line)
            for tok in toklist: #yield every number in that one line individually
                yield tok #so in terms of iteration, using get_next() will cause another line to be read

    print("reading 2D wall from GEQDSK")
  
    with open(filepath,'r') as file: #open file
 
        line = file.readline() #first line should be case, id number and dimensions
        if not line:
            raise IOError("ERROR: read_equilibrium_GEQDSK() cannot read from "+filepath)
         
        #extract case, id number and dimensions  
        conts = line.split()    #split by white space (no argument in .split())
        buffer_data['nZ_1D'] = np.asarray(int(conts[-1])) #same as nyefit or height dimension
        buffer_data['nR_1D'] = np.asarray(int(conts[-2])) #same as nxefit or width dimension
        buffer_data['idum'] = np.asarray(int(conts[-3]))
        #buffer_data['time'] = np.asarray(int(consts[-4])) #time slice time (in ms)
        #buffer_data['EFITshot'] = np.asarray(int(consts[-5])) #shot number

        flags = {} #read the case in (basically just ID data)
        flags['case'] = conts[0:-4] 
     
        #now use generator to read numbers from remaining lines in file
        token = file_numbers(file) #token now holds all lines in the file, containing all values. each processing.utils.get_next() call will grab the next number in a line (since processing.utils.get_next() returns the next value in the yield loop? check this plz)
         
        float_keys = [
        'rdim','zdim','rcentr','rleft','zmid',
        'rmaxis','zmaxis','simag','sibry','bcentr',
        'current','simag','xdum','rmaxis','xdum',
        'zmaxis','xdum','sibry','xdum','xdum']
     
        #read in all 0D floats
        for key in float_keys:                              
            buffer_data[key] = np.asarray(float(processing.utils.get_next(token))) #processing.utils.get_next(token) always yields just a single value, convert this to numpy array
     
        def read_1d(n): 
            """
            """
            buffer_data = np.zeros(n) #initialise blank lists
            for i in np.arange(n): #instead of using linspace or something makes a temporary numpy array of dimension n to iterate through
                buffer_data[i] = float(processing.utils.get_next(token))
            return buffer_data
     
        def read_2d(nx,ny):
            """
            notes:
                edited this function as it was not consistent with GEQDSK storage format
            """
            buffer_data = np.zeros((nx,ny))
            for j in np.arange(ny): #same technique here, iterate through one dimension and call read_1d along the other dimension
                buffer_data[:,j] = read_1d(nx)
            return buffer_data

        #read in the arrays
        buffer_data['fpol'] = read_1d(buffer_data['nR_1D']) #remember buffer_data['nR_1D'] holds width, buffer_data['nZ_1D'] holds height
        buffer_data['pres'] = read_1d(buffer_data['nR_1D'])
        buffer_data['ffprime'] = read_1d(buffer_data['nR_1D'])
        buffer_data['pprime'] = read_1d(buffer_data['nR_1D'])
        buffer_data['psirz'] = read_2d(buffer_data['nR_1D'],buffer_data['nZ_1D'])
        buffer_data['qpsi'] = read_1d(buffer_data['nR_1D'])
     
        #now deal with boundaries
        buffer_data['lcfs_n'] = np.array(int(processing.utils.get_next(token)))
        input_data['limitr'] = np.array(int(processing.utils.get_next(token)))
     
        def read_bndy(nb,nl): #number of boundaries and limiters
            if nb > 0: #read in boundaries
                rb = np.zeros(nb)
                zb = np.zeros(nb)
                for i in np.arange(nb):
                    rb[i] = float(processing.utils.get_next(token)) #read in R,Z pairs
                    zb[i] = float(processing.utils.get_next(token))
            else:
                rb = np.array(0)
                zb = np.array(0)
         
            if nl > 0: #read in limiters
                rl = np.zeros(nl)
                zl = np.zeros(nl)
                for i in np.arange(nl):
                    rl[i] = float(processing.utils.get_next(token))
                    zl[i] = float(processing.utils.get_next(token))
            else:
                rl = np.array(0)
                zl = np.array(0)
     
            return rb,zb,rl,zl
     
        buffer_data['lcfs_r'],buffer_data['lcfs_z'],input_data['rlim'],input_data['zlim'] = read_bndy(buffer_data['lcfs_n'],input_data['limitr'])

        print("finished reading 2D wall from GEQDSK")

    return input_data

def read_wall_ufile(filepath):
    """
    reads 2D wall limiter profile from TRANSP Ufile format

    notes:

    """

    print("reading 2D wall from TRANSP Ufile")

    with open(filepath,'r') as file:

        lines=file.readlines()
        for counter,line in enumerate(lines):
            if ';-# OF X0 PTS-' in line:
                number_points=int(line.split()[0])
                del(lines[0:counter+1])
                break

        input_data={}
        rlim=[]
        zlim=[]
        i=0
        for line in lines:
            split_line=line.split()
            for number in split_line:
                if i<number_points:
                    rlim+=[float(number)]
                    i+=1
                elif i<2.*number_points:
                    zlim+=[float(number)]
                    i+=1
                else:
                    break

        input_data['rlim']=np.asarray(rlim)
        input_data['zlim']=np.asarray(zlim)
        input_data['limitr']=np.asarray(len(input_data['rlim']))

    print("finished reading 2D wall from TRANSP Ufile")

    return input_data 

################################################################## Wall write functions

def dump_wall_LOCUST_2D(output_data,filepath):
    """
    dumps 2D wall limiter profile to LOCUST format

    notes:
        interpolates given R,Z wall points into 3600 equally-spaced poloidal points
    """

    print("writing 2D limiter wall to LOCUST format")

    with open(filepath,'w') as file:

        file.write("{} \n".format(processing.utils.fortran_string(3600,12))) #write the number of points (must be 3600 always)
        R_0=np.mean(output_data['rlim']) #write the major radius centre of the points
        file.write("{}\n".format(processing.utils.fortran_string(R_0,13,5,exponential=False)))

        #calculate angle and radius of limiter points with respect to geometric centre 
        Z_0=np.mean(output_data['zlim']) #should be roughly 0
        angles=np.arctan2(output_data['zlim']-Z_0,output_data['rlim']-R_0)
        angles[angles<0]=2.*pi+angles[angles<0]
        radii=(output_data['zlim']-Z_0)/np.sin(angles)
        angles,radii=processing.utils.sort_arrays(angles,radii) #now sort all arrays together in order of angle to stop any over wrapping

        angles=np.append(-1.*(2.*pi-angles[-1]),angles) #add last value to start to make periodic for interpolator 
        radii=np.append(radii[-1],radii)

        angles=np.append(angles,angles[1]+2.*pi) #add original first value to end
        radii=np.append(radii,radii[1])

        #now interpolate in polar coordinates onto 3600 points in theta, starting anti-clockwise from outboard side 
        angles_new=np.linspace(0.,(3599./3600.)*2.*pi,3600)

        radii_interpolator=processing.utils.interpolate_1D(angles,radii,function='linear',smooth=0.0001)
        radii_new=radii_interpolator(angles_new)

        for radius in radii_new: 
            file.write("{}\n".format(processing.utils.fortran_string(radius**2,13,7,exponential=False)))

    print("finished writing 2D limiter wall to LOCUST format")

def dump_wall_ASCOT(output_data,filepath):
    """
    dumps 2D wall outline to ASCOT input.wall_2d format
    
    notes:
        currently sets all divertor flags = 0 i.e. treats everything as 'wall'
    args:
        filename - output filename
        output_data - data structure holding wall with same wall variable names as GEQDSK equilibrium e.g. rlim,zlim - i.e. a GEQDSK equilibrium object 

    """

    print("dumping wall to ASCOT format")

    with open(filepath,'w') as file:

        file.write("{number_points} (R,z) wall points & divertor flag (1 = divertor, 0 = wall)\n".format(number_points=int(output_data['rlim'].size)))
        
        for r,z in zip(output_data['rlim'],output_data['zlim']):
            line=processing.utils.fortran_string(r,16,7)+processing.utils.fortran_string(z,16,7)+processing.utils.fortran_string(0.0,4,0,False)+"\n"
            file.write(line)

    print("finished dumping wall to ASCOT format")
 
################################################################## wall class
 
class Wall(classes.base_input.LOCUST_input):
    """
    class describing tokamak boundary wall for LOCUST
 
    inherited from LOCUST_input:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all input data in dictionary object
        self.LOCUST_input_type      string which holds this class' input type, this case = 'wall'
    class data
        self.data_format            data format of original data e.g. LOCUST
        self.filename               name of file in input_files folder
        self.filepath               full path of file in input_files folder  
        self.shot                   shot number
        self.run                    run number
        self.properties             data to hold additional class-specific information e.g. ion species
        key, value                  key for data dictionary to specify data entry holding value
        target                      external object to copy from
        filename                    name of file to write to
        filepath                    full path to output file in input_files folder
 
    notes:
    """
 
    LOCUST_input_type='wall'
 
    def read_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        read wall from file 
 
        notes:
        """
 
        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() - data_format required\n",data_format): #must always have data_format if reading in data
            pass
 
        elif data_format=='LOCUST_3D': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from LOCUST_3D - filename required\n",filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_wall_LOCUST_3D(self.filepath) #read the file

        elif data_format=='LOCUST_2D':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from LOCUST_2D - filename required\n",filename): #must check we have all info required for reading
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_wall_LOCUST_2D(self.filepath) #read the file

        elif data_format=='GEQDSK':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from GEQDSK - filename required\n",filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_wall_GEQDSK(self.filepath)

        elif data_format=='ufile':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot read_data() from Ufile - filename required\n",filename):
 
                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_input_files+filename
                self.properties={**properties}
                self.data=read_wall_ufile(self.filepath)

        else:
            print("ERROR: cannot read_data() - please specify a compatible data_format (LOCUST_3D/LOCUST_2D/GEQDSK/ufile)\n")            
 
    def dump_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        write wall to file
 
        notes: 
        """

        if not self.run_check():
            print("WARNING: run_check() returned false - insufficient data for LOCUST run:"+self.ID)

        if processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() - self.data and compatible data_format required\n",self.data,data_format):
            pass
         
        elif data_format=='LOCUST_2D':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to LOCUST_2D - filename required\n",filename):
                filepath=support.dir_input_files+filename
                dump_wall_LOCUST_2D(self.data,filepath)

        elif data_format=='ASCOT':
            if not processing.utils.none_check(self.ID,self.LOCUST_input_type,"ERROR: cannot dump_data() to ASCOT - filename required\n",filename):
                filepath=support.dir_input_files+filename
                dump_wall_ASCOT(self.data,filepath)                

        else:
            print("ERROR: cannot dump_data() - please specify a compatible data_format (LOCUST_2D/ASCOT)\n")


    def plot(self,some_equilibrium=False,key='Br',axes=['R','Z'],LCFS=False,limiters=False,real_scale=False,colmap=cmap_default,number_bins=20,fill=True,ax=False,fig=False): 
        """
        notes:
        args:
            some_equilibrium - corresponding equilibrium for plotting plasma boundary, scaled axes etc.
            key - select which data to plot
            axes - define plot axes in x,y order or as full list of indices/slices (see dfn_transform())
            LCFS - show plasma boundary outline (requires equilibrium arguement)
            limiters - toggles limiters on/off in 2D plots
            real_scale - plot to Tokamak scale
            colmap - set the colour map (use get_cmap names)
            transform - set to False if supplied dfn has already been cut down to correct dimensions
            number_bins - set number of bins or levels
            fill - toggle contour fill on 2D plots
            ax - take input axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        """

        import scipy
        import numpy as np
        import matplotlib
        from matplotlib import cm
        import matplotlib.pyplot as plt
        from mpl_toolkits import mplot3d #import 3D plotting axes
        from mpl_toolkits.mplot3d import Axes3D

        pass

#################################
 
##################################################################
 
###################################################################################################