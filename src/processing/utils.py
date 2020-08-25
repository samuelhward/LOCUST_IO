#processing.utils.py
 
'''
Samuel Ward
02/11/2017
----
supporting functions for LOCUST input/output classes
---
usage:
    see README.md for usage
 
notes:         
---
'''


###################################################################################################
#Preamble
 
import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway
try:
    import numpy as np
    import pathlib
    import copy
    import random
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1) 

try:
    import processing.process_input
except:
    raise ImportError("ERROR: LOCUST_IO/src/processing/process_input.py could not be imported!\nreturning\n")

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

################################################################################################### Supporting functions
 
def none_check(ID,LOCUST_input_type,error_message,*args):
    """
    generic function for checking if a None value appears in *args
 
    notes:
        message should be something specific to the section of code which called none_check
    """
    if all(arg is not None for arg in args):
        return False
    else:
        print("WARNING: none_check returned True (LOCUST_input_type={LOCUST_input_type}, ID={ID}): {message}".format(LOCUST_input_type=LOCUST_input_type,ID=ID,message=error_message))
        return True
          
def get_next(obj):
    """
    generic object iterator
 
    notes:
        every time get_next() is called, it will advance the counter one
    """
    pyVer = sys.version_info[0]
    if pyVer == 2:
        return obj.next() #NOTE how exactly does .next() work? will it accept tok in toklist values from file_numbers() and allow file_numbers to then read in the next bunch (since yield will pause the While True loop?)
    else:
        return next(obj)
     
def dict_set(data,**kwargs):
    """
    generalised function (based upon safe_set) to set values in dictionary 'data' after checking source data exists 
    
    notes:
    usage:
        dict_set(data,some_key=some_source,some_other_key=[1,2,3,4]) to set multiple values simultaneously
        dict_set(data,**{'some_key':100,'some_other_key':200}) equally
    """
    keys=kwargs.keys()
    values=kwargs.values()
    
    for key,value in zip(keys,values): #loop through kwargs
        if key is not None and value is not None: #use 'is' since we may want to set data to 0
            data[key]=value

def safe_set(target,source):
    """
    generalised function to set a target value to source if it exists
    
    notes:
    """
    if source is not None: #use 'is' since we may want to set data to 0
        target=source
 
def fortran_string(number_out,length,decimals=None,exponential=True):
    """
    produces a string stream of specified length, padded with space at the start

    notes:
        supposed to be a quick fix to fortran format descriptors
    args:
        number_out - variable holding number to write out 
        length - total number of characters for output string stream including white space 
        decimals - required for non-ints, how many decimal places the number is written to 
        exponential - true for output to exponential format, otherwise formatted like float
    """

    if decimals is not None:
        if exponential is True:
            output_string='{:.{decimals}e}'.format(number_out,decimals=decimals)        
        else:
            output_string='{:.{decimals}f}'.format(number_out,decimals=decimals)        
    else: #assume integer if decimals not specified
        output_string=str(number_out)

    number_spaces=length-len(output_string) #caclulate amount of padding needed from length of the number_out
    if number_spaces>=0:
        string_stream='{}{}'.format(' '*number_spaces,output_string) #construct the string stream
        return string_stream
    else:
        print('ERROR: fortran_string() too many decimal places for requested string length!')

def sort_arrays(main_array,*args):
    """
    sort an arbitrary number of arrays in parallel (*args) according to main_array

    notes:
        assumes all arrays are numpy arrays
    usage:
        to sort ascending according to array x: 
            x_sorted,y_sorted,z_sorted,...=sort_arrays(x,y,z,...)
    """

    sorted_indices=main_array.argsort() #get the order of sorted indices
    
    main_array=main_array[sorted_indices] #convert all the arrays
    returned_arrays=[main_array] #add arrays to list
    for array in args:
        array=array[sorted_indices]
        returned_arrays.append(array)

    return returned_arrays

def angle_pol(R_major,R,Z,Z_major=0.):
    """
    returns poloidal angle

    notes:
    args:
        R_major - major radius at geometric axis
        R - R coordinate at point of interest
        Z - Z coordinate at point of interest
        Z_major - value of Z at geometric axis
    """

    angle=np.arctan2(Z-Z_major,R-R_major)
    if angle<0: angle+=2.*np.pi
    return angle 

def minor_radius(R_major,R,Z,Z_major=0.):
    """
    returns minor radius of R Z points

    notes:
    args:
        R_major - major radius at geometric axis
        R - R coordinate at point of interest
        Z - Z coordinate at point of interest
        Z_major - value of Z at geometric axis
    """

    return np.sqrt((R-R_major)**2+(Z-Z_major)**2) 

def interpolate_2D(X_axis,Y_axis,Z_grid,function='multiquadric',type='RBS',smooth=0,rect_grid=True):
    """
    generate a 2D grid interpolator

    notes:
        keep as separate functions so can freely swap out interpolation method
        RBF - https://stackoverflow.com/questions/37872171/how-can-i-perform-two-dimensional-interpolation-using-scipy
            - high memory overhead, most accurate
        RBS - https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RectBivariateSpline.html#scipy.interpolate.RectBivariateSpline
              https://scipython.com/book/chapter-8-scipy/examples/two-dimensional-interpolation-with-scipyinterpolaterectbivariatespline/
    args:
        X_axis - 1D x-axis
        Y_axis - 1D y-axis
        Z_grid - 2D z-axis
        function - spline function
        type - name of interpolation function to use
        smooth - level of smoothing
        rect_grid - toggle whether X_axis, Y_axis are regular rectangular grid edges or arbitrary cordinates 
    usage:
        my_interpolator=interpolate_2D(X_axis,Y_axis,data)
        interpolated_value=my_interpolator(x,y) #may need to convert to float or take first array value here
    """

    try:
        import scipy.interpolate 
    except:
        raise ImportError("ERROR: interpolate_2D could not import scipy.interpolate module!\nreturning\n")
        return

    if type=='RBF':
        if rect_grid:
            Y_grid,X_grid=np.meshgrid(Y_axis,X_axis) #swap since things are defined r,z 
        else:
            Y_grid,X_grid=Y_axis,X_axis
        interpolator=scipy.interpolate.Rbf(X_grid,Y_grid,Z_grid,function=function,smooth=smooth)
    
    elif type=='RBS':
        interpolator=scipy.interpolate.RectBivariateSpline(X_axis,Y_axis,Z_grid) #normally order is other way in RBS but I have swapped my axes

    return interpolator

def interpolate_1D(X_axis,Y_axis,function='cubic',type='interp1d',smooth=0):
    """
    generate a 1D line interpolator

    notes:
        keep as separate functions so can freely swap out interpolation method
        uses Rbf - https://stackoverflow.com/questions/37872171/how-can-i-perform-two-dimensional-interpolation-using-scipy
    args:
        X_axis - 1D x-axis
        Y_axis - 1D y-axis
        function - spline function
        type - name of interpolation function to use
        smooth - level of smoothing
    usage:
        my_interpolator=interpolate_1D(X_axis,data)
        interpolated_value=my_interpolator(x)        
    """

    try:
        import scipy.interpolate 
    except:
        raise ImportError("ERROR: interpolate_1D could not import scipy.interpolate module!\nreturning\n")
        return

    if type=='RBF':
        interpolator=scipy.interpolate.Rbf(X_axis,Y_axis,function=function,smooth=smooth)
    elif type=='interp1d':
        interpolator=scipy.interpolate.interp1d(X_axis,Y_axis,kind=function)

    return interpolator

def RphiZ_to_XYZ(R,phi,RH=True):
    """
    converts R,phi positions to X,Y 

    args:
        R - R coordinates to convert
        phi - Toroidal angle coordinates to convert
        RH - toggle right-handed R phi Z coordinate system
    notes:
    """

    if RH is not True:
        phi=2.*constants.pi-phi

    X=np.array(R*np.cos(phi))
    Y=np.array(R*np.sin(phi))

    return X,Y

def XYZ_to_RphiZ(X,Y):
    """
    converts X,Y positions to R,phi

    notes:
    """

    phi=np.asarray(np.arctan2(Y,X))
    R=X*np.cos(phi)+Y*np.sin(phi)
    for counter,phi_ in enumerate(phi): 
        if phi_<0: phi[counter]+=2.*np.pi 
    return R,phi

def V_XYZ_to_V_RphiZ(X,Y,V_X,V_Y):
    """
    converts V_X,V_Y to V_R,V_phi given position X,Y

    notes:
    """

    R,phi=XYZ_to_RphiZ(X,Y)
    V_R=V_X*np.cos(phi)+V_Y*np.sin(phi)
    V_phi=-V_X*np.sin(phi)+V_Y*np.cos(phi)

    return V_R,V_phi

def V_RphiZ_to_V_XYZ(phi,V_R,V_phi):
    """
    converts V_R,V_phi to V_X,V_Y at given position phi

    notes:
    """

    V_X=V_R*np.cos(phi)-V_phi*np.sin(phi)
    V_Y=V_R*np.sin(phi)+V_phi*np.cos(phi) 

    return V_X,V_Y

def value_at_RZ(R,Z,quantity,grid):
    """
    generic function to interpolate value of 2D quantity at position R,Z

    args:
        R - list of R coordinates  
        Z - list of Z coordinates
        quantity - 2D quantity
        grid - quantity is stored on rectangular axes defined by grid['R_1D'] and grid['Z_1D']
    notes:
    """

    interpolator=interpolate_2D(grid['R_1D'],grid['Z_1D'],quantity,function='linear')
    value_at_coordinate=[]

    for r,z in zip(R,Z): #do element-wise to avoid implicitly interpolating 2x1D arrays onto a single 2D grid
        rz_point=interpolator(r,z)
        value_at_coordinate.append(np.squeeze(rz_point))

    return np.asarray(value_at_coordinate)

def flux_func_to_RZ(psi,quantity,equilibrium):
    """
    maps 1D flux function onto a 2D RZ grid equilibrium

    notes:
        assumes psi is consistent between psi, quantity and equilibrium i.e all against normalised poloidal flux, or Wb/rad etc.
        psi is measured against equilibrium['psirz']
    args:
        psi - 1D poloidal flux axis 
        quantity - 1D quantity mapped to psi
        equilibrium - equilibrium object with 2D psi grid
    """

    interpolator=interpolate_1D(psi,quantity,type='RBF')
    Z_2D,R_2D=np.meshgrid(equilibrium['Z_1D'],equilibrium['R_1D']) #this way around to order dimensions as [r,z]
    quantity_2D=copy.deepcopy(R_2D)

    for r_index in range(len(equilibrium['R_1D'])): #loop over grid, get psi at each point and feed to flux function interpolator
        for z_index in range(len(equilibrium['Z_1D'])):            
            quantity_2D[r_index,z_index]=interpolator(equilibrium['psirz'][r_index,z_index])

    return quantity_2D

def within_LCFS(R,Z,equilibrium):
    """
    determines whether R Z points are within LCFS

    args:
        R - array of R coordinates at points of interest
        Z - array of Z coordinates at points of interest
        equilibrium - equilibrium object with LCFS
    notes:
        assumes LCFS points are monotonic in poloidal angle
    returns:
        True if inside LCFS
        False if outside LCFS
    """
    distances_point_to_mag_axis=[]
    pol_angle_points=[]
    for r,z in zip(R,Z):
        distances_point_to_mag_axis.append(minor_radius(R=r,R_major=equilibrium['rmaxis'],Z=z,Z_major=equilibrium['zmaxis'])) #translate point to polar coordinates
        pol_angle_points.append(angle_pol(R=r,R_major=equilibrium['rmaxis'],Z=z,Z_major=equilibrium['zmaxis']))

    distances_point_to_mag_axis=np.asarray(distances_point_to_mag_axis)
    pol_angle_points=np.asarray(pol_angle_points)

    minor_rad_LCFS=[] #calculate translate LCFS to polar coordinates 
    pol_angle_LCFS=[] 
    for r,z in zip(equilibrium['lcfs_r'],equilibrium['lcfs_z']): 
        minor_rad_LCFS.append(minor_radius(R=r,R_major=equilibrium['rmaxis'],Z=z,Z_major=equilibrium['zmaxis']))
        pol_angle_LCFS.append(angle_pol(R=r,R_major=equilibrium['rmaxis'],Z=z,Z_major=equilibrium['zmaxis']))
    pol_angle_LCFS,minor_rad_LCFS=np.asarray(pol_angle_LCFS),np.asarray(minor_rad_LCFS)

    pol_angle_LCFS[pol_angle_LCFS<0]=2.*np.pi+pol_angle_LCFS[pol_angle_LCFS<0] #modulo angles
    pol_angle_LCFS,minor_rad_LCFS=sort_arrays(pol_angle_LCFS,minor_rad_LCFS)
    pol_angle_LCFS=np.concatenate(([pol_angle_LCFS[-1]-2.*np.pi],pol_angle_LCFS,[pol_angle_LCFS[0]+2.*np.pi])) #wrap around arrays for values that lie between last and first poloidal angle
    minor_rad_LCFS=np.concatenate(([minor_rad_LCFS[-1]],minor_rad_LCFS,[minor_rad_LCFS[0]]))

    interpolator_LCFS=interpolate_1D(pol_angle_LCFS,minor_rad_LCFS,function='linear',type='interp1d') #generate polar coordinate interpolator of LCFS
    
    distances_LCFS_to_mag_axis=[]
    for pol_angle_point in pol_angle_points:
        distances_LCFS_to_mag_axis.append(interpolator_LCFS(pol_angle_point)) #find radius of LCFS radius at that polar angle 

    return np.array([False if r_point>r_LCFS else True for r_point,r_LCFS in zip(distances_point_to_mag_axis,distances_LCFS_to_mag_axis)]) #is this larger than distance between point and mag axis?

def LCFS_crop(quantity,grid,equilibrium,crop_value=0.0,outside=True):
    """
    sets all values of quantity outside the LCFS = 0
    
    args:
        quantity - arbitrary 2D quantity
        grid - quantity is stored on rectangular axes defined by grid['R_1D'] and grid['Z_1D']
        equilibrium - equilibrium object with LCFS
        crop_value - regions outside LCFS will be set to this value
        outside - if outside then floor all values outside LCFS else floor all values inside
    notes:    
    """

    quantity_cropped=copy.deepcopy(quantity)

    Z_2D,R_2D=np.meshgrid(grid['Z_1D'],grid['R_1D'])    
    Z_2D_flat,R_2D_flat=Z_2D.flatten(),R_2D.flatten()
    nR_1D,nZ_1D=len(grid['R_1D']),len(grid['Z_1D'])
    within=np.array(within_LCFS(R_2D_flat,Z_2D_flat,equilibrium)).reshape(nR_1D,nZ_1D)

    for r_counter,r in enumerate(grid['R_1D']):
        for z_counter,z in enumerate(grid['Z_1D']):
            if within[r_counter,z_counter] and outside is False: #crop inside LCFS
                quantity_cropped[r_counter,z_counter]=crop_value
            elif not within[r_counter,z_counter] and outside is True: #crop outside LCFS
                quantity_cropped[r_counter,z_counter]=crop_value
    
    return quantity_cropped

def crop_2D(quantity,grid,reference,key,value,crop_value=0.0,over=True):
    """
    crops 

    args:
        quantity - 2D quantity to crop
        grid - quantity is stored on rectangular axes defined by grid['R_1D'] and grid['Z_1D']
        reference - reference is an object which holds its R_1D,Z_1D grid and the 2D quantity to crop according to e.g. poloidal flux, temperature map, pressure 
        key - name of quantity to crop in reference to
        value - reference[key] is compared to value to determine whether point is cropped
        crop_value - regions to be cropped will be set to this value
        over - if over then floor all reference[key] which are greater than value else floor all reference[key] under value
    notes:
    usage:
        cropped_quantity=crop_2D(quantity=pressure_map,grid=pressure_grid,reference=equilibrium,key='psirz',value=1.0,over=True) #this will remove all values of pressure which reside at a point where the poloidal flux is greater than 1.0
    """

    quantity_cropped=copy.deepcopy(quantity)

    Z_2D,R_2D=np.meshgrid(grid['Z_1D'],grid['R_1D'])
    Z_2D_flat,R_2D_flat=Z_2D.flatten(),R_2D.flatten()

    reference_values=value_at_RZ(R=R_2D_flat,Z=Z_2D_flat,quantity=reference[key],grid=reference) #calculate value of reference at the gridpoints where quantity is known
    reference_values=np.squeeze(reference_values)
    reference_values=reference_values.reshape(len(grid['R_1D']),len(grid['Z_1D']))

    for r_index,r in enumerate(grid['R_1D']): 
        for z_index,z in enumerate(grid['Z_1D']): 
            if reference_values[r_index,z_index]>value and over is True:
                quantity_cropped[r_index,z_index]=crop_value

    return quantity_cropped

def dot_product(vec1,vec2):
    """
    performs dot products of two vectors sharing same orthogonal basis
    
    notes:
        vec1,vec2 - lists of orthogonal vector components  
    """

    dot_product=0
    for vec_1_component,vec_2_component in zip(vec1,vec2):
        dot_product+=vec_1_component*vec_2_component

    return dot_product

def pitch_calc(particle_list,equilibria=None,perturbations=None,i3dr=-1,phase=0.):
    """
    calculates the pitch angles of a given set of particle positions for axisymmetric B_field 

    args:
        equilibrium - equilibrium object
        particle_list - generic particle list with positions particle_list['R'],particle_list['Z']
    notes:
        assumes RHS=R Phi Z (positive toroidal direction is counter-clockwise) 
    """

    if not equilibria and not perturbations:
        print("ERROR: pitch_calc needs equilibria or perturbations!\nreturning\n")
        return

    print("pitch_calc - start\n")

    Br_at_particles=np.zeros(len(particle_list['R']))
    Btor_at_particles=np.zeros(len(particle_list['R']))
    Bz_at_particles=np.zeros(len(particle_list['R']))

    if equilibria:
        for equilibrium in equilibria:
            B_r,B_tor,B_z=equilibrium.evaluate(particle_list['R'],particle_list['Z'])
            Br_at_particles+=B_r
            Btor_at_particles+=B_tor
            Bz_at_particles+=B_z

    if perturbations:
        for perturbation in perturbations:
            dB_r,dB_tor,dB_z=perturbation.evaluate(particle_list['R'],particle_list['phi'],particle_list['Z'],i3dr=i3dr,phase=phase)
            Br_at_particles+=dB_r
            Btor_at_particles+=dB_tor
            Bz_at_particles+=dB_z
    
    B_at_particles=np.sqrt(Br_at_particles**2+Btor_at_particles**2+Bz_at_particles**2)              
    
    Br_at_particles/=B_at_particles #calculate B_hat
    Btor_at_particles/=B_at_particles
    Bz_at_particles/=B_at_particles

    V_parallel=dot_product([particle_list['V_R'],particle_list['V_phi'],particle_list['V_Z']],[Br_at_particles,Btor_at_particles,Bz_at_particles])

    V=np.sqrt(particle_list['V_R']**2+particle_list['V_phi']**2+particle_list['V_Z']**2)
    V_pitch=V_parallel/V

    print("pitch_calc - finished\n")

    return V_pitch

def get_dfn_point(dfn,type='LOCUST',**kwargs):
    """
    returns magnitude of dfn at a point closest to that supplied

    args:
        dfn - distribution_function object
        type - type of distribution function
        kwargs - should define desired points in all possible dimensions to sample at
    usage:
        my_values=get_dfn_point(my_dfn,E=[1,2,3],V_pitch=[-1,0,1],R=[1,2,3],Z=[0,0,0],P=[-pi,-pi,-pi])
        my_values=get_dfn_point(my_dfn,type='TRANSP',E=[1,2,3],V_pitch=[-1,0,1],R=[1,2,3],Z=[0,0,0]) (for TRANSP)
    notes:
    """

    if type=='TRANSP': #TRANSP has non-standard way of defining dimensions

        dfn_values=[]

        for E,V_pitch,R,Z in zip(kwargs['E'],kwargs['V_pitch'],kwargs['R'],kwargs['Z']):

            diff_R=np.abs(dfn['R2D']-R)**2 #irregular grid so need to minimise distance to every point
            diff_Z=np.abs(dfn['Z2D']-Z)**2

            index_RZ=(diff_R+diff_Z).argmin()
            index_V_pitch=np.abs(dfn['V_pitch']-V_pitch).argmin()
            index_E=np.abs(dfn['E']-E).argmin()

            dfn_values.append(dfn['dfn'][index_RZ,index_V_pitch,index_E])

            print('R') #diagnostic printing out to show the point is in right ball park
            print(dfn['R2D'][index_RZ])
            print('Z') #diagnostic printing out to show the point is in right ball park
            print(dfn['Z2D'][index_RZ])
            print('V_pitch') #diagnostic printing out to show the point is in right ball park
            print(dfn['V_pitch'][index_V_pitch])  
            print('E') #diagnostic printing out to show the point is in right ball park
            print(dfn['E'][index_E])

    else:

        number_dimensions=len(kwargs.keys())
        for key in kwargs:
            number_points=len(kwargs[key])

        coordinate_array=np.ndarray(shape=(number_points,number_dimensions),dtype=int) #hold nD indices of all the N points we want to know in an NxnD array for an nD distribution function

        for key in kwargs: #go through supplied dimensions

            for dimension in dfn['dfn_index']: #find the dfn index of the dimension we're dealing with
                if key == dimension:
                    dimension_index=dfn['dfn_index'].tolist().index(key)

            for counter,point in enumerate(kwargs[key]):
                point_index=np.abs(dfn[key]-point).argmin() #find index of nearest grid point to the point we're dealing with along this dimension 
                coordinate_array[counter,dimension_index]=point_index 
                #print(key) #diagnostic printing out to show the point is in right ball park
                #print(dfn[key][point_index])

        dfn_values=[]
        for coordinate in coordinate_array:
            dfn_values.append(dfn['dfn'][tuple(coordinate)])


    dfn_values=np.asarray(dfn_values)
    return dfn_values

def Zeff_calc(density,charge):
    """
    calculates Zeff for given ion species
    
    args:
        density - list of relative densities given in arbitrary units
        charge - list of charges for given species in units of e
    usage:
        Zeff=Zeff_calc(density=[1,1,1],charge=[1,1,6]) #to calculate Zeff for equal parts Deuterium, Tritium and Carbon
    
    notes:

    """

    density=np.asarray(density)
    charge=np.asarray(charge)
    density_electron=np.sum(density*charge)

    Zeff=0
    for species_density,constants.species_charge in zip(density,charge):
        Zeff+=species_density*constants.species_charge**2/density_electron

    return Zeff

def Zeff_calc_density(Zeff,density,charge,fractions=False):
    """
    calculates missing impurity density fraction for given ion species and desired Zeff
    
    args:
        Zeff - desired Zeff
        density - list of relative densities given in arbitrary units
        charge - list of charges for given species in units of e
        fractions - toggle to return array of density fractions normalised to 1
    usage:
        relative_density=Zeff_calc_density(Zeff=2.5,density=[1,1],charge=[1,1,6]) #to calculate relative carbon density to achieve Zeff=2.5 in a hydrogen plasma
        density must be of length len(charge)-1 (always calculates density of final species in charge array)
    notes:

    """
    density=np.asarray(density)
    missing_charge=charge[-1]
    charge=np.asarray(charge[:-1])
    missing_density=np.sum(Zeff*density*charge-density*charge**2)/(missing_charge**2-Zeff*missing_charge)
    
    if fractions:
        density=np.append(density,missing_density)
        return density/np.sum(density)
    else:
        return missing_density

def knuth_shuffle(*args):
   """

   args:
      *args - arrays of same length to be shuffled in parallel
   notes:
      preserves index mapping between arrays
   """

   import random

   number_entries=len(args[0])
   range_number_entries=list(range(number_entries))[1:] #remove first entry since do not want entries swapping with themselves
   range_number_entries_reverse=range_number_entries.reverse()
   
   for i in range_number_entries:
      rand=random.uniform(0.,1.)
      rand_entry=int(rand*i+1)
      #iterate over all *args
      for arg in args:
         arg[rand_entry],arg[i]=arg[i],arg[rand_entry]
   return [arg for arg in args]


def sigmoid(x):
    """
    return sigmoid of x
    
    notes:
    """

    return 1 / (1 + np.exp(-x))

def sigmoid_derivative(x):
    """
    return derivative of sigmoid 

    notes:
    """

    return sigmoid(x)*(1-sigmoid(x))

def local_minimum_2D(x,y,X,Y,quantity,threshold=0.1):
    """
    args:
        x - x starting coordinate 
        y - y starting coordinate
        X - X grid dimension
        Y - Y grid dimension
        quantity - 2D quantity on rectilinear grid defined by X,Y 
        threshold - final result is accurate to threshold*grid_spacing
    returns:
        x - x minimum coordinate
        y - y minimum coordinate
        value - value at x,y
    notes:
    """

    if threshold>1.: 
        print("ERROR: local_minimum_2D() threshold must be <=1!\nreturning\n")
        return

    interpolator=interpolate_2D(X,Y,quantity,type='RBS')

    previous_value=interpolator(x,y)
    current_value=previous_value
    dx=X[1]-X[0]
    dy=Y[1]-Y[0]
    step=1.
    direction=1.

    def walk(x,y,dx,dy,current_value,axis):
        decreasing=True
        while decreasing:
            if axis is 'x':
                x+=dx*step*direction
            else:
                y+=dy*step*direction
            previous_value=current_value
            current_value=interpolator(x,y)
            if current_value>previous_value: decreasing=False
        return x,y,current_value

    while step>=threshold:
        x,y,current_value=walk(x,y,dx,dy,current_value,'x')
        x,y,current_value=walk(x,y,dx,dy,current_value,'y')
        step*=0.5
        direction*=-1.

    return x,y,previous_value

def extrapolate_kinetic_profiles(*kinetic_profiles,**kwargs):
    """
    extend kinetic profiles outside defined region

    args:
        kinetic_profiles - kinetic_profile objects e.g. temperatures, densities or rotations
        axis - extrapolate profile according to this independent variable
        start - value of axis to extrapolate from (deletes current values outside this)
        end - max value of axis to extrapolate to [axis units]
        decay_length - sets scale length of exponential decay from start [axis units]
        floor_distance - maintain constant profile once axis=start+floor_distance is reached [axis units]
        floor_value - maintain constant profile once kinetic_profile=floor_value is reached [kinetic_profile units]
        uniform_grid - toggle to dump extrapolated portion of profile on uniform grid (ensures all kinetic_profiles conform to same grid outside start)
        return_indices - toggle returning set of indices with each kinetic_profile respresenting region of extrapolation
    usage:
        extrapolated_profiles,indices=extrapolate_kinetic_profiles(temperature,density,axis='r_1D',start=0.4,end=3.,decay_length=1.,floor_distance=0.1,return_indices=True)
    notes:
        must supply either floor_distance or floor_value
        this routine may begin extrapolating before start if start lies between grid points - use return_indices argument for help 
        this routine may return profiles that do not extend exactly up to end value
    """

    axis=kwargs.get('axis','flux_pol_norm')
    start=kwargs.get('start',None)
    end=kwargs.get('end',2.)
    decay_length=kwargs.get('decay_length',.01)
    floor_distance=kwargs.get('floor_distance',None)
    floor_value=kwargs.get('floor_value',None)
    uniform_grid=kwargs.get('uniform_grid',True)
    return_indices=kwargs.get('return_indices',False)

    if sum([arg is not None for arg in [floor_distance,floor_value]])!=1:
        print("ERROR: extrapolate_kinetic_profiles() requires either 'floor_distance' or 'floor_value' args!\nreturning\n")
        return

    if decay_length<=0.: 
        print("ERROR: extrapolate_kinetic_profiles() requires decay_length>0!\nreturning\n")
        return

    quantity_to_extrapolate_dispatch={} #define map between type of kinetic profile and variable name of quantity we want to extrapolate
    quantity_to_extrapolate_dispatch['number_density']='n'
    quantity_to_extrapolate_dispatch['temperature']='T'
    quantity_to_extrapolate_dispatch['rotation']='rotation_ang'
    extrapolated_kinetic_profiles=[]
    extrapolated_indices=[]

    for kinetic_profile_ in kinetic_profiles:
        kinetic_profile=copy.deepcopy(kinetic_profile_)
        floor_value=kwargs.get('floor_value',None) #rest floor value since possibly edited below
    
        if kinetic_profile.LOCUST_input_type in quantity_to_extrapolate_dispatch.keys(): #figure out variable name of what we are extrapolating
            quantity_to_extrapolate=quantity_to_extrapolate_dispatch[kinetic_profile.LOCUST_input_type]
        else:
            print(f"ERROR: extrapolate_kinetic_profiles could not determine LOCUST_input_type of input (ID={kinetic_profile.ID})!\nskipping\n")

        #preprocess by cutting off outside where we want to extrapolate
        if start is None: 
            start=kinetic_profile[axis][-1]
        indices_to_extrapolate=np.where(kinetic_profile[axis]>start)[0] #remember indices of original data
        indices_not_extrapolated=np.where(kinetic_profile[axis]<=start)[0]
        kinetic_profile[axis]=np.delete(kinetic_profile[axis],indices_to_extrapolate)
        kinetic_profile[quantity_to_extrapolate]=np.delete(kinetic_profile[quantity_to_extrapolate],indices_to_extrapolate)
        axis_extrapolated=np.array([]) #store extrapolated profiles heres
        quantity_extrapolated=np.array([])
        
        #first stage is to smoothly transition to exponential decay via parabola 
        step_length=decay_length*0.01 #arbitrarily assign 100 steps per parabolic decay length
        smoothing_rate=1./(decay_length) #set parabola acceleration
        step_counter=0.
        starting_profile_value=kinetic_profile[quantity_to_extrapolate][-1]
        starting_axis_value=kinetic_profile[axis][-1]
        current_profile_value=starting_profile_value
        current_axis_value=kinetic_profile[axis][-1]
        current_gradient=(kinetic_profile[quantity_to_extrapolate][-1]-kinetic_profile[quantity_to_extrapolate][-2])/(kinetic_profile[axis][-1]-kinetic_profile[axis][-2])
        while current_gradient>-(1./decay_length)*np.exp((-1.*step_length)/decay_length)*starting_profile_value:
            step_counter+=1
            current_axis_value+=step_length
            current_gradient-=2.*step_counter*step_length*smoothing_rate*starting_profile_value
            current_profile_value+=current_gradient*step_length
            axis_extrapolated=np.concatenate((axis_extrapolated,[current_axis_value]))
            quantity_extrapolated=np.concatenate((quantity_extrapolated,[current_profile_value]))

        #next stage is to perform the exponential extrapolation
        starting_profile_value=current_profile_value
        step_counter=0
        step_length*=10. #100 steps per exponential decay length
        
        def end_condition(): #inject end condition logic - either reach distance from LCFS or profile value
            if floor_distance:
               return current_axis_value<starting_axis_value+floor_distance 
            elif floor_value:
               return current_profile_value>floor_value

        while end_condition() and current_axis_value<end:
            step_counter+=1
            current_axis_value+=step_length
            current_profile_value=starting_profile_value*np.exp((-1.*step_counter*step_length)/decay_length)
            axis_extrapolated=np.concatenate((axis_extrapolated,[current_axis_value]))
            quantity_extrapolated=np.concatenate((quantity_extrapolated,[current_profile_value]))

        floor_value=floor_value if floor_value else current_profile_value #XXX this needs fixing for floor distance! check density profiles
        if len(quantity_extrapolated)>0:
            quantity_extrapolated[-1]=floor_value #remove overshoot
            current_profile_value=floor_value

        #final stage is to extend as a flat profile to end value
        if current_axis_value>=end:
            pass
        else:
            distance_to_end=np.linspace(current_axis_value+step_length,end,100) #use largeish number here because cutting off just inside end later
            quantity_to_end=np.full(len(distance_to_end),floor_value)
            axis_extrapolated=np.concatenate((axis_extrapolated,distance_to_end))
            quantity_extrapolated=np.concatenate((quantity_extrapolated,quantity_to_end))
            
        #add flat section of profile
        kinetic_profile[axis]=np.concatenate((kinetic_profile[axis],axis_extrapolated))
        kinetic_profile[quantity_to_extrapolate]=np.concatenate((kinetic_profile[quantity_to_extrapolate],quantity_extrapolated))

        #cut off any data outside end
        inside_end=np.where(kinetic_profile[axis]<end)[0]
        last_index_before_end=inside_end[-1]
        kinetic_profile[axis]=kinetic_profile[axis][inside_end]
        kinetic_profile[quantity_to_extrapolate]=kinetic_profile[quantity_to_extrapolate][inside_end]
                
        #if requested, do some interpolation to sit all profiles on same grid outside of start
        if uniform_grid and len(quantity_extrapolated)>0:
            #often mag axis is not very accurate - so wise to only interpolate onto uniform grid for points outside LCFS
            #if you interpolate inside LCFS, then if rmaxis is inboard of true axis, you may observe rapid oscillations in LOCUST's profiles
            #at low flux_pol_norm, due to the fact it interpolates the kinetic profiles either side of the axis but sorts when dumping to LOCUST format  
            interpolator=interpolate_1D(kinetic_profile[axis],kinetic_profile[quantity_to_extrapolate],function='linear',type='interp1d') 
            index_before_extrap=indices_not_extrapolated[-1]
            index_after_extrap=index_before_extrap+1 
            axis_value_before_extrap=kinetic_profile[axis][index_before_extrap]
            axis_value_after_extrap=kinetic_profile[axis][index_after_extrap]
            kinetic_profile[axis]=np.concatenate((
                kinetic_profile[axis][indices_not_extrapolated],
                np.linspace(axis_value_before_extrap,kinetic_profile[axis][last_index_before_end],10000)[1:])) #use axis_value_before_extrap in linspace since axis_value_after_extrap will vary for different profiles
            indices_extrapolated=list(set(indices_not_extrapolated) ^ set(range(len(kinetic_profile[axis])))) #find indices of values within extrapolation zone
            kinetic_profile[quantity_to_extrapolate]=np.concatenate((
                kinetic_profile[quantity_to_extrapolate][indices_not_extrapolated],
                interpolator(kinetic_profile[axis][indices_extrapolated])))

        indices_extrapolated=list(set(indices_not_extrapolated) ^ set(range(len(kinetic_profile[axis])))) #find indices of values within extrapolation zone            

        extrapolated_kinetic_profiles.append(kinetic_profile) #add returns to results array
        extrapolated_indices.append(indices_extrapolated)

    if return_indices:
        return extrapolated_kinetic_profiles,extrapolated_indices
    else:
        return extrapolated_kinetic_profiles

def extrapolate_kinetic_profiles_ITER(equilibrium,*kinetic_profiles,**kwargs):
    """
    extend kinetic profiles outside defined region for ITER excel/IDS kinetic profiles

    args:
        equilibrium - equilibrium object with flux grids to extrapolate over
        kinetic_profiles - kinetic_profile objects e.g. temperatures, densities or rotations
        decay_length - sets scale length of exponential decay from start [metres]
        floor_distance - maintain constant profile once axis=start+floor_distance is reached [metres]
        floor_value - maintain constant profile once kinetic_profile=floor_value is reached [kinetic_profile units]
    usage:
        extrapolated_profiles=processing.utils.extrapolate_kinetic_profiles_ITER(equilibrium,*kinetic_profiles,decay_length=0.35)
    notes:
        XXX warning flux_tor and flux_tor_coord extrapolation will not work with ITER IDSs whose flux_tor do not extend past LCFS
        assumes kinetic_profiles run to flux_pol_norm>=1.
        assumes extrapolation region at flux_pol_norm>=1.
        extrapolates according to minor radius
        all kinetic profiles must be mapped to same supplied equilibrium
        will only extrapolate if all kinetic profiles contain necessary quantities
        dumps everything on same grid so can be dumped to IDS
        ITER rough numbers
            max distance between plasma + wall in ITER ~1.50m in divertor (outboard side ~0.35m, ~1m at about 11 o'clock)
            max outboard plasma radius in ITER ~8.20m
            mag axis in ITER R~6.412m Z~0.570m
            mag axis in ITER excel = 6.2m
            min inboard plasma radius in ITER ~4.20m
    """

    #get and check args
    decay_length=kwargs.get('decay_length',.03)
    floor_distance=kwargs.get('floor_distance',None)
    floor_value=kwargs.get('floor_value',None)
    if sum([arg is not None for arg in [floor_distance,floor_value]])!=1:
        print("ERROR: extrapolate_kinetic_profiles_ITER() requires either 'floor_distance' or 'floor_value' args!\nreturning\n")
        return

    #make sure all the profiles contain necessary data   
    if not all(['flux_pol_norm' in kinetic_profile.data for kinetic_profile in kinetic_profiles]) and all(['r_1D' in kinetic_profile.data for kinetic_profile in kinetic_profiles]): 
        print("ERROR: profiles supplied to extrapolate_kinetic_profiles_ITER() must all contain 'r_1D' and 'flux_pol_norm'!") 
        print("IDs of kinetic_profiles missing minor radius:\n{}".format([f'{kinetic_profile.ID}' for kinetic_profile in kinetic_profiles if 'r_1D' not in kinetic_profile.data and equilibrium is not None]))        
        print("IDs of kinetic_profiles missing flux:\n{}".format([f'{kinetic_profile.ID}' for kinetic_profile in kinetic_profiles if 'flux_pol_norm' not in kinetic_profile.data]))
        print("returning\n")
        return

    #first thing to do is accurately determine the magnetic axis location
    rmaxis_equilibrium,zmaxis_equilibrium,simag=equilibrium.mag_axis_calc(threshold=1.e-7)

    #create data along outboard midplane for interpolation 
    midplane_radius_major=np.linspace(rmaxis_equilibrium,equilibrium['R_1D'][-1],100)
    midplane_poloidal_flux=value_at_RZ(
    R=midplane_radius_major,
    Z=np.full(len(midplane_radius_major),
        zmaxis_equilibrium),
    quantity=equilibrium['psirz'],
    grid=equilibrium)
    midplane_poloidal_flux_norm=(midplane_poloidal_flux-simag)/(equilibrium['sibry']-simag)
    midplane_toroidal_flux=value_at_RZ( #XXX REMOVE THIS IF NOW NOT NEEDED
    R=midplane_radius_major,
    Z=np.full(len(midplane_radius_major),
        zmaxis_equilibrium),
    quantity=equilibrium['phirz'],
    grid=equilibrium)

    interpolator_r_min_of_psi_norm=interpolate_1D(midplane_poloidal_flux_norm,midplane_radius_major-rmaxis_equilibrium,function='linear',type='interp1d') #create r_minor(psi_norm) interpolator for later
    interpolator_psi_of_r_maj=interpolate_1D(midplane_radius_major,midplane_poloidal_flux,function='linear',type='interp1d') #create psi(r) interpolator for later
    interpolator_phi_of_r_maj=interpolate_1D(midplane_radius_major,midplane_toroidal_flux,function='linear',type='interp1d') #create phi(r) interpolator for later (just use equilibrium grid line closest to mag axis since edge is prone to getting noise)
    
    #XXX switching to just taking the closest grid line - delete above if keeping this 
    interpolator_phi_of_r_maj=interpolate_1D(equilibrium['R_1D'],equilibrium['phirz'][:,np.abs(equilibrium['Z_1D']-zmaxis_equilibrium).argmin()],function='linear',type='interp1d') #create phi(r) interpolator for later (just use equilibrium slice closest to mag axis since edge is prone to noise)

    '''#XXX
    import matplotlib.pyplot as plt #XXX
    plt.plot(midplane_radius_major,midplane_toroidal_flux)
    plt.plot(equilibrium['R_1D'],equilibrium['phirz'][:,np.abs(equilibrium['Z_1D']-zmaxis_equilibrium).argmin()],color='red')
    plt.show()
    '''#XXX

    #take flux_pol_norm stored in kinetic profiles as gospel and recalculate minor radius grid according to supplied equilibrium
    for kinetic_profile in kinetic_profiles: 
        if 'r_1D' in kinetic_profile.data: kinetic_profile.set(r_1D=np.asarray(interpolator_r_min_of_psi_norm(kinetic_profile['flux_pol_norm'])))

    #perform extrapolation
    extrapolated_kinetic_profiles,extrapolated_indices=extrapolate_kinetic_profiles(
    *kinetic_profiles,
    axis='r_1D',
    end=np.max(equilibrium['R_1D'])-rmaxis_equilibrium,
    decay_length=decay_length,
    floor_distance=floor_distance,
    floor_value=floor_value,
    return_indices=True,
    uniform_grid=True)

    r_LCFS=equilibrium.calc_r_LCFS() #find minor radius at outboard LCFS for normalisation later
    #determine toroidal and poloidal flux only at the new extrapolated positions
    for extrapolated_kinetic_profile,extrapolated_index in zip(extrapolated_kinetic_profiles,extrapolated_indices):

        index_not_extrapolated=list(set(extrapolated_index) ^ set(range(len(extrapolated_kinetic_profile['r_1D'])))) #find indices of values that were left alone
        
        #recalculate poloidal and toroidal flux axes        
        flux_pol=interpolator_psi_of_r_maj(extrapolated_kinetic_profile['r_1D']+rmaxis_equilibrium)
        flux_pol_norm=(flux_pol-simag)/(equilibrium['sibry']-simag)

        flux_tor=interpolator_phi_of_r_maj(extrapolated_kinetic_profile['r_1D']+rmaxis_equilibrium)
        flux_tor_LCFS=interpolator_phi_of_r_maj(r_LCFS+rmaxis_equilibrium)
        flux_tor_axis=interpolator_phi_of_r_maj(rmaxis_equilibrium)

        #calculate toroidal flux coordinate - rho_tor in IDS
        flux_tor_coord=np.sqrt(np.abs((flux_tor*2.*np.pi)/(np.pi*equilibrium['bcentr']))) 
        flux_tor_coord_LCFS=np.sqrt(np.abs((flux_tor_LCFS*2.*np.pi)/(np.pi*equilibrium['bcentr'])))
        flux_tor_coord_axis=np.sqrt(np.abs((flux_tor_axis*2.*np.pi)/(np.pi*equilibrium['bcentr'])))
        flux_tor_coord_norm=(flux_tor_coord-flux_tor_coord_axis)/(flux_tor_coord_LCFS-flux_tor_coord_axis)

        extrapolated_kinetic_profile.set(
        flux_pol=flux_pol,
        flux_pol_norm=flux_pol_norm,
        flux_tor=flux_tor,
        flux_tor_coord=flux_tor_coord,
        flux_tor_coord_norm=flux_tor_coord_norm)

        '''#XXX
        if extrapolated_kinetic_profile.properties['species']=='electrons' and 'n' in extrapolated_kinetic_profile.data:
            import matplotlib.pyplot as plt
            plt.plot(extrapolated_kinetic_profile['flux_tor'],extrapolated_kinetic_profile['n'])
            plt.plot(extrapolated_kinetic_profile['flux_tor'][extrapolated_index],extrapolated_kinetic_profile['n'][extrapolated_index],color='red')
            plt.show()
            plt.plot(extrapolated_kinetic_profile['r_1D'],extrapolated_kinetic_profile['flux_tor']) #XXX THIS SHOULD NOT BE FLAT EVER 
            plt.plot(extrapolated_kinetic_profile['r_1D'][extrapolated_index],extrapolated_kinetic_profile['flux_tor'][extrapolated_index],color='red')
            plt.show()
        '''#XXX

    return extrapolated_kinetic_profiles

def sample_1D(x,y,method='rejection',size=1000):
    """
    args: 
        x - x axis of pdf
        y - pdf
        method - rejection or inverse_transform (much quicker)
        size - number of samples to draw
    notes:
    """


    if method is 'rejection': 

        x_min,x_max=np.min(x),np.max(x)
        y_min,y_max=np.min(y),np.max(y)
        y_of_x=interpolate_1D(x,y)
        rng=np.random.default_rng()
        N=0
        samples=[]
        while N<size:
            rand_x=rng.uniform(x_min,x_max)
            rand_y=rng.uniform(y_min,y_max)
            if rand_y<=y_of_x(rand_x): #keep sample
                samples.append(rand_x)
                N+=1
        samples=np.asarray(samples)

    elif method is 'inverse_transform': 

        cdf=np.cumsum(y)
        cdf,x=sort_arrays(cdf,x)
        cdf,x=np.asarray(cdf),np.asarray(x)
        cdf_normalised=(cdf-np.min(cdf))/(np.max(cdf)-np.min(cdf))
        cdf_inverse_normalised=interpolate_1D(cdf_normalised,x,function='linear')
        rng=np.random.default_rng()
        rands=rng.random(size=size)
        samples=cdf_inverse_normalised(rands)

    return samples

def KS_test(n1,n2,data1,data2,alpha=None):
    """
    perform KS test on two empirical distribution functions data1 data2

    args:
        n1 - data1 sample size
        n2 - data2 sample size
        data1 - empirical distribution function 
        data2 - empirical distribution function 
        alpha - optional confidence interval for rejecting null hypothesis
    notes:
        data1, data2 assumed to be against same axis
    """

    #table of critical D coefficients given confidence interval alpha
    c_of_alpha={}
    c_of_alpha[0.1]=1.22
    c_of_alpha[0.05]=1.36
    c_of_alpha[0.025]=1.48
    c_of_alpha[0.01]=1.63
    c_of_alpha[0.005]=1.73
    c_of_alpha[0.001]=1.95
    c_of_alpha[None]=None

    if not any([alpha==alph for alph in c_of_alpha.keys()]):
        print(f"ERROR: KS_test given invalid alpha value - available options = {[a for a in c_of_alpha.keys()]}!\nreturning\n")
        return

    def normalise(array):
        """
        return normalised array
        """
        return (array-np.min(array))/(np.max(array)-np.min(array)) 

    def D_crit(n1,n2,alpha):
        """
        args:
            n1 - sample size
            n2 - sample size
            alpha - set probability of observing D>D_crit given null i.e. set confidence interval 
        notes:
            this is an approximation - can use binary search to iteratively calculate this value using p() below
            reject null if D>D_crit given alpha i.e. if D is in upper alpha% of D given null
        """
        return c_of_alpha[alpha]*np.sqrt((n1+n2)/(n1*n2))

    def p(n1,n2,D):
        """
        return probability that D_crit(alpha)>D_observed for n1,n2 where D_crit is KS statistic in top alpha% given null i.e. p = probability we measure this conclusion result given Null
        
        args:
            n1 - sample size
            n2 - sample size
            D - KS statistic
        notes:
            this parameter is irrespective of the shape of the empirical distribution functions
        """

        N=(n1*n2)/(n1+n2)
        x=D*(np.sqrt(N)+0.12+0.11/np.sqrt(N))
        result=0
        for j in range(1000):
            j+=1
            result+=(-1)**(j-1)*np.exp(-2.*j**2*x**2)
        return 2.*result

    data1_cum,data2_cum=np.cumsum(data1),np.cumsum(data2)
    data1_cum_norm=normalise(data1_cum)
    data2_cum_norm=normalise(data2_cum)

    D=np.max(np.abs(data1_cum_norm-data2_cum_norm))
    P=p(n1,n2,D)

    if alpha:
        reject=True if D>D_crit(n1,n2,alpha) else False

    return D,P,reject

#################################
 
##################################################################
 
###################################################################################################