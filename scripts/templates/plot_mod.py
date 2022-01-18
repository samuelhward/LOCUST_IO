#plot_mod.py
 
"""
Samuel Ward
12/08/20
----
plotting helpers
---
 
notes:
---
"""

###################################################################################################
#Preamble

import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway
try:
    import pathlib,matplotlib,multiprocessing,os,functools
    import numpy as np 
    import matplotlib.pyplot as plt 
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    from classes.input_classes.equilibrium import Equilibrium
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/equilibrium.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.perturbation import Perturbation
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/perturbation.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    from classes.input_classes.temperature import Temperature
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/temperature.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.number_density import Number_Density
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/number_density.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.rotation import Rotation
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/rotation.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.beam_deposition import Beam_Deposition
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/beam_deposition.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.input_classes.wall import Wall
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/input_classes/wall.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from classes.output_classes.distribution_function import Distribution_Function
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/output_classes/distribution_function.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    from classes.output_classes.particle_list import Final_Particle_List
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/output_classes/distribution_function.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    from classes.output_classes.rundata import Rundata
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/output_classes/distribution_function.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    from classes.output_classes.poincare import Poincare
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/output_classes/poincare.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    from classes.output_classes.orbits import Orbits
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/output_classes/orbits.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    from classes.output_classes.moments import Moments
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/output_classes/moments.py could not be imported!\nreturning\n")
    sys.exit(1)
try:
    from classes.output_classes.output_mesh import Output_Mesh
except:
    raise ImportError("ERROR: LOCUST_IO/src/classes/output_classes/output_mesh.py could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import processing.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/processing/utils.py could not be imported!\nreturning\n")
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
try:
    import support
except:
    raise ImportError("ERROR: LOCUST_IO/src/support.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    #cwd=pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    #sys.path.append(str(cwd.parents[1]))
    from .template_mod import *
except:
    raise ImportError("ERROR: templates/template_mod.py could not be imported!\nreturning\n") 
    sys.exit(1)

################################################################## 
#Main 


def read_locust_io_obj(filename,classtype='eq',*args,**kwargs):
    """read single locust_io object from filename and class type
    """

    data_format_dispatch={}
    data_format_dispatch['eq']='GEQDSK'
    data_format_dispatch['pert']='MARSF'
    data_format_dispatch['temp_e']='LOCUST'
    data_format_dispatch['temp_i']='LOCUST'
    data_format_dispatch['numden']='LOCUST'
    data_format_dispatch['rot']='LOCUST'
    data_format_dispatch['beamdepo']='LOCUST_FO_weighted'
    data_format_dispatch['wall']='LOCUST_3D'
    data_format_dispatch['dfn']='LOCUST'
    data_format_dispatch['fpl']='LOCUST'
    data_format_dispatch['rund']='LOCUST'
    data_format_dispatch['poinc1']='LOCUST'
    data_format_dispatch['poinc2']='LOCUST'
    data_format_dispatch['poinc3']='LOCUST'
    data_format_dispatch['orbit2D']='LOCUST'
    data_format_dispatch['orbit3D']='LOCUST'
    data_format_dispatch['pfc']='LOCUST'
    data_format_dispatch['mom']='LOCUST'

    classes_dispatch={}
    classes_dispatch['eq']=Equilibrium
    classes_dispatch['pert']=Perturbation
    classes_dispatch['temp_e']=Temperature
    classes_dispatch['temp_i']=Temperature
    classes_dispatch['numden']=Number_Density
    classes_dispatch['rot']=Rotation
    classes_dispatch['beamdepo']=Beam_Deposition
    classes_dispatch['wall']=Wall
    classes_dispatch['dfn']=Distribution_Function
    classes_dispatch['fpl']=Final_Particle_List
    classes_dispatch['rund']=Rundata
    classes_dispatch['poinc1']=Poincare
    classes_dispatch['poinc2']=Poincare
    classes_dispatch['poinc3']=Poincare
    classes_dispatch['orbit2D']=Orbits
    classes_dispatch['orbit3D']=Orbits
    classes_dispatch['pfc']=Output_Mesh
    classes_dispatch['mom']=Moments

    if filename:
        try:
            return classes_dispatch[classtype](ID=filename,data_format=data_format_dispatch[classtype],filename=filename,*args,**kwargs)
        except:
            return None
    else:
        print(f'could not open {filename}!')
        return None

def get_io_files(batch_data=None,classtype='eq',yield_filenames=False,return_lists=False,*args,**kwargs):
    """get all input/output objects/filenames for particular set of simulations
    """

    file_extension_dispatch={}
    file_extension_dispatch['dfn']='*.dfn'
    file_extension_dispatch['fpl']='ptcl_cache.dat'
    file_extension_dispatch['rund']='rundata*_1'
    file_extension_dispatch['poinc1']='Poincare_map_1.dat'
    file_extension_dispatch['poinc2']='Poincare_map_2.dat'
    file_extension_dispatch['poinc3']='Poincare_map_3.dat'
    file_extension_dispatch['orbit2D']='ORBIT_2D*'
    file_extension_dispatch['orbit3D']='ORBIT_3D*'
    file_extension_dispatch['pfc']='pfc_hitmap*'
    file_extension_dispatch['eq']='locust_eqm'
    file_extension_dispatch['pert']='BPLASMA_n?'
    file_extension_dispatch['temp_e']='profile_Te.dat'
    file_extension_dispatch['temp_i']='profile_Ti.dat'
    file_extension_dispatch['numden']='profile_ne.dat'
    file_extension_dispatch['rot']='profile_wT.dat'
    file_extension_dispatch['beamdepo']='ptcles.dat'
    file_extension_dispatch['wall']='locust_wall'
    file_extension_dispatch['mom']='*.h5'

    io_type_dispatch={}
    io_type_dispatch['dfn']='output'
    io_type_dispatch['fpl']='output'
    io_type_dispatch['rund']='output'
    io_type_dispatch['poinc1']='output'
    io_type_dispatch['poinc2']='output'
    io_type_dispatch['poinc3']='output'
    io_type_dispatch['orbit2D']='output'
    io_type_dispatch['orbit3D']='output'
    io_type_dispatch['pfc']='output'
    io_type_dispatch['mom']='output'
    io_type_dispatch['eq']='input'
    io_type_dispatch['pert']='input'
    io_type_dispatch['temp_e']='input'
    io_type_dispatch['temp_i']='input'
    io_type_dispatch['numden']='input'
    io_type_dispatch['rot']='input'
    io_type_dispatch['beamdepo']='input'
    io_type_dispatch['wall']='input'

    for parameter_string,dir in zip(batch_data.parameter_strings,batch_data.args_batch[f'LOCUST_run__dir_{io_type_dispatch[classtype]}']): #within each GPU folder the path to each file is the same

        dir=pathlib.Path(dir.strip("\'"))

        if classtype=='rund': #need to do something a bit special for rundata, since usually there's more than one but only one corresponds to the run where the dfn was calculated
            try:
                dfn_filename=list(dir.glob('F_*.dfn'))[0].parts[-1]
                date_ID=str(dfn_filename).strip('F_').strip('_TOTL.dfn')
                dir_filepaths=list(dir.glob(f'rundata_{date_ID}*_1'))
            except:
                yield None 
                continue
        else:
            dir_filepaths=list(dir.glob(file_extension_dispatch[classtype])) #get all filenames for runs corresponding to this choice of parameters    

        if dir_filepaths:
            if yield_filenames:
                if len(dir_filepaths)>1:
                    yield dir_filepaths
                else:
                    if return_lists:
                        yield [dir_filepaths[0]]
                    else:
                        yield dir_filepaths[0]
            else:
                if len(dir_filepaths)>1:
                    yield [read_locust_io_obj(filename=dir_filepath,classtype=classtype,*args,**kwargs) for dir_filepath in dir_filepaths]
                else:
                    if return_lists:
                        yield [read_locust_io_obj(filename=dir_filepaths[0],classtype=classtype,*args,**kwargs)]
                    else:
                        yield read_locust_io_obj(filename=dir_filepaths[0],classtype=classtype,*args,**kwargs)
        else:
            print(f"could not find {classtype} in {dir}")
            yield None

# for backwards compatability

def get_output_files(*args,**kwargs):    
    return get_io_files(*args,**kwargs)
def get_input_files(*args,**kwargs):    
    return get_io_files(*args,**kwargs)

def apply_func_parallel(func,classtype,batch_data,processes=4,chunksize=4,**func_args):
    """
    apply function func() in parallel across every piece of batch data 
    notes:
        func signature must take filename and classtype
    """
    filenames=list(get_io_files(batch_data=batch_data,classtype=classtype,yield_filenames=True))
    with multiprocessing.Pool(processes=processes) as pool:
        results=pool.map(functools.partial(func,classtype=classtype,**func_args),(filenames),chunksize=chunksize)
    return results

'''
XXX this decorator does not work since multiprocessing.Pool.map requires pickling
and decorated functions cannot be pickled

#decorator for applying a function to all locust_io objects in a given batch

#add extra top level because we need an arg (classtype) which is only 
#used by the decorator and not the function we are decorating 
#i.e. we want to change HOW we decorate without having classtype in the returned function signature
def apply_func(classtype,*args,**kwargs): #this is the signature of the decorate @apply_func_filename
    """function to decorate should operate on a locust_io object

    args:
        classtype - type of locust_io object e.g. fpl
    notes:
        decorated function signature is now the filename of the object to operate on
        
    """
    def decorator(func): #this defines the decorator 
        def inner(filename,*args,**kwargs): # this is now the new signature of the decorated function
            locust_io_obj=read_locust_io_obj(filename=filename,classtype=classtype)
            return func(locust_io_obj,*args,**kwargs)
        return inner
        #return functools.partial(inner,classtype=classtype)
    return decorator
#easiest if you start at the decorator() level, which always takes a func as its arg, which is always the proceeding "def func()" statement
#that function decorator() then defines an inner() function, with a given signature that does a certain thing, then returns inner()
#at this point it's crucial to remember that inner() is a function - so decorator is returning a function definition
#hence calling decorator() on func is actually now calling inner() and returning that result
#then at the top level there is another layer which simply takes an arg and changes the decorator() definition based on that arg
#i.e. the arg is changing how the behaviour of the decoration - calling the outermost level with its args returns the decorator, which then returns the next level etc.
#and this essentially means that @outermost_decorator(some_arg)=@final_decorator 
''' 

def get_divergence_files(batch_data):
    """
    notes:
        assume divergence is dumped on a rectangular grid
    """

    for run_number,(parameter_string,dir_output) in enumerate(zip(batch_data.parameter_strings,batch_data.args_batch['LOCUST_run__dir_output'])): #within each GPU folder the path to each output is the same
        #get the equilibrium to get an idea of domain to check over
        equilibrium=Equilibrium(ID='',data_format='GEQDSK',filename=batch_data.args_batch['RMP_study__filepath_equilibrium'][run_number],GEQDSKFIX1=True,GEQDSKFIX2=True)
        dir_output=pathlib.Path(str(dir_output).strip("\'"))
        dir_output_filepaths_RZ=list(dir_output.glob('field_data_divergence_check_RZ.out')) #get all filenames for runs corresponding to this choice of parameters    
        dir_output_filepaths_XY=list(dir_output.glob('field_data_divergence_check_XY.out')) #get all filenames for runs corresponding to this choice of parameters    
        if dir_output_filepaths_RZ and dir_output_filepaths_XY:
            try:

                field_data_RZ=Perturbation(ID=parameter_string,data_format='LOCUST_field_data',filename=dir_output_filepaths_RZ[0])
                field_data_XY=Perturbation(ID=parameter_string,data_format='LOCUST_field_data',filename=dir_output_filepaths_XY[0])

                field_data_RZ_nR=len(np.where(field_data_RZ['Z']==field_data_RZ['Z'][0])[0])
                field_data_RZ_nZ=np.where(field_data_RZ['Z']==field_data_RZ['Z'][0])[0][1]

                for quant in field_data_RZ.data:
                    field_data_RZ[quant]=field_data_RZ[quant].reshape(int(field_data_RZ_nR),int(field_data_RZ_nZ))
                field_data_RZ['R_1D']=field_data_RZ['R'][:,0]
                field_data_RZ['Z_1D']=field_data_RZ['Z'][0,:]

                field_data_XY_nphi=len(np.where(field_data_XY['R']==field_data_XY['R'][0])[0])
                field_data_XY_nR=np.where(field_data_XY['R']==field_data_XY['R'][0])[0][1]
                
                for quant in field_data_XY.data:
                    field_data_XY[quant]=field_data_XY[quant].reshape(int(field_data_XY_nR),int(field_data_XY_nphi)).T
                field_data_XY['phi']=field_data_XY['phi'][0,:]
                field_data_XY['R_1D']=field_data_XY['R'][:,0]

                yield field_data_RZ,field_data_XY

            except:
                yield None
        else:
            yield None
'''
def plot_kinetic_profiles(batch_data):
    """
    notes:
    """

    for parameters__database,parameters__sheet_name_kinetic_prof in zip(
            batch_data.parameters__databases,batch_data.parameters__sheet_names_kinetic_prof): 
        for parameters__kinetic_prof_tF_tE,parameters__kinetic_prof_Pr in zip(batch_data.parameters__kinetic_profs_tF_tE,batch_data.parameters__kinetic_profs_Pr):

            parameters__kinetic_prof_tF_tE_string=parameters__kinetic_profs_tF_tE__dispatch[parameters__kinetic_prof_tF_tE] #generate some variable string equivalents for later
            parameters__kinetic_prof_Pr_string=parameters__kinetic_profs_Pr__dispatch[parameters__kinetic_prof_Pr]

            parameters__sheet_name_rotation=target_IDS_dispatch[parameters__database][parameters__kinetic_prof_Pr_string][parameters__kinetic_prof_tF_tE_string]['sheet_name_rotation']
            parameters__var_name_rotation=target_IDS_dispatch[parameters__database][parameters__kinetic_prof_Pr_string][parameters__kinetic_prof_tF_tE_string]['var_name_rotation']
            RMP_study__filepath_kinetic_profiles = list((RMP_study__dir_input_database / parameters__database / folder_name_DataEq).glob('*.xlsx'))[0]

            #temperature=Temperature(ID='',data_format='EXCEL1',species='ions',filename=RMP_study__filepath_kinetic_profiles.relative_to(support.dir_input_files),sheet_name=parameters__sheet_name_kinetic_prof)
            #rotation=Rotation(ID='',data_format='EXCEL1',filename=RMP_study__filepath_kinetic_profiles.relative_to(support.dir_input_files),sheet_name=parameters__sheet_name_kinetic_prof,rotation_name=parameters__var_name_rotation,sheet_name_rotation=parameters__sheet_name_rotation)
            density=Number_Density(ID='',data_format='EXCEL1',species='deuterium',filename=RMP_study__filepath_kinetic_profiles.relative_to(support.dir_input_files),sheet_name=parameters__sheet_name_kinetic_prof)
'''

def plot_poincare_psi_theta(poincare_map,equilibrium,phi=0.,ax=False,fig=False):

    #interpolate poincare R Z to psi theta

    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if fig is False:
        fig_flag=False
    else:
        fig_flag=True

    if fig_flag is False:
        fig = plt.figure() #if user has not externally supplied figure, generate
    
    if ax_flag is False: #if user has not externally supplied axes, generate them
        ax = fig.add_subplot(111)

    equilibrium.set(psirz_norm=(equilibrium['psirz']-equilibrium['simag'])/(equilibrium['sibry']-equilibrium['simag']))

    R,Z=np.meshgrid(poincare_map['R'],poincare_map['Z'])
    R,Z=R.flatten(),Z.flatten()

    psi_norm_2D=processing.utils.value_at_RZ(
        R=R,
        Z=Z,
        quantity=equilibrium['psirz_norm'],
        grid=equilibrium
        ).reshape(
            len(poincare_map['Z']),
            len(poincare_map['R'])
        ).T

    theta_2D=processing.utils.angle_pol(
        R_major=equilibrium['rmaxis'],
        R=R,
        Z=Z,
        Z_major=equilibrium['zmaxis']
        ).reshape(
            len(poincare_map['Z']),
            len(poincare_map['R'])
        ).T*180./np.pi
    theta_2D[theta_2D>180]-=360

    phi_slice=np.argmin(np.abs(poincare_map['phi']-phi))
    inds=np.where((poincare_map['map'][:,:,phi_slice].flatten()==1))[0]

    ax.scatter(
        theta_2D.flatten()[inds][::1],
        psi_norm_2D.flatten()[inds][::1],
        s=0.05
    )
    ax.set_ylim([0.,1.])
    ax.set_xlim([-180.,180.])
    ax.set_xlabel(r'$\theta$ [deg]')
    ax.set_ylabel(r'$\psi_{\mathrm{n}}$')
    
    if ax_flag is False and fig_flag is False:
        plt.show()

def plot_poincare_q_theta(poincare_map,equilibrium,phi=0.,ax=False,fig=False):

    #interpolate poincare R Z to q theta

    if ax is False:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if fig is False:
        fig_flag=False
    else:
        fig_flag=True

    if fig_flag is False:
        fig = plt.figure() #if user has not externally supplied figure, generate
    
    if ax_flag is False: #if user has not externally supplied axes, generate them
        ax = fig.add_subplot(111)

    R,Z=np.meshgrid(poincare_map['R'],poincare_map['Z'])
    R,Z=R.flatten(),Z.flatten()

    q_rz=processing.utils.flux_func_to_RZ(psi=equilibrium['flux_pol'],quantity=equilibrium['qpsi'],equilibrium=equilibrium)
    q_2D=processing.utils.value_at_RZ(R=R,Z=Z,quantity=q_rz,grid=equilibrium).reshape(
        len(poincare_map['Z']),
        len(poincare_map['R'])).T

    theta_2D=processing.utils.angle_pol(
        R_major=equilibrium['rmaxis'],
        R=R,
        Z=Z,
        Z_major=equilibrium['zmaxis']
        ).reshape(
            len(poincare_map['Z']),
            len(poincare_map['R'])
        ).T*180./np.pi
    theta_2D[theta_2D>180]-=360

    phi_slice=np.argmin(np.abs(poincare_map['phi']-phi))
    inds=np.where((poincare_map['map'][:,:,phi_slice].flatten()==1))[0]

    ax.scatter(
        theta_2D.flatten()[inds][::1],
        q_2D.flatten()[inds][::1],
        s=0.05
    )
    #ax.set_ylim([0.,1.])
    ax.set_xlim([-180.,180.])
    ax.set_xlabel(r'$\theta$ [deg]')
    ax.set_ylabel(r'$q$')
    
    if ax_flag is False and fig_flag is False:
        plt.show()


def plot_perturbation_phi_theta(perturbations,equilibrium,psi,vminmax=None,number_bins=100,colmap=settings.cmap_default,ax=None,fig=None):

    # extract contour coordinates at psi in equilibrium
    X=equilibrium['R_1D'] #make a mesh
    Y=equilibrium['Z_1D']
    dx,dy=X[1]-X[0],Y[1]-Y[0]
    Y,X=np.meshgrid(Y-dy/2.,X-dx/2.) #offset ticks onto bin centres
    Z=(equilibrium['psirz']-equilibrium['simag'])/(equilibrium['sibry']-equilibrium['simag']) 
    levels=[
            psi
            ]
    fig_,ax_=plt.subplots(1)
    #ax.contourf(X,Y,Z,cmap='Greys')
    mesh=ax_.contour(X,Y,Z,levels=levels)

    if ax is None:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if fig is None:
        fig_flag=False
    else:
        fig_flag=True

    if fig_flag is False:
        fig = plt.figure() #if user has not externally supplied figure, generate

    if ax_flag is False: #if user has not externally supplied axes, generate them
        ax = fig.add_subplot(111)

    levels_coords=[]
    #total number of points to interpolate onto
    #this will be spread out over the contour's sub-paths (if any) in proportion to their length
    number_coords_new=400

    for contour in mesh.collections: 
        contour_x=[]
        contour_y=[]
        contour_coords_x=[]
        contour_coords_y=[]
        path_lengths=[]
        #scan each separate path to find longest
        #assume longest is closed flux surface and get rid of rest
        for path in contour.get_paths(): 
            v = path.vertices
            path_x=v[:,0]
            path_y=v[:,1]        
            coords=np.array([path_x,path_y])
            distance_total=np.cumsum(np.sqrt(np.sum(np.diff(coords,axis=1)**2,axis=0)))[-1]
            path_lengths.append(distance_total)
        for path_number,path in enumerate(contour.get_paths()): 
            if path_lengths[path_number]==np.max(np.array(path_lengths)):
                v = path.vertices
                contour_x.extend(v[:,0])
                contour_y.extend(v[:,1])
                path_x=v[:,0]
                path_y=v[:,1]        
                #interpolate along trajectory
                coords=np.array([path_x,path_y])
                distance=np.cumsum(np.sqrt(np.sum(np.diff(coords,axis=1)**2,axis=0)))
                distance=np.insert(distance,0,0.) #first distance is 0
                distance=(distance-distance[0])/(distance[-1]-distance[0]) #normalise
                interpolator=processing.utils.interpolate_1D(distance,coords,function='linear',type='interp1d',smooth=0)
                distance_coords_new=np.linspace(0,1,number_coords_new+1)[:-1] #this loops back onto itself, so add  extra point and remove the end
                coords_new=interpolator(distance_coords_new)
                contour_coords_x.extend(coords_new[0,:])
                contour_coords_y.extend(coords_new[1,:])
                #ax.scatter(*coords_new,color='y')
        levels_coords.append([contour_coords_x,contour_coords_y])
        #ax.scatter(contour_x,contour_y,color='r')
        # so now levels_coords[n,0,p] contains the x coordinate for point p along contour n
        # so now levels_coords[n,1,p] contains the y coordinate for point p along contour n
    #theta=np.linspace(0,2.*np.pi)
    #phi=np.linspace(0,2.*np.pi)
    #theta_2D,phi_2D=np.meshgrid(theta,phi)

    phi=np.linspace(-np.pi,np.pi,450)

    levels_coords=np.array(levels_coords)
    R_contour=levels_coords[0,0,:]
    Z_contour=levels_coords[0,1,:]
    R_flat=np.array(list(R_contour)*len(phi)).flatten()
    Z_flat=np.array(list(Z_contour)*len(phi)).flatten()
    phi_flat=np.array([np.full(len(R_contour),angle) for angle in phi]).flatten()

    dB_R=np.zeros(len(R_flat))
    dB_tor=np.zeros(len(R_flat))
    dB_Z=np.zeros(len(R_flat))
    for perturbation in perturbations:
        dB=perturbation.evaluate(R=R_flat,phi=phi_flat,Z=Z_flat,phase=0,i3dr=-1,mode_number=perturbation.mode_number) #evaluate poloidally
        dB_R+=dB[0]
        dB_tor+=dB[1]
        dB_Z+=dB[2]

    dB_R=dB_R.reshape(len(R_contour),len(phi),order='F')
    dB_tor=dB_tor.reshape(len(R_contour),len(phi),order='F')
    dB_Z=dB_Z.reshape(len(R_contour),len(phi),order='F')
    values=np.sqrt(dB_R**2+dB_tor**2+dB_Z**2)

    theta=processing.utils.angle_pol(R=R_contour,Z=Z_contour,R_major=0.639277243e1,Z_major=0.597384943)    
    theta[theta>np.pi]-=2.*np.pi

    if vminmax is None: vminmax=[np.min(values),np.max(values)]
    vmin,vmax=(_ for _ in vminmax)

    mesh=ax.contourf(phi*180/np.pi,theta*180/np.pi,values,levels=np.linspace(vmin,vmax,num=number_bins),colors=colmap(np.linspace(0.,1.,num=number_bins)),edgecolor='none',antialiased=True,vmin=vmin,vmax=vmax)

    if fig_flag is False:    
        fig.colorbar(mesh,ax=ax,orientation='horizontal')

    if ax_flag is True or fig_flag is True: #return the plot object
        return mesh

    if ax_flag is False and fig_flag is False:
        plt.show()

#functions to apply when reading files in parallel

def calc_PFC_power(filename,classtype='fpl'):
    io_obj=read_locust_io_obj(filename=filename,classtype=classtype)
    if io_obj:
        try:
            if classtype=='fpl':
                i=np.where(io_obj['status_flag']=='PFC_intercept_3D')[0]
                PFC_power=1.e6*io_obj['f']*np.sum((io_obj['V_R'][i]**2+io_obj['V_phi'][i]**2+io_obj['V_Z'][i]**2)*io_obj['FG'][i])*0.5*constants.mass_deuteron
            elif classtype=='rund':
                PFC_power=io_obj['PFC_power']['total']
        except:
            PFC_power=-1
    else:
        PFC_power=-1

    return PFC_power  

def calc_mean_loss_time(filename,classtype='fpl'):
    fpl=read_locust_io_obj(filename=filename,classtype=classtype)
    if fpl:
        i=np.where(fpl['status_flag']=='PFC_intercept_3D')[0]
        mean_loss_time=np.average(fpl['time'][i],weights=fpl['FG'][i])
    else:
        mean_loss_time=-1
    return mean_loss_time  

#################################
 
##################################################################

###################################################################################################