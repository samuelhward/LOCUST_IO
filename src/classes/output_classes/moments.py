#moments.py

"""
Samuel Ward
27/04/2018
----
class to handle LOCUST moments output data
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

def read_moments_LOCUST(filepath,**properties):
    """
    reads generic moments data output from LOCUST

    notes:

    """

    print("reading moments from LOCUST")

    input_data={} #initialise data dictionary

    try:
        import h5py
    except:
        raise ImportError("ERROR: read_moments_LOCUST could not import h5py!\nreturning\n")
        return

    with h5py.File(filepath,'r') as file:        

        input_data['flux_pol_norm_sqrt']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['sqrt(PSIn)']) 
        input_data['flux_pol_norm']=input_data['flux_pol_norm_sqrt']**2
        input_data['dVOL']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['dVOL'])
        input_data['beam_source']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Beam Source']['data'])
        input_data['beam_source_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Beam Source']['error'])
        input_data['density']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Density']['data'])
        input_data['density_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Density']['error'])
        input_data['energy_para']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['E (para)']['data'])
        input_data['energy_para_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['E (para)']['error'])
        input_data['energy_perp']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['E (perp)']['data'])
        input_data['energy_perp_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['E (perp)']['error'])
        input_data['energy']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Energy']['data'])
        input_data['energy_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Energy']['error'])
        input_data['J(NBCD)-raw']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['J (NBCD) (raw)']['data'])
        input_data['J(NBCD)-raw_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['J (NBCD) (raw)']['error'])
        input_data['NBI-heating-power(TOT)']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['NBI Heating Power (TOT)']['data'])
        input_data['NBI-heating-power(TOT)_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['NBI Heating Power (TOT)']['error'])
        input_data['NBI-heating-power(e-)']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['NBI Heating Power (e-)']['data'])
        input_data['NBI-heating-power(e-)_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['NBI Heating Power (e-)']['error'])
        input_data['NBI-heating-power(i1)']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['NBI Heating Power (i1)']['data'])
        input_data['NBI-heating-power(i1)_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['NBI Heating Power (i1)']['error'])
        input_data['residual-angular-momentum-density']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Residual Angular Momentum Density']['data'])
        input_data['residual-angular-momentum-density_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Residual Angular Momentum Density']['error'])
        input_data['torque-density(JxB-inst)']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Torque Density (JxB inst)']['data'])
        input_data['torque-density(JxB-inst)_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Torque Density (JxB inst)']['error'])
        input_data['torque-density(JxB-sweep)']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Torque Density (JxB sweep)']['data'])
        input_data['torque-density(JxB-sweep)_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Torque Density (JxB sweep)']['error'])
        input_data['torque-density(coll)']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Torque Density (coll)']['data'])
        input_data['torque-density(coll)_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Torque Density (coll)']['error'])
        
        ''' out of date or not needed for now

        input_data['temperature_e']=np.array(file['Input Data']['Kinetic Data']['Profiles (1D)']['Te']['data']) 
        input_data['temperature_i']=np.array(file['Input Data']['Kinetic Data']['Profiles (1D)']['Ti']['data']) 
        input_data['density_e']=np.array(file['Input Data']['Kinetic Data']['Profiles (1D)']['ne']['data']) 

        input_data['NBI-heating-power(i1)A']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['NBI-heating-power(i1)']['A'])
        input_data['NBI-heating-power(i1)Z']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['NBI-heating-power(i1)']['Z'])
        input_data['NBI-heating-power(i2)']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['NBI-heating-power(i2)']['data'])
        input_data['NBI-heating-power(i2)_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['NBI-heating-power(i2)']['error'])
        input_data['NBI-heating-power(i2)Z']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['NBI-heating-power(i2)']['Z'])
        input_data['NBI-heating-power(i2)A']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['NBI-heating-power(i2)']['A'])
        input_data['P(para)']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['P(para)']['data'])
        input_data['P(para)_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['P(para)']['error'])
        input_data['P(perp)']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['P(perp)']['data'])
        input_data['P(perp)_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['P(perp)']['error'])
        input_data['e-source']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['e-Source']['data'])
        input_data['e-source_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['e-Source']['error'])
        '''

    print("finished reading moments from LOCUST")

    return input_data

def read_moments_TRANSP(filepath,**properties):
    """
    reads moments from TRANSP <run_id>.CDF output file
    
    notes:
        time-resolved moments here are written out at corresponding TRANSP OUTTIMs
    """

    print("reading moments from TRANSP")

    try:
        from scipy.io import netcdf as ncdf
    except:
        raise ImportError("ERROR: read_moments_TRANSP could not import scipy.io.netcdf!\nreturning\n")
        return

    file=ncdf.netcdf_file(filepath,'r')
    input_data={}

    input_data['density']=np.array(file.variables['BDENS'].data)*1.e6 #convert to [m^-3]

    input_data['energy']=np.array(file.variables['UTHRM'].data)*1.e6 #convert to [m^-3]     thermal energy density
    input_data['energy_para']=np.array(file.variables['UFASTPA'].data)*1.e6 #convert to [m^-3]
    input_data['energy_perp']=np.array(file.variables['UFASTPP'].data)*1.e6 #convert to [m^-3]

    input_data['NBI-heating-power(i1)']=np.array(file.variables['PBI'].data)*1.e6 #convert to [m^-3]
    input_data['NBI-heating-power(e-)']=np.array(file.variables['PBE'].data)*1.e6 #convert to [m^-3]

    #input_data['torque-density(JxB-inst)']=np.array(file.variables['TQJJXBT'].data)*1.e6 #convert to [m^-3]
    input_data['torque-density(JxB-inst)_deposited']=np.array(file.variables['TQJBD'].data)*1.e6 #convert to [m^-3]

    input_data['torque-density(coll)']=np.array(file.variables['BPHCL'].data)

    input_data['residual-angular-momentum-density']=np.array(file.variables['BPHI'].data)

    input_data['beam_source']=np.array(file.variables['BDEP_D'].data)*1.e6 #convert to [m^-3]
    input_data['beam_source_captured']=np.array(file.variables['BPCAP'].data)[:,np.newaxis]
    input_data['PFC_power']=np.array(file.variables['BPLIM'].data)[:,np.newaxis]

    input_data['time']=np.array(file.variables['TIME3'].data)

    input_data['r/a ctr']=np.array(file.variables['X'].data)
    input_data['r/a bdy']=np.array(file.variables['XB'].data)
    input_data['flux_pol']=np.array(file.variables['PLFLX'].data) #Wb/rad
    input_data['flux_pol_norm']=(input_data['flux_pol']-input_data['flux_pol'][0])/(input_data['flux_pol'][-1]-input_data['flux_pol'][0])
    input_data['flux_pol_norm_sqrt']=np.sqrt(np.array(input_data['flux_pol_norm']))

    '''out of date or not needed for now
    input_data['torque-density(coll)']=file.variables['TQCOLNB'].data
    input_data['temperature_e']=file.variables['TE'].data
    input_data['temperature_i']=file.variables['TI'].data
    input_data['density_e']=file.variables['NE'].data
    '''

    file.close()
    file.close()
    
    print("finished reading moments from TRANSP")

    return input_data

def read_moments_ASCOT(filepath,**properties):
    """
    reads moments from ASCOT hdf5 output file

    notes:
        assumes si units
    """

    print("reading moments from ASCOT")

    try:
        import h5py
    except:
        raise ImportError("ERROR: read_moments_ASCOT could not import h5py!\nreturning\n")
    
    with h5py.File(filepath,'r') as file:
    
        input_data={} #initialise data dictionary

        input_data['flux_pol_norm_sqrt']=file['distributions/rhoDist/abscissae/dim1'].value
        input_data['flux_pol_norm']=input_data['flux_pol_norm_sqrt']**2

        #convert flux coordinates to bin centres
        for flux_coordinate in ['flux_pol_norm','flux_pol_norm_sqrt']:
            input_data[flux_coordinate]=(input_data[flux_coordinate][:-1]+input_data[flux_coordinate][1:])*0.5

        data=file['distributions/rhoDist/ordinate'].value
        data=np.swapaxes(data,-1,0)
        data=np.swapaxes(data,-2,1)
        data=np.squeeze(data) #so now we should have data[profile,flux_pol_norm,...]

        for profile in file['distributions/rhoDist/ordinates'].keys():    
            if profile[0:4]=='unit':
                pass
            elif profile[0:4]=='name':
                profile_path='distributions/rhoDist/ordinates/'+profile
                profile_number=int(profile[-1])
                profile_name=file[profile_path].value.decode('utf-8')
                input_data[profile_name]=data[profile_number,:]

        #rename variables to LOCUST_IO conventions
        ascot_names=['Power deposition to electrons','Power deposition to background species  1','parallel energy density','energy density','Pphi torque','total toroidal current density','jxB Torque']
        locust_io_names=['NBI-heating-power(e-)','NBI-heating-power(i1)','energy_para','energy','torque-density(JxB-inst)','J(NBCD)-raw','torque-density(JxB-sweep)']

        for ascot_name,locust_io_name in zip(ascot_names,locust_io_names):
            input_data[locust_io_name]=input_data.pop(ascot_name)

    print("finished reading moments from ASCOT")

    return input_data

def read_moments_IDS(shot,run,**properties):
    """
    reads LOCUST moments distributions IDS
 
    args: 
        shot - IDS shot number 
        run - IDS run number

    notes:
        assumes reading from first distribution
    """

    print("reading moments from IDS")

    try:
        import imas 
    except:
        raise ImportError("ERROR: read_moments_IDS could not import IMAS module!\nreturning\n")
        return
    
    input_IDS=imas.ids(int(shot),int(run)) #initialise new blank IDS
    input_IDS.open_env(settings.username,settings.imasdb,settings.imas_version)
    input_IDS.distributions.get() #open the file and get all the data from it

    input_data = {} #initialise blank dictionary to hold the data

    input_data['flux_pol_norm_sqrt']=np.sqrt(np.array(input_IDS.distributions.distribution[0].profiles_1d[0].grid.psi)) 
    input_data['flux_pol_norm']=np.array(input_IDS.distributions.distribution[0].profiles_1d[0].grid.psi)
    input_data['dVOL']=np.array(input_IDS.distributions.distribution[0].profiles_1d[0].grid.volume)
    input_data['beam_source']=np.array(input_IDS.distributions.distribution[0].profiles_1d[0].source[0].particles)
    #input_data['beam_source_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Beam Source']['error'])
    input_data['density']=np.array(input_IDS.distributions.distribution[0].profiles_1d[0].density_fast)
    #input_data['density_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Density']['error'])
    input_data['energy_para']=np.array(input_IDS.distributions.distribution[0].profiles_1d[0].pressure_fast_parallel)
    #input_data['energy_para_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['E (para)']['error'])
    input_data['energy_perp']=np.array(input_IDS.distributions.distribution[0].profiles_1d[0].pressure_fast)-input_data['energy_para']
    #input_data['energy_perp_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['E (perp)']['error'])
    input_data['energy']=np.array(input_IDS.distributions.distribution[0].profiles_1d[0].source[0].energy)
    #input_data['energy_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Energy']['error'])
    input_data['J(NBCD)-raw']=np.array(input_IDS.distributions.distribution[0].profiles_1d[0].current_fast_tor)
    #input_data['J(NBCD)-raw_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['J (NBCD) (raw)']['error'])
    input_data['NBI-heating-power(TOT)']=np.array(input_IDS.distributions.distribution[0].profiles_1d[0].thermalisation.energy)
    #input_data['NBI-heating-power(TOT)_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['NBI Heating Power (TOT)']['error'])
    input_data['NBI-heating-power(e-)']=np.array(input_IDS.distributions.distribution[0].profiles_1d[0].collisions.electrons.power_fast)
    #input_data['NBI-heating-power(e-)_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['NBI Heating Power (e-)']['error'])
    input_data['NBI-heating-power(i1)']=np.array(input_IDS.distributions.distribution[0].profiles_1d[0].collisions.ion[0].power_fast)
    #input_data['NBI-heating-power(i1)_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['NBI Heating Power (i1)']['error'])
    input_data['residual-angular-momentum-density']=np.array(input_IDS.distributions.distribution[0].profiles_1d[0].source[0].momentum_tor)
    #input_data['residual-angular-momentum-density_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Residual Angular Momentum Density']['error'])
    input_data['torque-density(JxB-inst)']=np.array(input_IDS.distributions.distribution[0].profiles_1d[0].torque_tor_j_radial)
    #input_data['torque-density(JxB-inst)_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Torque Density (JxB inst)']['error'])
    #input_data['torque-density(JxB-sweep)']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Torque Density (JxB sweep)']['data'])
    #input_data['torque-density(JxB-sweep)_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Torque Density (JxB sweep)']['error'])
    #input_data['torque-density(coll)']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Torque Density (coll)']['data'])
    #input_data['torque-density(coll)_err']=np.array(file['Output Data']['Fast Ions']['Profiles (1D)']['Torque Density (coll)']['error'])

    input_IDS.close()

    print("finished reading moments from IDS")
 
    return input_data


################################################################## Orbits class

class Moments(classes.base_output.LOCUST_output):
    """
    class describing a generic moments output for LOCUST
    
    inherited from LOCUST_output:
        self.ID                     unique object identifier, good convention to fill these for error handling etc
        self.data                   holds all output data in dictionary object
        self.LOCUST_output_type     string which holds this class' output type, this case = 'moments'
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

    LOCUST_output_type='moments'

    def read_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        read moments from file 

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
                self.data=read_moments_LOCUST(self.filepath,**properties)

        elif data_format=='TRANSP':
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot read_data() from TRANSP - filename required\n".format(self.ID),filename): #must check we have all info required for reading

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files / filename
                self.properties={**properties}
                self.data=read_moments_TRANSP(self.filepath,**properties)

        elif data_format=='ASCOT': #here are the blocks for various file types, they all follow the same pattern
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot read_data() from ASCOT - filename required\n".format(self.ID),filename): #must check we have all info required for reading

                self.data_format=data_format #add to the member data
                self.filename=filename
                self.filepath=support.dir_output_files / filename
                self.properties={**properties}
                self.data=read_moments_ASCOT(self.filepath,**properties)

        elif data_format=='IDS':
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot read_data() from distributions IDS - shot and run required\n".format(self.ID),shot,run):
 
                self.data_format=data_format
                self.shot=shot
                self.run=run
                self.properties={**properties}
                self.data=read_moments_IDS(self.shot,self.run,**properties)

        else:
            print("ERROR: {} cannot read_data() - please specify a compatible data_format (LOCUST/TRANSP/ASCOT/IDS)\n".format(self.ID))            

    def dump_data(self,data_format=None,filename=None,shot=None,run=None,**properties):
        """
        write moments to file

        notes: 
        """
        if processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot dump_data() - self.data and compatible data_format required\n".format(self.ID),self.data,data_format):
            pass
        
        elif data_format=='LOCUST':
            if not processing.utils.none_check(self.ID,self.LOCUST_output_type,"ERROR: {} cannot dump_data() to LOCUST - filename required\n".format(self.ID),filename):
                filepath=support.dir_output_files / filename
                dump_moments_LOCUST(self.data,filepath,**properties)
        else:
            print("ERROR: {} cannot dump_data() - please specify a compatible data_format (LOCUST)\n".format(self.ID))

    def plot(self,key,axis='flux_pol_norm',colmap=settings.cmap_default,colmap_val=np.random.uniform(),ax=False,fig=False):
        """
        plots moments

        notes:
            key - selects which data to plot
            axis - selects x axis of plot
            colmap - set the colour map (use get_cmap names)
            colmap_val - optional numerical value for defining single colour plots 
            ax - take input axes (can be used to stack plots)
            fig - take input fig (can be used to add colourbars etc)
        """

        import matplotlib
        from matplotlib import cm
        import matplotlib.pyplot as plt

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
        ax.set_title(self.ID)
       
        ax.plot(self[axis],self[key],color=colmap(colmap_val))
        ax.set_xlabel(axis)
        ax.set_ylabel(key)

        if ax_flag is False and fig_flag is False:
            plt.show()
 
#################################
 
##################################################################
 
###################################################################################################


#################################

##################################################################

###################################################################################################
