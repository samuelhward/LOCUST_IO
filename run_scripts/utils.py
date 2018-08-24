#utils.py
 
"""
Samuel Ward
24/08/2018
----
supporting functions for LOCUST run scripts
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
    import os
    import subprocess
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning")
    sys.exit(1)
try:
    from classes import support
except:
    raise ImportError("ERROR: support.py could not be imported!\nreturning") 
    sys.exit(1)

np.set_printoptions(precision=5,threshold=5) #set printing style of numpy arrays
cmap_default=matplotlib.cm.get_cmap('jet') #set default colourmap
plot_style_LCFS='m-' #set plot style for LCFS
plot_style_limiters='w-' #set plot style for limiters

pi=np.pi
e_charge=1.60217662e-19 #define electron charge
mass_neutron=1.674929e-27 #define mass of neutron
amu=1.66053904e-27
mass_deuterium=2.0141017781*amu




##################################################################
#Main Code


def TRANSP_get_fbm_FI_CDF(run_ID,number_files,particle_position=True,guiding_centre=True,device='d3d'):
    """
    notes:
        looks for files in output_files
    args:
        run_ID - full transp run number e.g. 157418S01
        number_files - total number of .DATA# files to extract CDF from 
        particle_position - toggle whether to generate set of CDFs at particle positions
        guiding_centre - toggle whether to generate set of CDFs at guiding centres
        device - device code for machine under study
    """

    project_dir=os.getcwd()
    os.chdir(support.dir_output_files) #change working directory to output files briefly due to bug in CCFE get_fbm installation
    
    for file_ID in range(number_files):
        file_ID+=1
        output_filename='{}{}{}{}'.format(run_ID,'_fi_',file_ID,'.cdf') #current output file
        output_filename_gc='{}{}{}{}'.format(run_ID,'_fi_',file_ID,'_gc.cdf')

        if guiding_centre:
            fbm_input = """
                        {run_ID}
                        {path}
                        {file_ID}
                        t
                        {device}
                        w
                        c
                        """.format(run_ID=run_ID,path='q',file_ID=file_ID,device=device)

            print("writing TRANSP FI netCDF file {}".format(output_filename_gc))
            proc = subprocess.Popen(['get_fbm'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
            out, err = proc.communicate(input=fbm_input.encode('utf-8'))
            os.rename(output_filename,output_filename_gc) #add '_gc' label onto file
            print("finished writing TRANSP FI netCDF file {}".format(output_filename_gc))

        if particle_position:
            fbm_input = """
                        {run_ID}
                        {path}
                        {file_ID}
                        t
                        {device}
                        w
                        p
                        """.format(run_ID=run_ID,path='q',file_ID=file_ID,device=device)

            print("writing TRANSP FI netCDF file {}".format(output_filename))
            proc = subprocess.Popen(['get_fbm'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
            out, err = proc.communicate(input=fbm_input.encode('utf-8'))
            print("finished writing TRANSP FI netCDF file {}".format(output_filename))

    os.chdir(project_dir) #change back to original working directory    


def TRANSP_get_fbm_FI_birth_deposition(run_ID,number_files,device='d3d'):
    """
    notes:
        looks for files in output_files
        assumes:
            all beams
            all energy components
            dumping at particle location
            dumping in x y z space
            uses Akima Hermite spline interpolation
            random seed = 1
            sample size = 1,000,000
    args:
        run_ID - full transp run number e.g. 157418S01
        number_files - total number of .DATA# files to extract CDF from 
        device - device code for machine under study
    """

    project_dir=os.getcwd()
    os.chdir(support.dir_output_files) #change working directory to output files briefly due to bug in CCFE get_fbm installation
    
    for file_ID in range(number_files):
        file_ID+=1
        output_filename='{}{}{}{}'.format(run_ID,'_fdep_nb_en_',file_ID,'.out') #current output file

        fbm_input = """
                    {run_ID}
                    {path}
                    {file_ID}
                    t
                    {device}
                    1
                    d
                    n
                    0
                    0
                    p
                    b
                    c
                    Y
                    1
                    1000000
                    {output_filename}
                    """.format(run_ID=run_ID,path='q',file_ID=file_ID,device=device,output_filename=output_filename)

        print("writing TRANSP FI random sample birth deposition {}".format(output_filename))
        proc = subprocess.Popen(['get_fbm'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        out, err = proc.communicate(input=fbm_input.encode('utf-8'))
        print("finished writing TRANSP FI random sample birth deposition {}".format(output_filename))

    os.chdir(project_dir) #change back to original working directory  


#################################
 
##################################################################
 
###################################################################################################