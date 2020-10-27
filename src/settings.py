#settings.py

'''
Samuel Ward
18/10/2018
----
File which holds LOCUST_IO settings
---
notes: 
---
'''


##################################################################
#Preamble


try:
    import subprocess
    import numpy as np
    import matplotlib
    from matplotlib import cm
    import getpass
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

##################################################################

#general
np.set_printoptions(precision=5,threshold=3) #set printing style of numpy arrays

#plotting
def cmap_custom(from_rgb,to_rgb):
    """
    generate custom colormaps
    args:
        from_rgb - list of starting r g b values  
        to_rgb - list of final r g b values
    notes:
    """
    from matplotlib.colors import LinearSegmentedColormap
    r1,g1,b1=from_rgb
    r2,g2,b2=to_rgb
    cdict={'red':((0,r1,r1),(1,r2,r2)),
            'green':((0,g1,g1),(1,g1,g1)),
            'blue':((0,b1,b1),(1,b1,b1))}
    cmap=LinearSegmentedColormap('custom cmap - {from_rgb}/{to_rgb}'.format(from_rgb=str(from_rgb),to_rgb=str(to_rgb)),cdict)
    return cmap

def colour_custom(rgba,N=256):
    """
    generate custom colours for colmaps

    args:
        rgba - array containing red, green, blue, alpha values 
        N - depth
    notes:
    """
    from matplotlib.colors import ListedColormap
    return ListedColormap(np.array(rgba)/np.array([N,N,N,1]))

def discrete_colmap(colmap_name,face_colour='white',number_bins=10):
    """
    return discretised version of matplotlib colour map

    args:
        colmap_name - name string of matplotlib colour map
        face_colour - name string of face colour
        number_bins - set number of bins or levels
    notes:
    """

    colmap=matplotlib.cm.get_cmap(colmap_name,number_bins)
    colmap.set_under(face_colour)
    colmap=list(colmap(np.arange(number_bins)))
    colmap[0]=face_colour
    colmap=matplotlib.colors.ListedColormap(colmap[:-1], "")
    return colmap

cmap_plasma=matplotlib.cm.get_cmap('plasma')
cmap_plasma_r=matplotlib.cm.get_cmap('plasma_r')
cmap_jet=matplotlib.cm.get_cmap('jet')
cmap_jet_r=matplotlib.cm.get_cmap('jet_r')
cmap_viridis=matplotlib.cm.get_cmap('viridis')
cmap_viridis_r=matplotlib.cm.get_cmap('viridis_r')
cmap_inferno=matplotlib.cm.get_cmap('inferno')
cmap_inferno_r=matplotlib.cm.get_cmap('inferno_r')
cmap_r=cmap_custom([1,0,0],[1,0,0]) #red
cmap_g=cmap_custom([0,1,0],[0,1,0]) #green
cmap_b=cmap_custom([0,0,1],[0,0,1]) #blue
cmap_y=cmap_custom([1,1,0],[1,1,0]) #yellow
cmap_m=cmap_custom([1,0,1],[1,0,1]) #magenta
cmap_c=cmap_custom([0,1,1],[0,1,1]) #cyan
cmap_w=cmap_custom([1,1,1],[1,1,1]) #white
cmap_k=cmap_custom([0,0,0],[0,0,0]) #black
cmap_default=discrete_colmap(colmap_name='inferno_r',face_colour='white',number_bins=10) #create default colourmap

plot_colour_LCFS=cmap_m(0.)
plot_line_style_LCFS='-' 
plot_colour_limiters=cmap_k(0.)
plot_line_style_limiters='-' 
plot_colour_gridlines=cmap_w(0.) 
plot_line_style_gridlines='-'
plot_line_style='solid' #set default plot line style
plot_contour_labels=False #toggle level labels for contour plots
colour_start_mark='red' #set plot starting markers default
marker_start_mark='o' 
markersize_start_mark=1
tick_frequency=5 #plot every Nth tick

matplotlib_rc={} #matplotlib RC settings
matplotlib_rc['font'] = {
                        'family'     : 'Bitstream Vera Sans', 
                        'serif'      : 'DejaVu Serif', 
                        'sans-serif' : 'Arial', 
                        'weight'     : 'normal',
                        'size'       :  20
                        }
matplotlib_rc['lines'] = {
                        'linewidth' : 2.,
                        'linestyle' : 'solid'
                        }
matplotlib_rc['axes'] = {    
                        'labelsize' : 30,
                        'titlesize' : 30
                        }
matplotlib_rc['figure'] = {    
                        'autolayout' : True,
                        'titlesize' : 30,
                        'figsize' : [10.,8.]
                        }
matplotlib_rc['legend'] = {
                        'loc' : 'best'    
                        }
for setting_type,setting in matplotlib_rc.items(): matplotlib.rc(setting_type, **setting) #enable settings

#system environment
try:
    username=getpass.getuser()
except:
    username='wards2'
imasdb='test'
system_default='TITAN'
imas_version='3'

#LOCUST git
repo_URL_LOCUST='ssh://git@git.iter.org/traj/locust.git'
branch_default_LOCUST='hot_fix/ITER'

#LOCUST_IO git
repo_URL_LOCUST_IO='https://github.com/armoured-moose/LOCUST_IO.git'
branch_default_LOCUST_IO='develop'

#NEMO git
repo_URL_NEMO='ssh://git@git.iter.org/heat/nemo.git'
branch_default_NEMO='develop'

#BBNBI git
repo_URL_BBNBI='ssh://git@git.iter.org/traj/ascot.git'
branch_default_BBNBI='develop'

commit_hash_default_LOCUST=None
commit_hash_default_LOCUST_IO=None
try: #get default commit hash
    git_ls_remote=subprocess.run(['git', 'ls-remote','-q',repo_URL_LOCUST],stdout=subprocess.PIPE).stdout.decode('utf-8').split()
    for counter,entry in enumerate(git_ls_remote):
        if branch_default_LOCUST in entry: #look for this branch in command output
            commit_hash_default_LOCUST=git_ls_remote[counter-1] #output of this command is in two columns (hash and branch)
except:
    commit_hash_default_LOCUST=None
try:
    git_ls_remote=subprocess.run(['git', 'ls-remote','-q',repo_URL_LOCUST_IO],stdout=subprocess.PIPE).stdout.decode('utf-8').split()
    for counter,entry in enumerate(git_ls_remote):
        if branch_default_LOCUST_IO in entry: #look for this branch in command output
            commit_hash_default_LOCUST_IO=git_ls_remote[counter-1] #output of this command is in two columns (hash and branch)
except:
    commit_hash_default_LOCUST_IO=None

#LOCUST source
LOCUST_dir_inputfiles_default='InputFiles' #default folder names according to LOCUST source code - src.run_scripts.LOCUST_run will create these directories during a run
LOCUST_dir_outputfiles_default='OutputFiles'
LOCUST_dir_cachefiles_default='CacheFiles'
LOCUST_dir_root_default='/tmp'

#################################

##################################################################

###################################################################################################
