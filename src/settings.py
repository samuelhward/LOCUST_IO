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
import numpy as np
import matplotlib
from matplotlib import cm

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

cmap_plasma=matplotlib.cm.get_cmap('plasma')
cmap_jet=matplotlib.cm.get_cmap('jet')
cmap_viridis=matplotlib.cm.get_cmap('viridis')
cmap_r=cmap_custom([1,0,0],[1,0,0]) #red
cmap_g=cmap_custom([0,1,0],[0,1,0]) #green
cmap_b=cmap_custom([0,0,1],[0,0,1]) #blue
cmap_y=cmap_custom([1,1,0],[1,1,0]) #yellow
cmap_m=cmap_custom([1,0,1],[1,0,1]) #magenta
cmap_c=cmap_custom([0,1,1],[0,1,1]) #cyan
cmap_w=cmap_custom([1,1,1],[1,1,1]) #white
cmap_k=cmap_custom([0,0,0],[0,0,0]) #black
cmap_default=cmap_plasma #set default colourmap

plot_style_LCFS='m-' #set plot style for LCFS
plot_style_limiters='k-' #set plot style for limiters
plot_style_gridlines='w'
plot_linewidth=0.5 #set plot line width
plot_contour_labels=False #toggle level labels for contour plots
font = {'family' : 'normal', #set figure font
        'weight' : 'bold',
        'size'   :  10}
matplotlib.rc('font', **font)
colour_start_mark='red' #set plot starting markers default
marker_start_mark='o' 
markersize_start_mark=1

#system environment
username='wards2'
imasdb='test'
system_default='TITAN'

#LOCUST git
repo_URL_LOCUST='ssh://git@git.iter.org/traj/locust.git'
repo_URL_LOCUST_IO='https://github.com/armoured-moose/LOCUST_IO.git'
branch_default_LOCUST='hot_fix/ITER'
branch_default_LOCUST_IO='develop'

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