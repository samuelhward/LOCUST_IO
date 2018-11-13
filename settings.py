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
import numpy
import matplotlib
from matplotlib import cm

##################################################################


#general
numpy.set_printoptions(precision=5,threshold=3) #set printing style of numpy arrays

#plotting
cmap_default=matplotlib.cm.get_cmap('jet') #set default colourmap
plot_style_LCFS='m-' #set plot style for LCFS
plot_style_limiters='w-' #set plot style for limiters
font = {'family' : 'normal', #set figure font
        'weight' : 'bold',
        'size'   :  20}
matplotlib.rc('font', **font)
colour_start_mark='red' #set colour of starting markers in plots

#################################

##################################################################

###################################################################################################