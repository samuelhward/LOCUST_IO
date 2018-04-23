#process_output.py

'''
Samuel Ward
25/1/2018
----
processing routines for LOCUST output data
---
notes:
	https://stackoverflow.com/questions/9455111/python-define-method-outside-of-class-definition
---
'''

##################################################################
#Preamble

try:
    import scipy.integrate
    import numpy as np
    import copy
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)
    
pi=np.pi

##################################################################
#Main Code

def dfn_integrate(some_dfn):
    """
    integrate a distribution function to get particles/cell

    notes:
        returns the Dfn array
    """

    dfn=copy.deepcopy(some_dfn['dfn'])

    for p in range(int(some_dfn['nP'])): #denormalise the velocity dimension
        for v in range(int(some_dfn['nV'])):
            dfn[p,v,:,:,:]*=some_dfn['v'][v]**2

    for r in range(int(some_dfn['nR'])): #convert from markers per s^3/m^6 to markers per cell
        if some_dfn['nP']>1:
            dfn[:,:,:,:,r]*=2.0*pi*some_dfn['R'][r]*(some_dfn['R'][1]-some_dfn['R'][0])*(some_dfn['Z'][1]-some_dfn['Z'][0])*(some_dfn['V'][1]-some_dfn['V'][0])*(some_dfn['L'][1]-some_dfn['L'][0])*(some_dfn['P'][1]-some_dfn['P'][0])
        else:
            dfn[:,:,:,:,r]*=2.0*pi*some_dfn['R'][r]*(some_dfn['R'][1]-some_dfn['R'][0])*(some_dfn['Z'][1]-some_dfn['Z'][0])*(some_dfn['V'][1]-some_dfn['V'][0])*(some_dfn['L'][1]-some_dfn['L'][0])*2.*pi

    return dfn

def dfn_collapse(some_dfn,dimensions=['R','Z']):
    """
    integrate and collapse the dfn to the supplied to dimensions

    notes:
        returns the Dfn array
    """

    dfn=copy.deepcopy(some_dfn['dfn'])

    dimensions_indices=[] #collapse dfn along specified axes
    if 'Z' not in dimensions:
        dimensions_indices.extend([4])
    if 'R' not in dimensions:
        dimensions_indices.extend([3])
    if 'L' not in dimensions:
        dimensions_indices.extend([2])
    if 'V' not in dimensions:
        dimensions_indices.extend([1])
    if 'P' not in dimensions:
        dimensions_indices.extend([0])

    for dimension in dimensions_indices: #axis denotes which dimension will be collapsed, so go in descending to get array shape correct
        dfn=np.sum(dfn,axis=dimension)
    
    return dfn

#################################

##################################################################

###################################################################################################