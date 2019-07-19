#plot_collision_operator.py

'''
Samuel Ward
13/05/2019
----
File which holds physical constants and other data
---
notes: 
---
'''


##################################################################
#Preamble

import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway

try:
    import numpy as np
    import matplotlib.pyplot as plt
    pi=np.pi
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import constants
except:
    raise ImportError("ERROR: LOCUST_IO/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)

##################################################################

def plot_collision_operator(Ti,Te,Ne,Einj,Pinj,ax,fig):
    """
    plot collision operator drift coefficients

    notes:
        adapted from coll IDL routine by Rob Akers
    args:
        Ti - ion temperature [eV]
        Te - electron temperature [eV]
        Ne - electron density [m**-3]
        Einj - injection energy [eV]
        Pinj - injection power [W]
    """

    try:
        from scipy.special import erf as erf
    except:
        raise ImportError("ERROR: LOCUST_IO/plot_scripts/plot_collision_operator.py could not import scipy.special.erf!\nreturning\n") 
        sys.exit(1)

    echg   = 1.602176487e-19            
    eps0   = 8.854187817e-12         
    mp     = 1.672621637e-27
    ma     = 2*mp
    mb     = 2*mp
    Ti     = Ti*echg                
    Te     = Te*echg
    n_e    = Ne
    n_i    = n_e                    
    Zi     = 1.0
    me     = 9.10938215e-31
    lnLe   = 17.0 #coulomb logs
    lnLi   = 21.0
    E      = np.linspace(1.5*Ti/echg,Einj,100)*echg  #in J  
    V      = np.sqrt(2.0*E/ma)
    Pinj  /= echg #in eV/s

    Vi     = np.sqrt(2.0*Ti/mb)
    Ve     = np.sqrt(2.0*Te/me)

    xe     = V/Ve
    xi     = V/Vi

    nui    = ((2*n_i*echg**4*lnLi/(8*pi*eps0**2*ma**2))*(mb/(2*Ti))**1.5)*2*((erf(xi) - 2*xi*np.exp(-xi**2)/np.sqrt(pi))*(ma/mb) - 2*xi*np.exp(-xi**2)/np.sqrt(pi))/xi**3

    nue    = ((2*n_e*echg**4*lnLe/(8*pi*eps0**2*ma**2))*(me/(2*Te))**1.5)*2*((erf(xe) - 2*xe*np.exp(-xe**2)/np.sqrt(pi))*(ma/me) - 2*xe*np.exp(-xe**2)/np.sqrt(pi))/xe**3


    ax.plot(E/(1000.*echg),-nue*E/1.0e-13,color='m')
    ax.plot(E/(1000.*echg),-nui*E/1.0e-13,color='g')
    ax.plot(E/(1000.*echg),-nui*E/1.0e-13-nue*E/1.0e-13,color='b')

    ax.set_xlabel('Energy [KeV]')
    ax.set_ylabel('Energy drift [10e-13 J/s]')
    plt.show()

    fig,ax=plt.subplots(1)
    ax.plot(E/(1000.*echg),-echg*Pinj/Einj/(-nui*E-nue*E),color='b')
    ax.set_xlabel('Energy [KeV]')
    ax.set_ylabel('# [eV^-1]')
    plt.show()



fig,ax=plt.subplots(1)
plot_collision_operator(Ti=1.0,Te=4119.0,Ne=5.96e19,Einj=80000,Pinj=1,fig=fig,ax=ax)

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import physical_constants as const
from scipy.special import erf

# Test particle parameters and energy grid
ma = const["deuteron mass"][0]
qa = const["elementary charge"][0]
Egrid = np.linspace(1, 120e3, 100)

# Plasma species parameters
Nion = 1
anum = np.array([2])
znum = np.array([1])
ne = 5.9567894e19
ni = np.array([5.9567894e19])
Te = 4.1185698e3
Ti = 9.4194805e3
va    = np.sqrt( const["elementary charge"][0]*2*Egrid/ma )

ax.set_xlabel("Energy [keV]")
ax.set_ylabel("Energy drift [J/s]")

ax.set_xlim(0, 120)

# Bacgkround quantities, the last item is the Coulomb logarithm
eback = [Te, -const["elementary charge"][0], const["electron mass"][0],
         17]
iback = [Ti,  const["elementary charge"][0], const["proton mass"][0],
         21]
for back in [eback, iback]:
    mb = back[2];
    qb = back[1];

    clog = back[3];
    cab  = qa*qa*qb*qb * clog \
           /( 8*np.pi* np.power(const["electric constant"][0],2) );

    vb = np.sqrt(2*back[0]*const["elementary charge"][0]/mb);
    x  = va/vb;

    h = 2*erf(x)/x;
    g = ( ( 4*np.exp(-x*x) / ( np.sqrt(np.pi)*x ) - 2*erf(x) / (x*x) ) ) / x;

    nu = (ne*cab/(ma*ma)) * (1 + ma/mb) * (1/np.power(vb,3))*g;
    D  = (ne*cab/(ma*ma)) * (1/vb) * h;

    ve = 2*nu*(0.5*va*va*ma) + ma*D;


    ax.plot(Egrid/1e3, ve/1e-13, marker='o', fillstyle='none');


plt.show()

#################################

##################################################################

###################################################################################################
