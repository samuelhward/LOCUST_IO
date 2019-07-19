#plot_collision_operator.py

'''
Samuel Ward
13/05/2019
----
plots the collision operator drag/drift terms as LOCUST sees it
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
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import constants
except:
    raise ImportError("ERROR: LOCUST_IO/constants.py could not be imported!\nreturning\n") 
    sys.exit(1)
try:
    from settings import *
except:
    raise ImportError("ERROR: LOCUST_IO/settings.py could not be imported!\nreturning\n") 
    sys.exit(1)


def calc_coulomb_logarithm(Et,n,T,A,At,Z,Zt,Bmod):
    """
    calculates the Coulomb logarithm at single test particle energy Et for arbitrary background species

    notes: 
        must pass all species information to this routine, since Rmax is sum over species
        returns lnL, an array with a Coulomb logarithm corresponding to each species
    args:
        Et - test particle energy [eV]
        n - array of densities for background species [m**-3]
        T - array of temperatures for background species [eV]
        A - array of atomic masses for background species [amu]
        At - atomic mass for test particle species [amu]
        Z - array of atomic charges for background species [int]
        Zt - atomic charge for test particle species [int]
        Bmod - magnetic field strength [T]
    """
    Et/=1000. #convert to keV
    T/=1000. #convert to KeV

    n=np.asarray(n)
    T=np.asarray(T)
    A=np.asarray(A)
    Z=np.asarray(Z)

    omega2 = 1.74*Z**2/A*n + 9.18e15*Z**2/A**2*Bmod**2
    vrel2  = 9.58e10*(T/A + 2.0*Et/At)
    rmx    = np.sqrt(1.0/np.sum(omega2/vrel2))

    vrel2  = 9.58e10*(3.0*T/A+2.0*Et/At)
    rmincl = 0.13793*np.abs(Z*Zt)*(A+At)/A/At/vrel2
    rminqu = 1.9121e-08*(A+At)/A/At/np.sqrt(vrel2)

    rmn=[]
    for rmincl_,rminqu_ in zip(rmincl,rminqu): 
        rmn.extend([np.max([rmincl_,rminqu_])])
    rmn=np.asarray(rmn)

    lnL=np.log(rmx/rmn)

    return np.full(len(n),17)

def calc_chandrasekhar_function(x,m,mt):
    """
    calculate Chandrasekhar function G(x)=(ERF(x) - x.d{ERF(x)}/dx)/2x**2

    notes:
        N.B. ERF(x) using function erf_7_1_26 breaks down at small x.
        Routine resorts to small x analytic approximation for x<0.05.
        Routine is only valid for x>=0.0
        Handbook of Mathematical Functions. Edited by Milton Abramowitz and 
        Irene A. Stegun, Dover Publications, Inc., New York, 1965.
        Error function and Fresnel Integrals, EQN. 7.1.26.
        Valid to |E(x)| <= 1.5e-7. Calculation in gpu precision
        f1 - = ERF(x) - (1+m/Mb).x.d{ERF(x)}/dx
        f2 - = ERF(x) - G(x)
    args:
        x - X
        m - background particle mass [kg]
        mt - test particle mass [kg]
    """

    a1 = 0.254829592
    a2 =-0.284496736
    a3 = 1.421413741
    a4 =-1.453152027
    a5 = 1.061405429
    p  = 0.327591100
          
    i  = np.where(x<0.0)[0]
    if len(i)>0:
        print("ERROR: calc_chandrasekhar_function invalid for x<0!\nreturning!\n")
        return

    t     = 1.0/(1.0+p*x)          
    exp_  = np.exp(-x**2)
    erf_s = 1.0 - ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*exp_
    G     = ( erf_s - x*2.0*exp_/np.sqrt(constants.pi) )/( 2.0*x**2)

    f1    = erf_s - (1.0+(m/mt))*x*2.0*exp_/np.sqrt(constants.pi)
    f2    = erf_s - G

    i     = np.where( x<0.005)[0]
    if len(i)!=0:     
        #eqn. for G breaks down at low x - revert to analytic approximation:
        x     = x[i]
        G[i]  = 2.0*x/(3.0*np.sqrt(constants.pi))
        f1[i] = (2.0/np.sqrt(constants.pi))*( ( (2.0/3.0) + (m/mt) )*x**3 - x*(m/mt) ) 
        f2[i] = (2.0/(3.0*np.sqrt(constants.pi)))*(2.0-x**2)*x

    return G,f1,f2


##################################################################

def plot_collision_operator(mass_t,mass_i,Zt,Zi,Ti,Ni,Einj,Pinj,Bmod,type='NRL',colmap=cmap_default,ax=False,fig=False):
    """
    plot collision operator drift coefficients for collisions of test particle against arbitrary background species

    notes:
        adapted from coll IDL routine by Rob Akers
        a lot of variable names are mapped to those in original IDL script
        XXX for now always uses LOCUST coulomb logarithm definition
        returns x,y of final plot in [J] [J/s]
    args:
        mass_t - test particle mass [amu]
        mass_i - array of background ion masses [amu]
        Zt - test particle charge [e]
        Zi - array of background ion charges [e]
        Ti - array of background ion temperatures [eV]
        Ni - array of background ion densities [m**-3]
        Einj - injection energy [eV]
        Pinj - injection power [W]
        Bmod - magnetic field strength [T]
        type - toggle various representations from different codes
        colmap - plotted line colour
        ax - ax object to add plot to
    """

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

    if type=='NRL': #plot the NRL formulary drag terms

        try:
            from scipy.special import erf as erf
        except:
            raise ImportError("ERROR: LOCzST_IO/plot_scripts/plot_collision_operator.py could not import scipy.special.erf!\nreturning\n") 
            sys.exit(1)

        echg   = constants.charge_e        
        eps0   = constants.epsilon_0         
        Mi     = mass_i*constants.amu #[kg]
        Mt     = mass_t*constants.amu #[kg]
        Ti    *= echg #[J]               
        n_i    = Ni                    
        E      = np.linspace(np.full(len(Ti),0.01*Einj),Einj,100)*echg #energy range to look at in [J] - 2D since energy range for each species
        V      = np.sqrt(2.0*E/Mt) #[m/s]
        Pinj  /= echg
        Vi     = np.sqrt(2.0*Ti/Mi) #[m/s]
        xi     = V/Vi

        lnL=[]
        for E_ in E: #calculate couloMi logarithms for each species at each energy (lnLe~17 lnL~21)
            lnL.append(calc_coulomb_logarithm(Et=E_/echg,n=n_i,T=Ti/echg,A=mass_i,At=mass_t,Z=Zi,Zt=Zt,Bmod=Bmod))
        lnL=np.array(lnL,ndmin=2) #now have 2D array, one axis for energy and one axis for species
        E=E.T[0] #collapse down to 1D array now
        nui=((2*n_i*echg**4*lnL/(8*constants.pi*eps0**2*Mt**2))*(Mi/(2*Ti))**1.5)*2*((erf(xi) - 2*xi*np.exp(-xi**2)/np.sqrt(constants.pi))*(Mt/Mi) - 2*xi*np.exp(-xi**2)/np.sqrt(constants.pi))/xi**3
        dE_dt=nui.T*E

    elif type=='LOCUST': #plot the LOCUST drag terms

        echg   = constants.charge_e        
        eps0   = constants.epsilon_0         
        Mi=mass_i*constants.amu #[kg]
        Mt=mass_t*constants.amu #[kg] 
        Ti*=echg #[J]
        n_i=Ni
        E= np.linspace(np.full(len(Ti),0.01*Einj),Einj,100)*echg #energy range to look at in [J]- 2D since energy range for each species
        V=np.sqrt(2.0*E/Mt) #test particle velocity [m/s]
        Vi=np.sqrt(2.0*Ti/Mi) #background ion thermal velocity [m/s]
        Xi    = V/Vi #currently nE by nion
        Xi=Xi.T
        gam_t = echg**4/(Mt**2*4.0*constants.pi*eps0**2)

        lnL=[]
        for E_ in E: #calculate coulomb logarithms for each species at each energy (lnLe~17 lnL~21)
            lnL.append(calc_coulomb_logarithm(Et=E_/echg,n=n_i,T=Ti/echg,A=mass_i,At=mass_t,Z=Zi,Zt=Zt,Bmod=Bmod))
        lnL=np.array(lnL,ndmin=2).T #now have 2D array, one axis for energy and one axis for species
        lnL=np.full(len(E)*len(n_i),17).reshape(len(E),len(n_i))

        dE_dt=[]
        V=V.T[0] #collapse down to 1D array now
        E=E.T[0]

        for Xi_,Mi_,Zi_,n_i_,lnL_ in zip(Xi,Mi,Zi,n_i,lnL):
            G_i,f1i,f2i=calc_chandrasekhar_function(Xi_,Mi_,Mt)
            H1i   =  -n_i_*Zi_**2*(Mt/Mi_)*f1i*lnL_*gam_t/V
            H2i   = 2*n_i_*Zi_**2*G_i*lnL_*gam_t/V
            H3i   = 2*n_i_*Zi_**2*f2i*lnL_*gam_t
            dE_dt.append(Mt*H1i)
        dE_dt=-1.*np.array(dE_dt,ndmin=2)
        
    for counter,drag in enumerate(dE_dt):
        ax.plot(E/(1000.*echg),(-drag/1.0e-13),color=colmap)#(counter/len(dE_dt))) #cycle through colours
        ax.set_xlabel('Energy [KeV]')
        ax.set_ylabel('Energy drift [10e-13 J/s]')
    ax.legend(tuple(['mass = {} amu'.format(str(mass)) for mass in mass_i]))

    if ax_flag is False and fig_flag is False:
        plt.show() 

    return E,-dE_dt #return x,y of plots in [J] [J/s]

  





#ASCOT equivalent

fig,ax=plt.subplots(1)

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
iback = [Ti,  const["elementary charge"][0], const["deuteron mass"][0],
         17]
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

#a,b=plot_collision_operator(mass_t=constants.mass_deuterium_amu,mass_i=np.array([constants.mass_neutron_amu,0.00054858]),Zt=1,Zi=np.array([1,1]),Ti=np.array([9419.,4119.]),Ni=np.array([5.96e19,5.96e19]),Einj=120000,Pinj=1,Bmod=1.5,fig=fig,ax=ax,type='LOCUST',colmap='g')
a,b=plot_collision_operator(mass_t=2,mass_i=np.array([2,0.00054858]),Zt=1,Zi=np.array([1,1]),Ti=np.array([9.4194805e3,4.1185698e3]),Ni=np.array([5.9567894e19,5.9567894e19]),Einj=120000,Pinj=1,Bmod=1.5,fig=fig,ax=ax,type='LOCUST',colmap='b')
ax.set_ylim([-3,0])
plt.show()
