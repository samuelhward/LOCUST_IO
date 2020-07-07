#plot_coulomb_logarithm.py

'''
Samuel Ward
13/05/2019
----
plots the coulomb logarithm as LOCUST sees it
---
notes: 
---
'''


##################################################################
#Preamble

import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway

if __name__=='__main__':

    try:
        import context
    except:
        raise ImportError("ERROR: context.py could not be imported!\nreturning\n")
        sys.exit(1)

try:
    import numpy as np
    import matplotlib.pyplot as plt
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

try:
    import run_scripts.utils
except:
    raise ImportError("ERROR: LOCUST_IO/src/run_scripts/utils.py could not be imported!\nreturning\n")
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


##################################################################

def plot_coulomb_logarithm(At,Ai,Zt,Zi,Ti,ni,Einj,Pdep,Bmod,code='LOCUST',colmap=settings.cmap_default,colmap_val=np.random.uniform(),line_style=settings.plot_line_style,ax=False,fig=False):
    """
    plot collision operator drift coefficients for collisions of test particle against arbitrary background species

    notes:
        adapted from coll IDL routine by Rob Akers
        a lot of variable names are mapped to those in original IDL script
        XXX for now always uses LOCUST coulomb logarithm definition
        returns x,y of final plot in [J] [J/s]
    args:
        At - test particle mass [amu]
        Ai - array of background ion masses [amu]
        Zt - test particle charge [e]
        Zi - array of background ion charges [e]
        Ti - array of background ion temperatures [eV]
        ni - array of background ion densities [m**-3]
        Einj - injection energy [eV]
        Pdep - injection power [W]
        Bmod - magnetic field strength [T]
        code - use coulomb logarithm from designated code (options=LOCUST,TRANSP,ASCOT)
        colmap - plotted line colour
        colmap_val - optional numerical value for defining single colour plots 
        line_style - set 1D line style
        ax - ax object to add plot to
        fig - take input fig (can be used to add colourbars etc)
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

    echg   = constants.charge_e        
    eps0   = constants.epsilon_0         
    Mi=Ai*constants.amu #[kg]
    Mt=At*constants.amu #[kg] 
    Ti*=echg #[J]
    n_i=ni
    E= np.linspace(1.5*Ti/echg,Einj,100)*echg #energy range to look at in [J]- 2D since energy range for each species
    V=np.sqrt(2.0*E/Mt) #test particle velocity [m/s]
    Vi=np.sqrt(2.0*Ti/Mi) #background ion thermal velocity [m/s]
    Xi    = V/Vi #currently nE by nion
    Xi=Xi.T
    gam_t = echg**4/(Mt**2*4.0*constants.pi*eps0**2)

    lnL=[]
    for E_ in E: #calculate coulomb logarithms for each species at each energy (lnLe~17 lnL~21)
        lnL.append(run_scripts.utils.calc_coulomb_logarithm(Et=E_/echg,n=n_i,T=Ti/echg,Ai=Ai,At=At,Z=Zi,Zt=Zt,Bmod=Bmod,code=code))
    lnL=np.array(lnL,ndmin=2).T #now have 2D array, one axis for energy and one axis for species
    E=E.T[0]

    for counter,lnL_ in enumerate(lnL):
        ax.plot(E/(1000.*echg),lnL_,color=colmap(colmap_val),linestyle=line_style)#(counter/len(lnL))) #cycle through colours
        ax.set_ylabel('Coulomb Logarithm')
        ax.set_xlabel('Energy [KeV]')
    ax.legend(tuple(['mass = {mass} [amu] temperature = {temp} [keV]'.format(mass=str(mass),temp=str((temp)/echg/1000)) for mass,temp in zip(Ai,Ti)]))

    if ax_flag is False and fig_flag is False:
        plt.show() 

    return E,lnL #return x,y of plots in [J] [J/s]


if __name__=='__main__': #plot collision operator and expected steady state distribution function (diffusion disabled)

    import argparse
    parser=argparse.ArgumentParser(description='plot coulomb logarithm')
    parser.add_argument('--Ti',nargs='+',type=float,action='store',default=[9419.,4119.],dest='Ti',help="temperature of each background species [eV] e.g. --Ti  9419 4119",required=False)
    parser.add_argument('--ni',nargs='+',type=float,action='store',default=[5.96e19,5.96e19],dest='ni',help="density of each background species [#/m^3] e.g. --ni  9419 4119",required=False)
    parser.add_argument('--Pdep',type=float,action='store',default=1.,dest='Pdep',help="power deposited to plasma [W]",required=False)
    parser.add_argument('--Einj',type=float,action='store',default=80000.,dest='Einj',help="test particle injection energy [eV]",required=False)
    parser.add_argument('--Bmod',type=float,action='store',default=1.5,dest='Bmod',help="|B| [T]",required=False)
    parser.add_argument('--Zt',type=int,action='store',default=1,dest='Zt',help="test particle charge [e]",required=False)
    parser.add_argument('--Zi',nargs='+',type=int,action='store',default=[1,1],dest='Zi',help="background species charges [e] e.g. --Zi 1 1",required=False)
    parser.add_argument('--At',type=float,action='store',default=constants.mass_deuteron_amu,dest='At',help="test particle mass [amu]",required=False)
    parser.add_argument('--Ai',nargs='+',type=float,action='store',default=[constants.mass_deuteron_amu,constants.mass_electron_amu],dest='Ai',help="background species masses [amu] e.g. --Ai 2.0141017781 0.00054858",required=False)
    parser.add_argument('--code',type=str,action='store',default='LOCUST',dest='code',help="specify coulomb logarithm used in this code (options=LOCUST,TRANSP,ASCOT)",required=False)

    args=parser.parse_args()

    fig,ax=plt.subplots(1)
    E,lnL=plot_coulomb_logarithm(At=args.At,
                                Ai=np.array(args.Ai),
                                Zt=args.Zt,
                                Zi=np.array(args.Zi),
                                Ti=np.array(args.Ti),
                                ni=np.array(args.ni),
                                Einj=args.Einj,
                                Pdep=args.Pdep,
                                Bmod=args.Bmod,
                                code=args.code,
                                fig=fig,
                                ax=ax,
                                colmap=settings.cmap_default)
    ax.set_title('coulomb logarithm for test particle mass={} amu'.format(args.At))
    plt.show()

#################################

##################################################################

###################################################################################################