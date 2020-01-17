#plot_tomographic.py

'''
Samuel Ward
17/01/2019
----
plots tomographic slices of combined quantities
---
usage:
    see README.md for usage
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

def plot_tomographic(colmaps=settings.cmap_default,colmap_val=np.random.uniform(),ax=False,fig=False,*args):
    """

    args:
    notes:
    """


    def update(time=0,harmonic_fundamental=2,omega=10,harmonic_additional=None):
        """
        
        args:
        notes:
        """

        ax1.cla() #clear axes
        ax2.cla() 
        for coil in range(number_coils): #update coil currents
            total_phase_covered=harmonic_fundamental*2.*np.pi
            phase_coil=total_phase_covered*(coil)/(number_coils)    
            phase_per_coil=2.*np.pi/number_coils
            first_point_to_plot=int(coil*(plot_points_per_coil))
            last_point_to_plot=int(first_point_to_plot+plot_points_per_coil)
            coil_current[first_point_to_plot:last_point_to_plot]=np.cos(phase_coil-omega*time+phase_per_coil)

        ax1.set_ylim(-3,3)
        ax1.plot(np.linspace(0.,np.pi*2.,len(coil_current)),coil_current,'k-')
        wave_fundamental=np.cos(harmonic_fundamental*np.linspace(0.,np.pi*2.,len(coil_current))-omega*time)
        ax1.plot(np.linspace(0.,np.pi*2.,len(coil_current)),wave_fundamental,'r-')
        ax2.plot(np.fft.fft(coil_current).real[0:15])

        if harmonic_additional:
            for harmonic in harmonic_additional:
                wave_harmonic=np.cos(harmonic*np.linspace(0.,np.pi*2.,len(coil_current))+omega*time+np.pi)
                ax1.plot(np.linspace(0.,np.pi*2.,len(coil_current)),wave_harmonic,'b-')

        if True:
            ax1.plot(np.linspace(0.,np.pi*2.,len(coil_current)),wave_harmonic+wave_fundamental,'g-')


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
        polar=True if axes==['X','Y'] else False
        ax = fig.add_subplot(111,polar=polar)


    fig,(ax1,ax2)=plt.subplots(2,1)
    time=np.linspace(0.,2.*np.pi/omega)
    harmonic_additional=[6]
    animation=FuncAnimation(fig,update,frames=time,fargs=[harmonic_fundamental,omega,harmonic_additional],repeat=True) #cycle through phases and make animation
    #plt.show()
    animation.save('RMP_rotation.gif',writer='pillow')

  
    
if __name__=='__main__': #plot collision operator and expected steady state distribution function (diffusion disabled)

    import argparse
    parser=argparse.ArgumentParser(description='plot drag coefficients and diffusionless steady-state distribution functions for given collision operator')
    parser.add_argument('--Ti',nargs='+',type=float,action='store',default=[9419.,4119.],dest='Ti',help="temperature of each background species [eV] e.g. --Ti  9419 4119",required=False)
    parser.add_argument('--ni',nargs='+',type=float,action='store',default=[5.96e19,5.96e19],dest='ni',help="density of each background species [#/m^3] e.g. --ni  1.6e19 3.3e19",required=False)
    parser.add_argument('--Pdep',type=float,action='store',default=1.,dest='Pdep',help="power deposited to plasma [W]",required=False)
    parser.add_argument('--Einj',type=float,action='store',default=80000.,dest='Einj',help="test particle injection energy [eV]",required=False)
    parser.add_argument('--type',type=str,action='store',default='LOCUST',dest='type',help="collision operator type [NRL/LOCUST]",required=False)
    parser.add_argument('--Bmod',type=float,action='store',default=1.5,dest='Bmod',help="|B| [T]",required=False)
    parser.add_argument('--Zt',type=int,action='store',default=1,dest='Zt',help="test particle charge [e]",required=False)
    parser.add_argument('--Zi',nargs='+',type=int,action='store',default=[1,1],dest='Zi',help="background species charges [e] e.g. --Zi 1 1",required=False)
    parser.add_argument('--At',type=float,action='store',default=constants.mass_deuteron_amu,dest='At',help="test particle mass [amu]",required=False)
    parser.add_argument('--Ai',nargs='+',type=float,action='store',default=[constants.mass_deuteron_amu,constants.mass_electron_amu],dest='Ai',help="background species masses [amu] e.g. --Ai 2.0141017781 0.00054858",required=False)

    args=parser.parse_args()

    fig,ax=plt.subplots(1)
    E,drag=plot_collision_operator(At=args.At,Ai=np.array(args.Ai),Zt=args.Zt,Zi=np.array(args.Zi),Ti=np.array(args.Ti),ni=np.array(args.ni),Einj=args.Einj,Pdep=args.Pdep,Bmod=args.Bmod,fig=fig,ax=ax,type=args.type,colmap=settings.cmap_default)
    ax.set_title('individual drag coefficients for test particle mass={} amu'.format(args.At))
    plt.show()

    drag_total=np.sum(drag,axis=0) #sum electron and ion drags and plot
    fig,ax=plt.subplots(1)
    ax.plot(E/(1000.*constants.charge_e),(drag_total/1.0e-13),color=settings.cmap_default(np.random.uniform()))
    ax.set_title('total drag for test particle mass={} amu'.format(args.At))
    ax.set_xlabel('Energy [KeV]')
    ax.set_ylabel('Energy drift [10e-13 J/s]')
    plt.show()
    
    #plot expected steady state distribution function (assuming no diffusion)
    fig,ax=plt.subplots(1)
    ax.plot(E/(1000.*constants.charge_e),-args.Pdep/args.Einj/(drag_total),color=settings.cmap_default(np.random.uniform()))
    ax.set_title('expected steady-state diffusionless distribution function')
    ax.set_xlabel('Energy [KeV]')
    ax.set_ylabel('# [eV**-1]')
    plt.show()

#################################

##################################################################

###################################################################################################



fig,(ax1,ax2)=plt.subplots(2,1)
time=np.linspace(0.,2.*np.pi/omega)
harmonic_additional=[6]
animation=FuncAnimation(fig,update,frames=time,fargs=[harmonic_fundamental,omega,harmonic_additional],repeat=True) #cycle through phases and make animation
#plt.show()
animation.save('RMP_rotation.gif',writer='pillow')

#################################

##################################################################

###################################################################################################
