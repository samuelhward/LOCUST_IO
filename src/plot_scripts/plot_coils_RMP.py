#plot_coils_RMP.py

'''
Samuel Ward
01/05/2020
----
plots RMP harmonics for given field according to MARS-F conventions
---
usage:
    see README.md for usage
notes:         
---
'''

###################################################################################################
#Preamble
 
import sys #have global imports --> makes less modular (no "from input_classes import x") but best practice to import whole input_classes module anyway

try:
    import numpy as np
    import scipy
    import cmath
    import matplotlib.pyplot as plt
    import matplotlib.patches
    from matplotlib.animation import FuncAnimation
except:
    raise ImportError("ERROR: initial modules could not be imported!\nreturning\n")
    sys.exit(1)

if __name__=='__main__':

    try:
        import context
    except:
        raise ImportError("ERROR: context.py could not be imported!\nreturning\n")
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

###################################################################################################
#Main

def plot_coils_RMP(phase_shift,n_0,n_range,coil_current,coil_rows=[1,2,3],tokamak='ITER',waveform='cos'):
    """
    plot generic RMP coils as well as associated theoretical harmonics

    args:
        phase_shift - phase shift of applied current waveform [degrees]
        n_0 - fundamental toroidal mode number [-]
        n_range - range of toroidal harmonics to include [-]
        coil_current - coil current [kAt] 
        coil_rows - list of coil rows to plot (1=uppermost) [int]
        tokamak - tokamak name [-]
        waveform - waveform type [str]
    notes:
    """


    for ax in axes: ax.cla() #clear axes

    rad_to_deg=360./2.*np.pi
    deg_to_rad=2.*np.pi/360.

    #define coil waveform table
    def current_waveform_cosine(n_0,phi,phase_shift):
        return coil_current*np.cos(n_0*(phi-phase_shift)*deg_to_rad) 

    waveforms={}
    waveforms['cos']=current_waveform_cosine
    waveform=waveforms[waveform]

    #define tokamak settings
    coil_data_options={} 
    coil_data_options['ITER']={} #ITER settings
    coil_data_options['ITER']['number_coils_per_row']=9
    coil_data_options['ITER']['number_rows']=3
    coil_data_options['ITER']['coil_offset']=np.array([30.,26.7,30.]) #Phi0 mapping from MARS-F origin to coil centres
    #coil_data_options['ITER']['coil_offset']=np.array([0.,0.,0]) #XXX see effect of changing Phi0 (just look at top and bottom coil rows)
    #XXXcoil_data_options['ITER']['coil_offset']=np.array([phase_shift,phase_shift,phase_shift]) #XXX scan Phi0
    coil_data_options['ITER']['coil_coverage']=np.array([29.4,20.9,30.5]) #toroidal coverage DELTAPhij between coils for each row - assuming equal coil spacing for a given row (DELTAPhij=DELTAPhi)
    coil_data_options['ITER']['phase_shift_coils']=np.array([0.,0.,0.])
    coil_data_options['ITER']['phase_shift_coils']=np.array([phase_shift,phase_shift,phase_shift]) #XXX
    #coil_data_options['ITER']['phase_shift_coils']=np.array([86.,0.,34.]) #XXX

    #derive extra information for settings dispatch
    for tokamak in coil_data_options.keys():
        coil_data_options[tokamak]['number_coils_total']=coil_data_options[tokamak]['number_coils_per_row']*coil_data_options[tokamak]['number_rows']
        coil_data_options[tokamak]['coil_spacing']=np.full(coil_data_options[tokamak]['number_rows'],360./coil_data_options[tokamak]['number_coils_per_row']) #space between coil centres deltaPhij per fow - assuming equal spacing for a given row
        coil_data_options[tokamak]['coil_locations']=np.array([phi_0+(j)*coil_spacing for coil_spacing,phi_0 in zip(coil_data_options[tokamak]['coil_spacing'],coil_data_options[tokamak]['coil_offset']) for j in np.arange(coil_data_options[tokamak]['number_coils_per_row'])]).reshape(coil_data_options[tokamak]['number_rows'],coil_data_options[tokamak]['number_coils_per_row']) #locations of coil centres phij
        coil_data_options[tokamak]['coil_index']=np.arange(coil_data_options[tokamak]['coil_locations'].size).reshape(coil_data_options[tokamak]['coil_locations'].shape)+1

    I_j=waveform(n_0=n_0,phi=coil_data_options[tokamak]['coil_locations'],phase_shift=coil_data_options[tokamak]['phase_shift_coils'][:,None])

    #manually fourier transform to n space (eq 4 and 7 in notes)
    n_axis=np.linspace(1,n_range,n_range,dtype=int)
    I_n_cos=np.zeros((coil_data_options[tokamak]['number_rows'],n_range))
    I_n_sin=np.zeros((coil_data_options[tokamak]['number_rows'],n_range))

    for coil_row in np.arange(coil_data_options[tokamak]['number_rows']):
        for n in n_axis: #calculate complex amplitude for each harmonic in each row
            I_n_cos[coil_row,n-1]=(2./(n*np.pi))*np.sin(n*coil_data_options[tokamak]['coil_coverage'][coil_row]*deg_to_rad/2.)*np.sum(I_j[coil_row]*np.cos(n*coil_data_options[tokamak]['coil_locations'][coil_row]*deg_to_rad))
            I_n_sin[coil_row,n-1]=(2./(n*np.pi))*np.sin(n*coil_data_options[tokamak]['coil_coverage'][coil_row]*deg_to_rad/2.)*np.sum(I_j[coil_row]*np.sin(n*coil_data_options[tokamak]['coil_locations'][coil_row]*deg_to_rad))

    I_n=np.array([complex(I_c,I_s) for coil_row in np.arange(coil_data_options[tokamak]['number_rows']) for I_c,I_s in zip(I_n_cos[coil_row],I_n_sin[coil_row])]).reshape(coil_data_options[tokamak]['number_rows'],n_range)

    #plotting options
    plot_fft_yueqiang=True
    plot_fft_numpy=False#True
    plot_fft_reconstruction=True
    plot_coils=True
    plot_perturbation_fundamental=True

    for coil_row,colour in zip(coil_rows,['blue','red','green']):
        #axes[0].plot(n_axis,I_n[coil_row].real,label=f'coil row = {coil_row}',color=colour,linestyle=':')
        #axes[0].plot(n_axisI_n[coil_row].imag,label=f'coil row = {coil_row}',color=colour,linestyle=':')
        
        coil_row-=1

        if plot_perturbation_fundamental:
            phi=np.linspace(0.,360.,100)
            I_j_fundamental=waveform(n_0=n_0,phi=phi,phase_shift=coil_data_options[tokamak]['phase_shift_coils'][coil_row]) #plot theoretical waveform
            axes[1].plot(phi,I_j_fundamental,color='magenta',linestyle='--',label='current profile')
        if plot_fft_yueqiang:
            axes[0].plot(n_axis,abs(I_n[coil_row]),label=f'coil row = {coil_row}',color=colour,marker='.',linestyle='-')
        if plot_fft_numpy: #calculate fft equivalent
            I_fft=np.abs(np.fft.fft(I_j)) 
            axes[0].plot(I_fft[coil_row,1:int(len(I_fft[coil_row])/2)],'k-',label=f'fft') #XXX
        if plot_fft_reconstruction:
            phi=np.linspace(0.,360.,360)
            waveform_reconstruction=np.zeros(len(phi))

            reconstruction_type='analytical'
            if reconstruction_type=='analytical':
                for n in n_axis:
                    #the following are two equivalent ways of thinking about/calculating this
                    #amplitude_this_n=abs(I_n[coil_row,n-1])*np.cos(n*phi*deg_to_rad-np.arctan2(I_n[coil_row,n-1].imag,I_n[coil_row,n-1].real))
                    amplitude_this_n=I_n.real[coil_row,n-1]*np.cos(n*phi*deg_to_rad)+I_n.imag[coil_row,n-1]*np.sin(n*phi*deg_to_rad)
                    plt.plot(phi,amplitude_this_n,color=settings.cmap_plasma(.3*n/n_axis[-1]+.4),label=f'n={n}')
                    waveform_reconstruction+=amplitude_this_n
            elif reconstruction_type=='individual':
                    pass 
            plt.plot(phi,waveform_reconstruction,color=settings.cmap_red_nice(0.),label='sum')

        #for n_highlight in n_highlights:
        #    axes[0].axvline(n_highlight,'m')
        if plot_coils:
            for coil in range(coil_data_options[tokamak]['number_coils_per_row']): #update coil currents
                rectangle=matplotlib.patches.Rectangle(
                    xy=(coil_data_options[tokamak]['coil_locations'][coil_row,coil]-coil_data_options[tokamak]['coil_coverage'][coil_row]/2.,0.0),
                    width=coil_data_options[tokamak]['coil_coverage'][coil_row],
                    height=I_j[coil_row,coil],
                    edgecolor='black',facecolor='none') #,label=f'coil row {coil_row}'
                axes[1].add_patch(rectangle)
        axes[1].scatter(coil_data_options[tokamak]['coil_locations'][coil_row],I_j[coil_row],color='m') #add on coil currents

    axes[0].set_xlabel('$n$')
    axes[0].set_xlim([0.,n_range])
    axes[0].set_ylim([-1.1*coil_current,1.1*coil_current])
    axes[0].set_title('Fourier transform amplitude')
    axes[1].set_xlabel('$\phi$')
    axes[1].legend(loc='right')
    axes[1].set_xlim([np.min(coil_data_options[tokamak]['coil_locations'])-np.max(coil_data_options[tokamak]['coil_coverage'])/2.,np.max(coil_data_options[tokamak]['coil_locations'])+np.max(coil_data_options[tokamak]['coil_coverage'])/2.])
    axes[1].set_xlim([0.,360])
    axes[1].set_ylim([-1.1*coil_current,1.1*coil_current])
    axes[1].set_title('Coil current [kAt]')

if __name__=='__main__':

    import argparse
    parser=argparse.ArgumentParser(description='plot RMP coils')

    parser.add_argument('--n_0',type=int,action='store',default=3,dest='n_0',help="fundamental toroidal mode number [-]",required=False)
    parser.add_argument('--n_range',type=int,action='store',default=20,dest='n_range',help="range of toroidal harmonics to include [-]",required=False)
    parser.add_argument('--coil_current',type=float,action='store',default=90.,dest='coil_current',help="coil current [kAt]",required=False)
    parser.add_argument('--coil_rows',nargs='+',type=int,action='store',default=[1],dest='coil_rows',help="list of coil rows to plot (1=uppermost) [int] e.g. 1 2 3",required=False)
    parser.add_argument('--tokamak',type=str,action='store',default='ITER',dest='tokamak',help="tokamak name [-]",required=False)
    parser.add_argument('--save',type=bool,action='store',default=False,dest='save',help="save gif? [-]",required=False)

    args=parser.parse_args()
    args={key:arg for key,arg in args._get_kwargs()}

    fig,axes=plt.subplots(2)
    phase_shift=np.linspace(0,360,360)
    animation=FuncAnimation(fig,plot_coils_RMP,frames=phase_shift,fargs=[args[arg] for arg in ['n_0','n_range','coil_current','coil_rows','tokamak']],repeat=True,interval=1) #cycle through phases and make animation
    plt.show()

    if args['save']: animation.save('plot_coils_RMP.gif',writer='pillow')

#################################

##################################################################

###################################################################################################