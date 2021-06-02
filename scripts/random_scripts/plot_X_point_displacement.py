'''
Python adaptation of Yueqiang's script for plotting X-point displacement from XPLASMA_SURFACE* files

also generates points along contour of fixed X-point displacement

expects input_files/X_point_displacement_data/* directory from Yueqiang
'''


def plot_X_point_displacement(case=5,n=3,LVV=90,fig=None,ax=None):
    """[summary]

    Args:
        case (int, optional): [description]. Defaults to 5.
        n (int, optional): [description]. Defaults to 3.
        LVV (int, optional): [description]. Defaults to 90.
        fig ([type], optional): [description]. Defaults to None.
        ax ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]
    """
    
    import context
    import numpy as np 
    import copy,support

    if ax is None:
        ax_flag=False #need to make extra ax_flag since ax state is overwritten before checking later
    else:
        ax_flag=True

    if fig is None:
        fig_flag=False
    else:
        fig_flag=True

    if fig_flag is False:
        fig = plt.figure() #if user has not externally supplied figure, generate

    KAPPROACH = 2
    PU10 = 0
    PL10=0
    ktorq = 0
    LSS = 'b-'
    ncase = case
    n1  = n

    if ncase == 1:
        SDIR_RMP     = 'case1_5MA_1d8T/n'+str(n1)+'_Pr03_tfte2/'
        B0EXP = 1.8               

    if ncase == 2:
        SDIR_RMP     = 'case2_7d5MA_2d65T/n'+str(n1)+'_Pr03_tfte2/' 
        B0EXP = 2.65               

    if ncase == 3:
        SDIR_RMP     = 'case3_7d5MA_5d3T/n'+str(n1)+'_Pr03_tfte2/'
        B0EXP = 5.3               

    if ncase == 6:
        SDIR_RMP     = 'case6_15MA_5d3T_Q5/n'+str(n1)+'_Pr03_tfte2/' 
        B0EXP = 5.3               

    if ncase == 5:
        SDIR_RMP     = 'case5_15MA_5d3T_Q10/n'+str(n1)+'_Pr03_tfte2/' 
        B0EXP = 5.3               

    if ncase == 8:
        SDIR_RMP     = 'case8_12d5MA_5d3T/n'+str(n1)+'_Pr03_tfte2/' 
        B0EXP = 5.3               

    kcasen = 'n='+str(n1)

    SLW     = 3 
    NR      = 101
    NZ      = 201
    Rplot   = 6.2   
    Zplot   = 0.0

    R0EXP = 6.2
    mu0   = 4e-7*np.pi
    Imid  = 100
    Ixpt  = 44
    sedge = 0.995
    kplot = 2

    if KAPPROACH==1 or KAPPROACH==2 or KAPPROACH==4: 
        KDIMENSION = 2
    if KAPPROACH==3:
        KDIMENSION = 1

    if KDIMENSION==2:
        DVV_WVa= np.zeros((1,4))  
    if KDIMENSION==1:
        DVV_WVa= np.zeros((1,3)) 

    DVV_CI = DVV_WVa # mid-plane plasma surface displacement 
    DVV_CX = DVV_WVa # X-point plasma surface displacement
    DVV_CR = DVV_WVa # ration of DVV_CX / DVV_CI

    DVV_JXB = DVV_WVa # mid-plane plasma surface displacement 
    DVV_REY = DVV_WVa # X-point plasma surface displacement
    DVV_NTV = DVV_WVa # ration of DVV_CX / DVV_CI
    DVV_TOT = DVV_WVa # ration of DVV_CX / DVV_CI

    # read field data 
    def load(filename):
        data=[]
        filename=support.dir_input_files / 'X_point_displacement_data' / filename
        with open(filename,'r') as file:
            lines=file.readlines()
            for line in lines:
                split_line=line.split()
                data.extend([float(number) for number in split_line])
        data=np.asarray(data)
        data=data.reshape(int(len(data)/6),6)
        return data
        
    if ktorq==0:
        SDIRE  = '_L.OUT'
        dataS = load(SDIR_RMP+'XPLASMA_SURFACE_FIO_n'+str(n1)+SDIRE)
        RS = dataS[:,1]
        ZS = dataS[:,3]
        XS_L = dataS[:,4]+1j*dataS[:,5]
        SDIRE  = '_M.OUT'
        dataS = load(SDIR_RMP+'XPLASMA_SURFACE_FIO_n'+str(n1)+SDIRE)
        XS_M = dataS[:,4]+1j*dataS[:,5]
        SDIRE  = '_U.OUT'
        dataS = load(SDIR_RMP+'XPLASMA_SURFACE_FIO_n'+str(n1)+SDIRE)
        XS_U = dataS[:,4]+1j*dataS[:,5]

    for KV in range(1):
        AU0 =LVV 
        AM0 = LVV 
        AL0 = LVV
        AU  = AU0 
        AM  = AM0 
        AL  = AL0

        PU = PU10*np.pi/180.
        PM = 0*np.pi/180.
        PL = PL10*np.pi/180.

    # now perform scan of ELM coil currents
    # assuming various criteria for the optimization
    if KAPPROACH==1: 
        Ica1 = np.linspace(0,90,41)
        Ica2 = np.linspace(0,90,41)
    if KAPPROACH==2:
        Ica1 = np.linspace(0,360,73)
        Ica2 = np.linspace(0,360,73)
    if KAPPROACH==3:
        Ica1 = np.linspace(0,360,361)
        Ica2 = 1
    if KAPPROACH==4: 
        Ica1 = np.linspace(0,1,21)
        Ica2 = np.linspace(0,360,201)

    #dispatch table for translations to data (to make it easier to interpolate) 
    translation_options={}
    translation_options[2]={}
    translation_options[3]={}
    translation_options[5]={}
    translation_options[2][3]=[0,0]
    translation_options[2][4]=[80,330]
    translation_options[3][3]=[150,200]
    translation_options[3][4]=[0,0]
    translation_options[5][3]=[0,0]
    translation_options[5][4]=[100,260]
    dx,dy=translation_options[ncase][n1]

    # note Yueqiang x-axis is upper coil row, so flip here since I prefer upper row being y-axis
    Ica2[Ica2>360-dx]-=360
    Ica1[Ica1>360-dy]-=360
    Ica2=np.sort(Ica2)
    Ica1=np.sort(Ica1)

    if KDIMENSION==2:
      Icc1,Icc2 = np.meshgrid(Ica1,Ica2)
    if KDIMENSION==1:
      Icc1 = Ica1

    Vnratio= copy.deepcopy(Icc1)
    VnMid  = copy.deepcopy(Icc1)
    VnXpt  = copy.deepcopy(Icc1)

    if KDIMENSION==1:
        Nk1 = 1
    if KDIMENSION==2:
        Nk1 = len(Ica2)

    for k1 in range(len(Ica1)):
        # if KDIMENSION>1:  print('k1='+str(k1))
        for k2 in range(Nk1):

            if KAPPROACH==1:
                Ic1 =Icc1[k1,k2]*np.exp(1j*PU)
                Ic2 =Icc2[k1,k2]*np.exp(1j*PM)
                Ic3 =Icc1[k1,k2]*np.exp(1j*PL)
            if KAPPROACH==2:
                Ic1 =AU*np.exp(1j*Icc1[k1,k2]*np.pi/180)
                Ic2 =AM*np.exp(1j*Icc1[k1,k2]*np.pi/180*0)
                Ic3 =AL*np.exp(1j*Icc2[k1,k2]*np.pi/180)
            if KAPPROACH==3:
                Ic1 =AU*np.exp(1j*Icc1[k1]*np.pi/180)
                Ic2 =AM*np.exp(1j*PM*np.pi/180)
                Ic3 =AL*np.exp(1j*PL*np.pi/180)
            if KAPPROACH==4:
                Ic1 =(1-Icc1[k1,k2])*AU*np.exp(1j*PU)
                Ic2 =(1-Icc1[k1,k2])*AM*np.exp(1j*PM)
                Ic3 =(1-Icc1[k1,k2])*AL*np.exp(1j*PL)
                Ic12 =Icc1[k1,k2]*AU2*np.exp(1j*PU2)*np.exp(1j*PM12*np.pi/180)
                Ic22 =Icc1[k1,k2]*AM2*np.exp(1j*PM2)*np.exp(1j*PM12*np.pi/180)
                Ic32 =Icc1[k1,k2]*AL2*np.exp(1j*PL2)*np.exp(1j*PM12*np.pi/180)


            if ktorq == 0:
                
                Vn     = XS_U*Ic1 + XS_M*Ic2 + XS_L*Ic3

                if KAPPROACH==4:
                    Vn2     = XS_U2*Ic12 + XS_M2*Ic22 + XS_L2*Ic32

                if KDIMENSION==2 and KAPPROACH != 4:
                    VnMid[k1,k2]  = np.abs(Vn[Imid])
                    VnXpt[k1,k2]  = np.abs(Vn[Ixpt])
                    Vnratio[k1,k2]= VnXpt[k1,k2]/VnMid[k1,k2]
                
                if KDIMENSION==2 and KAPPROACH == 4:
                    #   Vnn = real(exp(1j*n1*Icc2[k1,k2]*pi/180)*Vn + exp(1j*n2*Icc2[k1,k2]*pi/180)*Vn2)
                    Vnn = np.abs(np.exp(1j*n1*Icc2[k1,k2]*np.pi/180)*Vn + np.exp(1j*n2*Icc2[k1,k2]*np.pi/180)*Vn2)
                    VnMid[k1,k2]  = Vnn[Imid]
                    VnXpt[k1,k2]  = Vnn[Ixpt]
                    Vnratio[k1,k2]= VnXpt[k1,k2]/VnMid[k1,k2]
                

                if KDIMENSION==1:
                    VnMid[k1]  = np.abs(Vn[Imid])
                    VnXpt[k1]  = np.abs(Vn[Ixpt])
                    Vnratio[k1]= VnXpt[k1]/VnMid[k1]


    # now for plotting and calculating X-point displacement contour

    import matplotlib.pyplot as plt 
    import processing.utils 

    kgrid = KAPPROACH  
    if kgrid==1:
        SYLAB = '$|I^{\mathrm{U}}|=|I^{\mathrm{L}}|$ [kAt]|' 
        SXLAB = '$|I^{\mathrm{M}}|$ [kAt]'
    elif kgrid==2:
        SYLAB = '$\Phi_{\mathrm{U}}-\Phi_{\mathrm{M}}$ [deg]'
        SXLAB = '$\Phi_{\mathrm{L}}-\Phi_{\mathrm{M}}$ [deg]'
    elif kgrid==3:
        SYLAB = '$\Phi_{\mathrm{U}} - \Phi_{\mathrm{L}}$ [deg]'
    elif kgrid==4:
        SYLAB = '$\alpha_{2}$'
        SXLAB = '$\Phi$ [deg]'

    levels=np.linspace(np.min(VnXpt*1e+3),np.max(VnXpt*1e+3),7)
    mesh=ax.contourf(Icc2.T,Icc1.T,VnXpt.T*1e+3,levels=levels,cmap='Greys')
    cbar=fig.colorbar(mesh,ax=ax,orientation='vertical')
    # define contour levels we want to calculate points at
    levels=[
            2*np.max(VnXpt)*1e+3/3,
            ]
    levels_coords=[]
    mesh=ax.contour(Icc2.T,Icc1.T,VnXpt.T*1e+3,levels=levels,colors='m')
    cbar.add_lines(mesh)
    #total number of points to interpolate onto
    #this will be spread out over the contour's sub-paths (if any) in proportion to their length
    number_coords_new=8
    if case == 3 and n == 4: number_coords_new+=2 # strange topology/rounding errors means need to step in for this case

    for contour in mesh.collections: 
        contour_x=[]
        contour_y=[]
        contour_coords_x=[]
        contour_coords_y=[]
        path_lengths=[]
        # initial cycle through sub-paths to get relative lengths
        for path in contour.get_paths(): 
            v = path.vertices
            path_x=v[:,0]
            path_y=v[:,1]        
            coords=np.array([path_x,path_y])
            distance_total=np.cumsum(np.sqrt(np.sum(np.diff(coords,axis=1)**2,axis=0)))[-1]
            path_lengths.append(distance_total)
        for path_number,path in enumerate(contour.get_paths()): 
            v = path.vertices
            contour_x.extend(v[:,0])
            contour_y.extend(v[:,1])
            path_x=v[:,0]
            path_y=v[:,1]        
            #interpolate along trajectory
            coords=np.array([path_x,path_y])
            distance=np.cumsum(np.sqrt(np.sum(np.diff(coords,axis=1)**2,axis=0)))
            distance=np.insert(distance,0,0.) #first distance is 0
            distance=(distance-distance[0])/(distance[-1]-distance[0]) #normalise
            interpolator=processing.utils.interpolate_1D(distance,coords,function='cubic',type='interp1d',smooth=0)
            #allocate number of points based on relative length of sub-path
            number_coords_this_path=int(number_coords_new*path_lengths[path_number]/np.sum(path_lengths))
            distance_coords_new=np.linspace(0,1,number_coords_this_path+1)[:-1] #this loops back onto itself, so add  extra point and remove the end
            coords_new=interpolator(distance_coords_new)
            contour_coords_x.extend(coords_new[0,:])
            contour_coords_y.extend(coords_new[1,:])
            ax.scatter(*coords_new,color='y')
        levels_coords.append([contour_coords_x,contour_coords_y])
        #ax.scatter(contour_x,contour_y,color='r')

    ax.set_xlabel(SXLAB)
    ax.set_ylabel(SYLAB)

    levels_coords=np.array(levels_coords,ndmin=2).T
    for counter,(x,y) in enumerate(levels_coords): ax.text(x * (1 + 0.015), y * (1 + 0.015) , counter, fontsize=12)

    if case == 3 and n == 4: 
        rotation=45. 
        x=0.2
        y=0.05
    else:
        rotation=0. 
        x=0.03
        y=0.92

    ax.text(x=x,y=y,s=f'case = {case}, n={n}',horizontalalignment='left',rotation=rotation,
    fontsize=10,transform=ax.transAxes)#,color=settings.colour_custom(rgba=[100,100,100,1])(0.))

    if ax_flag is True or fig_flag is True: #return the plot object
        return mesh

    if ax_flag is False and fig_flag is False:
        plt.show() 

    # so now levels_coords[n][0] contains the x coordinates for points along contour n
    # so now levels_coords[n][1] contains the y coordinates for points along contour n

    
