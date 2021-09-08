### Data Dictionary

Since this package aims to bridge the gap between various file formats for different LOCUST inputs/outputs, I've included a description of the different variable names used within LOCUST_IO.

#### Equilibrium:

    0D data
        nR_1D                                                           #number of points in R (x or width)
        nZ_1D                                                           #number of points in Z (y or height)
        IDUM                                                            #dummy variable
        rdim                                                            #size of the R dimension in m
        zdim                                                            #size of the Z dimension in m
        rcentr                                                          #reference value of R at magnetic axis
        zcentr                                                          #reference value of Z at magnetic axis
        bcentr                                                          #vacuum toroidal magnetic field at rcentr
        rleft                                                           #R at left (inner) boundary
        zmid                                                            #Z at middle of domain (from origin)
        rmaxis                                                          #R at magnetic axis (O-point)
        zmaxis                                                          #Z at magnetic axis (O-point)
        simag                                                           #poloidal flux psi at magnetic axis [Weber / rad]
        sibry                                                           #poloidal flux psi at plasma boundary [Weber / rad]
        current                                                         #plasma current [Amps]   
        xdum                                                            #dummy variable - just contains int(zero)
        lcfs_n                                                          #number of points in the plasma boundary
        limitr                                                          #number of points in the wall boundary
    1D data
        fpol                                                            #poloidal flux function on uniform flux grid (1D array of f(psi)=R*B_toroidal [meter-Tesla]) (negative for positive plasma current in GEQDSK)
        pres                                                            #plasma pressure in nt/m^2 on uniform flux grid (1D array of p(psi) [Pascals])
        ffprime                                                         #F*d(F)/d(psi) where psi is poloidal flux per radian and  F=diamagnetic function=R*B_Phi 
        pprime                                                          #plasma pressure * d(plasma pressure)/d(psi) (check papers)
        qpsi                                                            #q profile (function of flux_pol)
        rlim                                                            #r coordinates of wall boundary
        zlim                                                            #z coordinates of wall boundary
        lcfs_r                                                          #r coordinates of plasma boundary
        lcfs_z                                                          #z coordinates of plasma boundary
        R_1D                                                            #R dimension [m]
        Z_1D                                                            #Z dimension [m]
        flux_pol                                                        #poloidal flux from magnetic axis up to the plasma boundary [Weber / rad]
        flux_tor                                                        #toroidal flux from magnetic axis up to the plasma boundary [Weber / rad]
        flux_tor_coord                                                  #toroidal flux coordinate sqrt(b_flux_tor/(pi*b0)) ~ sqrt(pi*r^2*b0/(pi*b0)) ~ r [m]
    2D data
        psirz                                                           #poloidal flux at coordinate [r,z] in [Weber / rad] 
        phirz                                                           #toroidal flux at coordinate [r,z] in [Weber / rad] 
        fpolrz                                                          #poloidal current function at coordinate [r,z] [m T] calculated by fpolrz_calc
        B_field_R                                                       #magnetic field Rcomponent at position [r,z] [T] calculate by B_calc
        B_field_tor                                                     #magnetic field phi component at position [r,z] [T] calculate by B_calc
        B_field_Z                                                       #magnetic field Z component at position [r,z] [T] calculate by B_calc

#### Beam Deposition:

    1D data
        R                                                               #R coordinates of particle p
        phi                                                             #toroidal angle of particle p [rad]
        Z                                                               #Z coordinates of particle p
        V_R                                                             #R component of v of particle p
        V_phi                                                           #toroidal component of v of particle p
        V_Z                                                             #Z component of v of particle p
        V_pitch                                                         #v_parallel/v
        E                                                               #energy of particle [eV]
        weight                                                          #Monte Carlo weights of each particle
        absorption_fraction                                             #absorption fraction
        absorption_scaling                                              #scaling to allow for missing particles
        number_particles                                                #number of particles in beam (if multiple numbers, .combine() may have been used)

#### Temperature:

    1D data
        T                                                               #temperature [eV]
        flux_pol                                                        #poloidal flux [Weber / rad]
        flux_pol_norm                                                   #normalised poloidal flux
        flux_pol_norm_sqrt                                              #sqrt(flux_pol_norm) 
        flux_tor_coord                                                  #toroidal flux coordinate sqrt(b_flux_tor/(pi*b0)) ~ sqrt(pi*r^2*b0/(pi*b0)) ~ r [m]
        flux_tor_coord_norm                                             #normalised toroidal flux coordinate
        flux_tor                                                        #toroidal flux [Weber / rad]
        q                                                               #safety factor
        r_1d                                                            #minor radius [metres]

#### Number Density:

    1D data
        n                                                               #number density [#/m^3]
        flux_pol                                                        #poloidal flux [Weber / rad]
        flux_pol_norm                                                   #normalised poloidal flux
        flux_pol_norm_sqrt                                              #sqrt(flux_pol_norm) 
        flux_tor_coord                                                  #toroidal flux coordinate sqrt(b_flux_tor/(pi*b0)) ~ sqrt(pi*r^2*b0/(pi*b0)) ~ r [m]
        flux_tor_coord_norm                                             #normalised toroidal flux coordinate
        flux_tor                                                        #toroidal flux [Weber / rad]
        q                                                               #safety factor
        r_1d                                                            #minor radius [metres]

#### Rotation

    1D data
        rotation_ang                                                    #rotation profile [rad/s]
        rotation_vel                                                    #rotation profile [m/s]
        flux_pol                                                        #poloidal flux [Weber / rad]
        flux_pol_norm                                                   #normalised poloidal flux
        flux_pol_norm_sqrt                                              #sqrt(flux_pol_norm) 
        flux_tor_coord                                                  #toroidal flux coordinate sqrt(b_flux_tor/(pi*b0)) ~ sqrt(pi*r^2*b0/(pi*b0)) ~ r [m]
        flux_tor_coord_norm                                             #normalised toroidal flux coordinate
        R_1D                                                            #major radius where corresponding rotation_vel is measured [m]

##### Perturbation

    1D data
        R_1D                                                            #R dimension [m]
        Z_1D                                                            #Z dimension [m]
        R_point_data                                                    #R coordinate of point_data.inp point for B field checking [m]
        phi_point_data                                                  #phi coordinate of point_data.inp point for B field checking [m]
        Z_point_data                                                    #Z coordinate of point_data.inp point for B field checking [m]
        time_point_data                                                 #time coordinate of point_data.inp point for B field checking [m]
        X_point_data                                                    #X coordinate of point_data.inp point for B field checking [m]
        Y_point_data                                                    #Y coordinate of point_data.inp point for B field checking [m]
    2D data
        R_2D                                                            #R coordinate of field grid [m]
        Z_2D                                                            #Z coordinate of field grid [m]
        dB_field_R_real                                                 #real R component of magnetic field [T] at point [r,z]
        dB_field_R_imag                                                 #imaginary R component of magnetic field [T] at point [r,z]
        dB_field_R                                                      #magnitude of R component of magnetic field [T] at point [r,z]
        dB_field_Z_real                                                 #real Z component of magnetic field [T] at point [r,z]
        dB_field_Z_imag                                                 #imaginary Z component of magnetic field [T] at point [r,z]
        dB_field_Z                                                      #magnitude of Z component of magnetic field [T] at point [r,z]
        dB_field_tor_real                                               #real toroidal component of magnetic field [T] at point [r,z]
        dB_field_tor_imag                                               #imaginary toroidal component of magnetic field [T] at point [r,z]
        dB_field_tor                                                    #magnitude of toroidal component of magnetic field [T] at point [r,z]
        B_field_R                                                       #magnetic field R component at position [r,z] [T] 
        B_field_tor                                                     #magnetic field phi component at position [r,z] [T]
        B_field_Z                                                       #magnetic field Z component at position [r,z] [T]
        B_field_mag                                                     #absolute magnitude of perturbation [T]

##### Wall

    0D data
        limitr                                                          #number of points in the wall boundary (2D wall)
    1D data
        rlim                                                            #r coordinates of wall boundary (2D wall)
        zlim                                                            #z coordinates of wall boundary (2D wall)

#### Orbits:

    0D data
        number_particles                                                #total number of particles
        number_timesteps                                                #total number of timesteps
    2D data
        R[t,p]                                                          #R for particle p at time step t
        phi[t,p]                                                        #phi for particle p at time step t
        Z[t,p]                                                          #Z for particle p at time step t

#### Final Particle List:

    0D data
        n                                                               #number of particles per GPU (blocks per grid * threads per block)
        ngpu                                                            #number of GPUs (OMP threads)
        niter                                                           #number simulation iterations
        npt_                                                            #number of particle info slots
        nphc                                                            #number particle cache levels (0th is final position, last is initial position)
        ntri                                                            #
        number_particles                                                #total number of particles = n*ngpu
    1D data
        R                                                               #R coordinate of particle  
        phi                                                             #phi coordinate of particle [rad]  
        Z                                                               #Z coordinate of particle  
        V_R                                                             #V_R coordinate of particle [m/s]   
        V_phi                                                           #V_phi coordinate of particle [m/s]  
        V_Z                                                             #V_Z coordinate of particle [m/s]  
        time                                                            #time coordinate of particle [s]
        dt                                                              #particle status flag  
        FG                                                              #particle weight [#/s]
        tet                                                             #tet ID holding particle
        V_R_next                                                        #V_R coordinate of particle at next timestep [m/s]  
        V_phi_next                                                      #V_phi coordinate of particle at next timestep [m/s]  
        V_Z_next                                                        #V_Z coordinate of particle at next timestep [m/s]  
        R_next                                                          #R coordinate of particle at next timestep   
        phi_next                                                        #phi coordinate of particle at next timestep   
        Z_next                                                          #Z coordinate of particle at next timestep   
        psi_initial                                                     #psi coordinate of particle at start of simulation   
        R_initial                                                       #R coordinate of particle at start of simulation   
        phi_initial                                                     #phi coordinate of particle at start of simulation   
        Z_initial                                                       #Z coordinate of particle at start of simulation   
        V_R_initial                                                     #V_R coordinate of particle at start of simulation   
        V_phi_initial                                                   #V_phi coordinate of particle at start of simulation   
        V_Z_initial                                                     #V_Z coordinate of particle at start of simulation   
        time_initial                                                    #time coordinate of particle at start of simulation   
        dt_initial                                                      #dt coordinate of particle at start of simulation   
        FG_initial                                                      #FG coordinate of particle at start of simulation   
        tet_initial                                                     #tet coordinate of particle at start of simulation   
        V_R_next_initial                                                #V_R_next coordinate of particle at start of simulation   
        V_phi_next_initial                                              #V_phi_next coordinate of particle at start of simulation   
        V_Z_next_initial                                                #V_Z_next coordinate of particle at start of simulation   
        R_next_initial                                                  #R_next coordinate of particle at start of simulation   
        phi_next_initial                                                #phi_next coordinate of particle at start of simulation   
        Z_next_initial                                                  #Z_next coordinate of particle at start of simulation   
        weight                                                          #particle weight [#/s]  
        f                                                               #power scaling; power of each marker = E*weight*f*1e6 
        E                                                               #marker energy [eV]
        status_flag                                                     #status value of particle at simulation end
        status_flags                                                    #numerical values corresponding to possible status_flags 
#### Distribution Function:

    0D data
        IDFTYP                                                          #dfn structure ID (= IDFTYP LOCUST flag)
        EQBM_MD5                                                        #checksum ID of the equilibrium
        nE                                                              #number of points in energy dimension
        dEh                                                             #energy bin width
        nV                                                              #number of points in velocity dimension
        dV                                                              #velocity bin width
        nR                                                              #number of R points on Dfn grid
        dR                                                              #R bin width
        nZ                                                              #number of Z points on Dfn grid
        dZ                                                              #Z bin width
        nV_pitch                                                        #number of points in V_phi/V dimension
        dV_pitch                                                        #pitch bin width
        dP_phi                                                          #toroidal canonical angular momentum bin width
        dmu                                                             #magnetic moment bin width
        nP                                                              #poloidal gyro-phase cell boundaries
        dP                                                              #special dimension bin width - simulation specific (e.g. gyrophase bin width)
        nPP                                                             #number of points in Pphi dimension
        nMU                                                             #number of points in Mu dimension
        nPSI                                                            #number of poloidal flux surface contours
        nPOL                                                            #number of poloidal cells for volume calculation
        dPPh                                                            #Pphi bin width
        dMUh                                                            #Mu bin width
        cpu_time                                                        #
        0_1                                                             #zero indicates fast ions only
        nc                                                              #number of slices that Dfn is written to file in
        nR_1D                                                           #number of R points on field grid
        nZ_1D                                                           #number of Z points on field grid
        Ai_1                                                            #first value of Ab
        Zi_1                                                            #first value of Zb
        Vsclh                                                           #velocity grid upper limit
        Vnrm                                                            #normalising velocity
        icoll                                                           #collisions (1=on 0=off)
        iscat                                                           #scattering (1=on 0=off)
        idiff                                                           #energy diffusion (1=on 0=off)
        iloss                                                           #charge exchange events (1=on 0=off)
        iterm                                                           #terminate if ptcl. leaves plasma
        niter                                                           #number of iterations for isym=1 simulation
        integrator                                                      #integrator type
        npnt                                                            #points per gyration
        0_2                                                             #the maximum integrator step size
        dt0                                                             #
        threads_per_block                                               #gpu threads per block
        blocks_per_grid                                                 #gpu blocks per grid
    1D data
        E                                                               #energy dimension of Dfn grid
        V                                                               #velocity dimension of Dfn grid
        PP                                                              #Pphi dimension of Dfn grid
        MU                                                              #Mu dimension of Dfn grid
        Jh                                                              #Jacobian
        Jh_s                                                            #Jacobian error
        R_1D                                                            #R dimension of field grid
        Z_1D                                                            #Z dimension of field grid
        npolh                                                           #
        R                                                               #R dimension of Dfn grid
        Z                                                               #Z dimension of Dfn grid
        V_pitch                                                         #pitch dimension of Dfn grid
        P_phi                                                           #toroidal canonical angular momentum
        dmu                                                             #magnetic moment
        P                                                               #special dimension - simulation specific (e.g. gyrophase)
        Ab                                                              #fast ion masses
        Zb                                                              #trace particle Z
        Pdep                                                            #normalised injected power
        tau_s                                                           #zeroth order slowing-down time
        E0                                                              #energy (plasma frame)
        EC                                                              #zeroth order critical energy
        rho                                                             #
        siglg                                                           #r.m.s. width of test src
        dfn_index                                                       #reference for names of each dfn dimension         
    2D data
        psirz                                                           #poloidal flux field [r,z] grid
        dVOL                                                            #
    3D data        
        dfn                                                             #Dfn grid (for IDFTYP=3)
        dfn_s                                                           #Dfn grid error (for IDFTYP=3)
    4D data
        csb                                                             #volume element data
    5D data
        dfn                                                             #Dfn grid (for IDFTYP!=3)
        dfn_s                                                           #Dfn grid error (for IDFTYP!=3)

LOCUST dumps distribution functions in unformatted binary format. Different run-time flag combinations will dictate whether the above fields are written to file. The Dfn grid is defined at the bin centres in SI units [s^3/m^6] (one must integrate to get particles per bin for plotting).
    
##### Moments

There are many moments of the distribution function which are provided by various codes, included below are those used by LOCUST_IO from LOCUST, however reading outputs from other codes may yield additional moments not yet shown here. Eventually all moments should be mapped and calculated to common variable names. 

    1D data
        flux_pol_norm                                                   #normalised poloidal flux
        flux_pol_norm_sqrt                                              #sqrt(flux_pol_norm) 
        dVOL                                                            #volume element size [m**3]
        beam_source                                                     #beam deposition density
        density                                                         #total beam density
        energy_para                                                     #parallel kinetic beam energy
        energy_perp                                                     #perpendicular kinetic beam energy
        energy                                                          #total beam energy
        J(NBCD)-raw                                                     #neutral beam current drive
        NBI-heating-power(TOT)                                          #total neutral beam heating power
        NBI-heating-power(e-)                                           #total neutral beam heating power to electrons
        NBI-heating-power(i1)                                           #total neutral beam heating power to ion species 1
        residual-angular-momentum-density                               #
        torque-density(JxB-inst)                                        #instantaneous JxB torque density
        torque-density(JxB-sweep)                                       #orbit-averaged JxB torque density
        torque-density(coll)                                            #collision torque density

##### Rundata

    0D data
        PFC_power                                                       
            component                                                   #PFC loss power resolved by component [W]
            total                                                       #total PFC loss power [W]
        time_total                                                      #kernel end time

##### Poincare

    0D data
        EQBM_MD5                                                        #checksum ID of the equilibrium
        nRP                                                             #number R zones
        nZP                                                             #number Z zones
        nTP                                                             #number toroidal zones
        npln                                                            #binning planes
        R0P                                                             #R domain minimum
        R1P                                                             #R domain maximum
        Z0P                                                             #Z domain minimum
        Z1P                                                             #Z domain maximum
        nEQ_R                                                           #number equilibrium R zones
        nEQ_Z                                                           #number equilibrium Z zones
        psi_equil_h                                                     #
    1D data
        R_equil_h                                                       #R equilibrium dimension [m]
        Z_equil_h                                                       #Z equilibrium dimension [m]
        R                                                               #R map dimension [m]
        Z                                                               #Z map dimension [m]
        phi                                                             #phi dimension [rad]
    3D data
        map                                                             #Poincare map [R,Z,phi]

##### Output Mesh

    1D data
        X                                                               #X coordinate of node [m]
        Y                                                               #Y coordinate of node [m]
        Z                                                               #Z coordinate of node [m]
        R                                                               #R coordinate of node [m]
        phi                                                             #phi coordinate of node [rad]
