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
        simag                                                           #poloidal flux psi at magnetic axis (Weber / rad)
        sibry                                                           #poloidal flux psi at plasma boundary (Weber / rad)
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
        R_1D                                                            #R dimension (m)
        Z_1D                                                            #Z dimension (m)
        flux_pol                                                        #poloidal flux from magnetic axis up to the plasma boundary (Weber / rad)
        flux_tor                                                        #toroidal flux from magnetic axis up to the plasma boundary (Weber / rad)
    2D data
        psirz                                                           #poloidal flux at coordinate [r,z] in (Weber / rad) 
        fpolrz                                                          #poloidal current function at coordinate [r,z] (m T) calculated by fpolrz_calc
    3D data
        B_field                                                         #magnetic field component [i] at position [r,z] (T) calculate by B_calc

#### Beam Deposition:

    1D data
        R                                                               #R coordinates of particle p
        phi                                                             #toroidal angle of particle p [rad]
        Z                                                               #Z coordinates of particle p
        V_R                                                             #R component of v of particle p
        V_tor                                                           #toroidal component of v of particle p
        V_Z                                                             #Z component of v of particle p
        V_pitch                                                         #v_parallel/v
        E                                                               #energy of particle (eV)
        weight                                                          #Monte Carlo weights of each particle
        absorption_fraction                                             #absorption fraction
        absorption_scaling                                              #scaling to allow for missing particles
        number_particles                                                #number of particles in beam (if multiple numbers, .combine() may have been used)

#### Temperature:

    1D data
        flux_pol                                                        #poloidal flux (Weber / rad)
        flux_pol_norm                                                   #normalised poloidal flux
        T                                                               #ion temperature (eV)
        T                                                               #electron temperature (eV)
        flux_tor_coord                                                  #toroidal flux coordinate
        flux_tor                                                        #toroidal flux (Weber / rad)
        q                                                               #safety factor

#### Number Density:

    1D data
        flux_pol                                                        #poloidal flux (Weber / rad)
        flux_pol_norm                                                   #normalised poloidal flux
        n                                                               #ion or electron number density (#/m^3)
        flux_tor_coord                                                  #toroidal flux coordinate
        flux_tor                                                        #toroidal flux (Weber / rad)
        q                                                               #safety factor

##### Perturbation

    2D data
        R_2D                                                            #R coordinate of field grid (m)
        Z_2D                                                            #Z coordinate of field grid (m)
        B_field_R_real                                                  #real R component of magnetic field (T)
        B_field_R_imag                                                  #imaginary R component of magnetic field (T)
        B_field_Z_real                                                  #real Z component of magnetic field (T)
        B_field_Z_imag                                                  #imaginary Z component of magnetic field (T)
        B_field_tor_real                                                #real toroidal component of magnetic field (T)
        B_field_tor_imag                                                #imaginary toroidal component of magnetic field (T)

#### Orbits:

    0D data
        number_particles                                                #total number of particles
        number_timesteps                                                #total number of timesteps
    3D data
        orbits[t,i,p]                                                   #spatial coordinate i for particle p at time step t

#### Final Particle List:

    0D data
        n                                                               #number of particles per GPU (blocks per grid * threads per block)
        ngpu                                                            #number of GPUs (OMP threads)
        niter                                                           #
        npt_                                                            #number of particle info slots
        nphc                                                            #
        ntri                                                            #
        number_particles                                                #total number of particles = n*ngpu
    1D data
        R                                                               #r coordinate of particle
        phi                                                             #phi coordinate of particle [rad]
        Z                                                               #z coordinate of particle
        V_R                                                             #v_r coordinate of particle
        V_tor                                                           #v_tor coordinate of particle
        V_Z                                                             #v_z coordinate of particle
        V_pitch                                                         #v_parallel/v
        energy                                                          #energy of particle
        time                                                            #time coordinate of particle
        status_flag                                                     #status value of particle at this time
        status_flags                                                    #possible status flags and their associated values

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
        nV_pitch                                                        #number of points in V_tor/V dimension
        dV_pitch                                                        #pitch bin width
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
        flux_pol_norm_sqrt                                              # 
        flux_pol_norm                                                   #
        dVOL=(dVOL)                                                     #
        beam_source                                                     #
        density                                                         #
        energy_para                                                     #
        energy_perp                                                     #
        energy                                                          #
        J(NBCD)-raw                                                     #
        NBI-heating-power(TOT)                                          #
        NBI-heating-power(e-)                                           #
        NBI-heating-power(i1)                                           #
        residual-angular-momentum-density                               #
        torque-density(JxB-inst)                                        #
        torque-density(JxB-sweep)                                       #
        torque-density(coll)                                            #