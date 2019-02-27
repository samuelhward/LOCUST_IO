### Processing Routines

LOCUST_IO contains a few simple physics routines to process data (please refer to docstrings for instructions): 


     process_input                               routines to manipulate inputs to LOCUST

        QTP_calc                                 calculates the missing quantity out of Q, toroidal or poloidal flux (given two)
        fpolrz_calc                              calculates the 1D flux function on the 2D computational grid
        B_calc                                   calculates the components of the axisymmetric magnetic field on the 2D computational grid
        mag_axis_calc                            calculate location of the magnetic axis
        transform_marker_velocities              transforms marker phase space velocities to LOCUST r,phi,z format

     process_output                              routines to manipulate outputs to LOCUST

        dfn_transform                            transform and integrate distribution function to coordinate system
        dfn_crop                                 crops dfn according to limits in any dimension
        particle_list_compression                opens and processes >>GB LOCUST particle lists in memory-efficient way

    run_scipts

        TRANSP_2_ASCOT                           convert full set of TRANSP inputs to ASCOT inputs
        TRANSP_2_LOCUST                          convert full set of TRANSP inputs to LOCUST inputs
        get_fbm_extract_fi                       extract fast ion distribution from TRANSP output using get_fbm