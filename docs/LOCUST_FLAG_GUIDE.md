## LOCUST Flag Guide

|      flag     |                                                                 description                                                                  |
|---------------|----------------------------------------------------------------------------------------------------------------------------------------------|
| `LEIID`=1     | McClements, Thyagaraja & Hamilton                                                                                                            |
| `LEIID`=2     | Runge-Kutta Fehlberg                                                                                                                         |
| `LEIID`=3     | Runge-Kutta Cash & Karp                                                                                                                      |
| `LEIID`=4     | Runge-Kutta Dormand & Prince                                                                                                                 |
| `LEIID`=5     | Runge-Kutta Goeken & Johnson                                                                                                                 |
| `LEIID`=6     | Boris                                                                                                                                        |
| `LEIID`=7     | Guiding-centre tracking                                                                                                                      |
| `LEIID`=8     | Strang-Splitting particle trajectory integrator                                                                                              |
| `TOKAMAK`=1   | set tokamak to ITER                                                                                                                          |
| `TOKAMAK`=2   | set tokamak to ASDEX-U                                                                                                                       |
| `TOKAMAK`=3   | set tokamak to MAST-U                                                                                                                        |
| `TOKAMAK`=4   | set tokamak to MAST                                                                                                                          |
| `TOKAMAK`=5   | set tokamak to JET                                                                                                                           |
| `TOKAMAK`=6   | set tokamak to W7-X                                                                                                                          |
| `TOKAMAK`=7   | set tokamak to CTF-York                                                                                                                      |
| `TOKAMAK`=8   | set tokamak to DIII-D                                                                                                                        |
| `TOKAMAK`=9   | set tokamak to ZOW-analytic                                                                                                                  |
| `TOKAMAK`=10  | set tokamak to STEP                                                                                                                          |
| `NOPFC`       | disable 3D CAD PFC model and hit checking                                                                                                    |
| `PFC2D`       | enable 2D limiter profile                                                                                                                    |
| `TOKHEAD`     | set data directory structure to contain tokamak name                                                                                         |
| `JXB2`        | something                                                                                                                                    |
| `PROV`        | use git provenance                                                                                                                           |
| `PITCHCUR`    | define pitch angle as parallel to plasma current instead of B                                                                                |
| `EBASE`       | dump distribution function against energy instead of velocity                                                                                |
| `EGSET`=x     | set maximum energy value of distribution function (overide Vsclh)                                                                            |
| `GEQDSKFIX1`  | flips sign of poloidal flux                                                                                                                  |
| `GEQDSKFIX2`  | flips sign of RB i.e. flips sign of toroidal flux                                                                                            |
| `CURANDOMP`   | Make CURAND wrappers openMP thread-safe                                                                                                      |
| `NTMAX`=N     | Limit iterations to nt=N                                                                                                                     |
| `NTMIN`=N     | Start Dfn. integration at N                                                                                                                  |
| `TIMAX`       | Max tracking time                                                                                                                            |
| `NO_RNG`      | Disable RNG (use with -Mcuda=emu)                                                                                                            |
| `CHEBYSHEV`   | Use Chebyshev polynomials for B                                                                                                              |
| `STATIC`      | Disable tracking                                                                                                                             |
| `TEST`        | Test mode (use tload for source)                                                                                                             |
| `SP`          | Terminate all but ptcl. #1 on each thread                                                                                                    |
| `NPLIM` =I    | Limit # ptcles to I                                                                                                                          |
| `RNGTYP`=1    | RNG generator = CURAND                                                                                                                       |
| `RNGTYP`=1    | RNG generator = SDK Mersenne Twister                                                                                                         |
| `RNGTYP`=1    | RNG generator = F90 intrinsic                                                                                                                |
| `IDFTYP`      | Dfn. storage type.                                                                                                                           |
| `FAVG`        | Used for checking orbit drift (~0.1)                                                                                                         |
| `NO_FLOW`     | Disable toroidal flow                                                                                                                        |
| `BP`          | Full double precision                                                                                                                        |
| `LOOSEBDRY`   | Allow a loose boundary for IDFTYP=3                                                                                                          |
| `FNLOAD`      | Load c.of.m DFn. for IHAG=2 test                                                                                                             |
| `POSMU`       | Force mu>=0 in DFn. (but not J)                                                                                                              |
| `BCOREC`      | Higher order Boris correction                                                                                                                |
| `MOVIE`       | Generate evolution movie data.                                                                                                               |
| `C1`          | Zeroth order collision operator                                                                                                              |
| `C2`          | 1st order collision operator                                                                                                                 |
| `CT`          | Disable thermal accumultion for -DC2                                                                                                         |
| `FO1`         | Trigger a re-goose every full orbit travelling upwards                                                                                       |
| `FO2`         | Trigger a re-goose every full orbit travelling downwards                                                                                     |
| `FO3`         | Trigger a re-goose upwards if V>=0, downwards otherwise                                                                                      |
| `FO4`         | Trigger a re-goose upwards if V<=0, downwards otherwise                                                                                      |
| `FO5`         | Re-goose as F03 for even ptcles, F04 for odd (compression will scramble)                                                                     |
| `BRIPPLE`     | Apply TF ripple (Boris only)                                                                                                                 |
| `ISEE`        | CURAND PRNG shift parameter (>=1)                                                                                                            |
| `VERBOSE`     | Turn on verbose output                                                                                                                       |
| `TXTMEM`      | Enable B-field texture memory (>=pgi12.8)                                                                                                    |
| `BCACHE`      | Turn "off" B-field caching                                                                                                                   |
| `TRILIN`      | Use tri-linear interpolation rather than sampling in dfn_bin_kernel                                                                          |
| `GTPB`        | Don't inhibit goosing for particles close to the trapped-passing boundary                                                                    |
| `NODUMP`      | Don't dump data files                                                                                                                        |
| `SMALLEQ`     | Allow 2D eqbm. to be smaller than void                                                                                                       |
| `OPENMESH`    | Allow open vacuum cells                                                                                                                      |
| `OPENTRACK`   | Allow ptcls born on mesh to soft leave                                                                                                       |
| `OPENTERM`    | Allow ptcls born off mesh to track                                                                                                           |
| `PFCL`        | Do not attempt rigourous PFC tracking                                                                                                        |
| `NOCOMPRESS`  | Don't apply compression every idsp steps                                                                                                     |
| `BICUB`       | Use method 2 bi-cubic interpolation                                                                                                          |
| `B3D`         | Enable 3D field                                                                                                                              |
| `BILIN`       | Use bilinear interpolation of 3D data                                                                                                        |
| `POINCARE`    | Generate Poincare map                                                                                                                        |
| `RKFIX`       | Use fixed step RK in Poincare map                                                                                                            |
| `RKBS23`      | Use BS23 Runge Kutta in GUIDE3                                                                                                               |
| `RKF45`       | Use Fehlberg Runge Kutta in GUIDE3                                                                                                           |
| `P5INT`       | Use 5th order solution to advance eqns.                                                                                                      |
| `HARD`        | Force hard failure of RK45                                                                                                                   |
| `GCFIX`       | Force G.C. to conserve CoM (2D B only) by spatially moving particle                                                                          |
| `CLM`         | Use OLD (incorrect) G.C. equations                                                                                                           |
| `NOGCCOR`     | Don't adjust ptcl. V to exactly match G.C. position                                                                                          |
| `GCCOL`       | Set to force scattering at G.C. (standard)                                                                                                    |
| `GCFG`        | If scattering at ptcl., set to force permanent switch to gyro tracking after V/V > 0.9, else G.C. <-> full gyro will toggle dependent on V/V |
| `G3IPFC0`     | In SPLIT mode, witch to full gyro for STAGE 2.                                                                                               |
| `ECUT`        | Energy termination cut [keV] (set = 10.0 to match ASCOT)                                                                                     |
| `PCYC`        | Allow ptcl. input file to wrap around                                                                                                        |
| `PSIREJ`      | Reject particles with PSIn<PSIREJ                                                                                                            |
| `SPLIT`       | SPLIT mode                                                                                                                                   |
| `SPALL`       | If not set, one chunk will be followed for stage 1, then splitting will occur up to nchunks for STAGE 2                                      |
| `BTASCOT`     | Para/diamagnetism is included in the 3D perturbation field                                                                                   |
| `SLOAD`       | STAGE 2 run                                                                                                                                  |
| `SLOAD2`      | ??? REMOVE ???                                                                                                                               |
| `PFCCON`      | Renormalize so that STAGE2 collected power matches STAGE 1                                                                                   |
| `PFCF`        | Only allow power to be collected to STAGE 1 matching component (-DPFCCON must be set)                                                        |
| `STELLARATOR` | Plasma device is a stellarator                                                                                                               |
| `SS`          | Apply stellarator symmetry rule                                                                                                              |
| `SSLIN`       | Reduce tri-cubic interpolation to tri-linear                                                                                                 |
| `TOK3D`       | Tokamak is 3D (bin toridal angle rather than gyro-phase)                                                                                     |
| `BPFLIP`      | Flip "2D" poloidal field (3D field perturbation is unaffected)                                                                               |
| `BTFLIP`      | Flip "2D" toroidal field                                                                                                                     |
| `GPUTIME`     | Enable GPU timing info                                                                                                                       |
| `WIPE`        | Store accumulated rather than scaled intermediate DFn.                                                                                       |
| `B3D_EX`      | 3D field perturbation is in MARS compressed format                                                                                           |
| `AVEC`        | 3D     -"-        in A-vector form                                                                                                           |
| `AR_GAUGE`    | Re-factor B3D_EX format B field into A-vector form (AR=0 gauge)                                                                              |
| `AT_GAUGE`    | -"-            (AT=0  -"- )                                                                                                                  |
| `AZ_GAUGE`    | -"-            (AZ=0  -"- )                                                                                                                  |
| `POINGR`      | Follow gyro-orbits when generating Poincare map                                                                                              |
| `ISEL`        | Force all ptcle states to ptcl. ISEL                                                                                                         |
| `VSET`        | Force all ptcl. velocities to VSET                                                                                                           |
| `ESET`        | Force all ptcl. energies   to ESET                                                                                                           |
| `UMP`         | Dump diagnostic data                                                                                                                         |
| `ISOL`        | Set magnetic field to Solov'ev analytic equilibrium                                                                                          |
| `BANAL`       | Set       -"-         infinite aspect ratio elliptical analytic equilibrium                                                                  |
| `HO`          | Re-goose every "half" poloidal orbit instead of "full"                                                                                       |
| `FILT`        | Apply median FILTER to CoM DFn.                                                                                                              |
| `BNCSOFT`     | Allow up to 100 bounces in Poincare plot tracking for STELLARATOR                                                                            |
| `GUIDE2`      | Higher order correction in G.C. tracking                                                                                                     |
| `MUITER`      | More corrections                                                                                                                             |
| `MU0`         | and even more.                                                                                                                               |
| `MUP2`        | And even more....this needs some work.                                                                                                       |
| `MUSEC`       | 2nd iteration in mu calculation this needs to be removed)                                                                                    |
| `NOGC`        | Do not re-position G.C. away from ptcl.                                                                                                      |
| `LEG`         | Determine divertor legs in B-field tuning process                                                                                            |
| `TETALL`      | Disable dual cycle tracking (always track in volume tets)                                                                                    |
| `TETCHECK`    | Extra level of tetrahedral cell conformity checking                                                                                          |
| `BCHECK`      | B-field diagnostic                                                                                                                           |
| `NCMACCESS`   | Overload fields with "constants" for memory access testing                                                                                   |
| `CNTSUP`      | Use with and w/o NCMACCESS for memory access testing                                                                                         |
| `PLOT` = N    | Plot N orbits                                                                                                                                |
| `PINFL`       | Pass Poincare data with n uncompressed                                                                                                       |
| `NOTUNE`      | Do not invoke B-field tuning                                                                                                                 |
| `WREAL`       | Pass back REAL weight to D.Fn. sum-redn.                                                                                                     |
| `INFLATE`     | Inflate markers in CYCLE II (w/o -DSPLIT)                                                                                                    |
| `FRQLOG`      | Generate procession and bounce Freq. data                                                                                                    |
| `CONLY`       | Delete PSIn>1 region of kinetic profiles                                                                                                     |
| `PROFPSI2`    | Kinetic profiles against sqrt(PSI)                                                                                                           |
| `TFILE`       | Convergence testing mode                                                                                                                     |
| `TTYPE1`      | Take <Rmin>, <Rmax> instead of Rmin,Rmax                                                                                                     |
| `RZGC`        | Fill DFn. at [R_gc,Z_gc]                                                                                                                     |
| `VPGC`        | -"-     [V_gc, V]                                                                                                                            |
| `MLOAD`       | Uniform source generator                                                                                                                     |
| `RANEN`       | Spread out -DMLOAD source in energy                                                                                                          |
| `ANALOGUE`    | Analogue mode (for GVR)                                                                                                                      |
| `LIMO` = Rlim | Artificial LFS limiter in POINCARE                                                                                                           |
| `LIMI` = Rlim | Artificial HFS limiter in POINCARE                                                                                                           |
| `VCCORE`      | Thermalize at lowest thermal ion velocity                                                                                                    |
| `OMEGAT`      | Enable 3D mode rotation                                                                                                                      |
| `EBASE`       | DFn. vs. energy rather than velocity                                                                                                         |
| `NEUT`        | Turn on CX + neutral tracking                                                                                                                |
| `CURT`        | Reduce binning frequency (test mode)                                                                                                         |
| `FINT`        | Store f(E) at every display period                                                                                                           |
| `FGFORCE`     | Apply "fixed Goosing" (test mode)                                                                                                            |
| `PINSEP`      | Set plasma constant outside LCFS                                                                                                             |
| `SOLCOL`      | Only apply collisions inside LCFS                                                                                                            |
| `GCBINB`      | Bin 1D profiles at the G.C. position                                                                                                         |
| `TREDG`       | Cut out ptcles. born at LCFS                                                                                                                 |
| `TRANBRTH`    | Ptcl. list is from TRANSP birth.cdf file                                                                                                     |
| `WLIST`       | Load in wght together with [x,v]                                                                                                             |
| `LNLNOB`      | Ignore B dependence in lnLambda                                                                                                              |
| `LNLNOE`      | Ignore E dependence in lnLambda                                                                                                              |
| `UHST`        | Allow Unresolved hits (soft termination)                                                                                                     |
| `LNLBT`       | Us BT in lnL instead of B                                                                                                                    |
| `BRELAX`      | Ignore markers that are off grid                                                                                                             |
| `PFCMOD`      | Output splot VTK files only if changed                                                                                                       |
| `TITHRM`      | Set thermalisation cut-off energy to TITHRM*Ti                                                                                               |
| `NCOILS`      | Specify total number of separate RMP coil rows                                                                                               |
| `VROT`        | Load rotation profile data                                                                                                                   |
|               |                                                                                                                                              |
