#------------------------------------------------------------------------------
#
# *** Description:
#
#     Makefile for LOCUST-gpu.
#
# *** Author:
#
#     Rob Akers, D3/1.36, Culham Centre for Fusion Energy, x6323
#
# *** Creation Date:
#
#     07/11/2012
#
#------------------------------------------------------------------------------

# Compilers:
FC      = pgfortran
CC      = nvcc
CUDALIB = cuda10.1

# Compile flags.  

# Algorithm options: (for POINCARE scan need -DUMP3D)

AOPTS   = $(FLAGS) -DBP -DORI $(UPHASE) $(MPHASE) $(LPHASE) $(NC)

# Choose between pgfortran and ifort:

ifeq ($(FC),pgfortran)
FCFLAGS = -Mpreprocess -mp -fast -Mlarge_arrays $(AOPTS)
else
FCFLAGS = -fpp -openmp -O3 -DINTEL -DCPU $(AOPTS) 
endif

# Target:

mars_build : prec_mod.o disp_mod.o v_spline.o splinck.o ibc_ck.o cspline.o \
	zonfind.o herm1ev.o evspline.o evspline_.o mkspline.o mars_read.o

	$(FC) $(FCFLAGS) -o mars_build \
	prec_mod.o disp_mod.o v_spline.o splinck.o ibc_ck.o cspline.o \
	zonfind.o herm1ev.o evspline.o evspline_.o mkspline.o mars_read.o



# Compile steps:

mars_read.o:mars_read.f90 prec_mod.o disp_mod.o
	$(FC) $(FCFLAGS) -c  mars_read.f90

prec_mod.o:./prec_mod.f90
	$(FC) $(FCFLAGS) -c  ./prec_mod.f90

disp_mod.o:./disp_mod.f90
	$(FC) $(FCFLAGS) -c  ./disp_mod.f90

v_spline.o:./v_spline.f90 prec_mod.o 
	$(FC) $(FCFLAGS) -c  ./v_spline.f90

splinck.o:./splinck.f90 prec_mod.o disp_mod.o
	$(FC) $(FCFLAGS) -c  ./splinck.f90

ibc_ck.o:./ibc_ck.f90 prec_mod.o disp_mod.o
	$(FC) $(FCFLAGS) -c  ./ibc_ck.f90

cspline.o:./cspline.f90 prec_mod.o disp_mod.o
	$(FC) $(FCFLAGS) -c  ./cspline.f90

zonfind.o:./zonfind.f90 prec_mod.o
	$(FC) $(FCFLAGS) -c  ./zonfind.f90

herm1ev.o:./herm1ev.f90 prec_mod.o
	$(FC) $(FCFLAGS) -c  ./herm1ev.f90

evspline.o:./evspline.f90 prec_mod.o
	$(FC) $(FCFLAGS) -c  ./evspline.f90

evspline_.o:./evspline_.f90 prec_mod.o
	$(FC) $(FCFLAGS) -c  ./evspline_.f90

mkspline.o:./mkspline.f90 prec_mod.o disp_mod.o
	$(FC) $(FCFLAGS) -c  ./mkspline.f90

# This entry allows you to type " make clean " to get rid of
# all object and module files
clean:
	 rm -f *.o *.mod *.ptx *.bin mars_build
