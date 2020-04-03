F90    = gfortran
MPIF90 = mpifort

FCOPTS = $(DFLAGS) -O3 -g
LDOPTS = $(DFLAGS) -O3 -g
LIBDIR = /usr/lib
LIBS   = -L$(LIBDIR) -llapack -lblas


