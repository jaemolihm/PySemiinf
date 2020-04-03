# Makefile for PySemiinf
include ../make.inc

default: all

all: kinds pysemiinf wanniertb

kinds:
	$(F90) ${FCOPTS} ${LIBS} -c kinds.f90

pysemiinf:
	$(F90) ${FCOPTS} ${LIBS} -c comms.f90
	f2py -L/usr/lib -llapack -lblas -c kinds.f90 comms.f90 semiinf.f90 -m PySemiinf

wanniertb:
	f2py -L/usr/lib -llapack -lblas -c kinds.f90  get_hk.f90 -m get_hk

clean:
	rm -f *.o *.so *.mod
