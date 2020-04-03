# Makefile for semiinf
sinclude make.inc

default: all

all: pysemiinf

pysemiinf:     
	( cd src ; $(MAKE) all || exit 1 )

clean : 
	( cd src ; $(MAKE) clean )
