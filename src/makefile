# Makefile for the two executables
# given upstream lb.dat, produce rb.dat, the approximate downstream concentrations
# creates executables CN, the traidtional central difference approximation, and 
# DISCUS, the semi-Lagrangian approach
# More info is available in README

CC=gcc
FCC=gfortran
SOLVERS=../lib
BINDIR=../bin
# -03 might cause problems in discus code
FFLAGS=-O 
CFLAGS=-O3
LDFLAGS=-lm

all: cnlib discuslib discus cn

# Object files used in the executables below
discuslib:	$(SOLVERS)/discus.o 
cnlib:	$(SOLVERS)/cn.o

# executables 
cn:	$(BINDIR)/CN
discus:	$(BINDIR)/DISCUS

$(SOLVERS)/discus.o: ./solve/discus_only_v2.f
	$(FCC) $(FFLAGS) -c ./solve/discus_only_v2.f -o $(SOLVERS)/discus.o

$(SOLVERS)/cn.o: ./solve/testtridJRM.c
	$(CC) $(CFLAGS) -c ./solve/testtridJRM.c -o $(SOLVERS)/cn.o

$(BINDIR)/CN: 	./solve/main.c cnlib
	$(CC) $(CFLAGS) $(LDFLAGS) ./solve/main.c $(SOLVERS)/cn.o -o $(BINDIR)/CN

$(BINDIR)/DISCUS:	./solve/main.c discuslib
	$(CC) $(CFLAGS) $(LDFLAGS) -lgfortran ./solve/main.c $(SOLVERS)/discus.o -o $(BINDIR)/DISCUS
clean:
	rm ../lib/*; rm ../bin/*

