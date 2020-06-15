######################################
#                                    #
#       Makefile for NadVibS         #
#                                    #
######################################

# Flags for user to tune
FFLAGS = -m64 -xCORE-AVX2 -mtune=core-avx2 -no-prec-div -fp-model fast=2 -O3
GA  = 
MPI = 
# On my laptop
# GA  = /home/yshen/Software/Programming/ga-5.4/openmpi-intel
# MPI = /home/yshen/Software/Programming/openmpi-3.1.4/intel
# On MARCC
# GA  = /home-4/yshen57@jhu.edu/data/yshen/Software/ga-5.4_i8
# MPI = /software/apps/mpi/openmpi/3.1/intel/17.0

# User does not have to take care of following variables
compiler = $(MPI)/bin/mpifort
flag = -mkl -I$(GA)/include -I$(MPI)/include -fpp \
-i8 -integer-size 64 -auto -assume byterecl $(FFLAGS)
lib = $(addprefix $(GA)/lib/, libcomex.a libga.a libarmci.a)
src = $(addprefix Source/, progdata.f90 filedata.f90 \
yarkony.f90 function.f90 secondary.f90 top.f90 main.f90)

# release
NadVibS.exe: $(src)
	$(compiler) $(flag) -ipo $^ $(lib) -o $@

# debug
debug.exe: main.o top.o secondary.o function.o yarkony.o progdata.o filedata.o
	$(compiler) $(flag) $^ $(lib) -o $@

%.o: Source/%.f90
	$(compiler) $(flag) -c $<
