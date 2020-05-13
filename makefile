######################################
#                                    #
#       Makefile for NadVibS         #
#                                    #
######################################

# Flags for user to tune
compiler = mpifort
FFLAGS = -m64 -xCORE-AVX2 -mtune=core-avx2 -no-prec-div -fp-model fast=2 -O3
GA  = /home-4/yshen57@jhu.edu/data/yshen/Software/ga-5.4_i8
MPI = /software/apps/mpi/openmpi/3.1/intel/18.0

# User does not have to take care of following variables
flag = -mkl -I$(GA)/include -I$(MPI)/include -fpp \
-i8 -integer-size 64 -auto -assume byterecl $(FFLAGS)
lib = $(addprefix $(GA)/lib/, libcomex.a libga.a libarmci.a)
src = $(addprefix source/, progdata.f90 filedata.f90 \
yarkony.f90 function.f90 secondary.f90 top.f90 main.f90)
# Faster compilation for debugging
obj = progdata.o filedata.o yarkony.o function.o secondary.o top.o main.o

# release
NadVibS.exe: $(src)
	$(compiler) $(flag) -ipo $^ $(lib) -o $@

# debug
debug.exe: $(obj)
	$(compiler) $(flag) $^ $(lib) -o $@

progdata.o: source/progdata.f90
	$(compiler) $(flag) -c $^

filedata.o: source/filedata.f90
	$(compiler) $(flag) -c $^

yarkony.o: source/yarkony.f90
	$(compiler) $(flag) -c $^

function.o: source/function.f90
	$(compiler) $(flag) -c $^

secondary.o: source/secondary.f90
	$(compiler) $(flag) -c $^

top.o: source/top.f90
	$(compiler) $(flag) -c $^

main.o: source/main.f90
	$(compiler) $(flag) -c $^