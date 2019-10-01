######################################
#                                    #
#       Makefile for NadVibS         #
#                                    #
######################################

compiler = mpifort
GAInc  = /home-4/yshen57@jhu.edu/data/yshen/Software/ga-5.4_i8/include
MPIInc = /software/apps/mpi/openmpi/3.1/intel/18.0/include
src = NadVibS.f90
exe = NadVibS.exe
flag = -m64 -xavx -mtune=corei7-avx -I$(GAInc) -I$(MPIInc) -fpp -i8 -integer-size 64 -auto -assume byterecl -ipo -O3 -no-prec-div -fp-model fast=2

$(exe): $(src)
	export BLAS_I8=yes
	$(compiler) $(flag) $^ libcomex.a libga.a libarmci.a -o $(exe)

clean:
	rm *.mod
