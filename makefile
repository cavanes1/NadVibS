######################################
#                                    #
#       Makefile for NadVibS         #
#                                    #
######################################

compiler = mpif90
GAInc  = /home-4/yshen57@jhu.edu/work/yshen/Software/ga-5-4/include
MPIInc = /software/apps/mpi/openmpi/3.1/intel/18.0/include
src = NadVibS.f90
exe = NadVibS.exe
flag = -I$(GAInc) -I$(MPIInc) -fpp -i8 -auto -assume byterecl -m64 -ipo -O3 -no-prec-div -fp-model fast=2 -march=core-avx2

$(exe): $(src)
	$(compiler) $(flag) $^ libcomex.a libga.a libarmci.a -o $(exe)

clean:
	rm $(exe)
	rm *.mod