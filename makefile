######################################
#                                    #
#       Makefile for NadVibS         #
#                                    #
######################################

# Flags for user to tune
GA  = /home/cavanes1/C/ga-5.8.2/GA
MPI = /data/apps/linux-centos8-cascadelake/intel-19.1.2.254/openmpi-3.1.6-d5iqzgx77lspks5kuha5bexs6ikhjrl7

OBJ = progdata.o filedata.o yarkony.o function.o secondary.o top.o main.o

# User does not have to take care of following variables
compiler = $(MPI)/bin/mpif90
FFLAGS = -O3 -fpp -i8 -auto -assume byterecl -I$(GA)/include -I$(MPI)/include
lib = $(addprefix $(GA)/lib/, libcomex.a libga.a libarmci.a)
LDFLAGS = -L/work1/soft/intel2018/mkl/lib/intel64 -lmkl_scalapack_ilp64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_ilp64 -lpthread -lm

NadVibS.exe: $(OBJ)
	$(compiler) $(FFLAGS) $(OBJ) -o NadVibS.exe $(lib) $(LDFLAGS)

%.o: %.f90
	$(compiler) $(FFLAGS) -c $< -o $@

clean:
	rm -f *.o *.mod a.out *.x *.a *.exe
