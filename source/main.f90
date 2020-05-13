!NadVibS: nonadiabatic vibronic spectrum simulation package
!Originate from NADVIBS.X by Michael Schuurman 2007
!
!Source code layout:
!    Modules: progdata, filedata
!        - contains the equivalent of global variables/common blocks
!        - progdata contains predominantly arithmetic variables for the running of the lanczos program
!        - filedata contains variables used for file I/O
!    program main
!        - NadVibS program loop
!        - this is the driver for the program. The lanczos iterations are run from this uppermost loop
!    Top level subroutines
!        - read_cmdline,read_basis,print_basis,read_lanczos,allocate_memory,intialize_elements,
!          update_lanczos,Hv,lanczos_step,orthog_vectors,compute_eigenvalues,print_footer,release_memory
!        - these subroutines are called by the program loop
!        - They are listed in the order they are called in the program loop
!    Hamiltonian matrix construction (Where is it?)
!        - Hvblock
!        - This is the master subroutine for construction of the H*PSI matrix-vector
!          product. This is, by far, the most time consuming step
!    Secondary subroutines
!        - All other subroutines are listed in alphabetical order
!    Functions
!        - All functions are listed in alphabetical order
program main
    use progdata,only:niter,nproc,orthog,orthogexact,myid,epsilon,nstates,dimen,maxstor,iiter,idroots,soroots,one,zero,maxdisk
    use filedata,only:outdir
    implicit none
#include "mpif.h"
#include "global.fh"
    logical::status
    integer::i,istat
    real*8::usermem,reqmem

!---------- Initialize ----------
    !Default orthogonalization approach:
    !    Approximate lanczos vector direct multiplies a cheap recursion relation
    !    Exact orthogonalization disabled
    !    Check orthogonality every 100 iterations
    call random_seed()
    call mpi_init(istat)
    call ga_initialize()
    usermem=30d0!Default memory = 30MB
    nproc=ga_nnodes()
    myid=ga_nodeid()
    if(myid==0) call read_cmdline(istat,usermem,outdir)
    call mpi_bcast(usermem,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(outdir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,istat)
    call ga_sync()
    call read_basis()
    call read_constants()
    if(myid==0) then
        call memory_test(usermem,reqmem)
        if(usermem<reqmem) write(*,'(1x,A71)')'Warning: estimated memory requirement is larger than user specification'
        write(*,*)'Scratch files written to: '//adjustl(outdir)
        call print_basis(usermem,reqmem)
    end if
!------------- End --------------

!----------- Run job ------------
    call initialize_elements()
    !Run lanczos diagonalization
    do i=iiter,niter
        call update_lanczos(i)
        call Hv(i)
        call lanczos_step(i)
        if(orthog.and.i.gt.1) then
            call set_omega(i)
            call partial_orthog(i)
        end if
    enddo
    call make_restart(niter)!Output
    if(myid.eq.0) then
      call compute_eigenvalues(niter)
      call print_footer()
    end if
    call GA_SYNC()!Clean up
    if(idroots.gt.0.or.soroots.gt.0) call identify_roots(niter)
    call release_memory()
!------------- End --------------

!---------- Clean up ------------
    call ga_terminate()
    call mpi_finalize(istat)
    if(myid.eq.0) write(*,*)'Mission complete'
!------------- End --------------

end program main