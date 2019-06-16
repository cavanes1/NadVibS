!NadVibS: nonadiabatic vibrational spectrum simulation package
!Originate from NADVIBS.X by Michael Schuurman 2007
!
!Source code layout:
!    Modules: progdata, filedata
!        - contains the equivalent of global variables/common blocks
!        - progdata contains predominantly arithmetic variables for the running of the lanczos program
!        - filedata contains variables used for file I/O
!    program main
!        - NADVIBS program loop
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

module progdata
    integer,parameter::zero=0,one=1,two=2,three=3
    real*8,parameter::zerodp=0d0,ELECMASSU=1822.888484929d0,&!Mass of a proton in atomic units
        MWAU2WAVE=5140.49227189d0,&!Conversion factor for mass weighted a.u. to cm-1, specifically, -> 219474.63*sqrt(1/ELECMASSU)
        AU2WAVE=219474.6313705d0,&!Conversion factor for a.u. to cm-1
        AU2EV=27.211384523d0!Conversion factor for a.u. to eV
    logical::neworigin,&!Shift reference geometry away from initial state
        orthog,&!TRUE if re-orthogonalization is requested
        orthogexact,&!TRUE if exact dot products are to be used to monitor lanczos vector orthogonality
        firststep,&!Global orthogonalization, when performed, must be done on subsequent steps
        restartrun,&!TRUE if current run is restarted from previous results
        savevecs!TRUE if lanczos vectors are to be saved
    integer::dimen,&!The product of the number of basis functions for each mode
        ordr,&!Order of the expansion
        nstblks,&!The number of unique state blocks
        nproc,&!Number of processors employed
        myid,&!Unique name of the current process
        nstates,&!Number of electronic states
        nmodes,&!Number of vibrational modes
        nmax,&!Maximum number of basis function in any given mode
        natoms,&!Number of atoms in molecule
        niter,&!Number of lanczos iterations to be performed
        iiter,&!Initial iteration, =1 if not a restart, =x otherwise
        totstart,totend,iterstart,iterend,&!Timing variable
        maxstor,&!The maximum number of lanczos vectors that can be stored on disk
        northog,&!The number of orthogonalization batches to perform
        chkorthog,&!If using recurrence, compute exact dot products every chkorthog iterations
        idroots,&!How many of the lead contributors to a vibronic level to identify
        soroots,&!How many roots to compute spin orbit parameters for
        nseg,&!Number of segements to divide lanczos vector into
        nirreps,&!Number of irreps in point group symmetry
        q1,q2,q3,&!Global array handles
        ndiag,&!Number of diagonal terms
        nmatel!Length of hamiltonian mat. el. array
    integer,allocatable,dimension(:)::nfunc,&!Array holding the number of harmonic oscillator functions for each node
        noterms,&!Number of unique coupling terms per order
        shft,&!Offset array for indexing basis sets
        loindex,&!Index to start orthogonalizing when doing partial orthogonalization
        roindex,&!Index to stop orthogonalizing lanczos vector when doing P.O.
        istate,&!Initial state of the reference (i.e. anion ground state)
        nalltrms,&!Total number of terms in a processor batch
        nztrms!Number of non-zero terms per order
    integer,allocatable,dimension(:,:)::nzindx,&!Index array of non-zero potential terms
        nzblks!Block locations of non-zero potential terms
    integer::ndblks!Number of blocks that have diagonal contributions, excluding (i,i), which do by definition
    integer,allocatable,dimension(:)::dblks!Identity of those blocks
    integer,dimension(:,:,:),allocatable::basind,bstart,rcshft!hterms: data for doing matrix vector product
    integer,dimension(:,:),allocatable::nelems!Number of non-zero elements in a particular section of Hamiltonian
    integer,dimension(:),allocatable::nsteps!Number of different matrix element types for a potential term
    integer,dimension(:,:),allocatable::basdif,stype,&!Delta values between harmonic oscillator matrix elements
        bmax,&!b value at which counter resets. Does that help?  
        iordr!Array of indicies for location of h.o. matrix elements in homatel
    integer,dimension(:),allocatable::dordr,&!Order of diagonal term i
        dmap,&!Location of term information for diagonal term i
        npirrep!Number of coordinates per irrep
    integer,dimension(:,:),allocatable::vbounds,&!Batch sizes
        cpindx!Coupling block indices
    integer,dimension(2)::lo,hi!Indices of the Lanczos vector a particular process is responsible for
    real*8,dimension(:),allocatable::concoef!dcoef: constant terms for each state block
    real*8,dimension(:,:),allocatable::nzcoef!Array of non-zero potential coefficients
    real*8::bjiconv,&!The degree of convergence before eigenvalues are printed -> 1e-(bijconv)
        epsilon,&!Machine precision of the architecture
        eta,&!Epsilon^(3/4)
        maxdisk,&!Maximum amount of disk available for use (in MB)
        ztoler!Any coefficient less than 1e-(ztoler) is set to 0
    real*8,dimension(:),allocatable::homatel,&!Harmonic oscillator matrix elements
        alpha,beta,&!Hold the alpha and beta constants used to construct tridiagonal matrix
        aomega,&!Value of omega for harmonic oscillator functions
        bomega!Alpha coefficients for determining overlap with reference state
    real*8,dimension(:,:),allocatable::omega!Monitors lanczos vector orthogonality, used only if re-orthog is requested
    real*8,dimension(:),allocatable::statew,&!Initial weights of the vibronic states
        cocoef,&!Constant coefficients
        dvec,&!d vector for determining overlap with reference state
        tmat!T matrix for determining overlap with reference state
    real*8,dimension(:,:),allocatable::scr1,scr2!Scratch buffers
end module progdata

module filedata!Contains info pertaining file I/O
    integer,parameter::BASISFILE=28,&!File containing information parameters required to run program BASISFILE = basis.in
        POTFILE=29,&!File containing potential information CONFILE = nadvibs.in
        OUTFILE=30,&!Majority of output information gets printed to OUTFILE = output.dat
        SCRFILE=31,&!A scratch file for printing out alpha and beta constants as the job progresses
        RESTARTFILE=32,&!File required to initiate a restart of lanczos algorithm
        ROOTINFO=33,&!Created only if idroots > 0
        ARCHIVE=34,&!DRA handle for archive file
        QRESTART=35,&!DRA handle for restart lanczos vectors
        SOINFO=36,&!File containing spinorbit parameters
        EIGVECS=37!File containing eigenvectors of H
    character*100::outdir!Director to output text files to
end module filedata

program main
    use progdata,only:niter,nproc,orthog,orthogexact,myid,epsilon,nstates,dimen,maxstor,iiter,idroots,soroots,one,zero,maxdisk
    use filedata,only:outdir
    implicit none
#include "mpif.h"
#include "global.fh"
    logical::status
    integer::i,istat,nconverged
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
        call memory_test(usermem,reqmem)
        if(myid==0) call print_basis(usermem,reqmem)
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
        write(*,*)'Mission complete'
    !------------- End --------------
end program main

!-------------- Top level subroutines: called by program block ---------------
!Reads the command line arguments.  Currently recognizes two:
!    -m $MEMORY : sets the memory avaiable to the program in MB
!    -o [dir]   : location of the scratch directories 
subroutine read_cmdline(ierr,memory,dir)
    integer,intent(inout)			:: ierr
    real*8, intent(inout)       :: memory
    CHARACTER(len=100), intent(inout)     :: dir
    CHARACTER(2),dimension(2)		:: args = (/'-m','-o'/)
    CHARACTER(len=100)                    :: argbuf
    integer				:: nargs,i,setmem,setout
    outdir = ""
    setmem = 0
    setout = 0
    nargs = iargc()
    do i = 1,nargs
        call getarg(i,argbuf)
        if(setmem.eq.1.and.(trim(argbuf).ne.trim(args(2)))) then
           read(argbuf,*)memory
           setmem = 0
        end if
        if(setout.eq.1.and.(trim(argbuf).ne.trim(args(1)))) then
           dir = trim(adjustl(argbuf))//'/'
           setout = 0
        end if
        if(trim(argbuf).eq.trim(args(1))) setmem = 1
        if(trim(argbuf).eq.trim(args(2))) setout = 1
    end do
end subroutine read_cmdline

subroutine read_basis()!Read the basic input from file basis.in
    use progdata
    use filedata, only: BASISFILE
    implicit none
#include 'mpif.h'
#include 'global.fh'
    integer                      :: i,j,istat
    integer                      :: numut
    if(myid.eq.0) call get_keywords()
    call ga_sync()
    call mpi_bcast(ordr,1,MPI_INTEGER8,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(idroots,1,MPI_INTEGER8,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(soroots,1,MPI_INTEGER8,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(niter,1,MPI_INTEGER8,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(iiter,1,MPI_INTEGER8,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(chkorthog,1,MPI_INTEGER8,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(nseg,1,MPI_INTEGER8,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(nirreps,1,MPI_INTEGER8,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(natoms,1,MPI_INTEGER8,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(nmodes,1,MPI_INTEGER8,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(nstates,1,MPI_INTEGER8,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(dimen,1,MPI_INTEGER8,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(ztoler,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(maxdisk,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(bjiconv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(restartrun,1,MPI_LOGICAL,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(neworigin,1,MPI_LOGICAL,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(savevecs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(orthog,1,MPI_LOGICAL,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(orthogexact,1,MPI_LOGICAL,0,MPI_COMM_WORLD,istat)
    if(myid.ne.0) then
        allocate(statew(nstates))
        allocate(istate(nmodes))
        allocate(nfunc(nmodes))
        allocate(npirrep(nirreps))
    end if
    call mpi_bcast(statew,nstates,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(istate,nmodes,MPI_INTEGER8,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(nfunc,nmodes,MPI_INTEGER8,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(npirrep,nirreps,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,istat)
    call ga_sync()
    call determine_epsilon()!Determine the machine precision for this architecture
    nstblks = numut(nstates,two)
    allocate(noterms(ordr))
    do i = 1,ordr
        noterms(i) = numut(nmodes,i)
    end do
    allocate(cpindx(nstblks,2))
    do i = 1,nstates
        do j = 1,i
            cpindx(i*(i-1)/2+j,1)=i
            cpindx(i*(i-1)/2+j,2)=j
        end do
    end do
end subroutine read_basis

subroutine read_constants()!Read in the potential term information from standard lanczos.in file
    use progdata, only: myid,ordr,nfunc,nmodes,nstates,neworigin,nstblks,nztrms,concoef,aomega,bomega, &
                        dvec,tmat,nzindx,nzblks,nzcoef,nztrms,noterms,ztoler,cpindx,AU2WAVE,zero,zerodp
    use filedata, only: POTFILE,OUTFILE
    implicit none
#include 'mpif.h'
#include 'global.fh'
    logical:: newterm,termchk
    integer:: i,j,k,l,m,n,p,cnt,loc,nblk,ioff
    integer:: istat
    real*8:: dval
    integer,dimension(ordr):: otab,uniq
    character*75:: commentline
    real*8,dimension(:,:,:),allocatable:: POTterms
    integer,dimension(:,:,:),allocatable:: nztemp1
    real*8,dimension(:,:,:),allocatable:: nztemp2
    ioff=0
    do i = 1,ordr
        if(noterms(i).gt.ioff)ioff=noterms(i)
    end do
    allocate(nztrms(ordr))
    allocate(concoef(nstblks))
    allocate(POTterms(ioff,ordr,nstblks))
    allocate(nztemp1(ordr,ordr*ioff,2*ordr+nstblks+1))
    allocate(nztemp2(ordr,ordr*ioff,nstblks))
    call setintarray(nztrms,ordr,zero)
    allocate(aomega(nmodes))
    if(neworigin) then
        allocate(bomega(nmodes))
        allocate(dvec(nmodes))
        allocate(tmat(nmodes**2))
    end if
    if(myid.eq.0) then!Only the root process will be able to see this file...
        open(POTFILE,file='nadvibs.in',access='sequential',form='formatted')
        read(POTFILE,fmt='(a75)') commentline
        read(POTFILE,*)(aomega(i),i=1,nmodes)
        do i = 1,nstblks
            read(unit=POTFILE,fmt='(a75)') commentline
            read(unit=POTFILE,fmt=*)concoef(i)
            do j = 1,ordr
                read(unit=POTFILE,fmt='(a75)') commentline
                read(unit=POTFILE,fmt=*)(POTterms(k,j,i),k=1,noterms(j))
            end do
        end do
        if(neworigin) then
            read(unit=POTFILE,fmt='(a75)') commentline
            read(unit=POTFILE,fmt=*)(bomega(i),i=1,nmodes)
            read(unit=POTFILE,fmt='(a75)') commentline
            read(unit=POTFILE,fmt=*)(dvec(i),i=1,nmodes)
            read(unit=POTFILE,fmt='(a75)') commentline
            read(unit=POTFILE,fmt=*)(tmat(i),i=1,nmodes**2)
        end if
        close(POTFILE)
    end if
    call ga_sync()
    do i = 1,nstblks
        do j = 1,ordr
            call mpi_bcast(POTterms(1:noterms(ordr),j,i),noterms(ordr),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,istat)
        end do
    end do
    call mpi_bcast(concoef,nstblks,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,istat)
    call mpi_bcast(aomega,nmodes,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,istat)
    if(neworigin) then
        call mpi_bcast(bomega,nmodes,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,istat)
        call mpi_bcast(dvec,nmodes,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,istat)
        call mpi_bcast(tmat,nmodes**2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,istat)
    end if
    do i = 1,nstblks
        do j = 1,ordr
            if(myid.eq.0)print *,'BLK=',i,' ORDR=',j
            m = 0
            call setintarray(otab,ordr,int(1))
    !Michael Schuurman's way of generating otab
            otab(j) = 0
            do k = 1,noterms(j)
                l = j
                do
                    if(l.eq.1) exit
                    if(otab(l).lt.otab(l-1)) exit
                    otab(l) = 1
                    l = l - 1
                end do
                otab(l) = otab(l) + 1
    !My preference is to use pseudo nmodes+1 counter satisfying former digit >= latter digit,
    !corresponding to the direct sum of an ordr-th order tensor's 1st dimension vector
            !otab(1)=0
            !do k = 1,noterms(j)
            !    otab(1)=otab(1)+1!Add 1 to the 1st digit
            !    do l=1,j-1!Carry to latter digits
            !        if(otab(l)>nmodes) then
            !            otab(l)=1
            !            otab(l+1)=otab(l+1)+1
            !        end if
            !    end do
            !    do l=j-1,1,-1!Modify to satisfy former digit >= latter digit
            !        if(otab(l)<otab(l+1)) otab(l)=otab(l+1)
            !    end do
    !Fuck, shitslide: potential term counting seems to suit only his specific definition
    !End of otab generation
                if(myid.eq.0)print *,'otab =',otab(1:j)
                m = m + 1
                call union(j,otab,cnt,uniq)
                if(abs(POTterms(m,j,i)).ge.ztoler)then
                    newterm = .true.
                    do n = 1,nztrms(cnt)
                        termchk = .true. 
                        do p = 1,cnt
                            if(uniq(p).ne.nztemp1(cnt,n,2*p-1).or.count(otab(1:j).eq.uniq(p)).ne.nztemp1(cnt,n,2*p))termchk=.false.
                        end do
                        if(termchk)then
                            newterm=.false.
                            loc = n
                            EXIT
                        end if
                    end do
                    if(newterm)then
                        nztrms(cnt) = nztrms(cnt) + 1
                        do n = 1,cnt
                            nztemp1(cnt,nztrms(cnt),2*n-1) = uniq(n)
                            nztemp1(cnt,nztrms(cnt),2*n)   = count(otab(1:j).eq.uniq(n))
                        end do
                        nztemp1(cnt,nztrms(cnt),2*cnt+1) = 1
                        nztemp1(cnt,nztrms(cnt),2*cnt+2) = i
                        nztemp2(cnt,nztrms(cnt),1) = POTterms(m,j,i)
                    else
                        nblk = nztemp1(cnt,loc,2*cnt+1)+1
                        nztemp1(cnt,loc,2*cnt+1) = nblk
                        nztemp1(cnt,loc,2*cnt+1+nblk) = i
                        nztemp2(cnt,loc,nblk) = POTterms(m,j,i)
                    end if
                end if
            end do
        end do
    end do
    cnt = sum(nztrms)
    allocate(nzindx(2*ordr,cnt))
    allocate(nzblks(1+2*nstates*nstates,cnt))
    allocate(nzcoef(nstates*nstates,cnt))
    call setintmatrix(nzindx,2*ordr,cnt,zero)
    call setintmatrix(nzblks,1+2*nstates*nstates,cnt,zero)
    call setmatrix(nzcoef,nstates*nstates,cnt,zerodp)
    ioff = 0
    do i = 1,ordr
        do j = 1,nztrms(i)
            dval = 1.
            do k = 1,i
                nzindx(2*k-1,ioff+j) = nztemp1(i,j,2*k-1)!Store the unique indices
                nzindx(2*k,ioff+j)   = nztemp1(i,j,2*k)!and how many times they occur
                dval = dval/Sqrt(aomega(nztemp1(i,j,2*k-1))**nztemp1(i,j,2*k))
            end do
            cnt=0!Store the number of blocks this term is non-zero
            do k = 1,nztemp1(i,j,2*i+1)
                cnt = cnt + 1
                nzblks(2*cnt,ioff+j)   = cpindx(nztemp1(i,j,2*i+1+k),1)
                nzblks(2*cnt+1,ioff+j) = cpindx(nztemp1(i,j,2*i+1+k),2)
                nzcoef(cnt,ioff+j)     = nztemp2(i,j,k)*dval
                if(nzblks(2*cnt,ioff+j).ne.nzblks(2*cnt+1,ioff+j))then
                    cnt = cnt + 1
                    nzblks(2*cnt,ioff+j)   = cpindx(nztemp1(i,j,2*i+1+k),2)
                    nzblks(2*cnt+1,ioff+j) = cpindx(nztemp1(i,j,2*i+1+k),1)
                    nzcoef(cnt,ioff+j)     = nztemp2(i,j,k)*dval
                end if
                nzblks(1,ioff+j) = cnt
            end do
        end do
        ioff = ioff + nztrms(i)
    end do
    print *,'NZTERMS: ',nztrms
    deallocate(POTterms)
    deallocate(nztemp1)
    deallocate(nztemp2)
end subroutine read_constants

subroutine print_basis(umem,rmem)!Print a summary of the basis set information from basis.in
    use progdata, only: ordr,natoms,nstates,nmodes,niter,dimen,orthog,orthogexact,maxstor,nfunc,ztoler, &
                        chkorthog,nproc,epsilon,maxdisk,statew,restartrun,iiter,soroots,bjiconv,nseg,   &
                        nirreps,npirrep,aomega,bomega,istate,nstblks,AU2WAVE,nzindx,nztrms,nzblks,zero,neworigin
    use filedata, only: outdir,OUTFILE
    implicit none
    real*8, intent(in)    :: umem,rmem
    integer                         :: i,j,k,l,pordr,ioff
    integer,dimension(ordr)         :: otab
    integer,dimension(nstblks,ordr) :: ntot
    real*8,dimension(nmodes) :: dpvec
    open(unit=OUTFILE,file='output.dat',status='replace')
    rewind(unit=OUTFILE)  
    write(unit=OUTFILE,fmt='(a)')   'PNADVIBS.X'
    write(unit=OUTFILE,fmt='(a)')   'VIBRONIC ENERGY LEVEL PROGRAM EMPLOYING A LANCZOS SOLVER' 
    write(unit=OUTFILE,fmt='(a)')   '--------------------------------------------------------'
    write(unit=OUTFILE,fmt='(a)')   'Michael Schuurman'
    write(unit=OUTFILE,fmt='(a)')   ' '
    write(unit=OUTFILE,fmt='(a)')	  ' '
    write(unit=OUTFILE,fmt=1007)adjustl(outdir)
    write(unit=OUTFILE,fmt='(a)')   ' '
    write(unit=OUTFILE,fmt='(a)')   'Input parameters read from BASIS.IN:'
    write(unit=OUTFILE,fmt='(a)')   '------------------------------------'
    write(unit=OUTFILE,fmt='(A48,I10)')'  Number of atoms:                              ',natoms
    write(unit=OUTFILE,fmt='(A48,I10)')'  Number of electronic states:                  ',nstates
    write(unit=OUTFILE,fmt='(A48,I10)')'  Order of the expansion:                       ',ordr
    write(unit=OUTFILE,fmt='(A48,I10)')'  Number of vibrational modes:                  ',nmodes
    write(unit=OUTFILE,fmt='(A48,I10)')'  Number of irreducible representations         ',nirreps
    write(unit=OUTFILE,fmt=1003)npirrep
    write(unit=OUTFILE,fmt='(A48,I10)')'  Number of Lanczos iterations:                 ',niter
    write(unit=OUTFILE,fmt='(A48,I10)')'  Number of processors used in execution:       ',nproc
    write(unit=OUTFILE,fmt='(A48,ES10.0)')'  Zero tolerance for potential coeffs.:         ',ztoler
    write(unit=OUTFILE,fmt='(A53,I5)')'  Compute spin-orbit parameters for first N roots:   ',soroots
    if(restartrun)write(unit=OUTFILE,fmt='(A48,I10)')'  Restarting Lanczos Procedure on step:         ',iiter
    write(unit=OUTFILE,fmt='(a53,l5)')'  Perform lanczos vector re-orthogonalization:       ',orthog
    if(orthog) then
        write(unit=OUTFILE,fmt='(a53,l5)')'    - exact dot products for vector orthog.:         ',orthogexact
        if(.not.orthogexact) write(unit=OUTFILE,fmt='(A48,I10)')'    - compute exact orthog. every X iterations: ',chkorthog
        write(unit=OUTFILE,fmt='(a50,f8.1)')'    - max. amount of disk available for storage:  ',maxdisk
        write(unit=OUTFILE,fmt='(A48,I10)')'    - max. number of lanczos vectors to store:  ',maxstor
    end if
    do i = 1,nstates
        write(unit=OUTFILE,fmt=1000)i,statew(i)
    end do
    write(unit=OUTFILE,fmt='(a50,es8.0)')'  Convergence criteria for eigenvalues (bji):     ',bjiconv
    write(unit=OUTFILE,fmt=1001)rmem
    write(unit=OUTFILE,fmt=1002)umem
    if(rmem>umem) write(*,'(1x,A71)')'Warning: estimated memory requirement is larger than user specification'
    write(unit=OUTFILE,fmt='(A48,I10)')'  Number of Segements per Lanczos Vector:        ',nseg
    write(unit=OUTFILE,fmt='(A42,I16)')'  Dimensionality of a single State Vector:',dimen
    write(unit=OUTFILE,fmt='(A42,I16)')'  Total Dimensionality of H matrix:       ',nstates*dimen
    write(unit=OUTFILE,fmt='(a42,es16.6)')'  Machine precision for this architecture:',epsilon
    ioff = 0
    call setintmatrix(ntot,nstblks,ordr,zero)
    do i = 1,ordr
        otab(i) = i
        do j = 1,nztrms(i)
            pordr = 0
            do k = 1,i
                pordr = pordr + nzindx(2*k,ioff+j)
            end do
            print *,'ordr=',i,' j=',j,' pordr=',pordr
            k = 1
            do
                if(k.gt.nzblks(1,ioff+j)) exit
                l = nzblks(2*k,ioff+j)*(nzblks(2*k,ioff+j)-1)/2+nzblks(2*k+1,ioff+j)
                ntot(l,pordr) = ntot(l,pordr) + 1
                if(nzblks(2*k,ioff+j).ne.nzblks(2*k+1,ioff+j))k = k + 1
                k = k + 1
            end do
        end do
        ioff = ioff + nztrms(i)
    end do
    write(unit=OUTFILE,fmt='(a)')''
    write(unit=OUTFILE,fmt='(a)')'  Number of Potential Terms per Block -----------'
    write(unit=OUTFILE,fmt='(a)')''
    write(unit=OUTFILE,fmt=1006)otab(1:ordr)
    do i = 1,nstblks
        write(unit=OUTFILE,fmt=1005)i,ntot(i,1:ordr)
    end do
    write(unit=OUTFILE,fmt='(a)')     ''
    write(unit=OUTFILE,fmt='(a)')     '  Basis set specification:'
    write(unit=OUTFILE,fmt='(a)')     '    mode  # func.    omega'
    do i=1,nmodes
        write(unit=OUTFILE,fmt=1004)i,nfunc(i),aomega(i)*AU2WAVE
    end do
    write(unit=OUTFILE,fmt='(a)') ''
    if(neworigin)then
        do i = 1,nmodes
            dpvec(i) = bomega(i)*AU2WAVE
        end do
    else
        do i = 1,nmodes
            dpvec(i) = aomega(i)*AU2WAVE
        end do
    end if
    !Temporary -- only one mode may be excited in istate
        j = 0
        do i = 1,nmodes
            if(istate(i).gt.0.and.j.eq.0)then
                j = 1
            else
                istate(i) = 0
            end if
        end do
    write(unit=OUTFILE,fmt='(a)')     ''
    write(unit=OUTFILE,fmt='(a)')     '  Initial state specification:'
    write(unit=OUTFILE,fmt='(a)')     '    mode  # quanta   omega'
    do i=1,nmodes
        write(unit=OUTFILE,fmt=1004)i,istate(i),dpvec(i)
    end do
    write(unit=OUTFILE,fmt='(a)') ''
    1000 format('  Initial weight of reference state ',i2,':              ',f5.2)
    1001 format('  Total Memory Required (GA+NADVIBS) (MB):',f16.0)
    1002 format('  Memory Available (MB):                  ',f16.0)
    1003 format('  Number of modes per irrep:                         ',8(i3))
    1004 format(4x,i3,i7,f12.3)
    1005 format('  Block ',i3,':',8(I6,'   '))
    1006 format('  ORDER:    ',8('    ',i2,'   '))
    1007 format('  Scratch files written to (blank if current directory): ',a100)
end subroutine print_basis

subroutine initialize_elements()
    use progdata, only: myid,nproc,nmodes,nstates,dimen,nfunc,q1,q2,q3,scr2,alpha,beta,totstart,totend,   &
                        aomega,iiter,niter,shft,orthog,epsilon,omega,loindex,roindex,zerodp,zero,one,two, &
                        firststep,restartrun,maxdisk,nseg,vbounds,savevecs,ordr,cpindx,nstblks,nztrms,    &
                        nzindx,stype,nsteps,nelems,nmax,homatel,iordr,nzblks,basdif,lo,hi,concoef,nmatel, &
                        ndiag,dmap,ndblks,dblks
    use filedata, only: outdir,OUTFILE,SCRFILE,QRESTART,ARCHIVE
    implicit none
#include 'mpif.h'
#include 'global.fh'
#include 'mafdecls.fh'
    integer		             :: i,j,k,l,m,n,ioff1,ioff2
    integer                            :: nload,seglen,stack,heap,global,nall
    integer                            :: istat,tcount,tmx,trm1,trm2,rlen
    integer,dimension(2)               :: qlo,qhi
    integer,dimension(nstblks,ordr)    :: ntcom,ntdet,ntsum
    integer,dimension(nstates,nstates) :: dchk
    integer,dimension(:),allocatable   :: ntot
    real*8                   :: mb2b,gauss_random,hoelem
    logical:: vecload,status,isdiag
    character*4:: procindex
    !Initialize and allocate local memory arrays
        call ga_sync()
        call system_clock(totstart,tcount,tmx)
        allocate(vbounds(nproc*nseg,2))
        allocate(shft(nmodes))
        allocate(alpha(niter))
        allocate(beta(0:niter)) 
        allocate(omega(2,0:niter))
        allocate(loindex(Ceiling(1.*niter/100.)))
        allocate(roindex(Ceiling(1.*niter/100.)))
        call setintmatrix(vbounds,nproc*nseg,2,zero)
        call setintmatrix(ntcom,nstblks,ordr,zero)
        call setintmatrix(ntdet,nstblks,ordr,zero)
        call setintmatrix(ntsum,nstblks,ordr,zero)
        call setintarray(shft,nmodes,one)
        nmax = 1
        do i = 1,nmodes
            if(nfunc(i).gt.nmax)nmax = nfunc(i)
            do j = i+1,nmodes
                shft(i) = shft(i)*nfunc(j)
            end do
        end do
        firststep = .true.
    !Initialize Global Arrays
        mb2b = 1024.*1024.
        stack  = Ceiling(mb2b)
        heap   = Ceiling(mb2b)
        global = Ceiling(3.*dimen*nstates/nproc)+100
        call ga_sync()
        if(ga_uses_ma()) then
            status = ma_init(MT_F_DBL,stack,heap+global)
            if(.not.status)print *,'MA_INIT ERROR.......'
        else
            call ga_set_memory_limit(ma_sizeof(MT_F_DBL,global,MT_F_BYTE))
        end if
        status = GA_CREATE(MT_F_DBL,dimen,nstates,"Q1",zero,nstates,q1)
        status = GA_CREATE(MT_F_DBL,dimen,nstates,"Q2",zero,nstates,q2)
        status = GA_CREATE(MT_F_DBL,dimen,nstates,"Q3",zero,nstates,q3)
        if(myid.eq.0)write(unit=OUTFILE,fmt=1000)
        do i = 1,nproc
            if(myid.eq.0)write(unit=OUTFILE,fmt=1001)i
            call NGA_DISTRIBUTION(q1,i-1,qlo,qhi)
            if(myid.eq.0)write(unit=OUTFILE,fmt=1002)qlo(1),qhi(1),qlo(2),qhi(2),(qhi(1)-qlo(1)+1)*qhi(2)
            call NGA_DISTRIBUTION(q2,i-1,qlo,qhi)
            if(myid.eq.0)write(unit=OUTFILE,fmt=1003)qlo(1),qhi(1),qlo(2),qhi(2),(qhi(1)-qlo(1)+1)*qhi(2)
            call NGA_DISTRIBUTION(q3,i-1,qlo,qhi)
            if(myid.eq.0)write(unit=OUTFILE,fmt=1004)qlo(1),qhi(1),qlo(2),qhi(2),(qhi(1)-qlo(1)+1)*qhi(2)
            seglen = nint(1.*(qhi(1) - qlo(1) + 1.)/nseg)
            do j = 1,nseg
                vbounds((i-1)*nseg+j,1) = qlo(1) + (j-1)*seglen 
                vbounds((i-1)*nseg+j,2) = qlo(1) + j*seglen - 1
            end do
            vbounds(i*nseg,2) = qhi(1)
        end do
        if(myid.eq.0) then
            write(unit=OUTFILE,fmt=1005)
            do i = 1,nproc
                write(unit=OUTFILE,fmt=1001)i
                do j = 1,nseg
                    write(unit=OUTFILE,fmt=1006)j,vbounds((i-1)*nseg+j,1),vbounds((i-1)*nseg+j,2),  &
                                                  vbounds((i-1)*nseg+j,2)-vbounds((i-1)*nseg+j,1)+1
                end do
            end do
            write(unit=OUTFILE,fmt=1007)
        end if
        lo = (/ vbounds(myid*nseg+1,1), one /)
        hi = (/ vbounds(myid*nseg+nseg,2), nstates /)
    !Setup the Hamiltonian
        !Allocate the buffer that can hold the local portion of a global array vector
            rlen = hi(1) - lo(1) + 1
            allocate(scr2(rlen,nstates))
            call ga_sync()
            call compute_allIndex()
        !Count the number of terms -- make sure predicted number agrees with the observed number of terms
            ioff1 = 0
            ioff2 = 0
            do i = 1,ordr
                do j = 1,nztrms(i)
                    trm1 = dimen
                    do l = 1,i
                       trm1 = trm1/nfunc(nzindx(2*l-1,ioff1+j))
                    end do
                    do k = 1,nsteps(ioff1+j) 
                        trm2 = 1
                        do l = 1,i
                            trm2 = trm2*(nfunc(nzindx(2*l-1,ioff1+j))-basdif(l,ioff2+k))
                        end do
                        l = 1
                        do
                            if(l.gt.nzblks(1,ioff1+j)) exit
                            m = nzblks(2*l,ioff1+j)*(nzblks(2*l,ioff1+j)-1)/2+nzblks(2*l+1,ioff1+j)
                            ntcom(m,i) = ntcom(m,i) + trm1*trm2
                            do n = 1,nproc*nseg
                             ntdet(m,i) = ntdet(m,i) + nelems(ioff2+k,n)
                            end do
                            if(nzblks(2*l,ioff1+j).ne.nzblks(2*l+1,ioff1+j))l = l + 1
                            l = l + 1
                        end do
                    end do
                    isdiag = .true.
                    do n = 1,i 
                        if(mod(nzindx(2*n,ioff1+j),2).ne.0)isdiag=.false.
                    end do
                    if(isdiag)then
                        l = 1
                        do 
                            if(l.gt.nzblks(1,ioff1+j)) exit
                            m = nzblks(2*l,ioff1+j)*(nzblks(2*l,ioff1+j)-1)/2+nzblks(2*l+1,ioff1+j)
                            ntdet(m,i) = ntdet(m,i) + hi(1) - lo(1) + 1
                            if(nzblks(2*l,ioff1+j).ne.nzblks(2*l+1,ioff1+j)) l = l + 1
                            l = l + 1
                        end do
                    end if
                    ioff2 = ioff2 + nsteps(ioff1+j)
                end do
                ioff1 = ioff1 + nztrms(i)
            end do
            call ga_sync()
            do i = 1,ordr
                call MPI_ALLREDUCE(ntdet(1:nstblks,i),ntsum(1:nstblks,i),nstblks,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,istat)
            end do
            if(myid.eq.0)then!Print out the results of the term counting
                write(OUTFILE,1008)
                trm1 = 0
                do i = 1,nstblks
                    do j = 1,ordr
                        trm1 = trm1 + ntsum(i,j)
                        write(OUTFILE,1009)j,i,ntsum(i,j)
                        if(ntsum(i,j).ne.ntcom(i,j)) write(OUTFILE,1010)j,i,ntcom(i,j),ntsum(i,j)
                    end do
                    write(OUTFILE,1011)
                end do
                write(OUTFILE,1007)
                write(OUTFILE,1012)trm1
            end if
            call ga_sync()
            nmatel = 0
            do i = 1,ordr
                j = mod(i,2)
                do
                    if(j.gt.ordr)EXIT
                    if(nmax-j.gt.0)nmatel = nmatel + nmax-j
                    j = j + 2
                end do
            end do
        !Construct array holding requisite h.o. matrix elements
            allocate(iordr(int(ordr/2.)+1,ordr))
            allocate(homatel(nmatel))
            ioff1 = 0
            do i = 1,ordr
                j = mod(i,2)
                k = 1
                do 
                    if(j.gt.ordr) exit
                    iordr(k,i) = ioff1
                    do l = 1,nmax-j
                        ioff1 = ioff1 + 1
                        homatel(ioff1) = hoelem(i,l,l+j)
                    end do
                    k = k + 1
                    j = j + 2
                end do
            end do
        call setintmatrix(dchk,nstates,nstates,zero)!Construct list of state blocks involving diagonal terms
        do i = 1,nstblks!Constant terms
            if(abs(concoef(i)).gt.1e-16)then
                dchk(cpindx(i,1),cpindx(i,2)) = 1
                dchk(cpindx(i,2),cpindx(i,1)) = 1
            end if
        end do
        !Potential terms
            do i = 1,ndiag
                do j = 1,nzblks(1,dmap(i))
                    k = nzblks(2*j,dmap(i))
                    l = nzblks(2*j+1,dmap(i))
                    dchk(k,l) = 1
                end do
            end do
            ndblks = 0
            allocate(dblks(2*nstates*nstates))
        do i = 1,nstates!Store the blocks that have non-zero diagonal terms
            do j = 1,i-1
                if(dchk(i,j).ne.0)then
                    ndblks = ndblks + 1
                    dblks(2*ndblks-1) = i
                    dblks(2*ndblks)   = j
                end if
            end do
        end do
    !Set restart dependent parameters
        write(procindex,'(i4)')myid
        procindex = adjustl(procindex)
        if(restartrun) then!If restarting from previous computation
            if(myid.eq.0)call load_restartinfo(vecload)
            nload = iiter - 1
            call ga_sync()
            call mpi_bcast(vecload,1,MPI_LOGICAL,0,MPI_COMM_WORLD,istat)
            call mpi_bcast(alpha,nload,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,istat)
            call mpi_bcast(beta(0:nload),nload+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,istat)
            if(orthog) then
                do i = 1,2
                    call mpi_bcast(omega(i,0:nload),nload+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,istat)
                end do
            end if
            if(savevecs) then
                open(unit=ARCHIVE,file=trim(adjustl(outdir))//'nadvibs.34.'//trim(procindex),access='direct', &
                     status='old',form='unformatted',recl=8*rlen)
                rewind(ARCHIVE)
            end if
            if(vecload) then
                open(unit=QRESTART,file=trim(adjustl(outdir))//'nadvibs.35.'//trim(procindex),access='direct', &
                     status='old',form='unformatted',recl=8*rlen)
                rewind(QRESTART)
                call read_ga(q1,QRESTART,1)
                call read_ga(q2,QRESTART,2)
                call read_ga(q3,QRESTART,3)
                close(QRESTART,status='delete')
            else
                call read_ga(q1,ARCHIVE,nload-2)
                call read_ga(q2,ARCHIVE,nload-1)
                call read_ga(q3,ARCHIVE,nload)
                call GA_SCALE(q3,beta(nload))
            end if
            call ga_sync()
        else
            omega(1,0) = 0.
            omega(1,1) = 1.
            omega(2,0) = 0.
            omega(2,1) = gauss_random(zerodp,dble(0.6))*epsilon*nstates*dimen
            omega(2,2) = 1.
            if(savevecs)then
                open(ARCHIVE,file=trim(adjustl(outdir))//'nadvibs.34.'//trim(procindex),access='direct', &
                     form='unformatted',status='replace',recl=8*rlen)
                rewind(ARCHIVE)
            end if
            call GA_ZERO(q2)
        end if
    call system_clock(totend,tcount,tmx)
    if(myid.eq.0) then
        open(SCRFILE,file='scratch.dat',status='replace')
        write(OUTFILE,'(a)') ' '
        write(OUTFILE,'(a)') '  Initialization complete.'
        write(OUTFILE,1013)1.*(totend-totstart)/tcount
    end if
    !Generate starting vector, if not restarting
        if(.not.restartrun) then
            call ga_sync()
            call generate_initialvec()
            beta(0) = Sqrt(GA_DDOT(q3,q3))
        end if
        if(myid.eq.0) write(OUTFILE,1014)
        call ga_sync()
        call system_clock(totstart,tcount,tmx)
    1020 format('x=',i3,' nelems(1,x)=',i6,' nelems(2,x)=',i6)
    1000 format(/,3x,'Distribution of Lanczos vectors -------------------')
    1001 format(/,3x,'PROCESS',I4)
    1002 format(3x,'Q1:',I15,' ->',I15,',',I15,' ->',I15,', TOTAL:',I19)
    1003 format(3x,'Q2:',I15,' ->',I15,',',I15,' ->',I15,', TOTAL:',I19)
    1004 format(3x,'Q3:',I15,' ->',I15,',',I15,' ->',I15,', TOTAL:',I19)
    1005 format(/,3x,'Segmentation Scheme -------------------------------')
    1006 format(3x,'Seg',I4,':',I20,' ->',I20,' LENGTH:',I20)
    1007 format(/,3x,'---------------------------------------------------')
    1008 format(' ------- TERM COUNTING FOR VIBRONIC HAMILTONIAN --------')
    1009 format('  TOT. NUM.',I4,'-INDEX TERMS IN BLK',I4,':',I20)
    1010 format('  !!!!!!!!!!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!!!!',/, &
                  '  DISAGREEMENT IN # OF ',i3,'-INDEX TERMS IN BLK ',i3,': predict=',i11,' != actual=',i11,/)
    1011 format('')
    1012 format('  TOTAL NUMBER OF TERMS IN ALL BLOCKS:',I20)
    1013 format('  Time Required:',F15.3,' secs.')
    1014 format(/,'  >>>>>>>> Running Lanczos Iterations <<<<<<<<< '/)
end subroutine initialize_elements

!The Lanczos alogrithm employed here has the following form:
!    r[1] != 0, q[0] = 0, beta[1] = ||r[1]||       -- intialize_elements
!    do j = 1, niter                               -- Main program loop
!        q2 = qnew[j]/beta[j]                      -- update_lanczos
!        u[j] = q2*A*q2 - beta[j]*q1               -- Hv and lanczos_step
!        alpha[j] = u[j] * q[j]                    -- lanczos_step
!        r[j+1] = u[j] - alpha[j]q[j]              -- lanczos_step
!        beta[j+1] = ||r[j+1]||                    -- lanczos_step
!        if(orthog)orthogonalize lanczos vectors   -- orthog_vectors
!    end do                                        -- Main program loop  
subroutine update_lanczos(iter)
    use progdata, only: dimen,beta,nstates,q1,q2,q3,myid,nproc,iterstart,iterend,one,savevecs
    use filedata, only: OUTFILE,ARCHIVE
    implicit none
#include 'global.fh'
#include 'mafdecls.fh'
    integer,intent(in)	:: iter
    integer               :: tcount,tmx,istat,ireq
    call system_clock(iterstart,tcount,tmx)
    call GA_COPY(q2,q1)
    call GA_SCALE(q3,1./beta(iter-1))
    call GA_COPY(q3,q2)
    call ga_sync()
    if(savevecs) then
        call write_ga(q2,ARCHIVE,iter)
        call ga_sync()
    end if
    call system_clock(iterend,tcount,tmx)
    if(myid.eq.0) then
        call write_scrfile(iter)
        write(unit=OUTFILE,fmt=1000)iter
        write(unit=OUTFILE,fmt=1001)1.*(iterend-iterstart)/tcount 
    end if
    1000 format(3x,'ITERATION: ',i5)
    1001 format(3x,' - Update Lanczos vector:                ',f14.3,' secs.')
end subroutine update_lanczos

subroutine Hv(iter)
    use progdata, only: myid,ordr,nproc,niter,dimen,nstates,nstblks,cpindx,nmodes,nfunc,zero,one,two,zerodp,nalltrms,  &
                        q2,q3,vbounds,scr1,scr2,nseg,concoef,shft,aomega,iterend,iordr,lo,hi,nmatel,                   &
                        stype,nztrms,basind,basdif,bstart,rcshft,nelems,nsteps,bmax,nzindx,nzblks,nzcoef,homatel,      &
                        ndiag,dordr,dmap,ndblks,dblks
    use filedata, only: OUTFILE
    implicit none
#include 'global.fh'
    integer,intent(in)                                       :: iter
    integer                                                  :: i,j,k,l,m,n,ioff1,ioff2
    integer                                                  :: ibat,nbat,rindx,cindx
    integer                                                  :: tend,tcount,tmx
    integer,dimension(1)                                     :: ld,qld
    integer,dimension(2)                                     :: qlo,qhi
    integer,dimension(nmodes)                                :: barray
    integer,dimension(ordr)                                  :: inds,icnt,bval
    integer,dimension(ordr+1)                                :: bcnt
    real*8                                         :: dkin,d1,dscale
    real*8,dimension(ordr)                         :: matel
    real*8,dimension(nstates,nstates)              :: dvals
    dscale = 1.
    qld(1) = 0
    ld(1) = hi(1) - lo(1) + 1
    nbat = nseg*nproc
    ibat = nseg*myid
    call setmatrix(scr2,ld(1),nstates,zerodp)
    do i = 1,nbat!Loop over the number of processors
        ibat = ibat + 1
        if(ibat.gt.nbat) ibat = ibat - nbat
        if(nalltrms(ibat).gt.0) then
            qlo = (/ vbounds(ibat,1),1 /)
            qhi = (/ vbounds(ibat,2),nstates /)
            if((qhi(1)-qlo(1)+1).ne.qld(1)) then
                qld(1) = qhi(1) - qlo(1) + 1
                if(allocated(scr1)) deallocate(scr1)
                allocate(scr1(qld(1),nstates))
            end if
            call NGA_GET(q2,qlo,qhi,scr1,qld)
            !If a 'diagonal/diagonal' block, do the constant, kinetic energy terms
            !Also, the <m|x^a|n> terms, where m==n for a=2,4,6, are done here as well
            if(ibat.ge.(nseg*myid+1).and.ibat.le.(nseg*myid+nseg)) then 
                call get_barray(qlo(1),barray)
                barray(nmodes) = barray(nmodes) - 1
                rindx = qlo(1) - lo(1)
                dkin = 0.
                do j = 1,nmodes    
                    dkin = dkin + (barray(j)-0.5)*aomega(j)
                end do 
                do j = 1,qld(1)
                    !Kinetic Energy
                    rindx = rindx + 1
                    k = nmodes
                    do
                        if(barray(k).lt.nfunc(k)) exit
                        dkin = dkin - (barray(k)-1)*aomega(k)
                        barray(k) = 1
                        k = k - 1
                    end do
                    barray(k) = barray(k) + 1
                    dkin = dkin + aomega(k)
                    !Constant TErm
                    do k = 1,nstblks
                        dvals(cpindx(k,1),cpindx(k,2)) = concoef(k)
                        dvals(cpindx(k,2),cpindx(k,1)) = concoef(k)
                    end do
                    !Diagonal terms       
                    do k = 1,ndiag
                        d1 = 1.
                        do l = 1,dordr(k)
                            d1 = d1*homatel(iordr(1,nzindx(2*l,dmap(k)))+barray(nzindx(2*l-1,dmap(k))))
                        end do
                        do l = 1,nzblks(1,dmap(k))
                            m = nzblks(2*l,dmap(k))
                            n = nzblks(2*l+1,dmap(k))
                            dvals(m,n) = dvals(m,n) + nzcoef(l,dmap(k))*d1
                        end do
                    end do
                    do k = 1,nstates
                        scr2(rindx,k) = scr2(rindx,k) + (dvals(k,k)+dkin)*scr1(j,k)
                    end do
                    do k = 1,ndblks
                        scr2(rindx,dblks(2*k-1)) = scr2(rindx,dblks(2*k-1)) + dvals(dblks(2*k-1),dblks(2*k))*scr1(j,dblks(2*k))
                    end do
                end do
            end if
            !Now do all remaining potential terms, <m|x^a|n> for m!=n, a=1,2,3,4,etc
            ioff1 = 0
            ioff2 = 0
            do j = 1,ordr
                do k = 1,nztrms(j)
                    do l = 1,j
                        inds(l) = nzindx(2*l-1,ioff1+k)
                        icnt(l) = nzindx(2*l,ioff1+k)
                    end do
                    do l = 1,nsteps(ioff1+k)
                        if(nelems(ioff2+l,ibat).gt.0) then
                            d1 = 1.
                            do m = 1,j
                                bval(m)  = basind(m,ioff2+l,ibat)
                                bcnt(m)  = bstart(m,ioff2+l,ibat)
                                n = iordr(stype(m,ioff2+l),icnt(m))+bval(m)-basdif(m,ioff2+l)
                                !if(n.lt.1.or.n.gt.nmatel)n=1
                                matel(m) = homatel(n)
                                d1 = d1*matel(m)
                            end do
                            cindx = rcshft(1,ioff2+l,ibat)
                            rindx = rcshft(2,ioff2+l,ibat)
                            bcnt(j+1) = 0
                        end if
                        do m = 1,nelems(ioff2+l,ibat)
                            do n = 1,nzblks(1,ioff1+k)
                                scr2(rindx,nzblks(2*n,ioff1+k)) = scr2(rindx,nzblks(2*n,ioff1+k))&
                                    +scr1(cindx,nzblks(2*n+1,ioff1+k))*nzcoef(n,ioff1+k)*d1
                            end do
                            cindx = cindx + 1
                            rindx = rindx + 1
                            bcnt(1) = bcnt(1) + 1
                            do n = 1,j
                                if(bcnt(n).le.bmax(n,ioff1+k)) exit
                                bcnt(n) = 1
                                bval(n) = bval(n) + 1
                                if(bval(n).gt.nfunc(inds(n)))then
                                    bcnt(n+1) = bcnt(n+1) + 1
                                    bval(n)   = 1 + basdif(n,ioff2+l)
                                    cindx     = cindx + basdif(n,ioff2+l)*shft(inds(n))
                                    rindx     = rindx + basdif(n,ioff2+l)*shft(inds(n))
                                end if
                                d1 = d1/matel(n)
                                matel(n) = homatel(iordr(stype(n,ioff2+l),icnt(n))+bval(n)-basdif(n,ioff2+l))
                                d1 = d1*matel(n)
                            end do
                        end do
                    end do
                    ioff2 = ioff2 + nsteps(ioff1+k)
                end do
                ioff1 = ioff1 + nztrms(j)
            end do
        end if
    end do
    call NGA_PUT(q3,lo,hi,scr2,ld,dscale)
    call ga_sync()
    deallocate(scr1)
    call system_clock(tend,tcount,tmx)
    if(myid.eq.0.and.iter.le.niter)write(unit=OUTFILE,fmt=1001)1.*(tend-iterend)/tcount
    iterend = tend
    1001 format(3x,' - Total matrix-vector product:          ',f14.3,' secs.')
end subroutine Hv

subroutine lanczos_step(iter)
    use progdata, only: myid,nstates,dimen,alpha,beta,iterstart,iterend,q1,q2,q3,orthog,one,two,three
    use filedata, only: OUTFILE
    implicit none
#include 'global.fh'
    integer, intent(in)		:: iter
    integer                       :: tend,tcount,tmx
    real*8	        :: onedp
    onedp = 1.
    alpha(iter) = GA_DDOT(q2,q3)
    call GA_ADD(onedp,q3,-alpha(iter),q2,q3)
    if(iter.gt.1)call GA_ADD(onedp,q3,-beta(iter-1),q1,q3)
    call ga_sync()
    beta(iter) = Sqrt(GA_DDOT(q3,q3))
    call system_clock(tend,tcount,tmx)
    if(myid.eq.0) then
        write(unit=OUTFILE,fmt=1001)1.*(tend-iterend)/tcount
        if(.not.orthog.or.iter.eq.1)write(unit=OUTFILE,fmt=1002)1.*(tend-iterstart)/tcount
    end if
    iterend = tend
    1001 format(3x,' - Take Lanczos step:                    ',f14.3,' secs.')
    1002 format(3x,'Time to complete iteration:              ',f14.3,' secs.'/) 	
end subroutine lanczos_step

subroutine set_omega(iter)!Sets the omega array & checks orthogonality of the lanczos vectors
    use progdata, only: myid,alpha,beta,dimen,nstates,epsilon,one,two,three,omega,orthogexact, &
                        zerodp,chkorthog,iterend,q1,q2,q3,niter
    use filedata, only: OUTFILE,ARCHIVE
    implicit none
#include 'global.fh'
    integer, intent(in)                  :: iter
    integer                              :: i,j,k,istat,ireq
    integer                              :: tend,tcount,tmx
    real*8                     :: randnorm,dp,gauss_random
    real*8,dimension(2,iter-1) :: ocheck
    if(orthogexact) then
        call setarray(omega(2,0:iter),iter+1,zerodp)
        do i = 1,iter
            call read_ga(q1,ARCHIVE,i)
            call ga_sync()
            omega(2,i) = GA_DDOT(q3,q1)
        end do
    else
        if(iter.lt.niter)omega(2,iter+1) = 1.0
        omega(1,iter) = gauss_random(zerodp,dble(0.6))*epsilon*nstates*dimen*beta(1)/beta(iter)
        do k = 1,iter-1
            omega(1,k) = (1/beta(iter))*                                             &
                         (beta(k)*omega(2,k+1) + (alpha(k)-alpha(iter))*omega(2,k) + &
                          beta(k-1)*omega(2,k-1) - beta(iter-1)*omega(1,k)) +        &
                          epsilon*(beta(k)+beta(iter))*gauss_random(zerodp,dble(0.3))
        end do
        do i = 1,iter
            dp = omega(2,i)
            omega(2,i) = omega(1,i)
            omega(1,i) = dp
        end do
        if(Mod(iter,chkorthog).eq.0) then
            call setmatrix(ocheck,two,iter-1,zerodp)
            do i = 1,iter-1
                call read_ga(q1,ARCHIVE,i)
                call ga_sync()
                ocheck(1,i) = GA_DDOT(q2,q1)
                ocheck(2,i) = GA_DDOT(q3,q1)
                ocheck(2,i) = ocheck(2,i)/beta(iter)
            end do
            do i = 1,iter-1
                do k = 1,2
                   if(abs(omega(k,i)).lt.abs(ocheck(k,i))) then
                      if(myid.eq.0) write(*,1001)iter,iter+k-1,i,omega(k,i),ocheck(k,i)
                      omega(k,i) = ocheck(k,i)
                   end if
                end do
            end do
        end if
    end if
    call system_clock(tend,tcount,tmx)
    if(myid.eq.0.and.(orthogexact.or.Mod(iter,chkorthog).eq.0)) then
        write(unit=OUTFILE,fmt=1002)1.*(tend-iterend)/tcount
    elseif(myid.eq.0) then
        write(unit=OUTFILE,fmt=1003)1.*(tend-iterend)/tcount
    end if
    iterend = tend
    1000 format('OMEGA(',i5,',',i5,')=',es18.10,' EXACT:',es18.10)
    1001 format('ITERATION ',i4,', OMEGA(',i4,',',i4,') -- ',es18.10,' --> ',es18.10)
    1002 format(3x,' - Determing Orthogonality, Exact:       ',f14.3,' secs.')
    1003 format(3x,' - Determing Orthogonality, Recurrence:  ',f14.3,' secs.')
    1004 format(3x,' OMEGA(',i3,',',i3,') -> recur=',es18.10,'  exact=',es18.10)
    1005 format(3x,' ITER=',i3,' J=',i3,' OMGEGA=',es18.10,' ACTUAL=',es18.10)
end subroutine set_omega

!This formulation is based on the Partial Orthogonalization procedure of Horst Simon, Mathematics of Computation, 42, 115-142, (1984)
subroutine partial_orthog(iter)
    use progdata, only: myid,alpha,beta,omega,dimen,nstates,niter,epsilon,eta,iterstart,iterend, &
                        one,two,zero,zerodp,q1,q2,q3,loindex,roindex,northog,firststep
    use filedata, only: OUTFILE,ARCHIVE
    implicit none
#include 'global.fh'
    integer,intent(in)		                        :: iter
    integer			                        :: i,j,noindices,istat,ireq
    integer                                               :: tend,tcount,tmx
    integer,dimension(Ceiling(1.*iter/10))                :: oindices,lindices,rindices
    real*8		                        :: dpval,onedp,gauss_random
    real*8,dimension(Ceiling(1.*iter/10.),iter) :: dpvec
    onedp=1d0
    if(firststep) then!Check if partial reorthogonalization is necessary
        noindices = 0
        northog = 0
        do i = 1,iter
            if(abs(omega(2,i)).ge.Sqrt(epsilon)) then
                noindices = noindices + 1
                oindices(noindices) = i
            end if
        end do
        if(noindices.gt.0) then
            do i = 1,noindices
                do j = oindices(i),1,-1
                    if(abs(omega(2,j)).lt.eta) exit
                end do
                if(j.lt.1)j=1
                lindices(i) = j
            end do
            do i = 1,noindices
                do j = oindices(i),iter
                    if(abs(omega(2,j)).lt.eta) exit
                end do
                if(j.gt.iter) j=iter
                rindices(i) = j
            end do
            northog = 1
            loindex(northog) = lindices(northog)
            roindex(northog) = rindices(northog)
            do i = 2,noindices
                if(lindices(i).ne.loindex(northog).and.rindices(i).ne.roindex(northog)) then
                    northog = northog + 1
                    loindex(northog) = lindices(i)
                    roindex(northog) = rindices(i)
                end if
            end do
            do i = 1,northog
                do j = 1,i-1
                    if(loindex(i).lt.loindex(j).and.roindex(i).gt.loindex(j).or. &
                       loindex(j).lt.loindex(i).and.roindex(j).gt.loindex(i))    &
                        write(OUTFILE,*)' -- ERROR: OVERLAPPING ORTHOGONALIZATION RANGES -- '
                end do
            end do
            if(minval(loindex,mask=loindex.gt.0).gt.1) then
                northog = northog + 1
                loindex(northog) = one
                roindex(northog) = one
            end if
        end if
    end if
    if(northog.gt.0) then!If partial orthogonalization is necessary...
        do i = 1,northog
            do j = loindex(i),roindex(i)
                call read_ga(q1,ARCHIVE,j)
                call ga_sync()
                dpvec(i,j-loindex(i)+1) = GA_DDOT(q3,q1)
            end do
        end do
        do i = 1,northog
            do j = loindex(i),roindex(i)
                call read_ga(q1,ARCHIVE,j)
                call ga_sync()
                call GA_ADD(onedp,q3,-dpvec(i,j-loindex(i)+1),q1,q3)
                omega(2,j) = epsilon*gauss_random(zerodp,dble(2.0))
            end do 
        end do
        firststep = .not.firststep
        call ga_sync()
        beta(iter) = Sqrt(GA_DDOT(q3,q3))
    end if
    call system_clock(tend,tcount,tmx)
    if(myid.eq.0) then
        if(northog.gt.0) then
            write(unit=OUTFILE,fmt=1000)loindex(1),roindex(1),1.*(tend-iterend)/tcount
        else
            write(unit=OUTFILE,fmt=1001)1.*(tend-iterend)/tcount
        end if
        write(unit=OUTFILE,fmt=1002)1.*(tend-iterstart)/tcount
    end if
    iterend = tend
    1000 format(3x,' - Perform Partial Orthog. ',i4,' - ',i4,':  ',f14.3,' secs.')
    1001 format(3x,' - No re-orthogonalization Required:     ',f14.3,' secs.')
    1002 format(3x,'Time to complete iteration:              ',f14.3,' secs.'/)
end subroutine partial_orthog		       

subroutine make_restart(iter)!Dumps alpha and beta coefficents to dump.dat
    use progdata
    use filedata, only: outdir,OUTFILE,RESTARTFILE,QRESTART
    implicit none
#include 'global.fh'
#include 'mafdecls.fh'
    integer,intent(in):: iter
    integer:: rlen,tcount,tmx
    integer:: istat,ireq
    character*4::procindex
    call system_clock(totend,tcount,tmx)
    write(procindex,'(i4)')myid
    procindex = adjustl(procindex)
    if(myid.eq.0) then
        write(unit=OUTFILE,fmt='(a)')'  - Lanczos Iterations Completed - '
        write(unit=OUTFILE,fmt=1003)1.*(totend-totstart)/tcount
        open(unit=RESTARTFILE,file='restart.log',status='replace')
        rewind(unit=RESTARTFILE)
        write(unit=RESTARTFILE,fmt='(a)')' TOTAL NUMBER OF ITERATIONS IN FILE'
        write(unit=RESTARTFILE,fmt='(i5)')iter
        write(unit=RESTARTFILE,fmt='(a)')' ITERATION NUMBER TO RESTART FROM'
        write(unit=RESTARTFILE,fmt='(i5)')iter
        write(unit=RESTARTFILE,fmt='(a)')' ALPHA COEFFICIENTS -- diagonal T matrix elements'
        write(unit=RESTARTFILE,fmt=*)alpha(1:iter)
        write(unit=RESTARTFILE,fmt='(a)')' BETA COEFFICIENTS -- off-diagonal T matrix elements'
        write(unit=RESTARTFILE,fmt=*)beta(0:iter)
        if(orthog) then
            write(unit=RESTARTFILE,fmt='(a)')' OMEGA(1) -- ORTHOGONALITY INFORMATION'
            write(unit=RESTARTFILE,fmt=*)omega(1,0:iter)
            write(unit=RESTARTFILE,fmt='(a)')' OMEGA(2) -- ORTHOGONALITY INFORMATION'
            write(unit=RESTARTFILE,fmt=*)omega(2,0:iter)
        end if
    end if
    call ga_sync()
    rlen = vbounds(myid*nseg+nseg,2) - vbounds(myid*nseg+1,1) + 1
    open(unit=QRESTART,file=trim(adjustl(outdir))//'nadvibs.35.'//trim(procindex),access='direct', &
         status='replace',form='unformatted',recl=8*rlen)
    rewind(QRESTART)
    call write_ga(q1,QRESTART,1)
    call write_ga(q2,QRESTART,2)
    call write_ga(q3,QRESTART,3)
    close(QRESTART)
    1003 format('  Time Required: ',f14.5,' secs.')
end subroutine make_restart

subroutine compute_eigenvalues(iter)!Computes the eigenvalues of the H
    use progdata, only: niter,dimen,nstates,alpha,beta,totstart,totend,aomega,nmodes,AU2WAVE,AU2EV,bjiconv
    use filedata, only: OUTFILE
    implicit none
    integer,intent(in)     :: iter
    integer                :: i,tcount,ngood,tmx,itotal
    integer,dimension(iter):: Tsign
    real*8                            :: dpval,thresh
    real*8, dimension(:,:),allocatable:: Teig
    real*8, dimension(:),allocatable  :: T
    real*8, dimension(iter)           :: Tvals,Tints,bji
    write(unit=OUTFILE,fmt='(a)') ''
    write(unit=OUTFILE,fmt='(a)') 'COMPUTED EIGENVALUES '
    write(unit=OUTFILE,fmt='(a)') '------------------------------'
    write(unit=OUTFILE,fmt='(a)') ''
    allocate(Teig(iter,iter))
    allocate(T(iter*(iter+1)/2))
    call make_tmatrix(T,iter)
    call compute_ritzsystem(T,iter,bji,Tvals,Teig,Tints,Tsign)
    thresh = 0.0000
    ngood = 1
    T(ngood) = 1
    do i = 2,iter
        if(bji(i).lt.bjiconv) then
            if(AU2WAVE*abs(Tvals(i)-Tvals(T(ngood))).lt.thresh) then
                Tints(T(ngood)) = Tints(T(ngood)) + Tints(i)
                if(bji(i).lt.bji(T(ngood))) then
                    Tints(i) = Tints(T(ngood))
                    T(ngood) = i
                end if
            else
                ngood = ngood + 1
                T(ngood) = i
            end if
        end if
    end do
    dpval = 0d0
    do i = 1,nmodes
     dpval = dpval + aomega(i)*0.5
    end do
    write(unit=OUTFILE,fmt=1001)dpval*AU2WAVE
    write(unit=OUTFILE,fmt=1002)Tvals(T(1))*AU2WAVE
    dpval = 0
    do i = 1,ngood
     if(Tints(T(i)).gt.dpval)dpval=Tints(T(i))
    end do
    write(unit=OUTFILE,fmt=1003)dpval
    write(unit=OUTFILE,fmt='(a)')''
    write(unit=OUTFILE,fmt=1004)
    do i = 1,ngood
        write(unit=OUTFILE,fmt=1005)i,int(T(i)),(Tvals(T(i))-Tvals(T(1)))*AU2WAVE,(Tvals(T(i))-Tvals(1))*AU2EV, &
                                    bji(T(i)),Tints(T(i))/dpval
    end do
    call system_clock(totstart,tcount,tmx)
    write(unit=OUTFILE,fmt='(a)') ' '
    write(unit=OUTFILE,fmt='(a)') '   All requested eigenvalues computed.'
    write(unit=OUTFILE,fmt=1000)1.*(totstart-totend)/tcount
    deallocate(T)
    deallocate(Teig)
    1000 format('  Time Required: ',f14.5,' secs.')
    1001 format('  ZPVE from nadvibs.in (reference state): ',f12.4,' cm-1')
    1002 format('  FIRST eigenvalue:                       ',f12.4,' cm-1')
    1003 format('  Intensity Factor (MAX Intensity):       ',f12.8)
    1004 format(' Root  Index     E(cm-1)       E(eV)     Convergence   Intensity')
    1005 format(i4,'  ',i4,'  ',f13.5,'  ',f10.5,'  ',es12.4,'  ',f8.4)
end subroutine compute_eigenvalues

subroutine print_footer()
    use filedata, only: OUTFILE
    write(unit=OUTFILE,fmt='(a)') ' '
    write(unit=OUTFILE,fmt='(a)') '------------ FINISHED EXECUTION ------------'   
    close(unit=OUTFILE)
end subroutine print_footer

!This subroutine computes the eigenvectors of the hamiltonian matrix.  It then determines diagnostic
!information about the roots, as well as dot products useful in spin-orbit computations
subroutine identify_roots(iter)
    use progdata, only: myid,dimen,nstates,nmodes,AU2WAVE,AU2EV,bjiconv,zerodp,idroots, &
                        q1,q2,q3,vbounds,nseg,soroots,lo,hi
    use filedata, only: outdir,SOINFO,ROOTINFO,ARCHIVE,EIGVECS
    implicit none
#include 'global.fh'
    integer,intent(in)          :: iter
    integer                     :: i,j,k,l,counter,nroots
    integer,dimension(1)        :: ld
    integer,dimension(2)        :: mxval,alo,ahi,blo,bhi
    integer,dimension(nmodes)   :: barray
    integer,dimension(2*idroots):: cindex
    integer,dimension(iter)     :: Tsign
    real*8,dimension(iter)           :: Tvals,Tints,bji
    real*8,dimension(2*nstates)      :: dp
    real*8,dimension(idroots)        :: cmag,cphase
    real*8,dimension(:),allocatable  :: T
    real*8,dimension(:,:),allocatable:: Tvecs
    real*8                           :: dpval,intfactor
    character*4:: procindex
    allocate(T(iter*(iter+1)/2))
    allocate(Tvecs(iter,iter))
    nroots = max(idroots,soroots)
    write(procindex,'(i4)')myid
    procindex = adjustl(procindex)
    ld(1) = nstates 
    if(nroots.gt.0.and.myid.eq.0)then
        open(unit=ROOTINFO,file='rootinfo.dat',status='replace')
        rewind(unit=ROOTINFO)
    end if
    if(soroots.gt.0)then
        if(myid.eq.0)then
            open(unit=SOINFO,file='sodata.in',status='replace')
            rewind(unit=SOINFO)
        end if
        open(unit=EIGVECS,file=trim(adjustl(outdir))//'nadvibs.37.'//trim(procindex),access='direct', &
            status='replace',form='unformatted',recl=8*(hi(1)-lo(1)+1))
        rewind(EIGVECS)
    end if
    call make_tmatrix(T,iter)
    call compute_ritzsystem(T,iter,bji,Tvals,Tvecs,Tints,Tsign)
    intfactor = 0
    do i = 1,nroots
        if(Tints(i).gt.intfactor)intfactor=Tints(i)
    end do
    if(myid.eq.0) then
        write(unit=ROOTINFO,fmt=1000)nroots
        if(soroots.gt.0) then
            write(unit=SOINFO,fmt=1001)soroots
            write(unit=SOINFO,fmt=1002)
            do i = 1,soroots
                write(unit=SOINFO,fmt=1003)i,Tvals(i)*AU2EV,Sqrt(Tints(i)/intfactor)*Tsign(i)
            end do
            if(nstates.eq.2)write(unit=SOINFO,fmt=1004)
            if(nstates.eq.3)write(unit=SOINFO,fmt=1005)
        end if
    end if
    call ga_sync()
    !Compute Eigenvectors
    do i = 1,nroots
        call GA_ZERO(q1)
        dpval = 1.D0
        do k = 1,iter
            call read_ga(q2,ARCHIVE,k)
            call GA_ADD(dpval,q1,Tvecs(k,i),q2,q1)
        end do
        if(i.le.soroots)call write_ga(q1,EIGVECS,i)
        call ga_sync()
        if(myid.eq.0)write(unit=ROOTINFO,fmt=1006)i
        !Analyze Eigenvectors
            if(myid.eq.0) then
                write(unit=ROOTINFO,fmt='(a)')' '
                write(unit=ROOTINFO,fmt=1007)i
                write(unit=ROOTINFO,fmt=1008)Tvals(i)*AU2WAVE,Tvals(i)*AU2EV,Tints(i)/intfactor
            end if
            call GA_ELEM_MULTIPLY(q1,q1,q2)
            do k = 1,5
                call NGA_SELECT_ELEM(q2,'max',dpval,mxval)     
                call NGA_GET(q1,mxval,mxval,cphase(k:k),ld)
                cmag(k)       = dpval
                cphase(k)     = sign(1.,cphase(k))
                cindex(2*k-1) = mxval(1)
                cindex(2*k)   = mxval(2)
                call GA_FILL_PATCH(q2,mxval(1),mxval(1),mxval(2),mxval(2),zerodp)
            end do
            dpval = 0.
            call ga_sync()
            do k = 1,nstates
                dp(k) = GA_DDOT_PATCH(q1,'N',1,dimen,k,k,q1,'N',1,dimen,k,k)
                dpval = dpval + dp(k)
            end do
            do k = 1,nstates
                dp(k) = 100.*abs(dp(k)/dpval)
            end do
            if(myid.eq.0) then
                do j = 1,5
                    call get_barray(cindex(2*j-1),barray)
                    do k = 1,nmodes
                        barray(k) = barray(k) - 1
                    end do
                    write(unit=ROOTINFO,fmt=1009)cindex(2*j),cindex(2*j-1),int(cphase(j)),cmag(j)
                    write(unit=ROOTINFO,fmt=1010)barray
                end do
                write(unit=ROOTINFO,fmt='(a)')' '
                if(nstates.eq.2)write(unit=ROOTINFO,fmt=1011)dp(1),dp(2)
                if(nstates.eq.3)write(unit=ROOTINFO,fmt=1012)dp(1),dp(2),dp(3)
                if(nstates.eq.4)write(unit=ROOTINFO,fmt=1016)dp(1),dp(2),dp(3),dp(4)
            end if
            if(myid.eq.0) write(unit=ROOTINFO,fmt='(a)')' --------------------------------------------------------------'
        !Compute spin-orbit parameters, if requested
            if(i.le.soroots)then
                do j = 1,i-1
                    counter = 0
                    call read_ga(q2,EIGVECS,j)
                    do k = 1,nstates
                        do l = 1,nstates
                            if(k.ne.l) then
                                counter = counter + 1
                                call GA_COPY_PATCH('N',q1,1,dimen,k,k,q3,1,dimen,l,l)
                                call ga_sync()
                                dp(counter) = GA_DDOT_PATCH(q3,'N',1,dimen,l,l,q2,'N',1,dimen,l,l)
                            end if
                        end do
                    end do
                    call ga_sync()
                    if(myid.eq.0) then
                        if(nstates.eq.2)write(unit=SOINFO,fmt=1013)i,j,dp(1),dp(2)
                        if(nstates.eq.3)write(unit=SOINFO,fmt=1014)i,j,dp(1),dp(2),dp(3),i,j,dp(4),dp(5),dp(6)
                    end if
                end do
            end if
    end do
    if(myid.eq.0) then
        if(soroots.gt.0) close(unit=SOINFO)
        write(unit=ROOTINFO,fmt='(a)') ' '
        write(unit=ROOTINFO,fmt='(a)') '------------ FINISHED EXECUTION ------------'   
        close(unit=ROOTINFO)
    end if
    if(soroots.gt.0)close(EIGVECS,status='delete')
    call ga_sync()
    deallocate(T)
    deallocate(Tvecs)
    1000 format(' Summary for computing ',i4,' eigenvectors corresponding to the lowest eigenvalues')
    1001 format(i5)
    1002 format(' Root   Energy (eV)      Sqrt(Int.)*sign')
    1003 format(i5,2x,es15.8,2x,es15.8)
    1004 format(' R(i)  R(j)   [S1(i)S2(j)]   [S2(i)S1(j)]')
    1005 format(' R(i)  R(j)   [S1(i)S2(j)]   [S1(i)S3(j)]  [S2(i)S1(j)]', &
               /' R(i)  R(j)   [S2(i)S3(j)]   [S3(i)S1(j)]  [S3(i)S2(j)]')
    1006 format(' Eigenvector ',i4,' has been computed...')
    1015 format('     --- Error in eigenvalue check: ',es14.6,' cm-1')
    1007 format(' ************ ROOT ',i4,' *************')
    1008 format(' Energy: ',f9.2,' cm-1, ',f9.5,' eV   Intensity: ',f6.4/)
    1009 format(' State vector: ',i2,' | Index: ',i12,' Phase: ',i2,' | Magnitude^2: ',f8.6)
    1010 format(' v(i): ',25(i3))
    1011 format(' Percent State (1,2): (',f6.2,' %, ',f6.2,' %)')
    1012 format(' Percent State (1,2,3): (',f6.2,' %, ',f6.2,' %, ',f6.2,' %)') 
    1016 format(' Percent State (1,2,3,4): (',f6.2,' %, ',f6.2,' %, ',f6.2,' %, ',f6.2,' %)') 
    1013 format(i5,1x,i5,2x,es14.6,1x,es14.6)
    1014 format(i5,1x,i5,2x,es14.6,1x,es14.6,1x,es14.6,/i5,1x,i5,2x,es14.6,1x,es14.6,1x,es14.6)
end subroutine identify_roots

subroutine release_memory()
    use progdata
    use filedata
    implicit none
#include 'global.fh'
    logical::lstat
    integer::istat
    if(allocated(scr1))deallocate(scr1)
    deallocate(scr2)
    deallocate(concoef)
    deallocate(nzindx)
    deallocate(nzblks)
    deallocate(nzcoef)
    deallocate(alpha)
    deallocate(beta)
    deallocate(aomega)
    deallocate(bomega)
    deallocate(omega)
    deallocate(nfunc)
    deallocate(shft)
    deallocate(statew)
    deallocate(loindex)
    deallocate(roindex)
    if(allocated(dmap))deallocate(dmap)
    if(allocated(dordr))deallocate(dordr)
    deallocate(bmax)
    deallocate(nsteps)
    deallocate(stype)
    deallocate(basind)
    deallocate(basdif)
    deallocate(bstart)
    deallocate(rcshft)
    deallocate(nelems)
    if(savevecs)close(ARCHIVE)
    if(myid.eq.0)then
        close(unit=SCRFILE,status='delete')
        close(unit=RESTARTFILE)
    end if
    lstat = GA_DESTROY(q1)
    if(.not.lstat)print *,'ERROR in GA_DESTROY(q1), ID=',myid
    lstat = GA_DESTROY(q2)
    if(.not.lstat)print *,'ERROR in GA_DESTROY(q2), ID=',myid
    lstat = GA_DESTROY(q3) 
    if(.not.lstat)print *,'ERROR in GA_DESTROY(q3), ID=',myid
    call ga_sync()
end subroutine release_memory
!------------------------------------ End ------------------------------------

!------------ Secondary subroutines listed in alphabetical order -------------
subroutine compute_allIndex()
    use progdata, only: myid,nseg,nproc,dimen,nmodes,nfunc,vbounds,ordr,nztrms,nzindx,lo,hi,ndiag,dordr,dmap, &
                        basind,basdif,bstart,rcshft,nelems,bmax,nsteps,zero,one,two,shft,nstblks,stype,nalltrms
    implicit none
    integer                             :: i,j,k,l,m,n,p,indx,nindx,ioff1,ioff2
    integer                             :: myrange,nbat,stride,ndir,ntot
    integer                             :: cindx,rindx,cshft,rshft
    integer                             :: ibat,ival
    integer                             :: get_batch
    integer,dimension(nmodes)           :: bc,br
    integer,dimension(ordr)             :: sdir,sdif,inds,icnt
    integer,dimension(ordr+1)           :: bcnt
    integer,dimension(ordr,sum(nztrms)) :: smax
    integer,dimension(5000)             :: buf1,buf2
    logical                             :: diagterm
    nbat = nproc*nseg
    ndiag = 0
    myrange = hi(1) - lo(1) + 1
    !Set the batch independent max parameters bmax,bshft
    allocate(nalltrms(nbat))
    allocate(nsteps(sum(nztrms)))
    allocate(bmax(ordr,sum(nztrms)))
    call setintmatrix(bmax,ordr,sum(nztrms),one)
    ioff1 = 0!Determine how many matrix element types there are per potential term
    do i = 1,ordr
        do j = 1,nztrms(i)
            ndir = 1
            do k = 1,i
                indx  = nzindx(2*k-1,ioff1+j)
                nindx = nzindx(2*k,ioff1+j)
                smax(k,ioff1+j) = nindx
                do 
                    if(smax(k,ioff1+j)<=nfunc(indx)-1) exit
                    smax(k,ioff1+j) = smax(k,ioff1+j) - 2
                end do
                if(smax(k,ioff1+j).lt.0)ndir = 0
                ndir = ndir*(smax(k,ioff1+j)+1)
                if(k==1)then
                    ival = nmodes
                else
                    ival = nzindx(2*k-3,ioff1+j)-1
                end if
                do l = indx+1,ival
                    bmax(k,ioff1+j) = bmax(k,ioff1+j)*nfunc(l)
                end do
            end do
            diagterm = .true.
            do k = 1,i
                if(mod(smax(k,ioff1+j),2).ne.0)then
                    diagterm = .false.
                    EXIT
                end if
            end do
            if(diagterm)then
                ndiag        = ndiag + 1
                buf1(ndiag)  = i
                buf2(ndiag)  = ioff1+j
            end if
            nsteps(ioff1+j) = ndir
        end do
        ioff1 = ioff1 + nztrms(i)
    end do
    ntot = sum(nsteps)
    if(ndiag>0)then!Allocate index array - for each iproc section of the Hamiltonian
        allocate(dordr(ndiag))
        allocate(dmap(ndiag))
    end if
    allocate(stype(ordr,ntot))
    allocate(basdif(ordr,ntot))
    allocate(nelems(ntot,nbat))
    allocate(basind(ordr,ntot,nbat))
    allocate(bstart(ordr,ntot,nbat))
    allocate(rcshft(two,ntot,nbat))
    call setintmatrix(nelems,ntot,nbat,zero)
    do i = 1,ndiag
        dordr(i) = buf1(i)
        dmap(i)  = buf2(i)
    end do
    ioff1 = 0
    ioff2 = 0
    do i = 1,ordr!Loop over n-index terms
        call setintarray(sdif,ordr,zero)
        call setintarray(sdir,ordr,zero)
        do j = 1,nztrms(i)!Loop over number of terms with i indices
            do k = 1,i
                sdir(k) = -smax(k,ioff1+j)
                inds(k) = nzindx(2*k-1,ioff1+j)
                icnt(k) = nzindx(2*k,ioff1+j)
            end do
            sdir(i) = sdir(i) - 2
            do k = 1,nsteps(ioff1+j)
                l = i
                do
                    if(sdir(l).lt.smax(l,ioff1+j)) exit
                    sdir(l) = -smax(l,ioff1+j)
                    l = l - 1
                end do
                sdir(l) = sdir(l) + 2
                stride = 0
                do l = 1,i
                    stype(l,ioff2+k) = int(abs(sdir(l)/2.))+1
                    basdif(l,ioff2+k) = abs(sdir(l))
                    stride = stride + sdir(l)*shft(inds(l))
                end do
                cindx = lo(1) + stride
                rindx = lo(1)
                rshft = 1
                if(cindx<=-myrange.or.cindx>dimen.or.stride==0)then
                    ibat=0
                else
                    if(cindx<=zero)then
                        rshft = rshft - (cindx-1)
                        rindx = rindx - (cindx-1)
                        cindx = 1
                    end if
                    ibat = get_batch(cindx)
                    cshft = cindx - vbounds(ibat,1) + 1           
                end if
                if(ibat.ne.0)then!If there are matrix elements handled by this processor
                    call get_barray(cindx,bc)
                    call get_barray(rindx,br)
                    do l = 1,i
                      sdif(l) = bc(inds(l)) - br(inds(l)) - sdir(l)
                    end do
                    ival = dot_product(sdif,sdif)
                    do       
                        if(ival==0.or.rshft>myrange.or.ibat>nbat) exit
                        cindx = cindx + 1
                        rindx = rindx + 1
                        rshft = rshft + 1
                        call get_barray(cindx,bc)
                        call get_barray(rindx,br)
                        do l = 1,i
                            sdif(l) = bc(inds(l)) - br(inds(l)) - sdir(l)
                        end do
                        ival = dot_product(sdif,sdif)
                        do
                            if(ibat>nbat) exit
                            if(cindx<=vbounds(ibat,2)) exit
                            ibat = ibat + 1
                        end do
                        if(ibat<=nbat) cshft = cindx - vbounds(ibat,1) + 1
                    end do
                    if(ibat<=nbat) then
                        do l = 1,i
                            bcnt(l) = 1
                            if(l==1)then
                                ival = nmodes
                            else
                                ival = inds(l-1)-1
                            end if
                            do m = inds(l)+1,ival
                                n = 1
                                do p = (m+1),ival 
                                    n = n*nfunc(p)
                                end do
                                bcnt(l) = bcnt(l) + (bc(m)-1)*n
                            end do
                        end do
                        do l = 1,i
                            basind(l,ioff2+k,ibat) = max(bc(inds(l)),br(inds(l)))
                            bstart(l,ioff2+k,ibat) = bcnt(l)
                        end do
                        rcshft(1,ioff2+k,ibat)   = cshft
                        rcshft(2,ioff2+k,ibat)   = rshft
                        do
                            if(rshft>myrange.or.ibat>nbat) exit
                            nelems(ioff2+k,ibat) = nelems(ioff2+k,ibat) + 1
                            ! begin constructing next term
                            cindx = cindx + 1
                            rindx = rindx + 1
                            cshft = cshft + 1
                            rshft = rshft + 1
                            bcnt(1) = bcnt(1) + 1
                            do l = 1,i
                                if(bcnt(l)<=bmax(l,ioff1+j)) exit
                                bcnt(l) = 1
                                bc(inds(l)) = bc(inds(l)) + 1
                                br(inds(l)) = br(inds(l)) + 1
                                ival = max(br(inds(l)),bc(inds(l)))
                                if(ival>nfunc(inds(l)))then
                                    bcnt(l+1) = bcnt(l+1) + 1
                                    ival = abs(sdir(l))
                                    bc(inds(l)) = 1 + ival*(1+sign(one,sdir(l)))/2
                                    br(inds(l)) = 1 + ival*(1-sign(one,sdir(l)))/2
                                    cindx = cindx + ival*shft(inds(l))
                                    rindx = rindx + ival*shft(inds(l))
                                    cshft = cshft + ival*shft(inds(l))
                                    rshft = rshft + ival*shft(inds(l))
                                end if
                            end do
                            if(cindx>vbounds(ibat,2)) then
                                do
                                    if(ibat>nbat) exit
                                    if(cindx<=vbounds(ibat,2)) exit
                                    ibat = ibat + 1
                                end do
                                if(ibat<=nbat) then
                                    cshft = cindx - vbounds(ibat,1) + 1
                                    do l = 1,i
                                        basind(l,ioff2+k,ibat) = max(bc(inds(l)),br(inds(l)))
                                        bstart(l,ioff2+k,ibat) = bcnt(l)
                                    end do
                                    rcshft(1,ioff2+k,ibat) = cshft
                                    rcshft(2,ioff2+k,ibat) = rshft
                                end if
                            end if
                        end do
                    end if
                end if
            end do
            ioff2 = ioff2 + nsteps(ioff1+j)
        end do
        ioff1 = ioff1 + nztrms(i)
    end do
    do i = 1,nbat
        ival = 0
        do j = 1,ntot
            ival = ival + nelems(j,i)
        end do
        if(i>=(myid*nseg+1).and.i<=(myid*nseg+nseg))ival = ival + vbounds(i,2)-vbounds(i,1)+1
        nalltrms(i) = ival
    end do
end subroutine compute_allIndex

!Given a tridiagonal matrix, assuming requisite lanczos vectors are on disk,
!compute the ritz vectors and ritz values.  The ritz vectors are writen to the
!file specified by RITZVECS, while the ritz values are returned in an array.  Will use the 
!1st n lanczos vectors from file.  The ritz vectors are written starting at vector roff+1 in the
!DA file RITZVECS
!    T:        tridiagonal array,                                                  intent(inout)
!    n:        dimension of T (n x n),                                             intent(in)
!    bji:      an array whose values indicated the convergence of an eigenvalue    intent(inout)
!    eigvals:  eigenvalues of T                                                    intent(inout) 
!    eigvecs:  eigenvectors of T                                                   intent(inout)
!    eigints:  intensities of eigenvalues of T                                     intent(inout)
!    eigsign:  sign of <0|Psi(i)>                                                  intent(inout)
subroutine compute_ritzsystem(T,n,bji,eigvals,eigvecs,eigints,eigsign)
    use progdata, only: one,beta
    implicit none
    integer,intent(in)                        :: n
    integer,dimension(n),intent(inout)        :: eigsign
    real*8,dimension(n*(n+1)/2),intent(inout) :: T
    real*8,dimension(n),intent(inout)         :: eigvals,eigints,bji
    real*8,dimension(n,n),intent(inout)       :: eigvecs
    real*8,dimension(5*n)                     :: tbuf
    integer                                   :: i
    call givens(n,n,n,T,tbuf,eigvals,eigvecs)
    call setintarray(eigsign,n,one)
    do i = 1,n!Test the beta(ji)
        bji(i) = abs(beta(n)*eigvecs(n,i))
        if(eigvecs(1,i).lt.0)eigsign(i)=-1
        eigints(i) = eigvecs(1,i)*eigvecs(1,i)
    end do
end subroutine compute_ritzsystem

!Determines the machine precision
subroutine determine_epsilon()
    use progdata, only: epsilon,eta
    implicit none
    real*8       :: ep1
    epsilon = dble(1)
    do
        epsilon = epsilon/dble(2)
        ep1 = epsilon + dble(1)
        if(ep1==dble(1)) exit
    end do
    eta = epsilon**(3d0/4d0)
end subroutine determine_epsilon

subroutine write_ga(ga,ga_file,recstart)
    use progdata, only: myid,nproc,nseg,nstates,vbounds,scr2
    implicit none
    integer, intent(in)   :: ga,ga_file,recstart
    integer               :: i,j,currec
    integer, dimension(2) :: lo,hi
    integer, dimension(1) :: ld
    lo = (/ vbounds(myid*nseg+1,1),    1 /)
    hi = (/ vbounds(myid*nseg+nseg,2), nstates /)
    ld = (/ hi(1) - lo(1) + 1 /)
    call nga_get(ga,lo,hi,scr2,ld)
    currec = (int(recstart)-1)*nstates + 1
    do i = 1,nstates
        write(ga_file,rec=currec)(scr2(j,i),j=1,ld(1))
        currec = currec + 1
    end do
end subroutine write_ga

subroutine read_ga(ga,ga_file,recstart)
    use progdata, only: myid,nproc,nseg,nstates,vbounds,scr2
    implicit none
    integer,intent(in)         :: ga,ga_file,recstart
    integer                    :: i,j,currec
    integer, dimension(2)      :: lo,hi
    integer, dimension(1)      :: ld
    real*8,parameter :: dscale = 1.
    lo = (/ vbounds(myid*nseg+1,1),    1 /)
    hi = (/ vbounds(myid*nseg+nseg,2), nstates /)
    ld = (/ hi(1) - lo(1) + 1 /)
    currec = (int(recstart)-1)*nstates + 1
    do i = 1,nstates
        read(ga_file,rec=currec)(scr2(j,i),j=1,ld(1))
        currec = currec + 1
    end do
    call NGA_PUT(ga,lo,hi,scr2,ld,dscale)
end subroutine read_ga

!Reduce Matrix to upper triangular form
subroutine gauss_elim(m,n,MATIN,MATOUT)
    implicit none
    integer,intent(in):: m,n
    real*8,dimension(m,n),intent(in)    :: MATIN
    real*8,dimension(m,n),intent(inout) :: MATOUT
    integer:: i,j,k,maxr
    real*8:: val
    do i = 1,m
        do j = 1,n
            MATOUT(i,j) = MATIN(i,j)
        end do
    end do
    do i = 1,m-1
        maxr = i
        val = abs(MATOUT(i,i))
        do j = i+1,m
            if(abs(MATOUT(j,i))>val) then
                val = abs(MATOUT(j,i))
                maxr = j
            end if
        end do
        if(maxr.ne.i)then
            do j = 1,n
                val = MATOUT(maxr,j)
                MATOUT(maxr,j) = MATOUT(i,j)
                MATOUT(i,j) = val
            end do
        end if
        do j = i+1,m
            if(MATOUT(j,i).ne.0) then
                val = -MATOUT(j,i)/MATOUT(i,i)
                do k = 1,n
                    MATOUT(j,k) = MATOUT(j,k) + val*MATOUT(i,k)
                end do
            end if
        end do
    end do
end subroutine gauss_elim

!Generate an inital lanczos vector, given a set of weights
subroutine generate_initialvec()
    use progdata, only: myid,nproc,nmodes,nstates,dimen,nfunc,shft,nseg,neworigin,vbounds,statew,scr2,q1,q3, &
                        aomega,bomega,dvec,tmat,nirreps,npirrep,totstart,totend,zero,zerodp,one,lo,hi,istate
    use filedata, only: OUTFILE
    implicit none
#include 'mpif.h'
#include 'global.fh'
    logical::fscaled
    integer:: i,j,k,l,m,n,ierr,tmx,tcount,lbuf,stateindex
    integer:: nlev,pstart,pindex,cstart,cindex,sindex,stindex
    integer:: i2,ndo,ndone,npack,ninit,imode,get_procid,get_batch
    integer,dimension(1)      :: ld,ldi
    integer,dimension(nmodes) :: bvec,irrep
    integer,dimension(nirreps):: irroff
    integer,dimension(:),allocatable:: labels,nlterms
    real*8::detA,detT,d1,d2,fac
    real*8,dimension(nmodes):: acoef,bcoef,A2d,A3d
    real*8,dimension(5*nmodes)::Ascr
    real*8,dimension(nmodes*(nmodes+1)/2)::A
    real*8,dimension(nmodes,nmodes)::Avec,Ainv,A1,A2,A3,tmatrix
    real*8,dimension(:),allocatable::Cm,Cmbuf,Cn,Cnbuf
    if(neworigin) then
        !Initialize variables and such
        if(myid==0) write(OUTFILE,1001)
        ld(1) = hi(1) - lo(1) + 1
        lbuf = int((vbounds(1,2) - vbounds(1,1)+1)*nstates/(3*nmodes+2))
        nlev = 0
        do i = 1,nmodes
            nlev = nlev + nfunc(i) - 1
        end do
        fscaled = .false.
        if(fscaled)then
            do i = 1,nmodes
                acoef(i) = 1.
                bcoef(i) = 1.
            end do
        else
            do i = 1,nmodes
                 acoef(i) = aomega(i)
                 bcoef(i) = bomega(i)
            end do
        end if
        do i = 1,nirreps
            if(i>1)then
                irroff(i) = irroff(i-1) + npirrep(i-1)
            else
                irroff(i) = 0
            end if
            do j = 1,npirrep(i)
                irrep(j+irroff(i)) = i
            end do
        end do
        !Construct Matrix Quantities
            do i = 1,nmodes
                do j = 1,nmodes
                    tmatrix(j,i) = tmat((i-1)*nmodes+j)
                end do
            end do
            !Get the determinant of the T matrix
            call gauss_elim(nmodes,nmodes,tmatrix,Avec)
            detT = 1.
            do i = 1,nmodes
                detT = detT*abs(Avec(i,i))
            end do
            if(myid==0) print *,'Det(T)=',detT
            !Construct A matrix
            call setarray(A,nmodes*(nmodes+1)/2,zerodp)
            do i = 1,nmodes
                do j = 1,i
                    i2 = i*(i-1)/2+j
                    do k = 1,nmodes
                        A(i2) = A(i2) + 0.5*bcoef(k)*tmatrix(k,i)*tmatrix(k,j)
                    end do
                    if(i==j)A(i2) = A(i2) + 0.5*acoef(i)
                end do
            end do
            !Construct A-1 matrix
            call GIVENS(nmodes,nmodes,nmodes,A,Ascr,A2d,Avec)
            detA = 1.
            do i = 1,nmodes
                detA = detA*abs(A2d(i))
            end do
            call setmatrix(Ainv,nmodes,nmodes,zerodp)
            do i = 1,nmodes
                Ainv(i,i) = 1/A2d(i)
            end do
            call EBC(A1,Avec,Ainv,nmodes,nmodes,nmodes)
            call EBCT(Ainv,A1,Avec,nmodes,nmodes,nmodes)
            !Construct A1, A2, and A3 matrices
            call setmatrix(A2,nmodes,nmodes,zerodp)
            call setmatrix(A3,nmodes,nmodes,zerodp)
            do i = 1,nmodes
                do j = 1,nmodes
                    A1(i,j) = Ainv(i,j)*Sqrt(acoef(i))*Sqrt(acoef(j))
                    do k = 1,nmodes
                        A2(i,j) = A2(i,j) + Ainv(i,k)*tmatrix(j,k)
                        do l = 1,nmodes
                            A3(i,j) = A3(i,j) + tmatrix(i,k)*Ainv(k,l)*tmatrix(j,l)
                        end do
                    end do
                    A2(i,j) = A2(i,j)*Sqrt(acoef(i))*Sqrt(bcoef(j))
                    A3(i,j) = A3(i,j)*Sqrt(bcoef(i))*Sqrt(bcoef(j))
                end do
            end do
            !Initialize terms used in recursion
            call setarray(A2d,nmodes,zerodp)
            call setarray(A3d,nmodes,zerodp)
            do i = 1,nmodes
                A1(i,i) = A1(i,i) - 1.
                do j = 1,nmodes
                    d1 = A3(j,i)
                    if(j==i) d1 = d1 - 2. 
                    A2d(i) = A2d(i) + A2(i,j)*dvec(j)*Sqrt(bcoef(j))
                    A3d(i) = A3d(i) + d1*dvec(j)*Sqrt(bcoef(j))
                end do
            end do
            !Initialize indices
            pstart = 1
            do
                if(nfunc(pstart)>1) exit
                pstart = pstart + 1
            end do
            cstart = pstart + 1
            !C00 = 1.
            if(myid==0)scr2(1,1)=1.
            call GA_SYNC()
        !Determine <0|n> for each n in the basis
            allocate(Cn(lbuf))
            allocate(Cnbuf(lbuf))
            allocate(labels(nmodes*lbuf))
            allocate(Cm(nmodes*lbuf))
            allocate(Cmbuf(nmodes*lbuf))
            allocate(nlterms(nlev+1))
            do i = 0,(nlev-1)
                call setintarray(bvec,nmodes,one)  
                nlterms(i+1) = 0
                ndone = 0
                pindex = pstart
                cindex = cstart + 1
                bvec(pindex) = bvec(pindex) + 1
                bvec(cindex) = bvec(cindex) - 1
                do
                    if(cindex>nmodes) exit
                    if(bvec(cindex)==nfunc(cindex).or.bvec(pindex)==1) then
                        cindex = cindex + 1
                        if(nfunc(cindex-1)>1)pindex = cindex - 1
                    else
                        bvec(pindex) = bvec(pindex) - 1
                        bvec(cindex) = bvec(cindex) + 1
                        if(cindex>cstart) then
                        do j = 1,cindex-1
                         bvec(j) = 1
                        end do
                        ndo = 0
                        do j = cindex,nmodes
                         ndo = ndo + bvec(j) - 1
                        end do
                        ndo = i - ndo
                        sindex = 1
                        npack = 0
                        do
                            if(npack==ndo) exit
                            if(bvec(sindex).lt.nfunc(sindex)) then
                                bvec(sindex) = bvec(sindex) + 1
                                npack = npack + 1
                            else
                                sindex = sindex + 1  
                            end if
                        end do
                        pindex = pstart
                        cindex = cstart 
                        end if
                        stindex = stateindex(nmodes,bvec,shft)
                        if(stindex>=lo(1).and.stindex<=hi(1))then
                            Cn(ndone+1) = scr2(stindex-lo(1)+1,1)
                        else
                            Cn(ndone+1) = 0d0
                        end if
                        do j = 1,nmodes
                            l = irrep(j)
                            d1 = 0.
                            do k = 1,npirrep(l)
                                i2 = stindex - shft(k+irroff(l))
                                if(i2>=lo(1).and.i2<=hi(1)) d1 = d1 + A1(j,k+irroff(l))*(bvec(k+irroff(l))-1)*scr2(i2-lo(1)+1,1)
                            end do
                            Cm(ndone*nmodes+j) = d1
                            labels(ndone*nmodes+j) = stindex + shft(j)
                            if(bvec(j)==nfunc(j)) labels(ndone*nmodes+j) = 0       
                        end do
                        ndone = ndone + 1
                    end if
                    if(ndone>=lbuf)then
                        call GA_SYNC()
                        call MPI_ALLREDUCE(Cn,Cnbuf,ndone,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
                        call MPI_ALLREDUCE(Cm,Cmbuf,ndone*nmodes,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
                        do j = 1,ndone
                            do k = 1,nmodes
                                i2 = labels((j-1)*nmodes+k)
                                if(i2>=lo(1).and.i2<=hi(1)) scr2(i2-lo(1)+1,1) =  2d0*Cmbuf((j-1)*nmodes+k) - A2d(k)*Cnbuf(j)
                            end do
                        end do
                        nlterms(i+1) = nlterms(i+1) + ndone
                        ndone = 0
                    end if
                end do
                call GA_SYNC()
                call MPI_ALLREDUCE(Cn,Cnbuf,ndone,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
                call MPI_ALLREDUCE(Cm,Cmbuf,ndone*nmodes,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
                do j = 1,ndone
                    do k = 1,nmodes
                        i2 = labels((j-1)*nmodes+k)
                        if(i2>=lo(1).and.i2<=hi(1)) scr2(i2-lo(1)+1,1) = 2.*Cmbuf((j-1)*nmodes+k) - A2d(k)*Cnbuf(j)
                    end do
                end do
                nlterms(i+1) = nlterms(i+1) + ndone
                if(myid==0)write(*,1003)i+1,nlev
            end do
            tcount = 0
            do j = 1,nlev
                tcount = tcount + nlterms(j)
                if(myid==0)write(OUTFILE,1004)j-1,nlterms(j)
            end do
            if(myid==0)write(OUTFILE,1004)nlev,1
            tcount = tcount + 1
            if(myid==0)write(OUTFILE,1005)tcount,dimen
            deallocate(labels)
            deallocate(nlterms)
            deallocate(Cn)
            deallocate(Cm)
            deallocate(Cnbuf)
            deallocate(Cmbuf)
        !Determine <m|n> by building up successive <m| states
        !    This assumes a single m in the initial basis is excited,
        !    more coding in necessary to make the recursions work for arbitrary <m|
            allocate(labels(lbuf*(2*nmodes+1)))
            allocate(Cn(lbuf))
            allocate(nlterms(nproc*nseg))
            imode = 1
            do
                if(istate(imode)>0.or.imode==nmodes) exit
                imode = imode + 1
            end do
            ninit = istate(imode)
            if(ninit>zero)then
                call NGA_PUT(q3,lo,hi,scr2,ld,dble(1))
                call setarray(scr2,ld(1),nstates,zerodp)
                call setintarray(nlterms,nproc*nseg,zero)
            end if  
            do i = 1,ninit ! do over levels, m=0, m=1, etc.
                if(myid==0) write(OUTFILE,1011)i
                call get_barray(lo(1),bvec)
                bvec(nmodes) = bvec(nmodes)-1
                ldi(1) = min(lbuf,ld(1))
                pstart = lo(1)
                pindex = lo(1)+ldi(1)-1
                call GA_GET(q3,pstart,pindex,1,1,Cn,ldi)
                ndone = 0
                j = 0
                do
                    if(j==lbuf.or.ndone==ld(1))then
                        do k = 1,nproc*nseg
                        if(nlterms(k)>0)then
                            npack = vbounds(k,2) - vbounds(k,1) + 1
                            ndo = 0
                            tcount = 0
                            do
                                if(tcount==nlterms(k)) exit
                                ldi(1) = min(lbuf,npack-ndo)
                                pstart = vbounds(k,1)+ndo
                                pindex = vbounds(k,1)+ndo+ldi(1)-1
                                call GA_GET(q3,pstart,pindex,1,1,Cn,ldi)
                                do l = 1,j
                                    m = (l-1)*(2*nmodes+1)
                                    do n = 1,nmodes
                                        if(labels(m+2*n)>=pstart.and.labels(m+2*n)<=pindex)then
                                            tcount = tcount + 1
                                            sindex = labels(m+2*n) - pstart + 1
                                            scr2(labels(m+1),1) = scr2(labels(m+1),1) + 2.*A2(n,imode)*labels(m+2*n+1)*Cn(sindex)
                                        end if
                                    end do
                                end do
                                ndo = ndo + ldi(1)
                            end do
                        end if
                        end do
                        ldi(1) = min(lbuf,ld(1)-ndone)
                        pstart = lo(1)+ndone
                        pindex = lo(1)+ndone+ldi(1)-1
                        if(ldi(1)>0)call GA_GET(q3,pstart,pindex,1,1,Cn,ldi)
                        call setintarray(nlterms,nproc*nseg,zero)
                        j = 0
                    end if
                    if(ndone==ld(1)) exit
                    k = nmodes
                    do
                        if(bvec(k).lt.nfunc(k)) exit
                        bvec(k) = one
                        k = k - 1
                    end do
                    bvec(k) = bvec(k) + 1
                    !Determine the C(m-1,n) contribution
                    scr2(ndone+1,1) = 2.*scr2(ndone+1,1)*(i-1)*(A3(imode,imode)-1.)
                    !Determine the C(m,n) contribution
                    scr2(ndone+1,1) = scr2(ndone+1,1) - A3d(imode)*Cn(j+1)
                    labels(j*(2*nmodes+1)+1) = ndone+1
                    !Determine where to find C(m,n-1) contributions
                    do k = 1,nmodes
                        if(bvec(k)>1)then
                            labels(j*(2*nmodes+1)+2*k)   = lo(1) + ndone - shft(k)
                            labels(j*(2*nmodes+1)+2*k+1) = bvec(k) - 1
                            i2 = get_batch(lo(1) + ndone - shft(k))
                            nlterms(i2) = nlterms(i2) + 1
                        else
                            labels(j*(2*nmodes+1)+2*k)   = -1
                            labels(j*(2*nmodes+1)+2*k+1) = bvec(k) - 1
                        end if       
                    end do
                    ndone = ndone + 1
                    j = j + 1
                end do
                if(i.lt.ninit)then!If(i.ne.ninit), swap the m+1,m levels
                    ndone = 0
                    do
                        if(ndone==ld(1)) exit
                        ldi(1) = min(lbuf,ld(1)-ndone)
                        pstart = lo(1)+ndone
                        pindex = lo(1)+ndone+ldi(1)-1
                        call GA_GET(q3,pstart,pindex,1,1,Cn,ldi)
                        call GA_PUT(q3,pstart,pindex,1,1,scr2(ndone+1:ndone+ldi(1),1:1),ldi)
                        do j = 1,ldi(1)
                            scr2(ndone+j,1) = Cn(j)
                        end do
                        ndone = ndone + ldi(1)
                    end do
                end if
            end do
        !Determine the appropriate scale factors to go from C(m,n) to <m|n>, and normalization of the Lanczos vector
            call get_barray(lo(1),bvec)
            bvec(nmodes) = bvec(nmodes)-1
            !Compute G
                d1 = 0.
                d2 = 1.
                do i = 1,nmodes
                    d1 = d1 - 0.5*bcoef(i)*(dvec(i)**2)
                end do
                do i = 1,nmodes
                    do j = 1,nmodes
                        d1 = d1 + 0.25*Sqrt(bcoef(i))*dvec(i)*A3(i,j)*dvec(j)*Sqrt(bcoef(j))
                    end do
                end do
                do i = 1,nmodes
                    d2 = d2*bcoef(i)*acoef(i)
                end do
                d2 = Sqrt(Sqrt(d2)*detT/detA)*exp(d1)
            !Determine state dependent scale factor
                do i = 1,ld(1)
                    j = nmodes
                    do
                        if(bvec(j).lt.nfunc(j)) exit
                        bvec(j) = one
                        j = j - 1
                    end do
                    bvec(j) = bvec(j) + 1
                    d1 = 1.
                    do j = 1,nmodes
                        d1 = d1*(2.**(bvec(j)-1))*fac(bvec(j)-1)*(2.**istate(j))*fac(istate(j))
                    end do
                    scr2(i,1) = d2*scr2(i,1)/Sqrt(d1)
                    do j = 2,nstates
                        scr2(i,j) = scr2(i,1)
                    end do
                end do
                call NGA_PUT(q3,lo,hi,scr2,ld,dble(1))
    else
        d1 = 1d0
        stindex = 1 
        if(myid==0)write(OUTFILE,1000)
        call GA_ZERO(q3)
        call GA_SYNC()
        call GA_FILL_PATCH(q3,stindex,stindex,1,nstates,d1)
    end if
    call GA_SYNC()
    d1 = Sqrt(GA_DDOT(q3,q3)/nstates)
    if(myid==0)then
        write(OUTFILE,1006)d1
        write(OUTFILE,1007)
    end if
    !Print out dominant contributions to initial vector
    call GA_ELEM_MULTIPLY(q3,q3,q1)
    sindex = 0
    do
        call NGA_SELECT_ELEM(q1,'max',d2,irrep)
        if(sindex==40.or.d2==0.) exit
        sindex = sindex + 1
        call GA_SYNC()
        call GA_FILL_PATCH(q1,irrep(1),irrep(1),1,nstates,zerodp)
        call get_barray(irrep(1),bvec)
        do i = 1,nmodes
            bvec(i) = bvec(i)-1
        end do
        if(myid==0) then
            write(OUTFILE,1008,advance='no')sindex
            write(OUTFILE,1009,advance='no')bvec
            write(OUTFILE,1010)d2
        end if
    end do
    !Scale each state vector according to the weights array
    d2=0.D0
    do i = 1,nstates
        d2 = d2 + statew(i)
    end do
    if(d2==0.)d2 = 1.D0
    do i = 1,nstates
        call GA_SCALE_PATCH(q3,1,dimen,i,i,statew(i)/(d2*d1))
    end do   
    !Do final normalization of initial Lanczos vector
    if(myid==0)write(OUTFILE,1022)
    d2 = 0.D0
    do i = 1,nstates
        d1 = Sqrt(GA_DDOT_PATCH(q3,'N',1,dimen,i,i,q3,'N',1,dimen,i,i))
        d2 = d2 + d1**2
        if(myid==0)write(OUTFILE,1023)i,d1,i,d1**2
    end do
    if(myid==0)then
        write(OUTFILE,1024)Sqrt(d2),d2
        write(OUTFILE,1025)
    end if
    call GA_SCALE(q3,1.D0/Sqrt(d2))
    d1 = GA_DDOT(q3,q3)
    if(myid==0)write(OUTFILE,1024)Sqrt(d1),d1
    call system_clock(totstart,tcount,tmx)
    if(myid==0)write(OUTFILE,1002)1.*(totstart-totend)/tcount
    !format
        1000 format(/,'  ',70('*'),/,'  Generating Seed Vector using ANION basis...')
        1001 format(/,'  ',70('*'),/,'  Generating Seed Vector by computing ANION/ORIGIN overlap...')
        1002 format(/,'  Time Required: ',f14.3,' secs.',/,'  ',70('*'),/)
        1003 format('  COMPUTED LEVEL ',i6,' of ',i6)
        1004 format('  LEVEL ',i5,', Number of terms:     ',i12)
        1005 format(/,'  Sum over levels:                  ',i12,/,      &
                      '  Size of direct product basis:     ',i12,/)
        1006 format(/,'  Overlap with the initial state:   ',f12.6)
        1007 format(/,'  Dominant Contributions to Overlap',/'  ---------------------------------')
        1008 format('  ',i2,'. ')
        1009 format(18(i2,','))
        1010 format(': ',f12.6)
        1011 format('  Building <m| = <',i2,'|')
        1020 format('ID=',i2,' INDX=',i6,' BVEC=',15(i3))
        1021 format('ID=',i2,' NORMINIT=',f15.10)
        1022 format(/,'  Normalization of Initial Lanczos Vector')
        1023 format('  ||STATE',i2,'|| = ',f6.4,' -- (||STATE',i2,'||)^2 = ',f6.4)
        1024 format('  || TOTAL || = ',f6.4,' -- (|| TOTAL ||)^2 = ',f6.4)
        1025 format('  Normalizing Initial Lanczos Vector...')
end subroutine generate_initialvec

!Returns the basis array corresponding to the state label, stateindex
subroutine get_barray(stateindex,barray)
    use progdata, only: nmodes,shft
    implicit none
    integer,intent(in)                     :: stateindex
    integer,dimension(nmodes),intent(inout):: barray
    integer                                :: i,stlabel,tmp
    stlabel = stateindex - 1
    if(stlabel.lt.0)stlabel=0
    do i = 1,nmodes
        barray(i) = int(stlabel/shft(i)) + 1
        stlabel = stlabel - (barray(i)-1)*shft(i)
    end do
end subroutine get_barray

subroutine get_keywords()
    use progdata
    use filedata, only: BASISFILE
    implicit none
    integer,dimension(10)          :: npirr
    integer,dimension(100)         :: basis,initstate
    integer                        :: i,itmp,restart,bconv,sodata,shiftref,get_restart_iter,reorthog
    real*8,dimension(10) :: weights
    NAMELIST /NADVIBS/         niter,natoms,nmodes,nstates,basis,restart,bconv,idroots,soroots,reorthog, &
                               chkorthog,nseg,ztoler,maxdisk,weights,shiftref,nirreps,npirr,ordr,initstate
    call setintarray(basis,int(100),one)
    call setintarray(initstate,int(100),zero)
    call setintarray(npirr,int(10),zero)
    call setarray(weights,int(10),zerodp)
    npirr(1) = nmodes
    ordr = 2
    nirreps = 1
    niter = 1
    natoms = 2
    nmodes = 1
    nstates = 1
    restart = 0
    bconv = 1
    idroots = 0
    soroots = 0
    reorthog = 0
    chkorthog = 100
    nseg = 1
    ztoler = 1.D-20
    maxdisk = 1000
    shiftref = 0
    dimen = 1
    neworigin  = .false.
    restartrun = .false.
    savevecs   = .false.
    orthog     = .false.
    orthogexact = .false.
    iiter = 1
    open(unit=BASISFILE,file='basis.in',access='sequential',form='formatted',status='old')
    read(unit=BASISFILE,NML=NADVIBS)
    close(unit=BASISFILE)
    if(sum(npirr).ne.nmodes)then
        write(*,1000)sum(npirr),nmodes
        nirreps = 1
        npirr(1) = nmodes
    end if
    allocate(npirrep(nirreps))
    do i = 1,nirreps
        npirrep(i) = npirr(i)
    end do
    if(reorthog>0) orthog = .true.
    if(reorthog>1) orthogexact = .true.
    idroots = abs(idroots)
    soroots = abs(soroots)
    if(restart.ne.0)then
     restartrun = .true.
     iiter = get_restart_iter()+1
     niter = niter + iiter - 1
    end if
    if(shiftref.ne.0)neworigin = .true.
    if(orthog.or.idroots>0.or.soroots>0)savevecs = .true.
    bjiconv = 10.**(-bconv)
    allocate(statew(nstates))
    allocate(nfunc(nmodes))
    allocate(istate(nmodes))
    do i = 1,nstates
        statew(i) = weights(i)
    end do
    do i = 1,nmodes
        istate(i) = initstate(i)
        nfunc(i)  = basis(i)
        dimen     = dimen*nfunc(i)
    end do
    !format
        1000 format(' SYMMETRY CONFLICT, ',i5,' != ',i5,' --> setting nirrep=1')
end subroutine get_keywords

!Extract all the information from a restart.log file
subroutine load_restartinfo(reloadall)
    use progdata, only: alpha,beta,omega,orthog
    use filedata, only: RESTARTFILE
    implicit none
    logical, intent(inout)::reloadall
    CHARACTER(75)::commentline
    CHARACTER(80)::command
    CHARACTER(20)::filename
    CHARACTER(8) ::searchterm,currterm
    integer::i,nload
    real*8::dpval
    open(unit=RESTARTFILE,file='restart.log',status='old')
    read(unit=RESTARTFILE,fmt='(a75)')commentline
    read(unit=RESTARTFILE,fmt='(i5)')i
    read(unit=RESTARTFILE,fmt='(a75)')commentline
    read(unit=RESTARTFILE,fmt='(i5)')nload
    reloadall = i==nload
    read(unit=RESTARTFILE,fmt='(a75)')commentline
    read(unit=RESTARTFILE,fmt=*)(alpha(i),i=1,nload)
    searchterm=' BETA CO'
    do 
        read(unit=RESTARTFILE,fmt=1001,end=10,ERR=10)currterm
        if(currterm==searchterm) exit
    end do
    read(unit=RESTARTFILE,fmt=*)(beta(i),i=0,nload)
    if(orthog) then
        searchterm=' OMEGA(1'
        do 
            read(unit=RESTARTFILE,fmt=1001,end=10,ERR=10)currterm
            if(currterm==searchterm) exit
        end do
        read(unit=RESTARTFILE,fmt=*)(omega(1,i),i=1,nload+2)
        searchterm=' OMEGA(2'
        do 
            read(unit=RESTARTFILE,fmt=1001,end=10,ERR=10)currterm
            if(currterm==searchterm) exit
        end do
        read(unit=RESTARTFILE,fmt=*)(omega(2,i),i=1,nload+2)
    end if   
    close(unit=RESTARTFILE)
    filename = 'restart.log'
    command = 'mv -f '//trim(filename)//' '//trim(filename)//'.old'
    call system(command)
    !goto and format
        10 STOP 'ERROR in load_restartinfo()'
        1001 format(a8)
end subroutine load_restartinfo

!Constructs the tridiagonal matrix T from alpha and beta. If n is less than the current iteration,
!make_tridiag will construct from the bottom right of the matrix up -- so the most current alpha and betas are always included
!Input:  n -- dimension of T to construct
!Output: T -- tridiagonal matrix, packed by columns in upper triangular form
subroutine make_tmatrix(T,n)
    use progdata, only: alpha,beta,maxstor,zerodp
    use filedata, only: OUTFILE
    implicit none
    integer,intent(in)::n
    real*8,intent(inout),dimension(n*(n+1)/2)::T
    integer::i,j,k,counter
    counter = 0
    do i = 1,n
        do j = 1,i
            counter = counter + 1
            if(j==i) then
                T(counter) = alpha(j)
            elseif(j==(i-1)) then
                T(counter) = beta(j)
            else
                T(counter) = 0.
            end if
        end do
    end do
    if(counter.ne.(n*(n+1)/2)) write(unit=OUTFILE,fmt=*)'POTENTIAL ERROR IN make_tmatrix - counter=',counter,', n=',n
end subroutine make_tmatrix

!Prints out a matrix to file
subroutine matrix_write(itape,n,m,a,wformat)
    integer, intent(in)                 	     :: n,itape,wformat
    real*8,intent(in),dimension(n,n) :: a
    integer                                    :: i,fstate
    assign 1000 to fstate
    if(wformat==1) assign 1001 to fstate
    do i=1,n
        write(unit=itape,fmt=fstate)i,(a(i,k),k=1,m)
    end do
    !format
        1000 format(1x,i3,6f10.5,/(4x,6f10.5))
        1001 format(1x,i3,6es14.5,/(4x,6es14.5))
end subroutine matrix_write

!Determines the total amount of memory the program will require in MB 
subroutine memory_test(umem,rmem)
    use progdata, only: dimen,nstates,nmodes,niter,nfunc,nstblks,orthog,maxstor,maxdisk,  &
                        restartrun,myid,nproc,nseg,ordr,noterms,nztrms,nsteps
    implicit none
    real*8,intent(in)   ::umem
    real*8,intent(inout)::rmem
    integer::i,i1,i2,vsize,lchunk
    real*8::ndoubles,nintegers
    ndoubles = 0d0
    nintegers = 0d0
    lchunk = Ceiling(1.*dimen/nproc)+100 ! To account for slight differences in vector lengths
    vsize = dimen*nstates
    maxstor = int(1024*1024*nproc*maxdisk/(8*vsize)) - 3
    !First we determine the memory required to hold the program data
    ndoubles = ndoubles + 2*niter       ! for alpha, beta
    ndoubles = ndoubles + niter*(niter+1)/2 !Tmat
    ndoubles = ndoubles + niter*niter   !Ymat
    ndoubles = ndoubles + 8*niter       !yval,bji,Tvals,escr
    ndoubles = ndoubles + 2*(niter+1)   !omega
    ndoubles = ndoubles + niter         !dpvec
    ndoubles = ndoubles + nmodes*nmodes !Tmat
    ndoubles = ndoubles + nmodes        !dvec
    ndoubles = ndoubles + 2*nmodes      !aomega,bomega
    !Storage required for lanczos vectors
    ndoubles = ndoubles + 3*nstates*lchunk + 100 ! global arrays memory
    !Storage required for scratch
    ndoubles = ndoubles + nstates*lchunk + nstates*Ceiling(1.*lchunk/nseg) ! scr1 and scr2
    nintegers = nintegers + nmodes !for shft
    nintegers = nintegers + nmodes !for nfunc
    nintegers = nintegers + Ceiling(2.*niter/100.) !loindex,roindex
    i1 = sum(nztrms)
    i2 = sum(nsteps)
    nintegers = nintegers + i1*nstates*nstates     ! nzcoef
    nintegers = nintegers + i1*2*nstates*nstates+1 ! nzblks
    nintegers = nintegers + i1*2*ordr              ! nzindx
    nintegers = nintegers + i1*ordr              ! bmax
    nintegers = nintegers + i1                   ! nsteps
    nintegers = nintegers + 2*i2*ordr            ! stype,basdif
    nintegers = nintegers + nproc*nseg*i2        ! nelems
    nintegers = nintegers + 2*nproc*nseg*i2*ordr ! basind, bstart
    nintegers = nintegers + 2*nproc*nseg         ! rcshft
    rmem = (ndoubles*8d0+nintegers*8d0)/(1024d0*1024d0) ! in MB
end subroutine memory_test

!Given a full symmetric matrix "a" packed by columns, return "m" in upper triangular form, packed by columns
subroutine totriang(n,a,m)
    integer,intent(in)					:: n
    integer						:: i,j,ij,mindex
    real*8,dimension(n*n),intent(in)		:: a
    real*8,dimension(n*(n+1)/2),intent(inout)	:: m
    mindex=1
    do i = 1,n
        do j = 1,i
            ij = i+n*(j-1)
            m(mindex) = a(ij)
            mindex =  mindex + 1
        end do
    end do
end subroutine totriang

subroutine write_scrfile(n)
    use progdata, only: alpha,beta
    use filedata, only: SCRFILE
    implicit none
    integer,intent(in)  :: n
    integer             :: i
    rewind(unit=SCRFILE)
    write(unit=SCRFILE,fmt=1000)n-1
    write(unit=SCRFILE,fmt='(a)')' ALPHA COEFFICIENTS'
    write(unit=SCRFILE,fmt=1001)alpha(1:n-1)
    write(unit=SCRFILE,fmt='(a)')' BETA COEFFICIENTS'
    write(unit=SCRFILE,fmt=1001)beta(0:n-1)
    !format
        1000 format('COEFFICIENT INFORMATION AFTER ',i4, ' ITERATIONS')
        1001 format(4(f18.14))
end subroutine write_scrfile

!Creates an n-length array, initializes elements to val
subroutine setarray(a,n,val)
    integer,intent(in)                          :: n
    real*8,intent(in)                 :: val
    real*8,intent(inout),dimension(n) :: a
    a=val
end subroutine setarray

!Creates an m x n matrix, intializes elements to val
subroutine setmatrix(a,m,n,val)
    integer, intent(in)                            :: m,n
    real*8, intent(in)                   :: val
    real*8, intent(inout),dimension(m,n) :: a
    integer 				         :: i,j
    a=val
end subroutine setmatrix

!Creates an n-length array, initializes elements to val
subroutine setintarray(a,n,val)
    integer, intent(in)                   :: n
    integer, intent(in)                   :: val
    integer, intent(inout),dimension(n)	:: a
    integer	                        :: i
    a=val
end subroutine setintarray

!Creates an m x n matrix, initializes elements to val
subroutine setintmatrix(a,m,n,val)
    integer, intent(in)                   :: m,n
    integer, intent(in)                   :: val
    integer, intent(inout),dimension(m,n)	:: a
    integer 			        :: i,j
    a=val
end subroutine setintmatrix

!Input:  nl order array list
!Output: nl order array uni, uni(1:nu) contains the unique elements in list
subroutine union(nl,list,nu,uni)
    integer, intent(in)                   :: nl
    integer, dimension(nl), intent(in)    :: list
    integer, intent(inout)                :: nu
    integer, dimension(nl), intent(inout) :: uni            
    integer                               :: i,j,k
    nu = 0
    do i = 1,nl
        k = 1
        do j = 1,nu
            if(list(i)==uni(j))then
                k=0
                exit
            end if
        end do
        if(k==1)then
            nu = nu + 1
            uni(nu) = list(i)
        end if
    end do
end subroutine union
!------------------------------------ End ------------------------------------

!---------------- All functions listed in alphabetical order -----------------
    !Exact factorial for N<=20
    !Double precision for 21<=N<=40
    !8 significant figures for N>=41
    real*8 function fac(N)
        implicit none
        integer,intent(in)::N
        integer::i
        real*8::temp
        select case(N)
            case(0)
                fac=1d0
            case(1)
                fac=1d0
            case(2)
                fac=2d0
            case(3)
                fac=6d0
            case(4)
                fac=24d0
            case(5)
                fac=120d0
            case(6)
                fac=720d0
            case(7)
                fac=5040d0
            case(8)
                fac=40320d0
            case(9)
                fac=362880d0
            case(10)
                fac=3628800d0
            case(11)
                fac=39916800d0
            case(12)
                fac=479001600d0
            case(13)
                fac=6227020800d0
            case(14)
                fac=87178291200d0
            case(15)
                fac=1307674368d3
            case(16)
                fac=20922789888d3
            case(17)
                fac=355687428096d3
            case(18)
                fac=6402373705728d3
            case(19)
                fac=121645100408832d3
            case(20)
                fac=243290200817664d4
            case(21)
                fac=5109094217170944d4
            case(22)
                fac=112400072777760768d4
            case(23)
                fac=2585201673888497664d4
            case(24)
                fac=62044840173323943936d4
            case(25)
                fac=15511210043330985984d6
            case(26)
                fac=403291461126605635584d6
            case(27)
                fac=10888869450418352160768d6
            case(28)
                fac=304888344611713860501504d6
            case(29)
                fac=8841761993739701954543616d6
            case(30)
                fac=26525285981219105863630848d7
            case(31)
                fac=822283865417792281772556288d7
            case(32)
                fac=26313083693369353016721801216d7
            case(33)
                fac=868331761881188649551819440128d7
            case(34)
                fac=29523279903960414084761860964352d7
            case(35)
                fac=103331479663861449296666513375232d8
            case(36)
                fac=3719933267899012174679994481508352d8
            case(37)
                fac=137637530912263450463159795815809024d8
            case(38)
                fac=5230226174666011117600072241000742912d8
            case(39)
                fac=203978820811974433586402817399028973568d8
            case(40)
                fac=815915283247897734345611269596115894272d9
            case default
                temp=(9d0*N)**3.141592653589793d0
                fac=Sqrt(6.283185307179586d0*N)*(N/2.718281828459045d0)**N*Exp(1d0/12d0/N-Log(9d0*N)/(temp-1d0/temp))
        end select
    end function fac
    
    integer function get_procid(index)
        use progdata,only:myid,vbounds,nproc,nseg
        implicit none
        integer,intent(inout)::index
        get_procid = 0
        do
            if(index.ge.vbounds(get_procid*nseg+1,1).and.index.le.vbounds(get_procid*nseg+nseg,2).or.get_procid.eq.nproc) exit
            get_procid = get_procid + 1
        end do
        if(get_procid.eq.nproc) then
            print *,'ID=',myid,' ERROR determining batch, get_procid=',get_procid,' INDEX=',index
            get_procid = 0
        end if
    end function get_procid
    
    integer function get_batch(index)
        use progdata,only:myid,vbounds,nproc,nseg
        implicit none
        integer,intent(inout)::index
        integer::batmax
        batmax=nproc*nseg
        get_batch=1
        do
            if((index.ge.vbounds(get_batch,1).and.index.le.vbounds(get_batch,2)).or.get_batch.gt.batmax)EXIT
            get_batch=get_batch+1
        end do
        if(get_batch.gt.nproc*nseg)then
            print *,'ID=',myid,' ERROR determining batch, get_batch=',get_batch,' INDEX=',index
            get_batch=0
        end if
    end function get_batch
    
    !Get the number of iterations run on previous from restart.log file
    integer function get_restart_iter()
        use filedata,only:RESTARTFILE  
        implicit none
        character*75::commentline
        integer::dummy
        open(unit=RESTARTFILE,file='restart.log',status='old')
            rewind(unit=RESTARTFILE)
            read(unit=RESTARTFILE,fmt='(a75)')commentline
            read(unit=RESTARTFILE,fmt='(i5)')dummy
            read(unit=RESTARTFILE,fmt='(a75)')commentline
            read(unit=RESTARTFILE,fmt='(i5)')get_restart_iter
        close(unit=RESTARTFILE)
    end function get_restart_iter
    
    !Returns the Nth Hermite polynomial evaluated at x
    real*8 function HermiteN(n,x)
        implicit none
        integer,intent(in)::n
        real*8,intent(in)::x
        integer::i
        real*8,dimension(3)::Hn
        Hn(1)=0d0
        Hn(2)=1d0
        do i = 0,n-1,1
            Hn(3)=2d0*x*Hn(2)-2d0*i*Hn(1)
            Hn(1)=Hn(2)
            Hn(2)=Hn(3)
        end do
        HermiteN=Hn(2)
    end function HermiteN
    
    !Returns the Nth Hermite polynomial evaluated at x, where x is imaginary
    real*8 function HermiteNI(n,x)
        implicit none
        integer,intent(in)::n
        real*8,intent(in)::x
        integer::k
        real*8::fac
        HermiteNI=0d0
        do k=0,n/2
            HermiteNI=HermiteNI+((-1)**(int((n-2.*k)/2.)))*((-1)**k)*fac(n)*((2*x)**(n-2*k))/(fac(k)*fac(n-2*k))
        end do
    end function HermiteNI
    
    !Return Sum[(ni+1/2)*wi,{i,NVIBS}]
    real*8 function href(nvals)
        use progdata, only: nmodes,aomega  
        implicit none
        integer,dimension(nmodes),intent(in)::nvals
        integer::i
        href=0d0
        do i=1,nmodes
            href=href+(1d0*nvals(i)-0.5d0)*aomega(i)
        end do
    end function href
    
    real*8 function hoelem(inttype,a,b)
        implicit none
        integer,intent(in)::inttype,a,b
        real*8::m,n
        !Basis functions are indexed 1 - nfunc
        !However, n/m runs from 0 - (nfunc-1)
        !Thus, the below recursion formulae employ (n/m - 1)
        m=a-1d0
        n=b-1d0
        hoelem=0d0
        select case(inttype)
        case(0)
            if(b.eq.a) then
                hoelem=1d0
            else
                hoelem=0d0
            end if
        case(1)
            if(b.eq.(a+1)) then
                hoelem=Sqrt((m+1)/2)
            elseif(b.eq.(a-1)) then
                hoelem=Sqrt(m/2)
            else
                hoelem=0d0
            end if
        case(2)
            if(b.eq.(a+2)) then
                hoelem=Sqrt((m+1)*(m+2)/4)
            elseif(b.eq.a) then
                hoelem=((2*m+1)/2)
            elseif(b.eq.(a-2)) then
                hoelem=Sqrt(m*(m-1)/4)
            else
                hoelem=0d0
            end if
        case(3)
            if(b.eq.(a+3)) then
                hoelem=Sqrt((m+1)*(m+2)*(m+3)/8)
            elseif(b.eq.(a+1)) then
                hoelem=Sqrt(((m+1)**3)/8)
            elseif(b.eq.(a-1)) then
                hoelem=Sqrt((m**3)/8)
            elseif(b.eq.(a-3)) then
                hoelem=Sqrt(m*(m-1)*(m-2)/8)
            else
                hoelem=0d0
            end if
        case(4)
            if(b.eq.(a+4)) then
                hoelem=Sqrt((m+1)*(m+2)*(m+3)*(m+4))/4
            elseif(b.eq.(a+2)) then
                hoelem=(2*m+3)*Sqrt((m+1)*(m+2))/2
            elseif(b.eq.a) then
                hoelem=3*(2*m*m+2*m+1)/4
            elseif(b.eq.(a-2)) then
                hoelem=(2*m-1)*Sqrt(m*(m-1))/2
            elseif(b.eq.(a-4)) then
                hoelem=Sqrt(m*(m-1)*(m-2)*(m-3))/4
            else
                hoelem=0d0
            end if
        case default
            hoelem=0d0
        end select
    end function hoelem
    
    !Return the current number of converged eigenvalues
    integer function nconverged(n)
        use progdata, only: bjiconv,beta
        implicit none
        integer,intent(in)::n
        integer::i
        real*8::eval,difeigs
        real*8,dimension(n*(n+1)/2)::T
        real*8,dimension(n,n)::eigvecs
        real*8,dimension(n)::eigvals
        real*8,dimension(5*n)::tbuf
        nconverged=0
        difeigs=1d-8
        eval=-1d0 
        call make_tmatrix(T,n)
        call givens(n,n,n,T,tbuf,eigvals,eigvecs)
        do i = 1,Ceiling(1.*n/2.)
            if((eigvals(i)-eval).gt.difeigs) then
                eval = eigvals(i)
                if(abs(beta(n)*eigvecs(n,i)).lt.1e-8) nconverged=nconverged+1
            end if
        end do
    end function nconverged
    
    !Returns the number of elements in the upper triangle of a n x n square matrix
    integer function numut(n,d)
        implicit none
        integer,intent(in)::n,d
        real*8::fac
        numut = anint(fac(n-1+d)/(fac(n-1)*fac(d)))
    end function numut
    
    real*8 function gauss_random(mean,sigma)
        implicit none
        real*8,intent(in)::mean,sigma
        real*8::r1,r2
        call random_number(r1)
        call random_number(r2)
        gauss_random=mean+sigma*dsqrt(-2d0*dlog(r1))*dcos(6.283185307179586d0*r2)
    end function gauss_random
    
    integer function stateindex(bl,barray,boff)
        integer,intent(inout)::bl
        integer,dimension(bl),intent(inout)::barray,boff
        integer::i
        stateindex=1
        do i=bl,1,-1
            stateindex=stateindex+(barray(i)-1)*boff(i)
        end do
    end function stateindex
!------------------------------------ End ------------------------------------

!---------- YARKONY ROUTINES: THESE MAY BE FOUND IN UTILITIES.COL.F ----------
    SUBROUTINE EVB(A,B,C,NB,NC)
    implicit none
    integer,intent(in)                          :: NB,NC
    real*8,dimension(NC),intent(inout):: A
    real*8,dimension(NB),intent(in)   :: B
    real*8,dimension(NB,NC),intent(in):: C
    integer                                     :: i,j
    
    do i = 1,NC
     A(i) = 0.
     do j = 1,NB
      A(i) = A(i) + B(j)*C(j,i)
     end do
    end do
    
    return
    end SUBROUTINE EVB
    !
    ! EBCT multiplies B and the transpose of C and puts the resulting matrix in A
    !
    SUBROUTINE EBCT(A,B,C,NI,NK,NJ)
    integer     			        :: NI,NK,NJ
    real*8,dimension(NI,NJ)	:: A
    real*8,dimension(NI,NK)	:: B
    real*8,dimension(NJ,NK)	:: C
    
    CALL setmatrix(A,NI,NJ,0)
    
    DO I=1,NI
       DO J=1,NJ
          T=0.0D+00
          DO K=1,NK
             T=T+B(I,K)*C(J,K)
          ENDDO
             A(I,J)=T
       ENDDO
    ENDDO
    
    RETURN
    END SUBROUTINE EBCT
    !
    ! EBTC multiplies the transpose of B by C and puts the reuslting matrix in A
    !
    SUBROUTINE EBTC(A,B,C,NI,NK,NJ)
      integer     			:: NI,NK,NJ
      real*8,dimension(NI,NJ)	:: A
      real*8,dimension(NK,NI)	:: B
      real*8,dimension(NK,NJ)	:: C
    
      CALL setmatrix(A,NI,NJ,0)
    
      DO I=1,NI
         DO J=1,NJ
            T=0.0D+00
            DO K=1,NK
               T=T+B(K,I)*C(K,J)
            ENDDO
            A(I,J)=T
         ENDDO
      ENDDO
    
      RETURN
    END SUBROUTINE EBTC
    !
    ! EBC multiplies B by C and puts the resulting matrix in A 
    !
    SUBROUTINE EBC(A,B,C,NI,NK,NJ)
      integer     			:: NI,NK,NJ
      real*8,dimension(NI,NJ)	:: A
      real*8,dimension(NI,NK)	:: B
      real*8,dimension(NK,NJ)	:: C
    
      CALL setmatrix(A,NI,NJ,0)
    
      DO I=1,NI
         DO J=1,NJ
            T=0.0D+00
            DO K=1,NK
               T=T+B(I,K)*C(K,J)
            ENDDO
            A(I,J)=T
         ENDDO
      ENDDO
      RETURN
    END SUBROUTINE EBC
    !
    !
    !
    SUBROUTINE GIVENS (NX,NROOTX,NJX,A,B,ROOT,VECT)
    !    real*8 VERSION BY MEL LEVY 8/72
    !62.3  GIVENS  -EIGENVALUES AND EIGENVECTORS BY THE GIVENS METHOD.
    !    BY FRANKLIN PROSSER, INDIANA UNIVERSITY.
    !    SEPTEMBER, 1967
    !    CALCULATES EIGENVALUES AND EIGENVECTORS OF REAL SYMMETRI!MATRIX
    !    STORED IN PACKED UPPER TRIANGULAR FORM.
    !
    !    THANKS ARE DUE TO F. E. HARRIS (STANFORD UNIVERSITY) AND H. H.
    !    MICHELS (UNITED AIRCRAFT RESEARCH LABORATORIES) FOR EXCELLENT
    !    WORK ON NUMERICAL DIFFICULTIES WITH EARLIER VERSIONS OF THIS
    !    PROGRAM.
    !
    !    THE PARAMETERS FOR THE ROUTINE ARE...
    !        NX     ORDER OF MATRIX
    !        NROOTX NUMBER OF ROOTS WANTED.  THE NROOTX SMALLEST (MOST
    !                NEGATIVE) ROOTS WILL BE CALCULATED.  IF NO VECTORS
    !                ARE WANTED, MAKE THIS NUMBER NEGATIVE.
    !        NJX    ROW DIMENSION OF VECT ARRAY.  SEE -VECT- BELOW.
    !                NJX MUST BE NOT LESS THAN NX.
    !        A      MATRIX STORED BY COLUMNS IN PACKED UPPER TRIANGULAR
    !               FORM, I.E. OCCUPYING NX*(NX+1)/2 CONSECUTIVE
    !               LOCATIONS.
    !        B      SCRATCH ARRAY USED BY GIVENS.  MUST BE AT LEAST
    !                NX*5 CELLS.
    !        ROOT   ARRAY TO HOLD THE EIGENVALUES.  MUST BE AT LEAST
    !               NROOTX CELLS LONG.  THE NROOTX SMALLEST ROOTS ARE
    !                ORDERED LARGEST FIRST IN THIS ARRAY.
    !        VECT   EIGENVECTOR ARRAY.  EACH COLUMN WILL HOLD AN
    !                EIGENVECTOR FOR THE CORRESPONDING ROOT.  MUST BE
    !                DIMENSIONED WITH -NJX- ROWS AND AT LEAST -NROOTX-
    !                COLUMNS, UNLESS NO VECTORS
    !                ARE REQUESTED (NEGATIVE NROOTX).  IN THIS LATTER
    !                CASE, THE ARGUMENT VECT IS JUST A DUMMY, AND THE
    !                STORAGE IS NOT USED.
    !                THE EIGENVECTORS ARE NORMALIZED TO UNITY.
    !
    !    THE ARRAYS A AND B ARE DESTROYED BY THE COMPUTATION.  THE RESULTS
    !    APPEAR IN ROOT AND VECT.
    !    FOR PROPER FUNCTIONING OF THIS ROUTINE, THE RESULT OF A FLOATING
    !    POINT UNDERFLOW SHOULD BE A ZERO.
    !    TO CONVERT THIS ROUTINE TO real*8 (E.G. ON IBM 360
    !    MACHINES), BE SURE THAT ALL REAL VARIABLES AND function
    !    REFERENCES ARE PROPERLY MADE real*8.
    !    THE VALUE OF -ETA- (SEE BELOW) SHOULD ALSO BE CHANGED, TO REFLECT
    !    THE INCREASED PRECISION.
    !
    !    THE ORIGINAL REFERENCE TO THE GIVENS TECHNIQUE IS IN OAK RIDGE
    !    REPORT NUMBER ORNL 1574 (PHYSICS), BY WALLACE GIVENS.
    !    THE METHOD AS PRESENTED IN THIS PROGRAM CONSISTS OF FOUR STEPS,
    !    ALL MODIFICATIONS OF THE ORIGINAL METHOD...
    !    FIRST, THE INPUT MATRIX IS REDUCED TO TRIDIAGONAL FORM BY THE
    !    HOUSEHOLDER TECHNIQUE (J. H. WILKINSON, COMP. J. 3, 23 (1960)).
    !    THE ROOTS ARE THEN LOCATED BY THE STURM SEQUENCE METHOD (J. M.
    !    ORTEGA (SEE REFERENCE BELOW).  THE VECTORS OF THE TRIDIAGONAL
    !    FORM ARE THEN EVALUATED (J. H. WILKINSON, COMP. J. 1, 90 (1958)),
    !    AND LAST THE TRIDIAGONAL VECTORS ARE ROTATED TO VECTORS OF THE
    !    ORIGINAL ARRAY (FIRST REFERENCE).
    !    VECTORS FOR DEGENERATE (OR NEAR-DEGENERATE) ROOTS ARE FORCED
    !    TO BE ORTHOGONAL, USING A METHOD SUGGESTED BY B. GARBOW, ARGONNE
    !    NATIONAL LABS (PRIVATE COMMUNICATION, 1964).  THE GRAM-SCHMIDT
    !    PROCESS IS USED FOR THE ORTHOGONALIZATION.
    !
    !    AN EXCELLENT PRESENTATION OF THE GIVENS TECHNIQUE IS FOUND IN
    !    J. M. ORTEGA-S ARTICLE IN -MATHEMATICS FOR DIGITAL COMPUTERS,-
    !    VOLUME 2, ED. BY RALSTON AND WILF, WILEY (1967), PAGE 94.
    !
        use progdata, only: epsilon
        IMPLICIT REAL*8 (A-H,O-Z)
        real*8 B(NX,5),A(1),ROOT(NROOTX),VECT(NJX,NROOTX)
        real*8 ETA,THETA,DEL1,DELTA,SMALL,DELBIG,THETA1,TOLER
        real*8 RPOWER,RPOW1,RAND1,FACTOR,ANORM,ALIMIT,SUM,TEMP
        real*8 AK,ROOTL,ROOTX,TRIAL,F0,SAVE,AROOT
        real*8 ELIM1,ELIM2
    !
    !** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !**  USERS PLEASE NOTE...
    !**  THE FOLLOWING TWO PARAMETERS, ETA AND THETA, SHOULD BE ADJUSTED
    !**  BY THE USER FOR HIS PARTICULAR MACHINE.
    !**  ETA IS AN INDICATION OF THE PRECISION OF THE FLOATING POINT
    !**  REPRESENTATION ON THE COMPUTER BEING USED (ROUGHLY 10**(-M),
    !**  WHERE M IS THE NUMBER OF DECIMALS OF PRECISION ).
    !**  THETA IS AN INDICATION OF THE RANGE OF NUMBERS THAT CAN BE
    !**  EXPRESSED IN THE FLOATING POINT REPRESENTATION (ROUGHLY THE
    !**  LARGEST NUMBER).
    !**  SOME RECOMMENDED VALUES FOLLOW.
    !**  FOR IBM 7094, UNIVA!1108, ETC. (27-BIT BINARY FRACTION, 8-BIT
    !**  BINARY EXPONENT), ETA=1.E-8, THETA=1.E37.
    !**  FOR CONTROL DATA 3600 (36-BIT BINARY FRACTION, 11-BIT BINARY
    !**  EXPONENT), ETA=1.E-11, THETA=1.E307.
    !**  FOR CONTROL DATA 6600 (48-BIT BINARY FRACTION, 11-BIT BINARY
    !**  EXPONENT), ETA=1.E-14, THETA=1.E307.
    !**  FOR IBM 360/50 AND 360/65 real*8 (56-BIT HEXADECIMAL
    !**  FRACTION, 7-BIT HEXADECIMAL EXPONENT), ETA=1.E-16, THETA=1.E75.
    !**
        THETA = huge(ELIM1)
        ETA = epsilon
    !** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
        DEL1=ETA/100.0D0
        DELTA=ETA**2*100.0D0
        SMALL=ETA**2/100.0D0
        DELBIG=THETA*DELTA/1000.0D0
        THETA1=1000.0D0/THETA
    !    TOLER  IS A FACTOR USED TO DETERMINE IF TWO ROOTS ARE CLOSE
    !    ENOUGH TO BE CONSIDERED DEGENERATE FOR PURPOSES OF ORTHOGONALI-
    !    ZING THEIR VECTORS.  FOR THE MATRIX NORMED TO UNITY, IF THE
    !    DIFFERENCE BETWEEN TWO ROOTS IS LESS THAN TOLER, THEN
    !    ORTHOGONALIZATION WILL OCCUR.
        TOLER = ETA*100.
    !
    !    INITIAL VALUE FOR PSEUDORANDOM NUMBER GENERATOR... (2**23)-3
        RPOWER=8388608.0D0
        RPOW1 = RPOWER/2.
        RAND1=RPOWER-3.0D0
    !
        N = NX
        NROOT = IABS(NROOTX)
        IF (NROOT.EQ.0) GO TO 1001
        IF (N-1) 1001,1003,105
    1003  ROOT(1) = A(1)
        IF(NROOTX .GT. 0)VECT(1,1)=1.0D0
        GO TO 1001
    105   CONTINUE
    !    NSIZE    NUMBER OF ELEMENTS IN THE PACKED ARRAY
        NSIZE = (N*(N+1))/2
        NM1 = N-1
        NM2 = N-2
    !
    !    SCALE MATRIX TO EUCLIDEAN NORM OF 1.  SCALE FACTOR IS ANORM.
        FACTOR=0.0D-0
        DO 70 I=1,NSIZE
    70   FACTOR=DMAX1(FACTOR,DABS(A(I)))
        IF(FACTOR .NE. 0.0D0)GO TO 72
    !    NULL MATRIX.  FIX UP ROOTS AND VECTORS, THEN EXIT.
        DO 78 I=1,NROOT
             IF (NROOTX.LT.0) GO TO 78
             DO 77 J=1,N
    77        VECT(J,I)=0.0D0
             VECT(I,I)=1.0D0
    78   ROOT(I)=0.0D0
        GO TO 1001
    !
    72   ANORM=0.0D0
        J = 1
        K = 1
        DO 80 I=1,NSIZE
             IF (I.NE.J) GO TO 81
             ANORM=ANORM+(A(I)/FACTOR)**2/2.0D0
             K = K+1
             J = J+K
             GO TO 80
    81         ANORM = ANORM + (A(I)/FACTOR)**2
    80    CONTINUE
        ANORM=DSQRT(ANORM*2.0D0)*FACTOR
        DO 91 I=1,NSIZE
    91    A(I) = A(I)/ANORM
        ALIMIT=1.0D-0
    !
    !    TRIDIA SECTION.
    !    TRIDIAGONALIZATION OF SYMMETRIC MATRIX
        ID = 0
        IA = 1
        IF (NM2.EQ.0) GO TO 201
        DO 200  J=1,NM2
    !    J       COUNTS ROW  OF A-MATRIX TO BE DIAGONALIZED
    !    IA      START OF NON-CODIAGONAL ELEMENTS IN THE ROW
    !    ID      INDEX OF CODIAGONAL ELEMENT ON ROW BEING CODIAGONALIZED.
             IA = IA+J+2
             ID = ID + J + 1
             JP2 = J+2
    !    SUM SQUARES OF NON-CODIAGONAL ELEMENTS IN ROW J
             II = IA
             SUM=0.0D-0
             DO 100 I=JP2,N
                  SUM = SUM + A(II)**2
    100        II = II + I
             TEMP = A(ID)
             IF (SUM.GT.SMALL) GO TO 110
    !    NO TRANSFORMATION NECESSARY IF ALL THE NON-CODIAGONAL
    !    ELEMENTS ARE TINY.
             B(J,1) = TEMP
             A(ID)=0.0D0
             GO TO 200
    !    NOW COMPLETE THE SUM OF OFF-DIAGONAL SQUARES
    110       SUM=DSQRT(SUM+TEMP**2)
    !    NEW CODIAGONAL ELEMENT
             B(J,1)=-DSIGN(SUM,TEMP)
    !    FIRST NON-ZERO ELEMENT OF THIS W-VECTOR
             B(J+1,2)=DSQRT((1.0D0+DABS(TEMP)/SUM)/2.0D0)
    !    FORM REST OF THE W-VECTOR ELEMENTS
             TEMP=DSIGN(0.5D0/(B(J+1,2)*SUM),TEMP)
             II = IA
             DO 130 I=JP2,N
                  B(I,2) = A(II)*TEMP
    130        II = II + I
    !    FORM P-VECTOR AND SCALAR.  P-VECTOR = A-MATRIX*W-VECTOR.
    !    SCALAR = W-VECTOR*P-VECTOR.
             AK=0.0D0
    !    IC      LOCATION OF NEXT DIAGONAL ELEMENT
             IC = ID + 1
             J1 = J + 1
             DO 190  I=J1,N
                  JJ = IC
                  TEMP=0.0D0
                  DO 180  II=J1,N
    !    I       RUNS OVER THE NON-ZERO P-ELEMENTS
    !    II      RUNS OVER ELEMENTS OF W-VECTOR
                       TEMP = TEMP + B(II,2)*A(JJ)
    !    CHANGE INCREMENTING MODE AT THE DIAGONAL ELEMENTS.
                       IF (II.LT.I) GO TO 210
                       JJ = JJ + II
                       GO TO 180
    210                  JJ = JJ + 1
    180             CONTINUE
    !    BUILD UP THE K-SCALAR (AK)
                  AK = AK + TEMP*B(I,2)
                  B(I,1) = TEMP
    !    MOVE I!TO TOP OF NEXT A-MATRIX -ROW-
    190        IC = IC + I
    !    FORM THE Q-VECTOR
             DO 150  I=J1,N
    150        B(I,1) = B(I,1) - AK*B(I,2)
    !    TRANSFORM THE REST OF THE A-MATRIX
    !    JJ      START-1 OF THE REST OF THE A-MATRIX
             JJ = ID
    !    MOVE W-VECTOR INTO THE OLD A-MATRIX LOCATIONS TO SAVE SPACE
    !    I       RUNS OVER THE SIGNIFICANT ELEMENTS OF THE W-VECTOR
             DO 160  I=J1,N
                  A(JJ) = B(I,2)
                  DO 170  II=J1,I
                       JJ = JJ + 1
    170            A(JJ)=A(JJ)-2.0D0*(B(I,1)*B(II,2)+B(I,2)*B(II,1))
    160        JJ = JJ + J
    200   CONTINUE
    !    MOVE LAST CODIAGONAL ELEMENT OUT INTO ITS PROPER PLACE
    201   CONTINUE
        B(NM1,1) = A(NSIZE-1)
        A(NSIZE-1)=0.0D-0
    !
    !    STURM SECTION.
    !    STURM SEQUENCE ITERATION TO OBTAIN ROOTS OF TRIDIAGONAL FORM.
    !    MOVE DIAGONAL ELEMENTS INTO SECOND N ELEMENTS OF B-VECTOR.
    !    THIS IS A MORE CONVENIENT INDEXING POSITION.
    !    ALSO, PUT SQUARE OF CODIAGONAL ELEMENTS IN THIRD N ELEMENTS.
        JUMP=1
        DO 320 J=1,N
             B(J,2)=A(JUMP)
             B(J,3) = B(J,1)**2
    320   JUMP = JUMP+J+1
        DO 310 I=1,NROOT
    310   ROOT(I) = +ALIMIT
        ROOTL = -ALIMIT
    !    ISOLATE THE ROOTS.  THE NROOT LOWEST ROOTS ARE FOUND, LOWEST FIRST
        DO 330 I=1,NROOT
    !    FIND CURRENT -BEST- UPPER BOUND
             ROOTX = +ALIMIT
             DO 335 J=I,NROOT
    335       ROOTX=DMIN1(ROOTX,ROOT(J))
             ROOT(I) = ROOTX
    !    GET IMPROVED TRIAL ROOT
    500       TRIAL=(ROOTL+ROOT(I))*0.5D0
             IF (TRIAL.EQ.ROOTL.OR.TRIAL.EQ.ROOT(I)) GO TO 330
    !    FORM STURM SEQUENCE RATIOS, USING ORTEGA-S ALGORITHM (MODIFIED).
    !    NOMTCH IS THE NUMBER OF ROOTS LESS THAN THE TRIAL VALUE.
             NOMTCH=N
             J=1
    360        F0 = B(J,2) - TRIAL
    370        CONTINUE
             IF(DABS(F0) .LT. THETA1)GO TO 380
             IF(F0 .GE. 0.0D0)NOMTCH=NOMTCH-1
             J = J + 1
             IF (J.GT.N) GO TO 390
    !    SINCE MATRIX IS NORMED TO UNITY, MAGNITUDE OF B(J,3) IS LESS THAN
    !    ONE, AND SINCE F0 IS GREATER THAN THETA1, OVERFLOW IS NOT POSSIBLE
    !    AT THE DIVISION STEP.
             F0 = B(J,2) - TRIAL - B(J-1,3)/F0
             GO TO 370
    380        J = J + 2
             NOMTCH = NOMTCH - 1
             IF (J.LE.N) GO TO 360
    390        CONTINUE
    !    FIX NEW BOUNDS ON ROOTS
             IF (NOMTCH.GE.I) GO TO 540
             ROOTL = TRIAL
             GO TO 500
    540        ROOT(I) = TRIAL
             NOM = MIN0(NROOT,NOMTCH)
             ROOT(NOM) = TRIAL
             GO TO 500
    330   CONTINUE
    !    REVERSE THE ORDER OF THE EIGENVALUES, SINCE CUSTOM DICTATES
    !    -LARGEST FIRST-.  THIS SECTION MAY BE REMOVED IF DESIRED WITHOUT
    !    AFFECTING THE REMAINDER OF THE ROUTINE.
    !    NRT = NROOT/2
    !    DO 10 I=1,NRT
    !    SAVE = ROOT(I)
    !    NMIP1 = NROOT - I + 1
    !!   ROOT(I) = ROOT(NMIP1)
    !10   ROOT(NMIP1) = SAVE
    !
    !    TRIVE!SECTION.
    !    EIGENVECTORS OF CODIAGONAL FORM
    !807  CONTINUE
    !    QUIT NOW IF NO VECTORS WERE REQUESTED.
        IF (NROOTX.LT.0) GO TO 1002
    !    INITIALIZE VECTOR ARRAY.
        DO 15 I=1,N
             DO 15 J=1,NROOT
    15   VECT(I,J)=1.0D-0
        DO 700 I=1,NROOT
             AROOT = ROOT(I)
    !    ORTHOGONALIZE IF ROOTS ARE CLOSE.
             IF (I.EQ.1) GO TO 710
    !    THE ABSOLUTE VALUE IN THE NEXT TEST IS TO ASSURE THAT THE TRIVEC
    !    SECTION IS INDEPENDENT OF THE ORDER OF THE EIGENVALUES.
             IF(DABS(ROOT(I-1)-AROOT) .LT. TOLER)GO TO 720
    710        IA = -1
    720        IA = IA + 1
             ELIM1 = A(1) - AROOT
             ELIM2 = B(1,1)
             JUMP = 1
             DO 750  J=1,NM1
                  JUMP = JUMP+J+1
    !    GET THE CORRECT PIVOT EQUATION FOR THIS STEP.
                  IF(DABS(ELIM1) .LE. DABS(B(J,1)))GO TO 760
    !    FIRST (ELIM1) EQUATION IS THE PIVOT THIS TIME.  CASE 1.
                  B(J,2) = ELIM1
                  B(J,3) = ELIM2
                  B(J,4)=0.0D0
                  TEMP = B(J,1)/ELIM1
                  ELIM1 = A(JUMP) - AROOT - TEMP*ELIM2
                  ELIM2 = B(J+1,1)
                  GO TO 755
    !    SECOND EQUATION IS THE PIVOT THIS TIME.  CASE 2.
    760             B(J,2) = B(J,1)
                  B(J,3) = A(JUMP) - AROOT
                  B(J,4) = B(J+1,1)
                  TEMP=1.0D0
                  IF(DABS(B(J,1)) .GT. THETA1)TEMP=ELIM1/B(J,1)
                  ELIM1 = ELIM2 - TEMP*B(J,3)
                  ELIM2 = -TEMP*B(J+1,1)
    !    SAVE FACTOR FOR THE SECOND ITERATION.
    755             B(J,5) = TEMP
    750        CONTINUE
             B(N,2) = ELIM1
             B(N,3)=0.0D-0
             B(N,4)=0.0D-0
             B(NM1,4)=0.0D-0
             ITER = 1
             IF (IA.NE.0) GO TO 801
    !    BACK SUBSTITUTE TO GET THIS VECTOR.
    790        L = N + 1
             DO 780 J=1,N
                  L = L - 1
    786             CONTINUE
                  ELIM1=VECT(L,I)-VECT(L+1,I)*B(L,3)-VECT(L+2,I)*B(L,4)
    !    IF OVERFLOW IS CONCEIVABLE, SCALE THE VECTOR DOWN.
    !    THIS APPROACH IS USED TO AVOID MACHINE-DEPENDENT AND SYSTEM-
    !    DEPENDENT CALLS TO OVERFLOW ROUTINES.
                  IF(DABS(ELIM1) .GT. DELBIG)GO TO 782
                  TEMP = B(L,2)
                  IF(DABS(B(L,2)) .LT. DELTA)TEMP=DELTA
                  VECT(L,I) = ELIM1/TEMP
                  GO TO 780
    !    VECTOR IS TOO BIG.  SCALE IT DOWN.
    782             DO 784 K=1,N
    784             VECT(K,I) = VECT(K,I)/DELBIG
                  GO TO 786
    780        CONTINUE
             GO TO (820,800), ITER
    !    SECOND ITERATION.  (BOTH ITERATIONS FOR REPEATED-ROOT VECTORS).
    820        ITER = ITER + 1
    890        ELIM1 = VECT(1,I)
             DO 830 J=1,NM1
                  IF (B(J,2).EQ.B(J,1)) GO TO 840
    !    CASE ONE.
                  VECT(J,I) = ELIM1
                  ELIM1 = VECT(J+1,I) - ELIM1*B(J,5)
                  GO TO 830
    !    CASE TWO.
    840             VECT(J,I) = VECT(J+1,I)
                  ELIM1 = ELIM1 - VECT(J+1,I)*TEMP
    830        CONTINUE
             VECT(N,I) = ELIM1
             GO TO 790
    !    PRODUCE A RANDOM VECTOR
    801        CONTINUE
             DO 802 J=1,N
    !    GENERATE PSEUDORANDOM NUMBERS WITH UNIFORM DISTRIBUTION IN (-1,1).
    !    THIS RANDOM NUMBER SCHEME IS OF THE FORM...
    !    RAND1 = AMOD((2**12+3)*RAND1,2**23)
    !    IT HAS A PERIOD OF 2**21 NUMBERS.
                  RAND1=DMOD(4099.0D0*RAND1,RPOWER)
    802       VECT(J,I)=RAND1/RPOW1-1.0D0
             GO TO 790
    !
    !    ORTHOGONALIZE THIS REPEATED-ROOT VECTOR TO OTHERS WITH THIS ROOT.
    800        IF (IA.EQ.0) GO TO 885
             DO 860 J1=1,IA
                  K = I - J1
                  TEMP=0.0D0
                  DO 870 J=1,N
    870             TEMP = TEMP + VECT(J,I)*VECT(J,K)
                  DO 880 J=1,N
    880             VECT(J,I) = VECT(J,I) - TEMP*VECT(J,K)
    860        CONTINUE
    885        GO TO (890,900), ITER
    !    NORMALIZE THE VECTOR
    900       ELIM1=0.0D0
             DO 904 J=1,N
    904       ELIM1=DMAX1(DABS(VECT(J,I)),ELIM1)
             TEMP=0.0D-0
             DO 910 J=1,N
                  ELIM2=VECT(J,I)/ELIM1
    910        TEMP = TEMP + ELIM2**2
             TEMP=1.0D0/(DSQRT(TEMP)*ELIM1)
             DO 920 J=1,N
                  VECT(J,I) = VECT(J,I)*TEMP
                  IF(DABS(VECT(J,I)) .LT. DEL1)VECT(J,I)=0.0D0
    920        CONTINUE
    700   CONTINUE
    !
    !    SIMVEC SECTION.
    !    ROTATE CODIAGONAL VECTORS INTO VECTORS OF ORIGINAL ARRAY
    !    LOOP OVER ALL THE TRANSFORMATION VECTORS
        IF (NM2.EQ.0) GO TO 1002
        JUMP = NSIZE - (N+1)
        IM = NM1
        DO 950  I=1,NM2
             J1 = JUMP
    !    MOVE A TRANSFORMATION VECTOR OUT INTO BETTER INDEXING POSITION.
             DO 955  J=IM,N
                  B(J,2) = A(J1)
    955        J1 = J1 + J
    !    MODIFY ALL REQUESTED VECTORS.
             DO 960  K=1,NROOT
                  TEMP=0.0D0
    !    FORM SCALAR PRODUCT OF TRANSFORMATION VECTOR WITH EIGENVECTOR
                  DO 970  J=IM,N
    970             TEMP = TEMP + B(J,2)*VECT(J,K)
                  TEMP = TEMP + TEMP
                  DO 980  J=IM,N
    980             VECT(J,K) = VECT(J,K) - TEMP*B(J,2)
    960        CONTINUE
             JUMP = JUMP - IM
    950   IM = IM - 1
    1002  CONTINUE
    !    RESTORE ROOTS TO THEIR PROPER SIZE.
        DO 95 I=1,NROOT
    95    ROOT(I) = ROOT(I)*ANORM
    1001  RETURN
        END SUBROUTINE GIVENS
!------------------------------------ End ------------------------------------