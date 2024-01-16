! Top level subroutines: called by program block

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
    dir = "./"
    setmem = 0
    setout = 0
    nargs = iargc()
    do i = 1,nargs
        call getarg(i,argbuf)
        if(setmem==1.and.(trim(argbuf).ne.trim(args(2)))) then
           read(argbuf,*)memory
           setmem = 0
        end if
        if(setout==1.and.(trim(argbuf).ne.trim(args(1)))) then
           dir = trim(adjustl(argbuf))//'/'
           setout = 0
        end if
        if(trim(argbuf)==trim(args(1))) setmem = 1
        if(trim(argbuf)==trim(args(2))) setout = 1
    end do
end subroutine read_cmdline

subroutine read_basis()
    use progdata
    use filedata, only: BASISFILE
    implicit none
#include 'mpif.h'
#include 'global.fh'
    integer                      :: i,j,istat
    integer                      :: numut
    if(myid==0) call get_keywords()
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

subroutine read_constants()!Read nadvibs.in and number the potential term (Hd expansion)
    use progdata, only: myid,ordr,nfunc,nmodes,nstates,neworigin,nstblks,nztrms,concoef,aomega,bomega, &
                        dvec,tmat,nzindx,nzblks,nzcoef,nztrms,noterms,ztoler,cpindx,AU2WAVE,zero,zerodp
    use filedata, only: POTFILE,OUTFILE
    implicit none
#include 'mpif.h'
#include 'global.fh'
    logical:: newterm,termchk
    integer:: i,j,k,l,n,p,cnt,loc,nblk,ioff
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
    if(myid==0) then!Only the root process will be able to see this file...
        write(*,*)'Reading diabatic Hamiltonian...'
        open(POTFILE,file='nadvibs.in',access='sequential',form='formatted')
            read(POTFILE,*); read(POTFILE,*)aomega(1:nmodes)
            do i = 1,nstblks
                read(POTFILE,*); read(POTFILE,*)concoef(i)
                do j = 1,ordr
                    read(POTFILE,*); read(POTFILE,*)POTterms(1:noterms(j),j,i)
                end do
            end do
            if(neworigin) then
                read(POTFILE,*); read(POTFILE,*)bomega(1:nmodes)
                read(POTFILE,*); read(POTFILE,*)dvec(1:nmodes)
                read(POTFILE,*); read(POTFILE,*)tmat(1:nmodes*nmodes)
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
            if(myid==0) write(*,'(1x,A7,I3,A8,I3)')'Block =',i,' order =',j
            call setintarray(otab,ordr,int(1))
            !My preference is to use pseudo nmodes+1 counter satisfying former digit >= latter digit,
            !corresponding to the direct sum of an ordr-th order tensor's 1st dimension vector
            otab(1)=0
            do k = 1,noterms(j)
                otab(1)=otab(1)+1!Add 1 to the 1st digit
                do l=1,j-1!Carry to latter digits
                    if(otab(l)>nmodes) then
                        otab(l)=1
                        otab(l+1)=otab(l+1)+1
                    end if
                end do
                do l=j-1,1,-1!Modify to satisfy former digit >= latter digit
                    if(otab(l)<otab(l+1)) otab(l)=otab(l+1)
                end do
                !End of otab generation
                call union(j,otab,cnt,uniq)
                if(abs(POTterms(k,j,i)).ge.ztoler)then
                    newterm = .true.
                    do n = 1,nztrms(cnt)
                        termchk = .true. 
                        do p = 1,cnt
                            if(uniq(p).ne.nztemp1(cnt,n,2*p-1).or.count(otab(1:j)==uniq(p)).ne.nztemp1(cnt,n,2*p))termchk=.false.
                        end do
                        if(termchk)then
                            newterm=.false.
                            loc = n
                            exit
                        end if
                    end do
                    if(newterm)then
                        nztrms(cnt) = nztrms(cnt) + 1
                        do n = 1,cnt
                            nztemp1(cnt,nztrms(cnt),2*n-1) = uniq(n)
                            nztemp1(cnt,nztrms(cnt),2*n)   = count(otab(1:j)==uniq(n))
                        end do
                        nztemp1(cnt,nztrms(cnt),2*cnt+1) = 1
                        nztemp1(cnt,nztrms(cnt),2*cnt+2) = i
                        nztemp2(cnt,nztrms(cnt),1) = POTterms(k,j,i)
                    else
                        nblk = nztemp1(cnt,loc,2*cnt+1)+1
                        nztemp1(cnt,loc,2*cnt+1) = nblk
                        nztemp1(cnt,loc,2*cnt+1+nblk) = i
                        nztemp2(cnt,loc,nblk) = POTterms(k,j,i)
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
    if (myid == 0) then
        write(*,'(1x,A15)',advance='no')'Non-zero terms:'
        do i=1,ordr-1
            write(*,'(I10)',advance='no')nztrms(i)
        end do
        write(*,'(I10)')nztrms(ordr)
        write(*,*)'Diabatic Hamiltonian has been read'
    end if
    deallocate(POTterms)
    deallocate(nztemp1)
    deallocate(nztemp2)
end subroutine read_constants

subroutine print_basis(umem,rmem)!Print a summary of job control information
    use progdata, only: ordr,natoms,nstates,nmodes,niter,dimen,orthog,orthogexact,maxstor,nfunc,ztoler, &
                        chkorthog,nproc,epsilon,maxdisk,statew,restartrun,iiter,soroots,bjiconv,nseg,   &
                        nirreps,npirrep,aomega,bomega,istate,nstblks,AU2WAVE,nzindx,nztrms,nzblks,zero,neworigin
    use filedata, only: outdir,OUTFILE,timestamp
    implicit none
    real*8, intent(in)    :: umem,rmem
    integer                         :: i,j,k,l,pordr,ioff
    integer,dimension(ordr)         :: otab
    integer,dimension(nstblks,ordr) :: ntot
    real*8,dimension(nmodes) :: dpvec
    open(unit=OUTFILE,file='output.dat',status='replace')
    write(OUTFILE,'(a)')   'NadVibS: nonadiabatic vibronic spectrum simulation package'
    write(OUTFILE,'(a)')   'Originate from NADVIBS.X by Michael Schuurman 2007'
    write(OUTFILE,'(a)')
    call timestamp()
    write(unit=OUTFILE,fmt='(a)') 'Calculation begins.'
    write(OUTFILE,'(a)')
    write(OUTFILE,'(a)')   'Input parameters read from basis.in:'
    write(OUTFILE,'(a)')   '------------------------------------'
    write(OUTFILE,'(A48,I10)')'  Number of atoms:                              ',natoms
    write(OUTFILE,'(A48,I10)')'  Number of electronic states:                  ',nstates
    write(OUTFILE,'(A48,I10)')'  Order of the expansion:                       ',ordr
    write(OUTFILE,'(A48,I10)')'  Number of vibrational modes:                  ',nmodes
    write(OUTFILE,'(A48,I10)')'  Number of irreducible representations         ',nirreps
    write(OUTFILE,'(A53,8I3)')'  Number of modes per irrep:                         ',npirrep
    write(OUTFILE,'(A48,I10)')'  Number of Lanczos iterations:                 ',niter
    write(OUTFILE,'(A48,I10)')'  Number of processors used in execution:       ',nproc
    write(OUTFILE,'(A48,ES10.0)')'  Zero tolerance for potential coeffs.:         ',ztoler
    write(OUTFILE,'(A53,I5)')'  Compute spin-orbit parameters for first N roots:   ',soroots
    if(restartrun) write(OUTFILE,'(A48,I10)')'  Restarting Lanczos Procedure on step:         ',iiter
    write(OUTFILE,'(A53,l5)')'  Perform lanczos vector re-orthogonalization:       ',orthog
    if(orthog) then
        write(OUTFILE,'(A53,l5)')'    - exact dot products for vector orthog.:         ',orthogexact
        if(.not.orthogexact) write(OUTFILE,'(A48,I10)')'    - compute exact orthog. every X iterations: ',chkorthog
        write(OUTFILE,'(A50,F8.1)')'    - max. amount of disk available for storage:  ',maxdisk
        write(OUTFILE,'(A48,I10)')'    - max. number of lanczos vectors to store:  ',maxstor
    end if
    do i = 1,nstates
        write(OUTFILE,'(A36,I2,A15,F5.2)')'  Initial weight of reference state ',i,':              ',statew(i)
    end do
    write(OUTFILE,'(a50,es8.0)')'  Convergence criteria for eigenvalues (bji):     ',bjiconv
    write(OUTFILE,'(A42,F16.0)')'  Total Memory Required (GA+NadVibS) (MB):',rmem
    write(OUTFILE,'(A42,F16.0)')'  Memory Available (MB):                  ',umem
    write(OUTFILE,'(A48,I10)')'  Number of Segements per Lanczos Vector:        ',nseg
    write(OUTFILE,'(A42,I16)')'  Dimensionality of a single State Vector:',dimen
    write(OUTFILE,'(A42,I16)')'  Total Dimensionality of H matrix:       ',nstates*dimen
    write(OUTFILE,'(a42,es16.6)')'  Machine precision for this architecture:',epsilon
    ioff = 0
    call setintmatrix(ntot,nstblks,ordr,zero)
    do i = 1,ordr
        otab(i) = i
        do j = 1,nztrms(i)
            pordr = 0
            do k = 1,i
                pordr = pordr + nzindx(2*k,ioff+j)
            end do
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
    write(OUTFILE,*); write(OUTFILE,'(a)')'  Number of Potential Terms per Block -----------'; write(OUTFILE,*)
    write(unit=OUTFILE,fmt=1006)otab(1:ordr)
    do i = 1,nstblks
        write(unit=OUTFILE,fmt=1005)i,ntot(i,1:ordr)
    end do
    write(OUTFILE,*); write(OUTFILE,'(a)')'  Basis set specification:'
    write(OUTFILE,'(a)')'    mode  # func.    omega'
    do i=1,nmodes; write(OUTFILE,'(4x,I3,I7,F12.3)')i,nfunc(i),aomega(i)*AU2WAVE; end do
    write(OUTFILE,*)
    if(neworigin)then
        do i=1,nmodes; dpvec(i)=bomega(i)*AU2WAVE; end do
    else
        do i=1,nmodes; dpvec(i)=aomega(i)*AU2WAVE; end do
    end if
    !Temporary -- only one mode may be excited in istate
        j=0
        do i = 1,nmodes
            if(istate(i).gt.0.and.j==0)then; j=1
            else; istate(i)=0; end if
        end do
    write(OUTFILE,*); write(OUTFILE,'(a)')'  Initial state specification:'
    write(OUTFILE,'(a)')'    mode  # quanta   omega'
    do i=1,nmodes; write(OUTFILE,'(4x,I3,I7,F12.3)')i,istate(i),dpvec(i); end do
    write(OUTFILE,*)
    1005 format('  Block ',i3,':',8(I6,'   '))
    1006 format('  ORDER:    ',8('    ',i2,'   '))
end subroutine print_basis

subroutine initialize_elements()
    use progdata, only: myid,nproc,nmodes,nstates,dimen,nfunc,q1,q2,q3,scr2,alpha,beta,totstart,totend,   &
                        aomega,iiter,niter,shft,orthog,epsilon,omega,loindex,roindex,zerodp,zero,one,two, &
                        firststep,restartrun,maxdisk,nseg,vbounds,savevecs,ordr,cpindx,nstblks,nztrms,    &
                        nzindx,stype,nsteps,nelems,nmax,homatel,iordr,nzblks,basdif,lo,hi,concoef,nmatel, &
                        ndiag,dmap,ndblks,dblks
    use filedata, only: outdir,OUTFILE,QRESTART,ARCHIVE
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
        if(myid==0)write(unit=OUTFILE,fmt=1000)
        do i = 1,nproc
            if(myid==0)write(unit=OUTFILE,fmt=1001)i
            call NGA_DISTRIBUTION(q1,i-1,qlo,qhi)
            if(myid==0)write(unit=OUTFILE,fmt=1002)qlo(1),qhi(1),qlo(2),qhi(2),(qhi(1)-qlo(1)+1)*qhi(2)
            call NGA_DISTRIBUTION(q2,i-1,qlo,qhi)
            if(myid==0)write(unit=OUTFILE,fmt=1003)qlo(1),qhi(1),qlo(2),qhi(2),(qhi(1)-qlo(1)+1)*qhi(2)
            call NGA_DISTRIBUTION(q3,i-1,qlo,qhi)
            if(myid==0)write(unit=OUTFILE,fmt=1004)qlo(1),qhi(1),qlo(2),qhi(2),(qhi(1)-qlo(1)+1)*qhi(2)
            seglen = nint(1.*(qhi(1) - qlo(1) + 1.)/nseg)
            do j = 1,nseg
                vbounds((i-1)*nseg+j,1) = qlo(1) + (j-1)*seglen 
                vbounds((i-1)*nseg+j,2) = qlo(1) + j*seglen - 1
            end do
            vbounds(i*nseg,2) = qhi(1)
        end do
        if(myid==0) then
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
            if(myid==0)then!Print out the results of the term counting
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
                    if(j.gt.ordr)exit
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
            if(myid==0) call load_restartinfo(vecload)
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
    if(myid==0) then
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
        if(myid==0) write(OUTFILE,1014)
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
    if(myid==0) then
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
    if(myid==0.and.iter.le.niter)write(unit=OUTFILE,fmt=1001)1.*(tend-iterend)/tcount
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
    if(myid==0) then
        write(unit=OUTFILE,fmt=1001)1.*(tend-iterend)/tcount
        if(.not.orthog.or.iter==1)write(unit=OUTFILE,fmt=1002)1.*(tend-iterstart)/tcount
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
        if(Mod(iter,chkorthog)==0) then
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
                      if(myid==0) write(*,1001)iter,iter+k-1,i,omega(k,i),ocheck(k,i)
                      omega(k,i) = ocheck(k,i)
                   end if
                end do
            end do
        end if
    end if
    call system_clock(tend,tcount,tmx)
    if(myid==0.and.(orthogexact.or.Mod(iter,chkorthog)==0)) then
        write(unit=OUTFILE,fmt=1002)1.*(tend-iterend)/tcount
    elseif(myid==0) then
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
    if(myid==0) then
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

subroutine make_restart(iter)!Save alpha and beta coefficents to restart.log
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
    if(myid==0) then
        write(unit=OUTFILE,fmt='(a)')'  - Lanczos Iterations Completed - '
        write(unit=OUTFILE,fmt=1003)1.*(totend-totstart)/tcount
        open(unit=RESTARTFILE,file='restart.log',status='replace')
        rewind(unit=RESTARTFILE)
        write(unit=RESTARTFILE,fmt='(a)')' TOTAL NUMBER OF ITERATIONS IN FILE'
        write(unit=RESTARTFILE,fmt=*)iter
        write(unit=RESTARTFILE,fmt='(a)')' ITERATION NUMBER TO RESTART FROM'
        write(unit=RESTARTFILE,fmt=*)iter
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
    open(unit=QRESTART,file=trim(adjustl(outdir))//'nadvibs.35.'//trim(procindex),&
    access='direct',status='replace',form='unformatted',recl=8*rlen)  
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
    write(OUTFILE,*); write(unit=OUTFILE,fmt=1004)
    do i = 1,ngood
        write(unit=OUTFILE,fmt=1005)i,int(T(i)),bji(T(i)),&
            (Tvals(T(i))-Tvals(T(1)))*AU2WAVE,(Tvals(T(i))-Tvals(1))*AU2EV,Tints(T(i))/dpval
    end do
    call system_clock(totstart,tcount,tmx)
    write(OUTFILE,*); write(OUTFILE,'(A)')'   All requested eigenvalues computed.'
    write(unit=OUTFILE,fmt=1000)1.*(totstart-totend)/tcount
    deallocate(T); deallocate(Teig)
    1000 format('  Time Required: ',f14.5,' secs.')
    1001 format('  ZPVE from nadvibs.in (reference state): ',f12.4,' cm-1')
    1002 format('  FIRST eigenvalue:                       ',f12.4,' cm-1')
    1003 format('  Intensity Factor (MAX Intensity):       ',f12.8)
    1004 format('  Root  Index  Convergence     E(cm-1)        E(eV)     Intensity')
    1005 format(     I6,    I7,      es13.4,         f15.5,      f12.5,     f12.7)
end subroutine compute_eigenvalues

subroutine print_footer()
    use filedata, only: OUTFILE,timestamp
    write(unit=OUTFILE,fmt='(a)') ' '
    call timestamp(); write(unit=OUTFILE,fmt='(a)') 'Calculation completed.'
    close(unit=OUTFILE)
end subroutine print_footer

!This subroutine computes the eigenvectors of the hamiltonian matrix,
!then determines diagnostic information about the roots & dot products useful in spin-orbit computation
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
    if(nroots.gt.0.and.myid==0)then
        open(unit=ROOTINFO,file='rootinfo.dat',status='replace')
    end if
    if(soroots.gt.0)then
        if(myid==0)then
            open(unit=SOINFO,file='sodata.in',status='replace')
            rewind(unit=SOINFO)
        end if
        open(unit=EIGVECS,file=trim(adjustl(outdir))//'nadvibs.37.'//trim(procindex),access='direct', &
            status='replace',form='unformatted',recl=8*(hi(1)-lo(1)+1))
    end if
    call make_tmatrix(T,iter)
    call compute_ritzsystem(T,iter,bji,Tvals,Tvecs,Tints,Tsign)
    intfactor = 0
    do i = 1,nroots
        if(Tints(i).gt.intfactor)intfactor=Tints(i)
    end do
    if(myid==0) then
        write(unit=ROOTINFO,fmt=1000)nroots
        if(soroots.gt.0) then
            write(unit=SOINFO,fmt=1001)soroots
            write(unit=SOINFO,fmt=1002)
            do i = 1,soroots
                write(unit=SOINFO,fmt=1003)i,Tvals(i)*AU2EV,Sqrt(Tints(i)/intfactor)*Tsign(i)
            end do
            if(nstates==2)write(unit=SOINFO,fmt=1004)
            if(nstates==3)write(unit=SOINFO,fmt=1005)
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
        if(i.le.soroots) call write_ga(q1,EIGVECS,i)
        call ga_sync()
        if(myid==0) write(unit=ROOTINFO,fmt=1006)i
        !Analyze Eigenvectors
            if(myid==0) then
                write(ROOTINFO,*); write(unit=ROOTINFO,fmt=1007)i
                write(unit=ROOTINFO,fmt=1008)Tvals(i)*AU2WAVE,Tvals(i)*AU2EV,Tints(i)/intfactor
            end if
            call GA_ELEM_MULTIPLY(q1,q1,q2)
            do k = 1,idroots
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
            if(myid==0) then
                do j = 1,idroots
                    call get_barray(cindex(2*j-1),barray)
                    do k = 1,nmodes
                        barray(k) = barray(k) - 1
                    end do
                    write(unit=ROOTINFO,fmt=1009)cindex(2*j),cindex(2*j-1),int(cphase(j)),cmag(j)
                    write(unit=ROOTINFO,fmt=1010)barray
                end do
                write(unit=ROOTINFO,fmt='(a)')' '
                if(nstates==2)write(unit=ROOTINFO,fmt=1011)dp(1),dp(2)
                if(nstates==3)write(unit=ROOTINFO,fmt=1012)dp(1),dp(2),dp(3)
                if(nstates==4)write(unit=ROOTINFO,fmt=1016)dp(1),dp(2),dp(3),dp(4)
                write(unit=ROOTINFO,fmt='(a)')' --------------------------------------------------------------'
            end if
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
                    if(myid==0) then
                        if(nstates==2)write(unit=SOINFO,fmt=1013)i,j,dp(1),dp(2)
                        if(nstates==3)write(unit=SOINFO,fmt=1014)i,j,dp(1),dp(2),dp(3),i,j,dp(4),dp(5),dp(6)
                    end if
                end do
            end if
    end do
    if(myid==0) then
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
    use progdata; use filedata; implicit none
#include 'global.fh'
    logical::lstat; integer::istat
    if(allocated(scr1)) deallocate(scr1)
    deallocate(scr2); deallocate(concoef)
    deallocate(nzindx); deallocate(nzblks); deallocate(nzcoef)
    deallocate(alpha); deallocate(beta)
    deallocate(aomega); deallocate(bomega); deallocate(omega)
    deallocate(nfunc); deallocate(shft); deallocate(statew)
    deallocate(loindex); deallocate(roindex)
    if(allocated(dmap)) deallocate(dmap)
    if(allocated(dordr)) deallocate(dordr)
    deallocate(bmax); deallocate(nsteps); deallocate(stype)
    deallocate(basind); deallocate(basdif)
    deallocate(bstart); deallocate(rcshft); deallocate(nelems)
    if(savevecs) close(ARCHIVE)
    if(myid==0) close(unit=RESTARTFILE)
    lstat = GA_DESTROY(q1)
    if(.not.lstat)print *,'ERROR in GA_DESTROY(q1), ID=',myid
    lstat = GA_DESTROY(q2)
    if(.not.lstat)print *,'ERROR in GA_DESTROY(q2), ID=',myid
    lstat = GA_DESTROY(q3) 
    if(.not.lstat)print *,'ERROR in GA_DESTROY(q3), ID=',myid
    call ga_sync()
end subroutine release_memory
