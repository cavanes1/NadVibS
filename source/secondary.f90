! Secondary subroutines listed in alphabetical order

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
                    exit
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

subroutine get_keywords()!Read main input file basis.in
    use progdata
    use filedata, only: BASISFILE
    implicit none
    integer,dimension(10)          :: npirr
    integer,dimension(100)         :: basis,initstate
    integer                        :: i,itmp,restart,bconv,sodata,shiftref,get_restart_iter,reorthog
    real*8,dimension(10) :: weights
    !Default values of input variables in basis.in
    restart=0; restartrun=.false.; niter=1; iiter=1
    call setintarray(basis,int(100),one); call setintarray(initstate,int(100),zero)
    nstates=1; ordr=2; natoms=2; nmodes=1
    nirreps=1; call setintarray(npirr,int(10),zero); npirr(1)=nmodes
    shiftref=0; neworigin=.false.
    bconv=0; idroots=0; soroots=0
    reorthog=0; chkorthog=100; orthog=.false.; orthogexact=.false.; savevecs=.false.
    nseg=1; ztoler=1d-20
    maxdisk=1000; call setarray(weights,int(10),1d0)
    NAMELIST /NADVIBS/  niter,natoms,nmodes,nstates,basis,restart,bconv,idroots,soroots,reorthog,&
                        chkorthog,nseg,ztoler,maxdisk,weights,shiftref,nirreps,npirr,ordr,initstate
    open(unit=BASISFILE,file='basis.in',access='sequential',form='formatted',status='old')
    read(unit=BASISFILE,NML=NADVIBS); close(unit=BASISFILE)
    if(sum(npirr).ne.nmodes)then
        write(*,'(1x,A51)')'Warning: symmetry is disabled due to conflict input'
        write(*,'(1x,A73)')'sum of all irreducible representation degrees != total degree of freedoms'
        nirreps=1; npirr(1)=nmodes
    end if
    allocate(npirrep(nirreps)); npirrep=npirr(1:nirreps)
    if(reorthog>0) orthog=.true.
    if(reorthog>1) orthogexact=.true.
    idroots=abs(idroots); soroots=abs(soroots)
    if(restart.ne.0)then
        restartrun = .true.
        iiter = get_restart_iter()+1
        niter = niter + iiter - 1
    end if
    if(shiftref.ne.0) neworigin=.true.
    if(orthog.or.idroots>0.or.soroots>0) savevecs=.true.
    bjiconv = 10d0**(-bconv)
    allocate(statew(nstates)); statew=weights(1:nstates)
    allocate(nfunc(nmodes)); nfunc=basis(1:nmodes)
    allocate(istate(nmodes)); istate=initstate(1:nmodes)
    dimen=1; do i=1,nmodes; dimen= dimen*nfunc(i); end do
    1000 format(' SYMMETRY CONFLICT, ',i5,' != ',i5,' --> setting nirrep=1')
end subroutine get_keywords

subroutine load_restartinfo(reloadall)!Extract all the information from restart.log file
    use progdata, only: alpha,beta,omega,orthog
    use filedata, only: RESTARTFILE
    implicit none
    logical, intent(inout)::reloadall
    CHARACTER(75)::commentline
    CHARACTER(8) ::searchterm,currterm
    integer::i,nload
    real*8::dpval
    open(unit=RESTARTFILE,file='restart.log',status='old')
        read(RESTARTFILE,*); read(RESTARTFILE,*)i
        read(RESTARTFILE,*); read(RESTARTFILE,*)nload; reloadall=i==nload
        read(RESTARTFILE,*); read(RESTARTFILE,*)alpha(1:nload)
        read(RESTARTFILE,*); read(RESTARTFILE,*)beta(0:nload)
        read(RESTARTFILE,*); read(RESTARTFILE,*)omega(1,0:nload)
        read(RESTARTFILE,*); read(RESTARTFILE,*)omega(2,0:nload)
    close(RESTARTFILE)
    call system('mv -f restart.log restart.log.old')
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