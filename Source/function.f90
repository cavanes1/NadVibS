! All functions listed in alphabetical order

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
        if((index.ge.vbounds(get_batch,1).and.index.le.vbounds(get_batch,2)).or.get_batch.gt.batmax)exit
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
    open(unit=RESTARTFILE,file='restart.log',status='old'); rewind(unit=RESTARTFILE)
        read(RESTARTFILE,*); read(RESTARTFILE,*)
        read(RESTARTFILE,*); read(RESTARTFILE,*)get_restart_iter
    close(RESTARTFILE)
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

!Returns the number of possibilities for picking out d balls of out n (with back)
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