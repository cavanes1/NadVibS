!=====================================================================================
program main
  implicit none
  real*8, parameter :: ev2cm=8065.51d0
  integer, parameter :: nmax=9999
  real*8 :: h(nmax), c(nmax)
  real*8 :: xstart, xend, dx, x, y, yt
  real*8 :: eshift
  real*8 :: fwhm
  integer :: nt
  integer :: i,j,ios,total

  open(unit=100,file='spec.line',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  do i=1,nmax
    read(100,*,iostat=ios) c(i),h(i)
    if(ios.ne.0) then
      total=i-1
      exit
    end if
  end do
  print*, total, 'valid points.'

  xstart=-0.1
  xend=0.4

  nt=2000
  fwhm=5.d-3

  eshift=2.7807
  !eshift=0.d0

  open(200,file='output')
  dx=(xend-xstart)/dble(nt+1)
  do i=0,nt+1
    x=xstart+dble(i)*dx
    y=0.d0
    do j=1,total
      call gaussian(x,yt,h(j),c(j),fwhm)
      y=y+yt
    end do
    write(200,"(2f15.8)") x+eshift,y
  end do

  stop
end
!=====================================================================================
