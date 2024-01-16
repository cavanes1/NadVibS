!======================================================================================
subroutine gaussian(x,y,h,c,fwhm)
  implicit none
  real*8, intent(in) :: x,h,c,fwhm
  real*8, intent(out) :: y
  real*8 :: sigma
  sigma=fwhm/(2.d0*sqrt(2.d0*log(2.d0)))
  y=h*exp(-(x-c)**2/(2.d0*sigma**2))
  return
end
!======================================================================================
