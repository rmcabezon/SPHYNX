!=====================================================!
!                                                     !
!     This work is distributed under CC_BY_NC_SA      !
!                                                     !
! Created by Ruben M. Cabezon and Domingo Garcia-Senz !
!               ruben.cabezon@unibas.ch               !
!               domingo.garcia@upc.edu                !
!                                                     !
!                        (2017)                       !
!                                                     !
!=====================================================!
!                                                     !
!               SPHYNX: cubicspline.f90               !
!                                                     !
!Coefficients for cubic spline kernl and gravitational!
!smoothing.                                           !
!=====================================================!

subroutine cubicspline()

  USE parameters

  IMPLICIT NONE

  ad1=-1.5d0
  ad2=0.75d0
  ad3=0.25d0
  pk4=4.d0
  adw1=2.d0*ad1
  adw2=3.d0*ad2
  adw3=-3.d0*ad3
  if(dim.eq.3)then
     pk=pim1
  else if(dim.eq.2) then
     pk=10.d0*pim1/7.d0
  else if(dim.eq.1) then
     pk=2.d0/3.d0
  else
     write(*,*)'Bad dim in parameters: ',dim
     stop
  endif
      
  return
end subroutine cubicspline
