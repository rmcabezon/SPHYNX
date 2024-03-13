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
!                 SPHYNX: profile.f90                 !
!                                                     !
! Calculates wall-clock time between calls.           !
!=====================================================!

  SUBROUTINE profile(mode)

    USE parameters,only:timeinit,timefin,timepass,iti

    IMPLICIT NONE

    INTEGER,intent(in)::mode
    DOUBLE PRECISION OMP_GET_WTIME

    if(mode.eq.0) then
       timeinit=OMP_GET_WTIME()

    else if(mode.eq.1)then
       timefin=OMP_GET_WTIME();iti=iti+1
       timepass(iti)=timefin-timeinit
    else if(mode.eq.2) then
       timefin=OMP_GET_WTIME()
       timepass(iti)=timefin-timeinit
    else
       print *,'Wrong mode value in profile!'
       stop
    endif

    RETURN
  END SUBROUTINE profile
