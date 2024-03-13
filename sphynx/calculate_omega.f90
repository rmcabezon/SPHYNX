!=====================================================!
!                                                     !
!     This work is distributed under CC_BY_NC_SA      !
!                                                     !
! Created by Ruben M. Cabezon and Domingo Garcia-Senz !
!               ruben.cabezon@unibas.ch               !
!               domingo.garcia@upc.edu                !
!                                                     !
!                        (2018)                       !
!                                                     !
!=====================================================!
!                                                     !
!              SPHYNX: calculate_omega.f90            !
!                                                     !
! Calculate grad_h terms.                             !
!=====================================================!

  SUBROUTINE calculate_omega

    USE parameters

    IMPLICIT NONE

    INTEGER i
    DOUBLE PRECISION dhdro


    if(flags)print *,'Calculate omega',id
    omega=0.d0
    omega(npini:npend)=1.d0
    vol=0.d0

    !$omp parallel private(i,dhdro)

    ! calculation of omega(i)
    if(l.ge.liniNR) then
      !$omp do schedule(static)
       do i=npini,npend
          dhdro=-h(i)/promro(i)/dim
          omega(i)=1.d0-dhdro*sumwh(i)
       enddo
       !$omp end do
    endif

    !$omp do schedule(static)
    do i=npini,npend
       vol(i)=xmass(i)/sumkx(i)
    enddo
    !$omp end do
    !$omp end parallel

  END SUBROUTINE calculate_omega
