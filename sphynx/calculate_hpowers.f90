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
!            SPHYNX: calculate_hpowers.f90            !
!                                                     !
! Calculates powers of smoothing length.              !
!=====================================================!

  SUBROUTINE calculate_hpowers

    USE parameters,only: h,hd,h2,h3,hm1,hm2,hm3,hm4,hm5,n

    IMPLICIT NONE
    
    INTEGER i

!$omp parallel private(i)
!$omp do schedule(static)
    do i=1,n
       hd(i)=2.d0*h(i)
       h2(i)=h(i)*h(i)
       h3(i)=h2(i)*h(i)
       hm1(i)=1.d0/h(i)
       hm2(i)=hm1(i)/h(i)
       hm3(i)=hm2(i)/h(i)
       hm4(i)=hm3(i)/h(i)
       hm5(i)=hm4(i)/h(i)
    enddo
!$omp end do
!$omp end parallel
    
    RETURN
  END SUBROUTINE calculate_hpowers
