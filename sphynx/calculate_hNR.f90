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
!             SPHYNX: calculate_hNR_eq.f90            !
!                                                     !
! Newton-Raphson for h.                               !
!=====================================================!

  SUBROUTINE calculate_hNR

    USE parameters

    IMPLICIT NONE

    INTEGER i
    DOUBLE PRECISION g1,g11,deltah

    print *,'calculate h Newton-Rapshon',id

!$OMP parallel private(i,g1,g11,deltah)
!$OMP do schedule(static)
    do i=npini,npend
       if(ready(i))cycle
       ! Find the zero of functions  g1 and g2
       g1=ballmass(i)*hm3(i)-promro(i)
       g11=-(dim*ballmass(i)*hm4(i)+sumwh(i))
       ! Solving the matrix
       deltah=-g1/g11
       if(l.le.liniNR) deltah=0.d0
       if(abs(deltah)/h(i).le.1.d-6) ready(i)=.true.
       if(abs(deltah/h(i)).ge.0.05d0)deltah=sign(1.d0,deltah)*0.05d0*h(i)
       if(isnan(deltah)) then
         print *,"NAN deltah",i,ballmass(i),hm4(i),sumwh(i)
         stop
       endif
       if(h(i)+deltah.le.0.d0) then
         print *,'negative h',i,h(i),deltah
         do while (h(i)+deltah.le.0.d0)
           deltah=deltah*0.9d0
         end do
         print *,'solved!',i,h(i),deltah
       endif
       h(i)=h(i)+deltah
    enddo
!$OMP end do
!$OMP end parallel
    localready(id+1)=.not.any(.not.ready(npini:npend))

  END SUBROUTINE calculate_hNR
