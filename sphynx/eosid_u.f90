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
!                 SPHYNX: eosid_u.f90                 !
!                                                     !
! Ideal equation of state (input: internal energy)    !
!=====================================================!

    subroutine eosid_u
      
      USE parameters,only:nmax,promro,p,u,c,rgasid,gamma,dudt,dpdt,temp,mui
      
      implicit none
      integer i
      double precision gammam1
      
      gammam1=gamma-1.d0

!$omp parallel private(i)
!$omp do schedule(static)
      do i=1,nmax
         p(i)=u(i)*promro(i)*gammam1
         dudt(i)=1.5d0*rgasid/mui(i)
         temp(i)=u(i)/dudt(i)
         dpdt(i)=rgasid*promro(i)/mui(i)
         c(i)=sqrt(gammam1*u(i))
      enddo
!$omp end do
!$omp end parallel
      return

    end subroutine eosid_u

