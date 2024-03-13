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
!                 SPHYNX: eosid_t.f90                 !
!                                                     !
! Ideal equation of state (input: temperature)        !
!=====================================================!

    subroutine eosid_t(ueos)
      
      USE parameters,only:nmax,promro,temp,mui,p,c,dudt,dpdt,rgasid
      
      implicit none
      integer i
      double precision,dimension(nmax),intent(out)::ueos
      double precision dpdr

!$omp parallel private(i,dpdr)
!$omp do schedule(static)
      do i=1,nmax
         p(i)=rgasid*promro(i)*temp(i)/mui(i)
         ueos(i)=3.d0/2.d0*rgasid*temp(i)/mui(i)
         dpdt(i)=rgasid*promro(i)/mui(i)
         dudt(i)=3.d0/2.d0*rgasid/mui(i)
         dpdr=p(i)/promro(i)
         c(i)=sqrt(dpdr)
      enddo
!$omp end do
!$omp end parallel
      
      return

    end subroutine eosid_t
