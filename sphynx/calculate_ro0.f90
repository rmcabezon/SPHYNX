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
!             SPHYNX: calculate_density.f90           !
!                                                     !
! Calculates density and omega                        !
!=====================================================!

  SUBROUTINE calculate_ro0

    USE parameters

    IMPLICIT NONE

    INTEGER i,j,k,ii,jj,iii

    DOUBLE PRECISION w1d,dw1,d05,v1,d1,d2,d3
    DOUBLE PRECISION dterh,norm,dnorm
    DOUBLE PRECISION jumpx,jumpy,jumpz

    if(flags)write(*,*) 'Rho0 calculation'

    !$omp parallel private(i,ii,k,j,jj,iii,d1,d2,d3,d05,v1,w1d,dw1,&
    !$omp                  dterh,jumpx,jumpy,jumpz,norm,dnorm)

    !Inicialization
    !$omp do schedule(guided)
    do i=npini,npend
       if(ready(i))cycle
       ro0(i)=0.d0
       sumwhro0(i)=0.d0
    enddo
    !$omp end do

    !$omp do schedule(guided)
    do i = npini,npend
       if(ready(i))cycle
       ii=1+dim*(i-1)
       iii=i-npini+1
       norm=pk(i)*hm3(i)
       dnorm=pk(i)*hm4(i)

       do k=1,nvi(i)-1
          j=neighbors(iii,k)
          jj=1+dim*(j-1)

          call apply_PBC(i,k,0,jumpx,jumpy,jumpz)

          d1=a(ii)-a(jj)-jumpx
          d2=a(ii+1)-a(jj+1)-jumpy
          d3=a(ii+2)-a(jj+2)-jumpz
          d05=sqrt(d1*d1+d2*d2+d3*d3)
          v1=d05/h(i)

          call Wkernel(v1,indice(i),w1d,dw1)

          dterh=-(dim*w1d+v1*dw1)
#ifdef EQMASS
          ro0(i)=ro0(i)+w1d
          sumwhro0(i)=sumwhro0(i)+dterh
#else
          ro0(i)=ro0(i)+mass(j)*w1d
          sumwhro0(i)=sumwhro0(i)+mass(j)*dterh
#endif
       enddo

       !Self-contribution
#ifdef EQMASS
       ro0(i)=(ro0(i)+1.d0)*norm*masspart
       sumwhro0(i)=(sumwhro0(i)-dim)*dnorm*masspart
       if(std_VE) then
         xmass(i)=masspart
       else
         xmass(i)=masspart/ro0(i)
       endif
#else
       ro0(i)=(ro0(i)+mass(i))*norm
       sumwhro0(i)=(sumwhro0(i)-mass(i)*dim)*dnorm
       if(std_VE) then
          xmass(i)=mass(i)
       else
          xmass(i)=mass(i)/ro0(i)
       endif
#endif
    enddo
    !$omp end do

    !$omp end parallel



    RETURN

  END SUBROUTINE calculate_ro0
