!=====================================================!
!                                                     !
!          This work is distributed under MIT         !
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

  SUBROUTINE calculate_density

    USE parameters

    IMPLICIT NONE

    INTEGER i,j,k,ii,jj,iii,ik1

    DOUBLE PRECISION w1d,dw1,w2d,dw2,d05,v1,v2,d1,d2,d3
    DOUBLE PRECISION dterh,norm,dnorm
    DOUBLE PRECISION dmy,jumpx,jumpy,jumpz


    if(flags)write(*,*) 'Density calculation'


    !Inicialization
    f=0.d0           !Gravitational force
    do i=npini,npend
       if(ready(i))cycle
       sumkx(i)=0.d0       !Volumen partition
       sumwh(i)=0.d0       !grad-h term
    enddo



    !$omp parallel private(i,ii,k,j,jj,iii,d1,d2,d3,d05,v1,w1d,dw1,&
    !$omp                  dterh,jumpx,jumpy,jumpz,norm,dnorm)

    !$omp do schedule(static)
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

          dterh=-(dim*w1d+v1*dw1)   !additional v1 for dw1

          sumkx(i)=sumkx(i)+xmass(j)*w1d
          sumwh(i)=sumwh(i)+xmass(j)*dterh
       enddo

       !adds self-contribution to density and omega
       sumkx(i)=(sumkx(i)+xmass(i))*norm
       sumwh(i)=(sumwh(i)-xmass(i)*dim)*dnorm

    enddo
    !$omp end do

    !$omp do schedule(static)
    do i=npini,npend
       if(ready(i))cycle
#ifdef EQMASS
          promro(i)=sumkx(i)*masspart/xmass(i)
#else
          promro(i)=sumkx(i)*mass(i)/xmass(i)
#endif
       if(std_VE) then
          sumwh(i)=sumwhro0(i)
       else
#ifdef EQMASS
          sumwh(i)=sumwh(i)*masspart/xmass(i)+&
            &   (promro(i)/ro0(i)-xmass(i)*pk(i)*hm3(i))*sumwhro0(i)
#else
          sumwh(i)=sumwh(i)*mass(i)/xmass(i)+&
            &   (promro(i)/ro0(i)-xmass(i)*pk(i)*hm3(i))*sumwhro0(i)
#endif
       endif
    enddo
    !$omp end do

    !$omp end parallel

    RETURN

  END SUBROUTINE calculate_density
