!=====================================================!
!                                                     !
!     This work is distributed under CC_BY_NC_SA      !
!                                                     !
! Created by Ruben M. Cabezon and Domingo Garcia-Senz !
!               ruben.cabezon@unibas.ch               !
!               domingo.garcia@upc.edu                !
!                                                     !
!                        (2020)                       !
!                                                     !
!=====================================================!
!                                                     !
!           SPHYNX: calculate_switches.f90            !
!                                                     !
! Calculates AV switches                              !
!=====================================================!

  SUBROUTINE calculate_switches

    USE parameters

    IMPLICIT NONE

    INTEGER i,ii,k,j,jj,iii

    DOUBLE PRECISION d1,d2,d3,d05,d02,jumpx,jumpy,jumpz
    DOUBLE PRECISION v1,w1d,norm,dji1,dji2,dji3
    DOUBLE PRECISION vij1,vij2,vij3,dmy,alfa0
    DOUBLE PRECISION vijsignal_ij,wij,vijrij,vijsignal
    DOUBLE PRECISION alfaloc,decay,alfadot
    DOUBLE PRECISION graddivvx,graddivvy,graddivvz,graddivv
    DOUBLE PRECISION kern11i,kern12i,kern13i
    DOUBLE PRECISION kern21i,kern22i,kern23i
    DOUBLE PRECISION kern31i,kern32i,kern33i
    DOUBLE PRECISION termA1,termA2,termA3

    call profile(0)
    if(flags)write(*,*) 'AV switches calculation'


    !Initialization

    !This rval is from Cullen & Dehnen 2010
    !vijsignal is from Read & Hayfield 2011
    !Follows pretty much everything form C&D (2010) and R&H (2011) for alfadot


    !$omp parallel private(i,ii,iii,k,j,jj,d1,d2,d3,d02,d05,v1,w1d,&
    !$omp                  norm,jumpx,jumpy,jumpz,alfa0,dmy,graddivv,&
    !$omp                  vij1,vij2,vij3,vijrij,wij,vijsignal_ij,vijsignal,&
    !$omp                  alfaloc,decay,alfadot,graddivvx,graddivvy,graddivvz,&
    !$omp                  kern11i,kern12i,kern13i,kern21i,kern22i,kern23i,&
    !$omp                  kern31i,kern32i,kern33i,termA1,termA2,termA3,&
    !$omp                  dji1,dji2,dji3)

    !$omp do schedule(static)
    do i=npini,npend
       ii=1+dim*(i-1)
       iii=i-npini+1
       vijsignal=-1.d60
       alfa0=alfa(i)
       graddivvx=0.d0
       graddivvy=0.d0
       graddivvz=0.d0
       graddivv=0.d0
       norm=pk(i)*hm3(i)
       decay=0.d0

       do k=1,nvi(i)-1
          j=neighbors(iii,k)
          jj=1+dim*(j-1)

          call apply_PBC(i,k,0,jumpx,jumpy,jumpz)


          d1=a(ii)-a(jj)-jumpx
          d2=a(ii+1)-a(jj+1)-jumpy
          d3=a(ii+2)-a(jj+2)-jumpz
          d02=d1*d1+d2*d2+d3*d3
          d05=sqrt(d02)

          v1=d05/h(i)

          call Wkernel_noderiv(v1,indice(i),w1d)

          vij1=v(ii)-v(jj)
          vij2=v(ii+1)-v(jj+1)
          vij3=v(ii+2)-v(jj+2)

          vijrij=vij1*d1+vij2*d2+vij3*d3
          if(vijrij.lt.0.d0) then
             wij=vijrij/d05
             vijsignal_ij=c(i)+c(j)-3.d0*wij
          else
             vijsignal_ij=1.d-40*c(i)
          endif
          vijsignal=max(vijsignal,vijsignal_ij)

          dji1=-d1
          dji2=-d2
          dji3=-d3

          kern11i=c11(i)*dji1
          kern12i=c12(i)*dji2
          kern13i=c13(i)*dji3

          kern21i=c12(i)*dji1
          kern22i=c22(i)*dji2
          kern23i=c23(i)*dji3

          kern31i=c13(i)*dji1
          kern32i=c23(i)*dji2
          kern33i=c33(i)*dji3

          termA1=(kern11i+kern12i+kern13i)*w1d
          termA2=(kern21i+kern22i+kern23i)*w1d
          termA3=(kern31i+kern32i+kern33i)*w1d

          dmy=vol(j)*(divv(i)-divv(j))
          graddivvx=graddivvx+dmy*termA1
          graddivvy=graddivvy+dmy*termA2
          graddivvz=graddivvz+dmy*termA3
       enddo

       graddivv=norm*sqrt(graddivvx**2+graddivvy**2+graddivvz**2)

       if(divv(i).lt.0.d0) then
          dmy=h2(i)*graddivv
          alfaloc=alfamax*dmy/(dmy+h(i)*abs(divv(i))+0.05d0*c(i))
       else
          alfaloc=0.d0
       endif
       if(alfaloc.ge.alfa(i))then
          alfa(i)=alfaloc
       else
          decay=h(i)/(decay_constant*vijsignal)
          if(decay.eq.0.d0.or.isnan(decay)) then
            print *,'Zero or isnan decay in calculate_switches'
            print *,'i,h(i),vijsignal'
            print *,i,h(i),vijsignal
            print *,'decay:',decay
            stop
          endif

          if(alfaloc.ge.alfamin) then
             alfadot=(alfaloc-alfa0)/decay
          else
             alfadot=(alfamin-alfa0)/decay
          endif
          alfa(i)=alfa0+alfadot*dt
       endif
    enddo
    !$omp end do
    !$omp end parallel

    call profile(1)

    RETURN
  END SUBROUTINE calculate_switches
