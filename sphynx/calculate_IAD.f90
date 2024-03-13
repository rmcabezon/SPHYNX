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
!              SPHYNX: calculate_IAD.f90              !
!                                                     !
! Calculates IAD terms as explained in                !
! Garcia-Senz, Cabezon, Escartin, A&A, 538, A9 (2012) !
!=====================================================!

  SUBROUTINE calculate_IAD

    USE parameters

    IMPLICIT NONE

    INTEGER i,j,k,iii,ii,jj

    DOUBLE PRECISION w1d,d05,v1,d1,d2,d3
    DOUBLE PRECISION norm
    DOUBLE PRECISION dmy,det,jumpx,jumpy,jumpz

    REAL*16::tau11,tau12,tau13,tau22,tau23,tau33

    call profile(0)
    if(flags)write(*,*) 'Calculating IAD terms'

    !Inicialization

    det=0.d0
    c11=0.d0
    c12=0.d0
    c13=0.d0
    c22=0.d0
    c23=0.d0
    c33=0.d0
    checkInorm=0.d0
    checkdeltax=0.d0
    checkdeltay=0.d0
    checkdeltaz=0.d0

    !$omp parallel private(i,ii,k,j,jj,iii,d1,d2,d3,d05,v1,w1d,norm,&
    !$omp                  dmy,jumpx,jumpy,jumpz,det,tau11,tau12,&
    !$omp                  tau13,tau22,tau23,tau33)
    !$omp do schedule(static)
    do i = npini,npend
       ii=1+dim*(i-1)
       iii=i-npini+1
       tau11=0.d0
       tau12=0.d0
       tau13=0.d0
       tau22=0.d0
       tau23=0.d0
       tau33=0.d0
       norm=pk(i)*hm3(i)

       do k=1,nvi(i)-1
          j=neighbors(iii,k)
          jj=1+dim*(j-1)

          call apply_PBC(i,k,0,jumpx,jumpy,jumpz)

          d1=a(ii)-a(jj)-jumpx
          d2=a(ii+1)-a(jj+1)-jumpy
          d3=a(ii+2)-a(jj+2)-jumpz
          d05=sqrt(d1*d1+d2*d2+d3*d3)
          v1=d05/h(i)

          call Wkernel_noderiv(v1,indice(i),w1d)

          !If vol(j) is smaller than 1.d-33 we could have a denormal number
          !here that gives 0 in det. Check this.
          dmy=w1d*vol(j)

          tau11=tau11+d1*d1*dmy
          tau12=tau12+d1*d2*dmy
          tau13=tau13+d1*d3*dmy
          tau22=tau22+d2*d2*dmy
          tau23=tau23+d2*d3*dmy
          tau33=tau33+d3*d3*dmy

          checkInorm(i)=checkInorm(i)+dmy
          checkdeltax(i)=checkdeltax(i)+d1*dmy
          checkdeltay(i)=checkdeltay(i)+d2*dmy
          checkdeltaz(i)=checkdeltaz(i)+d3*dmy

       enddo

       !taus should be multiplied by norm, so det should be norm**3
       !but we are interested on cXX vectors which are tau**2/det
       !That is 1/norm, so it is enough to multiply det by norm and
       !that gives cxx/norm.
       det=(tau11*tau22*tau33+2.d0*tau12*tau23*tau13-&
            &tau11*tau23*tau23-tau22*tau13*tau13-&
            &tau33*tau12*tau12)*norm

       c11(i)=(tau22*tau33-tau23*tau23)/det
       c12(i)=(tau13*tau23-tau33*tau12)/det
       c13(i)=(tau12*tau23-tau13*tau22)/det
       c22(i)=(tau11*tau33-tau13*tau13)/det
       c23(i)=(tau13*tau12-tau11*tau23)/det
       c33(i)=(tau11*tau22-tau12*tau12)/det

       checkInorm(i)=(checkInorm(i)+vol(i))*norm !self-contribution
       checkdeltax(i)=checkdeltax(i)*norm
       checkdeltay(i)=checkdeltay(i)*norm
       checkdeltaz(i)=checkdeltaz(i)*norm
    enddo
    !$omp end do
    !$omp end parallel

! If masspart~1.d30 then det~1.d100 and it can lead to float point errors.
! So we work with tau per mass unit, therefore det is per mass unit^3, so
! c terms should be 1./mass

    call profile(1)

    RETURN
  END SUBROUTINE calculate_IAD
