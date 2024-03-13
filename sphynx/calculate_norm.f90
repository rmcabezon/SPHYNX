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
!              SPHYNX: calculate_norm.f90             !
!                                                     !
!Normalization factor for sinc kernels as explained in!
!         Garcia-Senz, et al. A&A, 570 (2014)         !
!=====================================================!

   SUBROUTINE calculate_norm(n,dim,norm,dnorm)
     IMPLICIT NONE
     DOUBLE PRECISION, intent(out):: norm,dnorm
     DOUBLE PRECISION, intent(in):: n
     INTEGER, intent(in):: dim
     DOUBLE PRECISION n05,a0,a1,a2,a3

     n05=sqrt(n)

     if(dim.eq.3) then
        a0= 2.7012593d-2
        a1= 2.0410827d-2
        a2= 3.7451957d-3
        a3= 4.7013839d-2
        norm= a0+a1*n05+a2*n+a3*n*n05
        dnorm= 0.5d0*a1/n05+a2+1.5d0*a3*n05
     else if (dim.eq.2) then
        a0= 5.2245027d-2
        a1= 1.3090245d-1
        a2= 1.9358485d-2
        a3=-6.1642906d-3
        norm= a0+a1*n+a2/n+a3/n/n
        dnorm= a1-a2/n/n-2.d0*a3/n/n/n
     else if (dim.eq.1) then
        a0=-1.5404568d-2
        a1= 3.6632876d-1
        a2=-4.6519576d-4
        a3= 7.3658324d-2
        norm= a0+a1*n05+a2*n+a3/n05
        dnorm= 0.5d0*a1/n05+a2-0.5d0*a3/n/n05
     else
        write(*,*) 'Wrong value of dim!: ', dim
        stop
     endif

     RETURN
   END SUBROUTINE calculate_norm
