!=====================================================!
!                                                     !
!          This work is distributed under MIT         !
!                                                     !
! Created by Ruben M. Cabezon and Domingo Garcia-Senz !
!               ruben.cabezon@unibas.ch               !
!               domingo.garcia@upc.edu                !
!                                                     !
!                        (2017)                       !
!                                                     !
!=====================================================!
!                                                     !
!              SPHYNX: init_scenario.f90              !
!                                                     !
! Initial values for several model variables.         !
!=====================================================!

SUBROUTINE init_scenario

  USE parameters

  IMPLICIT NONE
  INTEGER i,j,ii
  DOUBLE PRECISION dmax,energytot,width,ener0

  write(*,*) 'start initial scenario'
!--------------------  User can change this  --------------------
  if (iterini.eq.1) then
    energytot=1.d0     !Sedov test
    width=0.1d0
    ener0=energytot/pi**(3.d0/2.d0)/1.d0/width**3
    do i=1,n
       u(i)=ener0*dexp(-(radius(i)**2/width**2))+1.d-08
    enddo
    temp(:)=u(:)/dudt(:)
  endif

!----------------------------------------------------------------

  call calculate_hpowers

  !Desplacement and inicialization
  dmax=maxval(radius)                 !Estimation system size.
  rad=2.d0*sqrt(7.d0)*dmax            !Size of the tree root.
  despl=rad

  print *,'Rad:',rad

  a(:)=a(:)+despl

  !call indexx(n,radius,indx)       !Creates sorted vector for first iteration

  !Initial velocity is zero?
  if(inivel0)v=0.d0

  !Volume Elements

#ifdef EQMASS
  !Set to the STD until 5 iteration before NR update starts
!  if(iterini.le.liniNR-5) then
!     xmass(:)=masspart
!  else
     xmass(:)=masspart/promro(:)
!  endif
#else
  !Set to the STD until 5 iteration before NR update starts
!  if(iterini.le.liniNR-5) then
!     xmass(:)=mass(:)
!  else
!     xmass(:)=mass(:)/promro(:)
!  endif
#endif

  write(*,*) 'End initial scenario'

  RETURN
END SUBROUTINE init_scenario
