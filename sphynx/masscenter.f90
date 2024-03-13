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
!                SPHYNX: masscenter.f90               !
!                                                     !
! Calculates the center of mass.                      !
!=====================================================!

subroutine masscenter(dataini,datafin)

#ifdef EQMASS
  USE parameters,only:a,dim,masspart,cm
#else
  USE parameters,only:a,dim,mass,cm
#endif

  implicit none
  integer i,ii
  DOUBLE PRECISION npart,mtotal
  integer,intent(in)::dataini,datafin

  cm=0.d0
#ifdef EQMASS
  mtotal=dble(datafin-dataini+1)
  do i=dataini,datafin
     ii=1+dim*(i-1)
     cm(1)=cm(1)+a(ii)
     cm(2)=cm(2)+a(ii+1)
     cm(3)=cm(3)+a(ii+2)
  enddo
#else
  mtotal=0.d0
  do i=dataini,datafin
     ii=1+dim*(i-1)
     cm(1)=cm(1)+a(ii)*mass(i)
     cm(2)=cm(2)+a(ii+1)*mass(i)
     cm(3)=cm(3)+a(ii+2)*mass(i)
  enddo
  do i=dataini,datafin
     mtotal=mtotal+mass(i)
  enddo
#endif
  do i=1,dim
     cm(i)=cm(i)/mtotal
  enddo

  return
end subroutine masscenter
