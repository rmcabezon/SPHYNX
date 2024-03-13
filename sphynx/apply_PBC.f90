!=====================================================!
!                                                     !
!          This work is distributed under MIT         !
!                                                     !
! Created by Ruben M. Cabezon and Domingo Garcia-Senz !
!               ruben.cabezon@unibas.ch               !
!               domingo.garcia@upc.edu                !
!                                                     !
!                       (2018)                        !
!                                                     !
!=====================================================!
!                                                     !
!               SPHYNX: apply_PBC.f90                 !
!                                                     !
! Calculates box displacements to apply PBC.          !
!=====================================================!


  subroutine apply_PBC(i,k,ik1,jumpx,jumpy,jumpz)
    
    use parameters,only: cube,dim,xbox,ybox,zbox,ncubes,npini
    
    implicit none
    
    integer, intent(in)::i,k
    integer, intent(in)::ik1
    double precision, intent(out)::jumpx,jumpy,jumpz

    integer ik11,iii

    jumpx=0.d0
    jumpy=0.d0
    jumpz=0.d0

    !This is for finding corresponding cube to neighbor.
    !Otherwise, user provides cube number
    if(ncubes.eq.1) return

    iii=i-npini+1
    
    if(ik1.eq.0) then
       ik11=cube(iii,k)
    else
       ik11=ik1
    endif

    if(dim.eq.2) then
       if(ik11.gt.1.and.ik11.le.3) then
          jumpx=xbox*(-1.d0)**ik11
       else if(ik11.gt.3.and.ik11.le.5) then
          jumpy=ybox*(-1.d0)**ik11
       else if(ik11.eq.6) then
          jumpx=xbox
          jumpy=ybox
       else if(ik11.eq.7) then
          jumpx=xbox
          jumpy=-ybox
       else if(ik11.eq.8) then
          jumpx=-xbox
          jumpy=ybox
       else if(ik11.eq.9) then
          jumpx=-xbox
          jumpy=-ybox
       endif
    else if(dim.eq.3) then
       if(ik11.gt.1.and.ik11.le.3) then
          jumpx=xbox*(-1.d0)**ik11
       else if(ik11.gt.3.and.ik11.le.5) then
          jumpy=ybox*(-1.d0)**ik11
       else if(ik11.eq.6) then
          jumpx=xbox
          jumpy=ybox
       else if(ik11.eq.7) then
          jumpx=xbox
          jumpy=-ybox
       else if(ik11.eq.8) then
          jumpx=-xbox
          jumpy=ybox
       else if(ik11.eq.9) then
          jumpx=-xbox
          jumpy=-ybox
       else if(ik11.eq.10) then
          jumpz=zbox
       else if(ik11.gt.10.and.ik11.le.12) then
          jumpx=xbox*(-1.d0)**ik11
          jumpz=zbox
       else if(ik11.gt.12.and.ik11.le.14) then
          jumpy=ybox*(-1.d0)**ik11
          jumpz=zbox
       else if(ik11.eq.15) then
          jumpx=xbox
          jumpy=ybox
          jumpz=zbox
       else if(ik11.eq.16) then
          jumpx=xbox
          jumpy=-ybox
          jumpz=zbox
       else if(ik11.eq.17) then
          jumpx=-xbox
          jumpy=ybox
          jumpz=zbox
       else if(ik11.eq.18) then
          jumpx=-xbox
          jumpy=-ybox
          jumpz=zbox
       else if(ik11.eq.19) then
          jumpz=-zbox
       else if(ik11.gt.19.and.ik11.le.21) then
          jumpx=xbox*(-1.d0)**ik11
          jumpz=-zbox
       else if(ik11.gt.21.and.ik11.le.23) then
          jumpy=ybox*(-1.d0)**ik11
          jumpz=-zbox
       else if(ik11.eq.24) then
          jumpx=xbox
          jumpy=ybox
          jumpz=-zbox
       else if(ik11.eq.25) then
          jumpx=xbox
          jumpy=-ybox
          jumpz=-zbox
       else if(ik11.eq.26) then
          jumpx=-xbox
          jumpy=ybox
          jumpz=-zbox
       else if(ik11.eq.27) then
          jumpx=-xbox
          jumpy=-ybox
          jumpz=-zbox
       endif
    endif

  END subroutine apply_PBC
