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
!          SPHYNX: treewalkmod_grav_mefec.f90         !
!                                                     !
! Calculates gravitational force and potential.       !
!=====================================================!

    SUBROUTINE treewalk

      USE parameters

      IMPLICIT NONE
      DOUBLE PRECISION d1,d2,d3,dc,ai1,ai2,ai3
      DOUBLE PRECISION d02,d05,d32,g0,r5,r7,rqr,c1,c2,c3
      DOUBLE PRECISION mefec,vgr,hij,v1,v2
      DOUBLE PRECISION ugraq,gt1,gt2,gt3


      INTEGER nodo,nodo0,multip,i,j,ii
      LOGICAL jump

      DOUBLE PRECISION, DIMENSION(3) ::  r12,qr,fq

      if(flags)write(*,*) 'Gravity calculation'
      call profile(0)

      multip=0

      ! Tree-walk calculating gravitational force and potential.

      !Inicialization
      f=0.d0           !Gravitational force
      ugrav=0.d0       !Gravitational potential

!$omp parallel private(i,ii,ai1,ai2,ai3,nodo,nodo0,d1,d2,d3,dc,&
!$omp                  r12,d02,jump,d05,d32,v1,v2,g0,hij,vgr,mefec,&
!$omp                  multip,r5,r7,qr,rqr,c1,c2,c3,fq,ugraq,j,gt1,gt2,gt3)


!$omp do schedule(guided)
      do 206 i = npini,npend

         ii=1+3*(i-1)

!     inclusion of quadrupolar term in gravitational potential.
         ai1=a(ii)
         ai2=a(ii+1)
         ai3=a(ii+2)
         nodo=1
 31      continue
         d1=dabs(ai1-xce(nodo))
         d2=dabs(ai2-yce(nodo))
         d3=dabs(ai3-zce(nodo))

         dc=4.d0*h(i)+dxx(nodo)/2.d0
         nodo0=ntotal+2

!     checking the intersection of the hierarchical cube against the
!     'oversized' yardstick dc of particle i.
         if((d1.le.dc) .and. (d2.le.dc) .and. (d3 .le.dc)) goto 35

!     nodo0 is the node having no branches intersecting h(i).
         nodo0=nodo
 36      continue
         r12(1)=ai1-xmc(nodo)
         r12(2)=ai2-ymc(nodo)
         r12(3)=ai3-zmc(nodo)
         d02=r12(1)**2+r12(2)**2+r12(3)**2

!     check aplicability of the multipolar expansion.
         jump=.false.
 322     if((jump).or.(dxx2(nodo) .le. tol*d02)) then


!        Gravity g is smoothed when the node is just a close particle.
            if(jump) then
               jump=.false.
               goto 32
            endif

            d05=dsqrt(d02)
            d32=1.d0/d05/d02
            if(dxx(nodo) .eq. 0.d0) then
               j=nnod(nodo)
               v1=d05/h(i)
               v2=d05/h(j)
               if(v1 .gt. 2.d0 .and. v2.gt.2.d0) then
                  g0=mcm(nodo)*d32
                  goto 148
               endif
               hij=h(i)+h(j)
               vgr=d05/hij
               mefec=vgr**3
               if(mefec.gt.1.) mefec=1.d0
#ifdef EQMASS
               g0=mefec*d32
#else
               g0=mefec*mass(j)*d32
#endif

 148           continue
               f(ii)=f(ii)-g0*r12(1)
               f(ii+1)=f(ii+1)-g0*r12(2)
               f(ii+2)=f(ii+2)-g0*r12(3)

               ugrav(i)=ugrav(i)+g0*d02
               multip=multip+1
            else

!     Direct calculation of g if the node is a particle.
               g0=mcm(nodo)*d32

               f(ii)=f(ii)-g0*r12(1)
               f(ii+1)=f(ii+1)-g0*r12(2)
               f(ii+2)=f(ii+2)-g0*r12(3)

               ugrav(i)=ugrav(i)+g0*d02

!     Quadrupolar contribution.
               r5=d02*d02*d05
               r7=r5*d02
               qr(1)=r12(1)*qxx(nodo)+r12(2)*qxy(nodo)+r12(3)*qxz(nodo)
               qr(2)=r12(1)*qxy(nodo)+r12(2)*qyy(nodo)+r12(3)*qyz(nodo)
               qr(3)=r12(1)*qxz(nodo)+r12(2)*qyz(nodo)+r12(3)*qzz(nodo)
               rqr=r12(1)*qr(1)+r12(2)*qr(2)+r12(3)*qr(3)
               c1=-7.5d0/r7*rqr
               c2=3.d0/r5
               c3=0.5d0*trq(nodo)
               fq(1)=c1*r12(1)+c2*(qr(1)+c3*r12(1))
               fq(2)=c1*r12(2)+c2*(qr(2)+c3*r12(2))
               fq(3)=c1*r12(3)+c2*(qr(3)+c3*r12(3))
               f(ii)=f(ii)+fq(1)
               f(ii+1)=f(ii+1)+fq(2)
               f(ii+2)=f(ii+2)+fq(3)

               ugraq=-1.5d0/r5*rqr+c3*d32
               ugrav(i)=ugrav(i)+ugraq

               multip=multip+1
            endif
 32         continue
            if(nodo .eq. nh2(npa(nodo))) then
               nodo=npa(nodo)
               if(nodo.eq.0) goto 301
               goto 32
            else
               if(nodo .le. nodo0) then
                  nodo=nodo+1
                  goto 31
               else
                  nodo=nodo+1
                  goto 36
               endif
            endif
         else
            nodo=nh1(nodo)
            goto 36
         endif
 35      continue

         if(dxx(nodo) .eq. 0.d0) then
            j=nnod(nodo)
            if(i.eq.j) then
               jump=.true.
               goto 322
            endif
            d02=d1*d1+d2*d2+d3*d3
            d05=sqrt(d02)
            if(d05 .gt. hd(i) .and. d05 .gt. hd(j)) then
               g0=1.d0/d05/d02
               goto 128
            endif
            v1=d05/h(i)
            v2=d05/h(j)
            hij=h(i)+h(j)
            vgr=d05/hij
            mefec=vgr**3
            if(mefec.gt.1.d0) mefec=1.d0
            g0=mefec/d05/d02

 128        continue

!     ------------------------------------------------------------
!     ------------------------------------------------------------
!
!     GRAVITY
!
!     ------------------------------------------------------------
!     ------------------------------------------------------------

            gt1=g0*(ai1-xmc(nodo))
            gt2=g0*(ai2-ymc(nodo))
            gt3=g0*(ai3-zmc(nodo))

#ifdef EQMASS
            f(ii)=f(ii)-gt1
            f(ii+1)=f(ii+1)-gt2
            f(ii+2)=f(ii+2)-gt3
            ugrav(i)=ugrav(i)+g0*d02
#else
            f(ii)=f(ii)-gt1*mass(j)
            f(ii+1)=f(ii+1)-gt2*mass(j)
            f(ii+2)=f(ii+2)-gt3*mass(j)
            ugrav(i)=ugrav(i)+g0*d02*mass(j)
#endif

            jump=.true.
            goto 322
!     closes the 'if' coming fromthe condition  if(dxx .eq. 0. above)
         endif
         nodo=nh1(nodo)
         goto 31
 301     continue

206   enddo
!$omp end do
!$omp end parallel

#ifdef EQMASS
      f(:)=f(:)*masspart
      ugrav(:)=ugrav(:)*masspart
#endif

      if(flags)write(*,*) 'Gravity calculation done!',id

      call profile(1)

!     ------------------------------------------------------------
!     ------------------------------------------------------------
!
!     ==>    End of gravity calculation   <==
!
!     ------------------------------------------------------------
!     ------------------------------------------------------------
!
      RETURN
    END SUBROUTINE treewalk
