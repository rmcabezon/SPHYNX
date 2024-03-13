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
!                  SPHYNX: update.f90                 !
!                                                     !
! Updates position, velocity and internal energy.     !
!=====================================================!

!----------------------------------------------------------------------------!
! This subroutine updates the position, velocities, internal energy and h of !
! all the particles. Uses a 2nd order centered scheme for integrating the    !
! equations.                                                                 !
!                                                                            !
! Velocities are updated in two ways depending on the variable 'relax'.      !
! If relax=.true. we are performing a relaxation procedure to let the        !
! particles to find their equilibrium position, erasing artificial pressure  !
! gradients derived from the initial random distribution of the particles.   !
! Thus particles are allowed to move over a sphere with constant radius      !
! and damping in velocity is applied to avoid unphysically large velocities. !
! In this case is better to calculate without gravity.                       !
!                                                                            !
! If relax=.false. we are in a normal SPH simulation.                        !
!............................................................................!

    SUBROUTINE update

      USE parameters

      IMPLICIT NONE
      DOUBLE PRECISION ax,ay,az
      DOUBLE PRECISION ux,uy,uz,vrad,vradx,vrady,vradz
      DOUBLE PRECISION aux,aux2,ax0,ay0,az0,val1,val2,val3
      DOUBLE PRECISION,DIMENSION(dim)::vcm

      INTEGER,DIMENSION(n-1)::indx3
      DOUBLE PRECISION,DIMENSION(n-1)::rd,yed,tempd,rdaux,yedaux,tempdaux
      DOUBLE PRECISION,DIMENSION(n)::radiusold

      INTEGER count,i,ii,j,k,idmy
      LOGICAL endhere,corrected


      call profile(0)
      if(flags) write(*,*) 'Update'

!************************************************************************
!.. Standard SPH update.
!************************************************************************

      do i=1,n
         ii=1+dim*(i-1)
         ac(ii)=-(gradp(ii)-g*f(ii))
         ac(ii+1)=-(gradp(ii+1)-g*f(ii+1))
         ac(ii+2)=-(gradp(ii+2)-g*f(ii+2))
      enddo

!.. No movement for the first DELAY iterations to allow h stabilization.
      if(l.lt.iterini+delay) ac=0.d0

      !accel modulus
      acmod=0.d0
      do i=1,n
         ii=1+dim*(i-1)
         acmod(i)=acmod(i)+ac(ii)*ac(ii)+ac(ii+1)*ac(ii+1)+ac(ii+2)*ac(ii+2)
      enddo
      do i=1,n
         acmod(i)=sqrt(acmod(i))
      enddo


      if(l.eq.iterini) then
         dtold=dt
         energyold(:)=energy(:)
         aviscold(:)=avisc(:)
         aviscuold(:)=aviscu(:)
         do i=1,n
            ii=1+dim*(i-1)
            a0(ii)=a(ii)-v(ii)*dt
            a0(ii+1)=a(ii+1)-v(ii+1)*dt
            a0(ii+2)=a(ii+2)-v(ii+2)*dt
         enddo
      endif

      aux=dt+dtold*.5d0
      aux2=(dt+dtold)*.5d0

      v0=v !save v_i to v_(i-1) before updating to v_(i+1)

!$omp parallel private(i,ii,ax,ay,az,ax0,ay0,az0,val1,val2,val3,&
!$omp                  ux,uy,uz,vrad,vradx,vrady,vradz)
!$omp do schedule(static)
      do i=1,n
         ii=1+dim*(i-1)
         ax=a(ii)-despl                  !a_i
         ay=a(ii+1)-despl
         az=a(ii+2)-despl

         ax0=a0(ii)-despl                !a_(i-1)
         ay0=a0(ii+1)-despl
         az0=a0(ii+2)-despl

         val1=(ax-ax0)/dtold
         val2=(ay-ay0)/dtold
         val3=(az-az0)/dtold

         v(ii)=val1+ac(ii)*aux            !v_(i+1)
         v(ii+1)=val2+ac(ii+1)*aux
         v(ii+2)=val3+ac(ii+2)*aux

         !a0 has despl included. It has no effect in applying PBC to them.
         a0(ii)=a(ii)     !save a_i to a_(i-1) before updating to a_(i+1)
         a0(ii+1)=a(ii+1)
         a0(ii+2)=a(ii+2)

         a(ii)=ax+val1*dt+(v(ii)-val1)*dt*aux2/aux
         a(ii+1)=ay+val2*dt+(v(ii+1)-val2)*dt*aux2/aux
         a(ii+2)=az+val3*dt+(v(ii+2)-val3)*dt*aux2/aux
      enddo
!$omp end do



      if(ncubes.gt.1.and.l.ne.iterini) then
         !Assumes that everything is NOT CM-centered for PBC!!!

!$omp do schedule(static)
         do i=1,n                   !Check limits and applies PBC
            ii=1+dim*(i-1)          !a0 must have PBC applied when a has it
            if(a(ii).gt.xbox) then  !otherwise the distance with a0 for a particle
               a(ii)=a(ii)-xbox     !that appears on the other side of the box
               a0(ii)=a0(ii)-xbox   !would be unphysical (~xbox)
            endif
            if(a(ii+1).gt.ybox) then
               a(ii+1)=a(ii+1)-ybox
               a0(ii+1)=a0(ii+1)-ybox
            endif
            if(a(ii+2).gt.zbox) then
               a(ii+2)=a(ii+2)-zbox
               a0(ii+2)=a0(ii+2)-zbox
            endif
            if(a(ii).lt.0.d0)then
               a(ii)=a(ii)+xbox
               a0(ii)=a0(ii)+xbox
            endif
            if(a(ii+1).lt.0.d0)then
               a(ii+1)=a(ii+1)+ybox
               a0(ii+1)=a0(ii+1)+ybox
            endif
            if(a(ii+2).lt.0.d0)then
               a(ii+2)=a(ii+2)+zbox
               a0(ii+2)=a0(ii+2)+zbox
            endif
         enddo
!$omp end do
      endif

!$omp end parallel

      if(l.lt.liniNR) then
          !h(:)=1.35d0*(masspart/promro(:))**third
          h(:)=h(:)*.5d0*(1.d0+hfac*ni0/dble(nvi(:)))**hexp
          call calculate_hpowers
      endif

!      call masscenter(1,n)

!            do j=0,dim-1
!               do i=1,n
!                  a(i+j*n)=a(i+j*n)-cm(j+1)
!               enddo
!            enddo


      val1=dt*dt*.5d0/dtold
      val2=dt+val1
      if(l.lt.iterini+delay) then
         val1=0.d0
         val2=0.d0
      endif

      if(.not.relax) then
         !$omp parallel private(i)
         !$omp do schedule(static)
         do i=1,n
            !         u(i)=u(i)+(.5d0*(energy(i)+avisc(i))+aviscu(i))*val2-&  !U_(i+1)
            !              &(.5d0*(energyold(i)+aviscold(i))+aviscuold(i))*val1
            temp(i)=temp(i)+((.5d0*(energyT(i)+avisc(i))+aviscu(i))*val2-&  !U_(i+1)
                 &(.5d0*(energyTold(i)+aviscold(i))+aviscuold(i))*val1)/dudt(i)
         enddo
         !$omp end do
         !$omp end parallel
      endif

      !do i=npini,npend
      !   if(temp(i).le.0.d0.or.isnan(temp(i))) then
      !      ii=1+dim*(i-1)
      !     write(68,formatout) i,u(i),dble(l),a(ii),a(ii+1),h(i),dble(nvi(i)),&
      !      &       temp(i),energyT(i),avisc(i),aviscu(i),energyTold(i),aviscuold(i),&
      !      &       val1,val2,dt,dudt(i)
      !   endif
      !enddo



      if(l.ge.iterini+delay) then
            !For the next timestep:           !a_(i-1) is already done before.
         energyold(:)=energy(:)         !du/dt_(i-1)
         energyTold(:)=energyT(:)
         aviscold(:)=avisc(:)           !d(AV)/dt_(i-1)
         aviscuold(:)=aviscu(:)
         dtold=dt                       !dt_(i-1)
      endif

!      call masscenter(1,n)

      radius(:)=0.d0
      !$omp parallel private(i,ii)
      if(ncubes.eq.1) then
         !$omp do schedule(static)
         do i=1,n
            ii=1+dim*(i-1)
            radius(i)=sqrt(a(ii)*a(ii)+a(ii+1)*a(ii+1)+a(ii+2)*a(ii+2))
         enddo
         !$omp end do
      else
         !$omp do schedule(static)
         do i=1,n
            ii=1+dim*(i-1)
            radius(i)=sqrt((a(ii)-cm(1))**2+(a(ii+1)-cm(2))**2+(a(ii+2)-cm(3))**2)
         enddo
         !$omp end do
      endif

      !$omp do schedule(static)
      do i=1,n*dim
         a(i)=a(i)+despl
      enddo
      !$omp end do
      !$omp end parallel


      if(estabil)call estabilmod()!Creates radius(i) sorted with central dens

      !Smoothing length

! !      if(l.gt.liniNR) then
!
!       !$omp parallel private(i)
!       if(l.le.(liniNR-5)) then
! #ifdef EQMASS
!          !$omp do schedule(static)
!          do i=1,n
!             xmass(i)=masspart
!          enddo
!          !$omp end do
! #else
!          !$omp do schedule(static)
!          do i=1,n
!             xmass(i)=mass(i)
!          enddo
!          !$omp end do
! #endif
!       endif
!       !$omp end parallel

      endhere=.false.
      do i=npini,npend
         ii=1+dim*(i-1)
         if(isnan(a(ii)).or.isnan(a(ii+1))) then
            endhere=.true.
            write(*,*) i,'ISNAN!',l
            write(42,*) i,'ISNAN!',l
            write(42,'(40(1x,es12.5))') a(ii),a(ii+1),a(ii+2),v(ii),v(ii+1),v(ii+2),&
                 &f(ii),f(ii+1),f(ii+2),gradp(ii),gradp(ii+1),gradp(ii+2),&
                 &v0(ii),v0(ii+1),v0(ii+2),ac(ii),ac(ii+1),ac(ii+2),&
                 &a0(ii),a0(ii+1),a0(ii+2),(cm(j),j=1,dim),&
                 &radius(i),h(i),dble(nvi(i)),p(i),promro(i),u(i),ye(i),&
                 &c(i),dble(id),dble(l)
         endif
      enddo
      if(endhere)stop

 12   format(1x,i6,30(1x,es12.5))
      call profile(1)

      RETURN
    END SUBROUTINE update
