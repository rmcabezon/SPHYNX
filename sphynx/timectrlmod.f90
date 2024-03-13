!=====================================================!
!                                                     !
!         This work is distributed under MIT          !
!                                                     !
! Created by Ruben M. Cabezon and Domingo Garcia-Senz !
!               ruben.cabezon@unibas.ch               !
!               domingo.garcia@upc.edu                !
!                                                     !
!                        (2018)                       !
!                                                     !
!=====================================================!
!                                                     !
!                SPHYNX: timectrlmod.f90              !
!                                                     !
! NCalculates the next time-step.                     !
!=====================================================!

      SUBROUTINE timectrl

      USE parameters

      IMPLICIT NONE

      INTEGER i,icons

      DOUBLE PRECISION dtmin,dtcour,dtro,dtacc
      DOUBLE PRECISION vmod,cabmu,dmy

      call profile(0)
      if(flags)write(*,*) 'Timestep control'
      dtmin=1.d60



!     ------ COURANT ------
      dtcour=dtmin
      dtro=dtmin
      dtacc=dtmin

      if(timecour) then
         do i=1,n
            if(maxvsignal(i).ne.0.d0) then
               dmy=dkcour*h(i)/maxvsignal(i)
            else
               dmy=dkcour*h(i)/c(i)
            endif
            dtcour=min(dmy,dtcour)
         enddo
      endif


!     ------ IN DENSITY ------

      if(timero) then
         do i=1,n
            if(divv(i).eq.0.d0)cycle
            dmy=dkro/abs(divv(i))
            dtro=min(dmy,dtro)
         enddo
      endif

!     ------ IN ACCELERATION  ------

      if(timeacc) then
         do i=1,n
            if(acmod(i).eq.0.d0)cycle
            dmy=dkacc*sqrt(h(i)/acmod(i))
            dtacc=min(dmy,dtacc)
         enddo
      endif

      if(l.ge.iterini+delay) then
         tt=tt+dtnma
         dmy=dtnma
         dtnma=dmin1(dtro,dtcour,dtacc)
         if(dtnma.ge.1.1d0*dmy)dtnma=1.1d0*dmy

         if(dtnma.gt.maxtstepaux)dtnma=maxtstepaux

      endif


      if(timesave.and.id.eq.0) then
 1       format(10(1x,es12.5))
         open(1,file='timectrl.d',position='append')
         write(1,1) dtcour,dtro,dtacc,dtnma,tt
         close(1)
      endif

      call profile(1)

      RETURN
    END SUBROUTINE timectrl
