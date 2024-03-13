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
!                   SPHYNX: init.f90                  !
!                                                     !
! Initial values for several code variables.          !
!=====================================================!

    SUBROUTINE init

      USE parameters

      IMPLICIT NONE
      INTEGER i

      !Initialize number of particles
      n=nmax
      if(iterini.eq.1)nvi(:)=int(ni0)

      !Initialize times and counters
      dtnma=deltaini             !initial timestep
      tt=timeini                 !total time
      writeout=.false.
      maxtstepaux=maxtstep       !initial value of auxiliar maxtstep
      if(iterini.eq.1) then
        if(switches) then
          alfa(:)=alfamin        !initial AV coeficient
        else
          alfa(:)=alfamax
        endif
      endif

      if(delay.eq.0) then
         findballmass=.false.
      else
         findballmass=.true.
      endif

      !ALWAYS Initialize Cubic Spline kernel, needed for multipolar terms
      !ag1,ag2,etc...
      call cubicspline()
      if(kernel.eq.1) then   !Do nothing (already initialized)
         continue

      !If kernel is harmonic overwrite pk
      !initial harmonic kernel index (or indexes if more are defined)
      else if(kernel.eq.2) then
         indice(:)=nbaseline
         do i=1,n
            call calculate_norm(indice(i),dim,pk(i),dpk(i))
         enddo
      else if (kernel.eq.3) then
         indice(:)=nbaseline
         do i=1,n
            call calculate_norm(indice(i),dim,pk(i),dpk(i))
         enddo
         open(1,file='../sphynx/sinx_x.dat')
         do i=1,20001
            read(1,*) uk(i),fk(i),fdk1(i),fdk2(i)
         enddo
         close(1)
      else
         write(*,*) 'Not defined kernel! Kernel=',kernel
         stop
      endif

      if(eos.eq.4) then
        call read_helm_table
      endif


      ALLOCATE (neighbors(npend-npini+1,int(nvmax)))
      ALLOCATE (cube(npend-npini+1,int(nvmax)))


      RETURN
    END SUBROUTINE init
