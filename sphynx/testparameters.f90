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
!              SPHYNX: testparameters.f90             !
!                                                     !
! Checks correctness of parameters.                   !
!=====================================================!

 SUBROUTINE testparameters

    USE parameters

    IMPLICIT NONE

    integer status

    status=0
    if(iterini.ne.1.and.timeini.eq.0.d0) then
       status=1
       write(*,*) 'WARNING: This is a restart (iterini!=1) with timeini!=0'
    endif

    if(.not.timecour.and..not.timero.and..not.timeacc) then
       status=2
       write(*,*) 'ERROR: No time-step criteria selected'
    endif

    if(dim.ne.3) then
       status=2
       write(*,*) 'ERROR: wrong dimension selected (dim)'
    endif

    if(ncubes.ne.1.and.ncubes.ne.9.and.ncubes.ne.27) then
       status=2
       write(*,*) 'ERROR: wrong number of cubes for PBC (ncubes)'
    endif

    if(switches.and.balsara) then
       status=2
       write(*,*) 'ERROR: AV switches and Balsara cannot be true at the same time'
    endif


    if(status.eq.0) then
       write(*,*) 'parameters testing .................. OK'
    else if(status.eq.1) then
       write(*,*) 'parameters testing .................. WARNING(s)!'
    else
       write(*,*) 'parameters testing .................. FAILED!'
       stop
    endif


    RETURN
  END SUBROUTINE testparameters
