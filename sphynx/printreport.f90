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
!                SPHYNX: printreport.f90              !
!                                                     !
! Writes a report with the most important parameters  !
! of the simulation.                                  !
!=====================================================!

  SUBROUTINE printreport

    USE parameters

    IMPLICIT NONE
    character stat*8,statl*30,date*8,time*10

    open(1,file='REPORT')
1   format(1x,a20,5x,i10)
2   format(1x,a20,5x,a8)
3   format(1x,a20,5x,f5.1,1x,a1,1x,a6)
4   format(1x,a20,5x,a30)
5   format(1x,a20,5x,es8.1)
6   format(1x,a20,5x,f8.4)

    call DATE_AND_TIME(date, time)

    write(1,*) '*********************************************************'
    write(1,*) '             SIMULATION ON ',date
    write(1,*) '*********************************************************'
    write(1,*)
    write(1,*)
    write(1,*) '  SPH PARAMETERS'
    write(1,*) '------------------'
    write(1,1) 'Number of particles:',nmax
    write(1,1) 'Max. iterations:',nnl
    write(1,1) 'Initial iteration:',iterini
    if(iterini.ne.1) write(1,*)'This is a RESTART!'
!---------------------------------
    write(1,3) 'Neighbors:', ni0
!---------------------------------
    if(kernel.eq.1) then
       stat='CS'
    else if (kernel.eq.2) then
       stat='Sinc'
    else if (kernel.eq.3) then
       stat='Sinc_int'
    endif
    write(1,2) 'Kernel:',stat
    if(kernel.eq.2.or.kernel.eq.3)write(1,3)'Sinc exponent:',nbaseline
!---------------------------------
    if (gradh) then
       stat='YES'
    else
       stat='NO'
    endif
    write(1,2) 'grad-h:',stat
!---------------------------------
    if (balsara) then
       stat='YES'
    else
       stat='NO'
    endif
    write(1,2) 'Balsara correction:',stat
    if (switches) then
       stat='YES'
    else
       stat='NO'
    endif
    write(1,2) 'AV switches corr.:',stat
    if (switches) then
       write(1,3) 'AV alphamax:',alfamax
       write(1,3) 'AV alphamin:',alfamin
    else
       write(1,3) 'AV alpha:',alfamax
    endif
    write(1,3) 'AV beta:',betaAV
    write(1,3) 'AV thermal cond.:',alfau
!---------------------------------
    if (std_VE) then
       stat='YES'
    else
       stat='NO'
    endif
    write(1,2) 'STD Vol.El.:',stat

    if (Atmin.lt.0.d0.and.Atmax.lt.0.d0) then
       stat='YES'
    else if (Atmin.ge.1.d50) then
       stat='NO'
    else
       stat='Mixed'
    endif
    write(1,2) 'Crossed Xmass:',stat
    write(1,5) 'Atmin:',Atmin
    write(1,5) 'Atmax:',Atmax
    if(stat.eq.'Mixed') write(1,3) 'ramp:',ramp
!---------------------------------
    if (NRitermax.gt.0) then
       stat='YES'
    else
       stat='NO'
    endif
    write(1,2) 'NR for h update:',stat
    if (NRitermax.gt.0) then
       write(1,1) 'NR starts in:',liniNR
       write(1,1) 'NR max. iterations:',NRitermax
    endif
!---------------------------------
    if (ncubes.gt.1) then
       stat='YES'
    else
       stat='NO'
    endif
    write(1,2) 'PBC:',stat
    if (ncubes.gt.1) then
       write(1,1) 'Number of cubes:',ncubes
       write(1,6) 'xbox:',xbox
       write(1,6) 'ybox:',ybox
       write(1,6) 'zbox:',zbox
    endif
!---------------------------------
    write(1,*)
    write(1,*)
    write(1,*) '     PHYSICS'
    write(1,*) '------------------'
!---------------------------------
    write(1,5) 'Initial time:',timeini
    write(1,5) 'Initial timestep:',deltaini
    write(1,5) 'Max. timestep:',maxtstep
!---------------------------------
    if(gravity) then
       stat= 'Oct-tree'
    else
       stat='NONE'
    endif
    write(1,2) 'Gravity calc.:',stat
    if(gravity) write(1,3) 'Grav. tol:',tol
!---------------------------------
    if (inivel0) then
       stat='YES'
    else
       stat='NO'
    endif
    write(1,2) 'Initial vel. = 0:',stat
!---------------------------------
    if(eos.eq.1) then
       statl='Id.gas gamma=5/3 (T input)'
    else if(eos.eq.2) then
       statl='Id.gas gamma=5/3 (u input)'
    else if(eos.eq.3) then
      statl='Polytropic'
    else if(eos.eq.4) then
      statl='Helmholtz'
    endif
    write(1,4) 'EOS:',statl
!---------------------------------
    write(1,*)
    write(1,*)
    write(1,*) '     COMMENTS'
    write(1,*) '-------------------'
    write(1,'(a15,a20)') 'Initial model: ',inputfile
    write(1,'(a100)') comments
!--------------------------------
    write(1,*)
    write(1,*)
#ifdef EQMASS
    stat='YES'
#else
    stat='NO'
#endif
    write(1,2) 'Same-mass particles:',stat


    close(1)
  end SUBROUTINE printreport
