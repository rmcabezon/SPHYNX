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
!               SPHYNX: sphynx_hybrid.f90             !
!                                                     !
! Main code of SPHYNX.                                !
!=====================================================!

    program sphynx_hybrid

      use mpi

      use parameters

      implicit none

      INTEGER j,i,k,NRiter,idmy
      INTEGER ierr,auxi(nmax)
      INTEGER, DIMENSION(MPI_STATUS_SIZE):: sts
      DOUBLE PRECISION aux4,dmy1,dmy2
      DOUBLE PRECISION aux(nmax),aux1(nmax),aux2(n3)
      DOUBLE PRECISION statev(n4+2)

      LOGICAL auxl(nodes)

!.....MPI initialization
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      new_comm=MPI_COMM_WORLD

      if(nproc.ne.nodes) then
         write(*,*) 'NPROC = ',nproc,', should be the same as NODES = ',nodes
         write(*,*) 'Correct this in parameters.f90'
         call MPI_FINALIZE(ierr)
         stop
      endif

      if(id.lt.nproc-1) then
         npini=id*(nmax/nproc)+1
         npend=(id+1)*(nmax/nproc)
         print '(a21,i2,1x,a3,i2,1x,a33,i9,1x,a12,i9,a27,i2,a1)',&
              &'Hi! I am node number ', id,'of ',nproc,&
              &'nodes. Calculating from particle ',npini,'to particle ',&
              &npend,'.'
      else if(id.eq.nproc-1) then
         npini=id*(nmax/nproc)+1
         npend=nmax
         print '(a21,i2,1x,a3,i2,1x,a33,i9,1x,a12,i9,a27,i2,a1)',&
              &'Hi! I am node number ', id,'of ',nproc,&
              &'nodes. Calculating from particle ',npini,'to particle ',&
              &npend,'.'
      endif

!********************************************************
!.....Tests consistency of given parameters and report.
!********************************************************
      call testparameters
      call printreport

!********************************************************
!.....Initialization of Cubic Spline if necessary.
!********************************************************
      call init

!********************************************************
!.....Read initial data and initialization of scenario.
!********************************************************

      call readdata
      call init_scenario
      if(id.eq.0)open(5,file='timing.d',position='append')

      call buildtree
      call first_findneighbors


      timepass=0.d0

!.....Start iterations.
      do 100 l=iterini,nnl

         iti=0

         write(*,*) l,id
         dt=dtnma

!********************************************************
!1.....Builds oct-tree and calculates grav. moments.
!********************************************************
         call buildtree    !1 - Buildtree

!********************************************************
!2.....Finds neighbors.
!3.....Neighbors communications
!4.....Calculates density and omega (grad-h terms) iteratively.
!5.....Density communications
!********************************************************
         call findneighbors        !2 - Neighbors

         call profile(0)           !3 - Neighbors Comms
         aux1=0.d0;aux1(npini:npend)=ballmass(npini:npend)
         call MPI_ALLREDUCE(aux1,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);ballmass=aux;aux=0.d0
         auxi=0
         call MPI_ALLREDUCE(nvi,auxi,nmax,MPI_INTEGER,MPI_SUM,new_comm,ierr);nvi=auxi;auxi=0
         call profile(1)


         call profile(0)           !4 - Density + Omega + avglogrho + Equalization
         if(l.le.liniNR) then
            NRiter=1
         else
            NRiter=NRitermax
         endif

         ready=.false.
         localready=.false.
         i=0
         do while(any(.not.localready).and.i.lt.NRiter)
            i=i+1
            if(.not.localready(id+1)) call calculate_ro0
            aux=0.d0;aux(npini:npend)=ro0(npini:npend);ro0=0.d0
            call MPI_ALLREDUCE(aux,ro0,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr)
            aux=0.d0;aux(npini:npend)=xmass(npini:npend);xmass=0.d0
            call MPI_ALLREDUCE(aux,xmass,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr)
            if(.not.localready(id+1)) call calculate_density !Density
            if(NRiter.gt.1) then
               if(.not.localready(id+1)) call calculate_hNR
               aux1=0.d0;aux1(npini:npend)=h(npini:npend);aux=0.d0
               call MPI_ALLREDUCE(aux1,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);h=aux;aux=0.d0
               call calculate_hpowers
               auxl=.false.;auxl(id+1)=localready(id+1)
               call MPI_ALLREDUCE(auxl,localready,nodes,MPI_LOGICAL,MPI_LOR,new_comm,ierr)
            endif
         enddo

         print *,'NR iterations',i

         call calculate_omega      !Omega
         call profile(1)

         call profile(0)           !5 - Density comm
         aux1=0.d0;aux1(npini:npend)=sumwh(npini:npend)
         call MPI_ALLREDUCE(aux1,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);sumwh=aux;aux=0.d0
         aux1=0.d0;aux1(npini:npend)=promro(npini:npend)
         call MPI_ALLREDUCE(aux1,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);promro=aux;aux=0.d0
         aux1=0.d0;aux1(npini:npend)=omega(npini:npend)
         call MPI_ALLREDUCE(aux1,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);omega=aux;aux=0.d0
         aux1=0.d0;aux1(npini:npend)=sumkx(npini:npend)
         call MPI_ALLREDUCE(aux1,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);sumkx=aux;aux=0.d0
         aux1=0.d0;aux1(npini:npend)=indice(npini:npend)
         call MPI_ALLREDUCE(aux1,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);indice=aux;aux=0.d0
         aux1=0.d0;aux1(npini:npend)=vol(npini:npend)
         call MPI_ALLREDUCE(aux1,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);vol=aux;aux=0.d0
         ro2=promro*promro
         call profile(1)


!********************************************************
!6.....Calculates IAD terms.
!7....IAD communications
!********************************************************

         call calculate_IAD       !6 - IAD terms

         call profile(0)          !7 - IAD terms comm
         aux=0.d0
         call MPI_ALLREDUCE(c11,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);c11=aux;aux=0.d0
         call MPI_ALLREDUCE(c12,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);c12=aux;aux=0.d0
         call MPI_ALLREDUCE(c13,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);c13=aux;aux=0.d0
         call MPI_ALLREDUCE(c22,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);c22=aux;aux=0.d0
         call MPI_ALLREDUCE(c23,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);c23=aux;aux=0.d0
         call MPI_ALLREDUCE(c33,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);c33=aux;aux=0.d0
         call MPI_ALLREDUCE(checkInorm,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);checkInorm=aux;aux=0.d0
         call MPI_ALLREDUCE(checkdeltax,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);checkdeltax=aux;aux=0.d0
         call MPI_ALLREDUCE(checkdeltay,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);checkdeltay=aux;aux=0.d0
         call MPI_ALLREDUCE(checkdeltaz,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);checkdeltaz=aux;aux=0.d0

         call profile(1)

!********************************************************
!8....EOS and thermodynamic variables.
!********************************************************
         call eostot            !8 - EOS

!********************************************************
!9....Calculates Div v term for AV triggers.
!10....Div v communications.
!********************************************************
         call calculate_divv    !9 - Div-v

         call profile(0)        !10 - Div-v comm
         call MPI_ALLREDUCE(divv,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);divv=aux;aux=0.d0
         call MPI_ALLREDUCE(curlv,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);curlv=aux;aux=0.d0
         call profile(1)

!********************************************************
!11....Calculates AV switches
!12....AV switches communication.
!********************************************************
         if(switches) then
            call calculate_switches    !11 - AV switches

            call profile(0)            !12 - AV switches comm
            aux1=0.d0;aux1(npini:npend)=alfa(npini:npend);aux=0.d0
            call MPI_ALLREDUCE(aux1,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);alfa=aux;aux=0.d0
            call profile(1)
         else
            iti=iti+2
         endif

!********************************************************
!13....Calculates momentum and energy equations.
!14....Momentum and energy communications.
!********************************************************
         call momeqn            !13 - Mom + Energy eqs

         call profile(0)        !14 - Mom + Energy eqs comm
         aux=0.d0
         aux2=0.d0
         aux4=0.d0
         call MPI_ALLREDUCE(gradp,aux2,n3,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);gradp=aux2
         call MPI_ALLREDUCE(energy,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);energy=aux;aux=0.d0
         call MPI_ALLREDUCE(energyT,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);energyT=aux;aux=0.d0
         call MPI_ALLREDUCE(avisc,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);avisc=aux;aux=0.d0
         call MPI_ALLREDUCE(aviscu,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);aviscu=aux;aux=-1.d90
         call MPI_ALLREDUCE(maxvsignal,aux,nmax,MPI_DOUBLE_PRECISION,MPI_MAX,new_comm,ierr);maxvsignal=aux;aux=0.d0
         call MPI_ALLREDUCE(mark_ramp,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);mark_ramp=aux;aux=0.d0
         call profile(1)

!********************************************************
!15.....Calculates gravitational force
!16.....Gravity communications.
!********************************************************
         if (gravity) then

            call treewalk          !15 - Gravity calculation


            call profile(0)        !16 - Gravity comm.
            call MPI_ALLREDUCE(f,aux2,n3,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);f=aux2
            call MPI_ALLREDUCE(ugrav,aux,n,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);ugrav=aux
            call profile(1)
         else
            iti=iti+2
         endif

!********************************************************
!17....Integration of velocities, position and energy.
!********************************************************

         call update            !17 - Update

!********************************************************
!18....Timestep evaluation.
!********************************************************
         call timectrl          !18 - Timestep

!********************************************************
!19.....Conservation laws.
!********************************************************
         if(conserv)call conservation      !19 - Conservation laws


!********************************************************
!20.....Write output.
!********************************************************
         if(id.eq.0)call output            !20 - Output

!********************************************************
!.......Write timing.
!********************************************************
         if(id.eq.0)write(5,formatin) &
              (timepass(j),j=1,iti),&  !1-20 Timing of individual sections
              &sum(timepass)-timepass(20),&         !21 - Total time - output time
              &timepass(3)+timepass(5)+timepass(7)+timepass(10)+timepass(12)+timepass(14) !22 - Comms total time

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

100   enddo
      if(id.eq.0) close(5)

      DEALLOCATE (neighbors)
      DEALLOCATE (cube)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
    END PROGRAM sphynx_hybrid
