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
!                 SPHYNX: readdata.f90                !
!                                                     !
! Input of data.                                      !
!=====================================================!

    SUBROUTINE readdata

      USE parameters

      IMPLICIT NONE

      DOUBLE PRECISION dmy                       !dummy real variables
      INTEGER idmy                               !dummy integer variable
      INTEGER i,ii,j,k,ierr
      LOGICAL endhere
      DOUBLE PRECISION, DIMENSION(n3):: aloc


      if(flags) write(*,*) 'Reading data...'

!--------------------  User can change this  --------------------
      if(iterini.eq.1) then
         open(1,file='../initialmodels/'//inputfile)!,form='unformatted')
         do i=1,n
            read(1,*) a(i),a(i+n),a(i+n2),h(i),promro(i)
            radius(i)=sqrt(a(i)**2+a(i+n)**2+a(i+n2)**2)
         enddo
         a(:)=a(:) + 0.5d0
         mui(:)=10.d0
         dudt(:)=3.d0/2.d0*rgasid/mui(:)
         masspart=1.d0/dble(n)
         print *,masspart
         if(id.eq.0) then
            open(22,file='REPORT',position='append')
            write(22,'(1x,a20,5x,es23.16)') 'masspart:',masspart
            close(22)
         endif
         !Reshaping
         !$omp parallel private(i,ii)
         !$omp do schedule(static)
         do i=1,n
            ii=1+dim*(i-1)
            aloc(ii)=a(i)
            aloc(ii+1)=a(i+n)
            aloc(ii+2)=a(i+n2)
         enddo
         !$omp end do
         !$omp end parallel
         a(:)=aloc(:)
         ballmass(:)=promro(:)*h(:)**3                                                                                                                                                                             

! RESTART -------------------------------------------------------
      else

         !open file
         
         !read data
         
         !initialize masspart
         
      endif


!----------------------------------------------------------------

      close(1)

      !Tests validity of initial data
      endhere=.false.
      do i=1,n*dim
         if(isnan(a(i))) then
            write(*,*) 'NAN position at: (',i,',',i/n,')',a(i)
            endhere=.true.
         else if(isnan(v(i))) then
            write(*,*) 'NAN velocity at: (',i,',',i/n,')',v(i)
            endhere=.true.
         endif
         if(endhere)stop
      enddo

#ifdef EQMASS
      if(isnan(masspart).or.masspart.le.0.d0) then
          write(*,*) 'Incorrect masspart',i,masspart
          endhere=.true.
      endif
#endif

      do i=1,n
         if(isnan(u(i))) then
            write(*,*) 'Incorrect U at i=',i,u(i)
            endhere=.true.
         endif
#ifdef EQMASS
#else
         if(isnan(mass(i)).or.mass(i).le.0.d0) then
            write(*,*) 'Incorrect mass at i=',i,mass(i)
            endhere=.true.
         endif
#endif
         if(isnan(h(i)).or.h(i).le.0.d0) then
            write(*,*) 'Incorrect smoothing-length at i=',i,h(i)
            endhere=.true.
         endif
         if(endhere)stop
      enddo

      if(flags) write(*,*) 'Reading data... Done!'


      RETURN
    END SUBROUTINE readdata
