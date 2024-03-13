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
!                SPHYNX: outputmod.f90                !
!                                                     !
! Output of data.                                     !
!=====================================================!

    SUBROUTINE output

      USE parameters

      IMPLICIT NONE
      INTEGER i,j,ii
      DOUBLE PRECISION macum
      DOUBLE PRECISION,DIMENSION(dim)::ps
      DOUBLE PRECISION,DIMENSION(nmax)::vrad,gpr,fgr
      DOUBLE PRECISION,DIMENSION(nmax)::xpos,ypos,zpos
      CHARACTER*14 nomfs1
      CHARACTER*6 iteration
      CHARACTER*8 prefix

      call profile(0)


      !user-defined conditions for writing
      if(mod(l,iodelay).eq.0.or.(l.eq.iterini.and.outini))writeout=.true.
      if(checkdens)writeout=.true.
      !if(l.eq.36)writeout=.true.

      if(flags.and.writeout)write(*,*)'Writing data'

      if(writeout) then
         open(10,file='outputtimes.d',position='append')
         write(10,'(1x,i6,1x,a17,2(1x,es17.10))') l,'                 ',tt
         close(10)
      endif

      if(writeout) then
         writeout=.false.
         prefix='data/s1.'
         write(iteration,'(i6.6)') l
         nomfs1=prefix//iteration
         write(*,*) 'Output file: ',nomfs1
         open(10,file=nomfs1)!,form='unformatted')
!--------------------  User can change this  --------------------
         do i=1,n
            ii=1+dim*(i-1)
            ps(1)=a(ii)-despl-0.5d0
            ps(2)=a(ii+1)-despl-0.5d0
            ps(3)=a(ii+2)-despl-0.5d0
            radius(i)=sqrt(ps(1)**2+ps(2)**2+ps(3)**2)
            vrad(i)=(ps(1)*v(ii)+ps(2)*v(ii+1)+ps(3)*v(ii+2))/radius(i)
            write(10,formatout) ind(i),(ps(j),j=1,dim),h(i),u(i),promro(i),&
                 &    v(ii),v(ii+1),v(ii+2),radius(i),curlv(i),p(i),&
                 &    divv(i),omega(i),f(ii),f(ii+1),f(ii+2),gradp(ii),gradp(ii+1),&
                 &    gradp(ii+2),temp(i),avisc(i)*0.5d0,energy(i)*.5d0,&
                 &    ballmass(i),xmass(i),vrad(i),dble(nvi(i)),alfa(i),&
                 &    mark_ramp(i),sumwh(i),c(i),sumkx(i),checkInorm(i)
         enddo
!----------------------------------------------------------------

         close(10)

      endif

      call profile(1)

      RETURN
    END SUBROUTINE output
