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
!           SPHYNX: first_findneighbors.f90           !
!                                                     !
! Initial h adjustment to get reasonable nvi(i)       !
!=====================================================!

  SUBROUTINE first_findneighbors

    USE mpi
    USE parameters

    IMPLICIT NONE

    INTEGER nodo,count,i,j,ii,jj,ik1,ierr
    DOUBLE PRECISION axce,ayce,azce,d1,d2,d3,dc,d05
    DOUBLE PRECISION v1,v2,nij,dmy,ah,jumpx,jumpy,jumpz
    DOUBLE PRECISION, DIMENSION(nt)::xcet,ycet,zcet
    DOUBLE PRECISION, DIMENSION(n)::aux,aux1
    LOGICAL recalc,jump
    LOGICAL, DIMENSION(nmax)::find_neighbors

    recalc=.true.
    find_neighbors=.true.
    count=0
    nvi(:)=0
    xcet(:)=0.d0
    ycet(:)=0.d0
    zcet(:)=0.d0

    !We don't exclude here i because we are only interested in getting nvi


    do while(recalc)
       recalc=.false.
       count=count+1
       if(flags)write(*,*) 'Readjusting for first neighbors',count,id
       if(count.eq.200)then
          print *,'Too many iterations!'
          stop
       endif

       do ik1=1,ncubes

          call apply_PBC(0,0,ik1,jumpx,jumpy,jumpz)

          xcet(1:ntotal)=xce(1:ntotal)+jumpx
          ycet(1:ntotal)=yce(1:ntotal)+jumpy
          zcet(1:ntotal)=zce(1:ntotal)+jumpz

          !$omp parallel private(i,ii,nodo,axce,ayce,azce,d1,d2,d3,dc,j,&
          !$omp                  jj,d05,v1,v2)
          !$omp do schedule(guided)
global:   do i = npini,npend
             if(find_neighbors(i)) then
                ii=1+dim*(i-1)
                nodo=0
51              nodo=nodo+1
                if(nodo.gt.ntotal) then
                   print *,'Error',id,nodo,ntotal
                   print *, i
                   print *,a(ii),a(ii+1),a(ii+2)
                   print *,xcet(nodo),xcet(nodo-1)
                   print *,jumpx,jumpy,jumpz
                   stop "endhere"
                endif
                axce=a(ii)-xcet(nodo)
                ayce=a(ii+1)-ycet(nodo)
                azce=a(ii+2)-zcet(nodo)
                d1=dabs(axce)
                d2=dabs(ayce)
                d3=dabs(azce)
                dc=4.d0*h(i)+dxx(nodo)/2.d0
                if((d1.le.dc).and.(d2.le.dc).and.(d3.le.dc)) goto 55
52              if (nodo .eq. nh2(npa(nodo))) then
                   nodo=npa(nodo)
                   if(nodo .eq. 0) then
                      if(nvi(i).le.int(nvmax).and.nvi(i).ge.int(nvmin).and.ik1.eq.ncubes)find_neighbors(i)=.false.
                      cycle global
                   endif
                   goto 52
                else
                   goto 51
                endif

!     If this condition is fulfilled means that it has found a particle.
!     If not, there is no need to continue with the operations.


55              if (dxx(nodo) .eq. 0.d0) then
                   j=nnod(nodo)
                   d05=sqrt(d1*d1+d2*d2+d3*d3)
                   v1=d05/h(i)
                   v2=d05/h(j)
                   if(v1.le.2.d0) then
!                   if(v1.le.2.d0.or.v2.le.2.d0) then
                      nvi(i)=nvi(i)+1
                   endif
                   goto 52
                endif
                nodo=nh1(nodo)-1
                goto 51
             endif
          enddo global
          !$omp end do
          !$omp end parallel
       enddo

       recalc=any(find_neighbors(npini:npend))

       if(recalc) then
         !$omp parallel private(i)
         !$omp do schedule(static)
          do i=npini,npend
             if(find_neighbors(i))then
                if(nvi(i).gt.int(nvmax).or.nvi(i).lt.int(nvmin)) then
                   print *,i,ind(i),nvi(i)
                   !h(i)=h(i)*.5d0*(1.d0+hfac*ni0/dble(nvi(i)))**hexp
                   h(i)=min(1.01d0*h(i),h(i)*.5d0*(1.d0+(ni0/dble(nvi(i)))**hexp))
                endif
                nvi(i)=1
             endif
          enddo
          !$omp end do
          !$omp end parallel
       endif
    enddo

    aux1=0.d0;aux1(npini:npend)=h(npini:npend);aux=0.d0
    call MPI_ALLREDUCE(aux1,aux,nmax,MPI_DOUBLE_PRECISION,MPI_SUM,new_comm,ierr);h=aux
    call calculate_hpowers

    RETURN
  END SUBROUTINE first_findneighbors
