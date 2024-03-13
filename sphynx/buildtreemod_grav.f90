!=====================================================!
!                                                     !
!          This work is distributed under MIT         !
!                                                     !
! Created by Ruben M. Cabezon and Domingo Garcia-Senz !
!               ruben.cabezon@unibas.ch               !
!               domingo.garcia@upc.edu                !
!                                                     !
!                      (Jul 2017)                     !
!                                                     !
!=====================================================!
!                                                     !
!            SPHYNX: builtreemod_grav.f90             !
!                                                     !
! Tree construction and gravitational force evaluation!
!=====================================================!



!    ************************************************************
!    ------------------------------------------------------------
!                        SUBROUTINE  buildtree
!    ------------------------------------------------------------
!    ************************************************************
!
!         calculation of the gravitational acceleration  f(i)
!
!     ----------------------------------------------------------
!                               GRAV3D
!     ----------------------------------------------------------
!
!     Subroutine that builds an oc-tree to facilitate the calculation
!     of the self-gravity of a particle system. At each hierarchical level
!     (given by the value of n) the algorithm divides any subcube containing
!     two or more particles in 8 subcubes with 1/2 side, until reaching a
!     subcube with one or none particles in it.
!
!     The tree is built horizontally from the top (root) until the last
!     branch. Each hierarchical level n is divided in nm sublevels (cubes)
!     with an a half size; if the sublevel is empty it is ignored but if the
!     sublevel contains just one particle their coordinates are stored.
!     If the cube has more than one particle we store the mass inside the cube,
!     the coordinates of the center of mass and the quadrupolar moment.
!     Those nodes having more than one particle are afterwards re-ordered
!     (parameter nm) and we move to the next hierarchical level containing the
!     descendants (childrens) of the re-ordered nodes (the fathers).
!
!     Simultaneously, the index of the particles hosted by a subcube is also
!     re-ordered so that when descending to the next level each subcube is
!     looped only by the particles of the node father and not by all paticles
!     of the system.
!
!     The output of the algorithm provides the characteristic size of the node,
!     the coordinates and mass of the center of mass, and the quadripolar
!     moment of the different hierarchical levels.
!
!
!     Definition of variables:
!
!     x(i),y(i),z(i)    coordinates of the particles
!
!     ntotal            sequential number of the tree's nodes,
!                       at least equal to the number of particles
!
!     n                 number of particles
!
!     npn,np            gives the order of the first and last particle of
!                       each subcube at level n
!
!     nje               jerarchy level or height of the tree
!
!     nm                number of subcubes at level n with 2 or more particles
!
!     dx,dy,dz          characteristic size of the subcube at level n
!
!     xi,xd, etc.       left and right coordinates of the subdividing segment
!
!     npp               acumulates the renumbering of all particles
!
!     npp0              acumulates the renumbering of the particles within
!                       each one of the son subcubes, in order to improve the
!                       computational time in the next jerarchic level
!
!     xmct,ymct,zmct    center of mass coordinates of each cube
!
!     xcet,ycet,zcet    center coordinates of each cube
!
!     xia(nm,1), etc.   origen coordinates of the subcubes sorted by nm and
!                       that will be used in the next level
!
!     nh1(j)            order by ntotal of the j node's first son
!
!     nh2(j)            order by ntotal of the j node's last son
!
!     npa(j)            order by ntotal of j's father
!
!     dxx               size of subcube
!
!     mcm               subcube mass
!
!-----------------------------------------------------------------------------


    SUBROUTINE buildtree

      USE parameters

      IMPLICIT NONE

      DOUBLE PRECISION xt,dx,dy,dz
      DOUBLE PRECISION xcma,ycma,zcma,pmt,xec,yec,zec
      DOUBLE PRECISION qxxa,qxya,qxza,qyya,qyza,qzza
      DOUBLE PRECISION x1,y1,z1
      INTEGER nje,nm,nnma,npp,nnm,nmmm,i1,i2,i3,npp0
      INTEGER nodo,i,ii,j,npointer1,idmy,nsatura,nnmax

      DOUBLE PRECISION,DIMENSION(nmax) :: x,y,z
      DOUBLE PRECISION,DIMENSION(ntree,1) :: xia,yia,zia
      DOUBLE PRECISION,DIMENSION(ntree,2) :: xi,yi,zi,xd,yd,zd
      DOUBLE PRECISION,DIMENSION(dim) :: r12
      INTEGER,DIMENSION(nt) :: npn,np,npointer
      INTEGER, DIMENSION(ntree) :: nauxa,naux


      if(flags)write(*,*) 'Tree building'
      call profile(0)

      do 34 i=1,n
         ii=1+dim*(i-1)
         x(i)=a(ii)
         y(i)=a(ii+1)
         z(i)=a(ii+2)
         npointer(i)=i
34    enddo

      idmy=2**dim
      do 33 i=1,idmy
         nauxa(i)=0
33    enddo
      nsatura=0
      nnmax=0
      nnma=0
      ntotal=0
      npn(1)=1
      npn(2)=n
      xt=2.d0*despl
      nje=0
      nm=1
      xia(1,1)=0.d0
      yia(1,1)=0.d0
      zia(1,1)=0.d0
      xce(:)=0.d0
      yce(:)=0.d0
      zce(:)=0.d0
      do while(nm.ne.0)
         nje=nje+1
         npp=0
         dx=xt/2.d0**nje
         dy=dx
         dz=dx
         do 70 j=1,nm
            np(2*j-1)=npn(2*j-1)
            np(2*j)=npn(2*j)
            xi(j,1)=xia(j,1)
            yi(j,1)=yia(j,1)
            zi(j,1)=zia(j,1)
            naux(j)=nauxa(j)
70       enddo
         nnm=nm
         nnmax=max(nnm,nnmax)
         if(nnm .ge. nnma) nnma=nnm
         if(nnm .gt. ntree-10) nsatura=1
         nm=0

         do 5 j=1,nnm
            nmmm=0
            do 12 i1=1,2
               xi(j,i1)=xi(j,1)+dble(i1-1)*dx
               xd(j,i1)=xi(j,i1)+dx
               do 13 i2=1,2
                  yi(j,i2)=yi(j,1)+dble(i2-1)*dy
                  yd(j,i2)=yi(j,i2)+dy
                  do 15 i3=1,2
                     zi(j,i3)=zi(j,1)+dble(i3-1)*dz
                     zd(j,i3)=zi(j,i3)+dz
                     npp0=0
                     xcma=0.d0
                     ycma=0.d0
                     zcma=0.d0
                     pmt=0.d0

                     xec=(xi(j,i1)+xd(j,i1))/2.d0
                     yec=(yi(j,i2)+yd(j,i2))/2.d0
                     zec=(zi(j,i3)+zd(j,i3))/2.d0

                     if(gravity) then
                        qxxa=0.d0
                        qxya=0.d0
                        qxza=0.d0
                        qyya=0.d0
                        qyza=0.d0
                        qzza=0.d0
                     endif

                     do 101 i=np(2*j-1),np(2*j)
                        if((xi(j,i1).le.x(i)).and.(x(i).lt.xd(j,i1)))then
                           if((yi(j,i2).le.y(i)).and.(y(i).lt.yd(j,i2)))then
                              if((zi(j,i3).le.z(i)).and.(z(i).lt.zd(j,i3)))then
                                 npp=npp+1
                                 npp0=npp0+1
#ifdef EQMASS
                                 xcma=xcma+x(i)
                                 ycma=ycma+y(i)
                                 zcma=zcma+z(i)
                                 pmt=pmt+1.d0
#else
                                 xcma=xcma+x(i)*mass(i)
                                 ycma=ycma+y(i)*mass(i)
                                 zcma=zcma+z(i)*mass(i)
                                 pmt=pmt+mass(i)
#endif
!     calculates cuadrupolar term of gravitational potential

                                 if(gravity)then
                                    r12(1)=x(i)-xec
                                    r12(2)=y(i)-yec
                                    r12(3)=z(i)-zec
#ifdef EQMASS
                                    qxxa=qxxa+r12(1)*r12(1)
                                    qxya=qxya+r12(1)*r12(2)
                                    qxza=qxza+r12(1)*r12(3)
                                    qyya=qyya+r12(2)*r12(2)
                                    qyza=qyza+r12(2)*r12(3)
                                    qzza=qzza+r12(3)*r12(3)
#else
                                    qxxa=qxxa+r12(1)*r12(1)*mass(i)
                                    qxya=qxya+r12(1)*r12(2)*mass(i)
                                    qxza=qxza+r12(1)*r12(3)*mass(i)
                                    qyya=qyya+r12(2)*r12(2)*mass(i)
                                    qyza=qyza+r12(2)*r12(3)*mass(i)
                                    qzza=qzza+r12(3)*r12(3)*mass(i)
#endif
                                 endif

                                 x1=x(npp)
                                 y1=y(npp)
                                 z1=z(npp)
                                 npointer1=npointer(npp)
                                 x(npp)=x(i)
                                 y(npp)=y(i)
                                 z(npp)=z(i)
                                 npointer(npp)=npointer(i)
                                 x(i)=x1
                                 y(i)=y1
                                 z(i)=z1
                                 npointer(i)=npointer1

                              endif
                           endif
                        endif
101                  enddo

                     if(npp0-1.ge.0) then !For npp0-1<0 does nothing
                        if(npp0-1.gt.0)then
                           ntotal=ntotal+1
                           if(ntotal .gt. nt-50) nsatura=1
                           nm=nm+1
                           nmmm=nmmm+1
                           xce(ntotal)=(xi(j,i1)+xd(j,i1))/2.d0
                           yce(ntotal)=(yi(j,i2)+yd(j,i2))/2.d0
                           zce(ntotal)=(zi(j,i3)+zd(j,i3))/2.d0
                           npa(ntotal)=naux(j)
                           nauxa(nm)=ntotal
                           xmc(ntotal)=xcma/pmt
                           ymc(ntotal)=ycma/pmt
                           zmc(ntotal)=zcma/pmt

!     quadrupolar moment of an agregate of particles respect its CM

                           if(gravity) then
                              r12(1)=xec-xmc(ntotal)
                              r12(2)=yec-ymc(ntotal)
                              r12(3)=zec-zmc(ntotal)
                              qxx(ntotal)=qxxa-pmt*r12(1)*r12(1)
                              qxy(ntotal)=qxya-pmt*r12(1)*r12(2)
                              qxz(ntotal)=qxza-pmt*r12(1)*r12(3)
                              qyy(ntotal)=qyya-pmt*r12(2)*r12(2)
                              qyz(ntotal)=qyza-pmt*r12(2)*r12(3)
                              qzz(ntotal)=qzza-pmt*r12(3)*r12(3)
                              trq(ntotal)=qxx(ntotal)+qyy(ntotal)+qzz(ntotal)
                           endif

                           mcm(ntotal)=pmt

                           dxx(ntotal)=dx
                           dxx2(ntotal)=dx*dx
                           npn(2*nm-1)=npp-npp0+1
                           npn(2*nm)=npp
                           xia(nm,1)=xi(j,i1)
                           yia(nm,1)=yi(j,i2)
                           zia(nm,1)=zi(j,i3)

                        else       !Case npp0-1=0

                           ntotal=ntotal+1
                           if(ntotal .gt. nt-50) nsatura=1
                           npa(ntotal)=naux(j)
                           nmmm=nmmm+1
                           xce(ntotal)=x(npp)
                           yce(ntotal)=y(npp)
                           zce(ntotal)=z(npp)

                           xmc(ntotal)=x(npp)
                           ymc(ntotal)=y(npp)
                           zmc(ntotal)=z(npp)
#ifdef EQMASS
                           mcm(ntotal)=1.d0
#else
                           mcm(ntotal)=mass(npp)
#endif

                           if(gravity) then
                              qxx(ntotal)=0.d0
                              qxy(ntotal)=0.d0
                              qxz(ntotal)=0.d0
                              qyy(ntotal)=0.d0
                              qyz(ntotal)=0.d0
                              qzz(ntotal)=0.d0
                              trq(ntotal)=0.d0
                           endif

                           dxx(ntotal)=0.d0
                           dxx2(ntotal)=0.d0

                           nnod(ntotal)=npointer(npp)
                        endif
                     endif
                     if(ntotal.ge.nt) then
                        write(*,*) 'Tree out of range!!',ntotal,nt
                        write(*,*) i,j,' - 2'
                        write(*,*) x(i),y(i),z(i)
                        write(*,*) x(j),y(j),z(j)
                        write(*,*) promro(i),promro(j)
                        write(*,*) h(i),h(j)
                        write(*,*) nvi(i),nvi(j)
                        stop
                     endif
15                enddo
13             enddo
12          enddo
            nh2(naux(j))=ntotal
            nh1(naux(j))=ntotal-nmmm+1
5        enddo
      enddo

      if(nsatura.ne.0) then
          print *,'Saturated tree',ntotal,nt-50
          print *,'Too depth tree',nnmax,ntree-10
          stop
      endif

      call profile(1)
      RETURN
    END SUBROUTINE buildtree
