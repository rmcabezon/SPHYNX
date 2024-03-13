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
!                  SPHYNX: eosmod.f90                 !
!                                                     !
! Calls different EOS.                                !
!=====================================================!

      SUBROUTINE eostot

      USE parameters

      IMPLICIT NONE
      INTEGER i
      DOUBLE PRECISION, DIMENSION(nmax):: ueos
      INTEGER ii

      call profile(0)
      if(flags)print *,'Calculating EOS'

      p=0.d0
      c=0.d0
      pro=0.d0
      dpdt=0.d0
      dudt=0.d0
      dudv=0.d0
      ueos=0.d0
      TdPdTro=0.d0

      if (eos.eq.1) then

         call eosid_t(ueos)
         u(:)=ueos(:)
         !input: temp,promro,mui
         !This populates p,u,c,dudt,dpdt,dpdr

      else if (eos.eq.2) then
         if (l.eq.iterini.and.inienergy) then
            call eosid_t(ueos)
         else
            call eosid_u
         endif
         !input: promro,u,gamma
         !This populates p,temp,c,cv,dpdt

       else if (eos.eq.3) then

          p(:)=kpol*promro(:)**gammapol
          c(:)=sqrt(p(:)/promro(:))


       else if (eos.eq.4) then

         do i=1,n
            if(temp(i).ge.1.d3) then
                call helmeos(temp(i),promro(i),abar_eos(i),zbar_eos(i),&
                            & p(i),dpdt(i),u(i),dudt(i),c(i))
            else   !use Id.Gas id T is out of limits for Helmholtz.
                p(i)=rgasid*promro(i)*temp(i)/mui(i)
                u(i)=3.d0/2.d0*rgasid*temp(i)/mui(i)
                dudt(i)=u(i)/temp(i)
                dpdt(i)=rgasid*promro(i)/mui(i)
                c(i)=sqrt(p(i)/promro(i))
            endif
         end do
       else

         write(*,*) 'Not defined EOS!'
         stop

      endif

      cmax=maxval(c)
      if(l.eq.iterini.and.inienergy)u=ueos

      if(std_VE) then
         pro(:)=p(:)!/(sumkx(:))**2
         TdPdTro(:)=temp(:)*dpdt(:)
      else
         pro(:)=p(:)/(sumkx(:))
         TdPdTro(:)=temp(:)*dpdt(:)/(sumkx(:))
      endif

      call profile(1)

      RETURN
      END SUBROUTINE eostot
