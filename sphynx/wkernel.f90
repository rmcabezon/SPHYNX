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
!                 SPHYNX: wkernel.f90                 !
!                                                     !
! Interpolation kernel calculation.                   !
!=====================================================!


    SUBROUTINE Wkernel(vloc,nloc,wloc,dwloc)

      USE parameters

      IMPLICIT NONE
      INTEGER j
      DOUBLE PRECISION,INTENT(in)::vloc,nloc
      DOUBLE PRECISION,INTENT(out)::wloc,dwloc
      DOUBLE PRECISION v2,v3,pi2v,dmy,dmy1,dmy2

! WARNING!!!
! All dW/dv should be divided by an extra v which comes from the term dr/dx_i:
! Ex: dW/dx=dW/dv*dv/dr*dr/dx=W'*(1/h)*(x/r)=W'/(v*h^2)
! In this way we avoid a sqrt to calculate r in each kernel call.
! h^2 is taken into account in the main program as h5m(i) (for 3D), hence dwloc
! should always be returned divided by vloc.
! Take this into account when using dw in the main program (f.e: grad-h terms)
! where dw should be multiplied by an extra v.

      if(kernel.eq.1) then          !Cubic Spline

         if(vloc.le.1.d0) then
            v2=vloc*vloc
            v3=v2*vloc
            wloc=1.d0+ad1*v2+ad2*v3
            dwloc=adw1+adw2*vloc
         else if(vloc.gt.1.d0 .and. vloc.le.2.d0) then
            wloc=ad3*(2.d0-vloc)**3
            dwloc=adw3*(2.d0-vloc)**2
         else
            wloc=0.d0
            dwloc=0.d0
         endif

      else if (kernel.eq.2) then   !Harmonic kernel

         if(vloc.eq.0.d0) then
            wloc=1.d0
            dwloc=0.d0
         else if(vloc.le.2.d0) then
            pi2v=pihalf*vloc
            wloc=((sin(pi2v))/(pi2v))**nloc
            dwloc=nloc*pihalf*wloc*(1.d0/tan(pi2v)-1.d0/(pi2v))
         else
            wloc=0.d0
            dwloc=0.d0
         endif

      else if (kernel.eq.3) then   !Interpolated harmonic kernel

         if(vloc.eq.0.d0) then
            wloc=1.d0
            dwloc=0.d0
         else if (vloc.lt.2.d0) then
            j=1+int(vloc*10000)
            wloc=fk(j)+fdk1(j)*(vloc-uk(j))
            if(wloc.lt.0.d0) wloc=0.d0
            if(j.eq.20001) then
               dwloc=0.d0
            else
               dwloc=fdk2(j)
            endif
         else
            wloc=0.d0
            dwloc=0.d0
         endif

         wloc=wloc**nloc
         dwloc=nloc*wloc*dwloc

      else if (kernel.eq.4) then   !Wendland 3,3

         v2=vloc*vloc
         v3=v2*vloc
         dmy=max(1.d0-vloc/2.d0,0.d0)
         dmy1=dmy**7
         dmy2=4.d0*v3+6.25d0*v2+4.d0*vloc+1.d0
         wloc=dmy*dmy1*dmy2
         dwloc=-4.d0*dmy2*dmy1+dmy*dmy1*(12.d0*v2+12.5d0*vloc+4.d0)

      endif

      RETURN
    END SUBROUTINE Wkernel
