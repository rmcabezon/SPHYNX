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
!                SPHYNX: parameters.f90               !
!                                                     !
! Main module with all parameters and global variables!
!                                                     !
!=====================================================!

MODULE parameters

  IMPLICIT NONE
  DOUBLE PRECISION rad,rocentral,despl
  DOUBLE PRECISION dt,dtnma,tt,dtold
  DOUBLE PRECISION mtotal,cmax
  DOUBLE PRECISION ad1,ad2,ad3,adw1,adw2,adw3,pk4
  DOUBLE PRECISION refrad,refax,refay,refaz
  DOUBLE PRECISION maxtstepaux,timeinit,timefin,masspart

  INTEGER l,n,npini,npend,id,iti,nproc,ijk
  INTEGER ntotal,lastcoldim,new_comm

  LOGICAL writeout,findballmass
  LOGICAL checkdens,conduction

  INTEGER,PARAMETER::nmax = 1000000  !Maximum number of particles

  INTEGER,PARAMETER::n2=2*nmax       !Multiples of nmax
  INTEGER,PARAMETER::n3=3*nmax
  INTEGER,PARAMETER::n4=4*nmax
  INTEGER,PARAMETER::n5=5*nmax
  INTEGER,PARAMETER::n7=7*nmax

  INTEGER,PARAMETER::nodes=4         !Amount of nodes.
  INTEGER,PARAMETER::dim = 3         !Spatial imensions
  INTEGER,PARAMETER::nnl = 480       !Iterations
  INTEGER,PARAMETER::iodelay = 120   !Steps until output
  INTEGER,PARAMETER::nt = 3*nmax     !Max # of nodes in tree
  INTEGER,PARAMETER::ntree = nmax/2  !Max depth in tree
  INTEGER,PARAMETER::ncubes = 27     !Cubes for periodic BC
                                     !1 - No periodic BC
                                     !9 - Full periodic BC 2D
                                     !27- Full periodic BC 3D
  INTEGER,PARAMETER::kernel = 2      !Choosing kernel
                                     !1 - Cubic Spline
                                     !2 - Harmonic
                                     !3 - Harmonic interpolated
  LOGICAL,PARAMETER::oldnorm=.false. ! Choose old interpolation for norm
  INTEGER,PARAMETER::NRitermax = 10  !Maximum NR iterations for h
  INTEGER,PARAMETER::liniNR = nnl+1     !Iteration to start NR for h

  DOUBLE PRECISION,PARAMETER::scale=1.d0      !spatial scale (cm)
  DOUBLE PRECISION,PARAMETER::xbox=1.d0       !Box size x
  DOUBLE PRECISION,PARAMETER::ybox=1.d0       !Box size y
  DOUBLE PRECISION,PARAMETER::zbox=1.d0       !Box size z
  DOUBLE PRECISION,PARAMETER::ni0=100.d0      !Number of neighbors
  DOUBLE PRECISION,PARAMETER::nvmax=175.d0    !Maximum number of neighbors(Int.)
  DOUBLE PRECISION,PARAMETER::nvmin=30.d0     !Minimum number of neighbors(Int.)
  DOUBLE PRECISION,PARAMETER::fixfactor=3.d0  !Increased % for the neigh. array
  DOUBLE PRECISION,PARAMETER::hfac=1023.d0      !Factor for h update
  DOUBLE PRECISION,PARAMETER::hexp=1.d0/10.d0  !Exponent for h update
  DOUBLE PRECISION,PARAMETER::tol = 0.5d0     !Tolerance for opening nodes
  DOUBLE PRECISION,PARAMETER::dkcour = 0.2d0  !Factor for courant timestep
  DOUBLE PRECISION,PARAMETER::dkro = .06d0    !Factor for density timestep
  DOUBLE PRECISION,PARAMETER::dkacc = .02d0   !Factor for acceleration timestep
  DOUBLE PRECISION,PARAMETER::deltahindex=0.d0!Maximum index change
  DOUBLE PRECISION,PARAMETER::lambdac=0.5d0   !Parameter for index change
  DOUBLE PRECISION,PARAMETER::nbaseline=6.d0  !Minimum index
  DOUBLE PRECISION,PARAMETER::alfamax=1.d0    !Artificial viscosity constant
  DOUBLE PRECISION,PARAMETER::alfau=0.05d0    !Artificial viscosity thermal. Put it to zero if not wanted
  DOUBLE PRECISION,PARAMETER::alfamin=0.05d0  !AV floor (switches)
  DOUBLE PRECISION,PARAMETER::betaAV=2.d0     !Artificial viscosity constant beta
  DOUBLE PRECISION,PARAMETER::decay_constant=0.2     !Alfa decay (switches)
  DOUBLE PRECISION,PARAMETER::p0=1.d0         !Initial constant pressure
  DOUBLE PRECISION,PARAMETER::pext=2.d14 !2.7d13     !External pressure corr.
  DOUBLE PRECISION,PARAMETER::pmin=1.d2       !Minimum pressure corr.

  !Min. Atwood number in ramp function in momentum equation (crossed/uncrossed selection)
  !Complete uncrossed option (Atmin>=1.d50, Atmax it doesn't matter).
  !Complete crossed (Atmin and Atmax negative)
  DOUBLE PRECISION,PARAMETER::Atmin=0.1d0
  DOUBLE PRECISION,PARAMETER::Atmax=0.2d0
  DOUBLE PRECISION,PARAMETER::ramp=1.d0/(Atmax-Atmin)

  DOUBLE PRECISION,PARAMETER::maxtstep=1.d0   !Biggest allowed timestep

  LOGICAL,PARAMETER::flags =.true.       !Activate verbosity
  LOGICAL,PARAMETER::relax = .false.      !Activate angular relaxation
  LOGICAL,PARAMETER::gravity =.false.     !Calculate gravity
  LOGICAL,PARAMETER::estabil = .false.    !Write central density and radius
  LOGICAL,PARAMETER::conserv =.true.     !Tracking of conservation laws
  LOGICAL,PARAMETER::timesave =.true.    !Tracking of timesteps
  LOGICAL,PARAMETER::balsara = .false.    !Balsara corrections
  LOGICAL,PARAMETER::switches =.true.    !AV switches
  LOGICAL,PARAMETER::gradh = .true.       !Includes grad-h terms
  LOGICAL,PARAMETER::conductionini= .false.!Thermal conduction (not in connection with AV)
  LOGICAL,PARAMETER::std_VE= .false.       ! use xmass=mass(i) or new  xmass=mass(i)/ro0(i)
  LOGICAL,PARAMETER::outini=.true.        !Write output file at first iteration
  LOGICAL,PARAMETER::timecour =.true.    !Courant timestep
  LOGICAL,PARAMETER::timero =.true.      !Densidad timestep
  LOGICAL,PARAMETER::timeacc = .false.    !Acceleration timestep
  LOGICAL,PARAMETER::KH= .false.           !set KH test in EOS
  LOGICAL,PARAMETER::usepext= .false.      !Add external pressure
  LOGICAL,PARAMETER::noisetr = .false.    !Simplest version of noise trigger by Rosswog

  CHARACTER*43,PARAMETER::formatin='(37(1x,es23.16))' !Input format
  CHARACTER*43,PARAMETER::formatout='(1x,i7,45(1x,es12.5),1x,i2,1x,l1)'  !Output format

  !RESTART parameters
  INTEGER,PARAMETER::iterini = 1           !Begining of first iteration
  INTEGER,PARAMETER::delay = 0             !timesteps of delay to h adjustment
  DOUBLE PRECISION,PARAMETER::timeini = 0.d0      !Origin of time
  DOUBLE PRECISION,PARAMETER::deltaini = 1.d-6    !Initial timestep
  LOGICAL,PARAMETER::inivel0 = .false.            !Initial velocity=0
  LOGICAL,PARAMETER::inienergy = .false.          !Initialize U with EOS
  CHARACTER*20,PARAMETER::inputfile='sedov_glass.txt'!'sun_5Mv2' !Input file
  CHARACTER*100,PARAMETER::comments=&
       &'Sedov test in 3D'
  !EOS
  DOUBLE PRECISION, PARAMETER::gamma=5.d0/3.d0
  INTEGER,PARAMETER::eos = 1  !Choose EOS
                              !1 - Ideal gas with gamma with Temp
                              !2 - Ideal gas with gamma with Eint
                              !3 - Polytropic
                              !4 - Helmholtz

  CHARACTER*24,PARAMETER::EOStab='../eos/eosLS220'      !Tabulated EOS
  DOUBLE PRECISION,PARAMETER::kpol = 42479.0909309992d0 !Polytropic constant
  DOUBLE PRECISION,PARAMETER::gammapol = 2.d0           !Polytropic exponent

  !Nuclear species for EOS
  DOUBLE PRECISION,DIMENSION(nmax)::abar_eos,zbar_eos

  !Physical and math constants
  DOUBLE PRECISION,PARAMETER::g=6.6726d-8     !Gravitation constant
!  DOUBLE PRECISION,PARAMETER::g=1.d0     !Gravitation constant
  DOUBLE PRECISION,PARAMETER::cvel=2.997924562d10       !Speed of light
  DOUBLE PRECISION,PARAMETER::pi=3.14159265358979323846d0
  DOUBLE PRECISION,PARAMETER::pihalf=pi/2.d0
  DOUBLE PRECISION,PARAMETER::pim1=1.d0/pi
  DOUBLE PRECISION,PARAMETER::third=1.d0/3.d0
  DOUBLE PRECISION,PARAMETER::twothird=2.d0/3.d0
  DOUBLE PRECISION,PARAMETER::msun=1.989d33
  DOUBLE PRECISION,PARAMETER::rsun=6.9598d10
  DOUBLE PRECISION,PARAMETER::MeV=1.602192d-6
  DOUBLE PRECISION,PARAMETER::mb=1.674d-24    !Baryon mass
  DOUBLE PRECISION,PARAMETER::K2MeV=8.617342857d-11 !Conversion from K to MeV
  DOUBLE PRECISION,PARAMETER::rgasid=8.317d7


  !Arrays
  INTEGER,DIMENSION(nmax)::nvi,indx,ind
  INTEGER,DIMENSION(0:nt)::nh1,nh2
  INTEGER,DIMENSION(nt)::nnod,npa
  INTEGER,ALLOCATABLE,DIMENSION(:,:)::neighbors,cube

#ifdef EQMASS
#else
  DOUBLE PRECISION,DIMENSION(nmax)::mass
#endif
  DOUBLE PRECISION,DIMENSION(nmax)::temp,promro,ro2,u,p,c,pro,promro0
  DOUBLE PRECISION,DIMENSION(nmax)::mue,mui,dudt,dudv,dpdt,ro0,sumwhro0
  DOUBLE PRECISION,DIMENSION(nmax)::s,alfa,TdPdTro,ye
  DOUBLE PRECISION,DIMENSION(nmax)::mindist,avgro,xmass,sumkx
  DOUBLE PRECISION,DIMENSION(nmax)::h,hd,h2,h3,hm1,hm2,hm3,hm4,hm5,omega
  DOUBLE PRECISION,DIMENSION(nmax)::energy,avisc,radius,sumwh,ballmass
  DOUBLE PRECISION,DIMENSION(nmax)::vol,energyT,energyTold
  DOUBLE PRECISION,DIMENSION(nmax)::energyold,aviscold,aviscu,aviscuold
  DOUBLE PRECISION,DIMENSION(nmax)::acmod,x3li
  DOUBLE PRECISION,DIMENSION(nmax)::c11,c12,c13,c22,c23,c33
  DOUBLE PRECISION,DIMENSION(nmax)::magcorr,gradcorr,maxvsignal
  DOUBLE PRECISION,DIMENSION(nmax)::kappa,kappa0,vturb
  DOUBLE PRECISION,DIMENSION(nmax)::indice,pk,dpk,mark_ramp
  DOUBLE PRECISION,DIMENSION(nmax)::checkInorm,checkdeltax,checkdeltay,checkdeltaz


  DOUBLE PRECISION,DIMENSION(20001)::uk,fk,fdk1,fdk2
  DOUBLE PRECISION,dimension(25)::timepass

  DOUBLE PRECISION,DIMENSION(nmax)::diff

  DOUBLE PRECISION,DIMENSION(nt)::xce,yce,zce,dxx,dxx2,xmc,ymc,zmc
  DOUBLE PRECISION,DIMENSION(nt)::qxx,qxy,qxz,qyy,qyz,qzz,trq,mcm

  DOUBLE PRECISION,DIMENSION(nmax*dim)::a,a0,v,v0,f,ac,gradp

  DOUBLE PRECISION,DIMENSION(dim)::cm

  DOUBLE PRECISION,DIMENSION(nmax)::ugrav,divv,curlv
  DOUBLE PRECISION,DIMENSION(nmax)::dvxdx,dvxdy,dvxdz,dvydx,dvydy,dvydz,dvzdx,dvzdy,dvzdz

  LOGICAL,DIMENSION(nmax) ::ready
  LOGICAL,DIMENSION(nodes)::localready

END MODULE parameters
