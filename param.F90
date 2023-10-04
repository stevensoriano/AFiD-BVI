!===========================================================
! Declaration of global variables
!***********************************************************      
      module param
        implicit none
!==========================================================
!       read from input file bou.in
!==========================================================
        integer   :: n2, n3,n1
        integer   :: nsst, nwrit, nread, ntst, ireset
        real      :: tframe,tpin,tmax,walltimemax
        real      :: xlen, ylen, zlen
        real      :: pra,dt,resid,cflmax,tsta
        integer   :: starea
        real      :: dtmax,cfllim 
        real      :: tl,epsstar,kfmax
        integer   :: nson,idtv,imovie
        integer   :: imlsfor,imlsstr,imlsref,pread
        integer   :: wcheck,fcheck,wcub,wexp,nel
        real      :: wcon, wscl, sclf, rhop,tmovie
        character*50  gtsfx,datfx


!=================================================
!       end of input file
!=================================================
        real :: time
!******* Grid parameters**************************
        real :: dx2,dx3,dx1
        real :: dx2q,dx3q,dx1q
         
        real, allocatable, dimension(:) :: xc,xm
        real, allocatable, dimension(:) :: yc,ym
        real, allocatable, dimension(:) :: zc,zm
!==========================================================
!******* Grid indices**************************************
        integer, allocatable, dimension(:) :: jmv,jpv
        integer, allocatable, dimension(:) :: imv,ipv
        integer, allocatable, dimension(:) :: jmhv
        integer, allocatable, dimension(:) :: kmv,kpv
!============================================================
!******* Variables for FFTW and Poisson solver****************
        real, dimension(13) :: ifx1
        integer*8 :: fwd_plan,bck_plan
        integer*8 :: fwdplan_1d,bckplan_1d
        real, allocatable, dimension(:) :: ao,ap,af
        real, allocatable, dimension(:) :: ak1,ak2,ak3
        
!===========================================================
!******* Other variables ***********************************
        integer  :: n2m, n3m, n1m,n1mh,n2mh
        integer  :: iaxsy
        real :: keta
        real :: ren, pec
        real :: pi
        real :: al,ga,ro
        real :: beta
        real :: qqmax,qqtot
        real :: re
        real :: vortexrad, wirev
        integer :: ntime
        real, dimension(1:3) :: vmax
        real, dimension(1:3) :: gam,rom,alm
        complex,allocatable,dimension(:,:) :: term1a,term1b
        complex,allocatable,dimension(:,:) :: term2a,term2b
        complex,allocatable,dimension(:,:) :: term3a,term3b
        real, dimension(6,7,7,7) :: bcoefs
        logical :: ismaster = .false.
        logical :: solvestructure = .false.
        logical :: mlsforcing = .false.
      end module param
      
!************* End of param module******************************
!===============================================================
!******* 2D arrays, dynamically allocated by each process*******
      module local_arrays
      use param
        implicit none
        real,allocatable,dimension(:,:,:) :: vx,vy,vz
        real,allocatable,dimension(:,:,:) :: qbuf,forcx,forcy,forcz
        real,allocatable,dimension(:,:,:) :: pr,rhs
        real,allocatable,dimension(:,:,:) :: ru1,ru2,ru3
        real,allocatable,dimension(:,:,:) :: qcap
        real,allocatable,dimension(:,:,:) :: dph,dq
      end module local_arrays

!===============================================================
      module stat_arrays
       implicit none
       real,allocatable, dimension(:,:,:) :: vx_me,vx_rms,pr_me,pr_rms
       real,allocatable, dimension(:,:,:) :: vy_me,vz_me,vy_rms,vz_rms 
       integer :: timeint_cdsp
        real :: vxvyvz_rms_vol
      end module stat_arrays
!=====================================================       
      module mpih
        implicit none
        include 'mpif.h'
        integer :: myid, numtasks, ierr
        integer, parameter :: master=0
        integer, parameter :: lvlhalo=3
        integer :: MDP = MPI_DOUBLE_PRECISION
        integer :: MCP = MPI_DOUBLE_COMPLEX
        integer :: STATUS(MPI_STATUS_SIZE,4)
        integer :: req(1:4)
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
      end module mpih
      
      module mpi_param
        implicit none
        integer :: istart,iend,jstart,jend, kstart,kend
        integer :: jstartp,jendp
        integer :: dj,dk,mydata,mydatam
        integer :: djp
        integer, allocatable, dimension(:) :: offsetj,offsetk
        integer, allocatable, dimension(:) :: offsetjp
        integer, allocatable, dimension(:) :: countj,countk
        integer, allocatable, dimension(:) :: countjp
        integer, allocatable, dimension(:) :: countf
        integer(8), allocatable, dimension(:) :: offsetf 
      end module mpi_param
!====================================================
      module mls_param
        use param
        implicit none
        integer :: Nparticle
        parameter (Nparticle=1)
!     --------variables for structural solver------------------------
      integer, parameter :: max_n_edge_of_vert=50

      integer, dimension(:,:),allocatable :: n_edge_of_vert

      integer, dimension(:,:,:), allocatable :: v_cell
      integer, dimension(:,:,:), allocatable :: vert_of_edge
      integer, dimension(:,:,:), allocatable :: vert_of_face
      integer, dimension(:,:,:), allocatable :: edge_of_face
      integer, dimension(:,:,:), allocatable :: vert_of_vert
      integer, dimension(:,:,:), allocatable :: edge_of_vert
      integer, dimension(:,:,:), allocatable :: face_of_edge
      integer, dimension(:,:,:), allocatable :: v1234

      integer, dimension(:,:),allocatable::dum_for
      integer, dimension(:,:,:),allocatable::pind
      integer,dimension(:,:,:),allocatable::bboxind
 
      real,dimension(:,:),allocatable::theta0,theta
      real,dimension(:,:),allocatable::dist0,dist
      real,dimension(:,:),allocatable::sur0,sur
      real,dimension(:,:),allocatable::dismax
      real,dimension(:,:),allocatable::mvol,press_face,vforc_face

      real,dimension(:),allocatable::Volume0,Volume
      real,dimension(:),allocatable::Surface0,Surface
      real,dimension(:),allocatable::shwtx,shwty,shwtz1,shwtz2

      real, dimension(:,:,:), allocatable :: xyz,xyz0
      real, dimension(:,:,:), allocatable :: xyzv,xyzv1,xyza
      real, dimension(:,:,:), allocatable :: xyza0,dpdn

      real,dimension(:,:,:),allocatable::tri_ver,vel_tri,acc_tri
      real,dimension(:,:,:),allocatable::tri_bar,tri_nor
      
      real,dimension(:,:,:),allocatable::fpxyz,fexyz,fbxyz,fvxyz
      real,dimension(:,:,:),allocatable::fatxyz,falxyz,fxyz
      real,dimension(:,:,:),allocatable::fsxyz

      real, dimension(:,:),allocatable :: sca_nod1,sca_nod2,sca_nod3

      integer :: maxnv,maxne,maxnf,nsstep
      real    :: ke,kb,kv,kat,kal,kvis,kcol,Cv,uspm
      real    :: thickness,rhos,usVolele
      real    :: celvol,Hboxx

      end module mls_param

!====================================================
      module mls_local
        use param
        implicit none
        real, allocatable, dimension(:,:,:) :: for_xc, for_yc, for_zc
      end module mls_local
!====================================================
      module pointparticle
       use param
       implicit none
       logical :: withppart = .false.
       integer :: Npointpart,Onparticle
       integer :: iresetp
       real  :: timeONp,p1p2,toutpp
       real :: usfroude,dbd,dbd2,rhohat1,rhohat2,stokes1,stokes2

       real, allocatable, dimension(:,:,:) :: for_x_part,for_y_part, &
                                             for_z_part
       real, allocatable, dimension(:,:,:) :: for_xc_part,for_yc_part,&
                                             for_zc_part
       real, dimension(:), allocatable :: qVal1,qVal2,qVal3      ! for timestepping
       real, dimension(:), allocatable :: qValo1,qValo2,qValo3   ! for timestepping
       real, dimension(:), allocatable :: kalb1,kalb2,kalb3
       real, dimension(:), allocatable :: kalbo1,kalbo2,kalbo3
       real, dimension(:), allocatable :: vort1,vort2,vort3
       real, dimension(:), allocatable :: vorto1,vorto2,vorto3!
       real, dimension(:,:), allocatable :: xp,xpo
       real, dimension(:), allocatable :: q1p,q2p,q3p    ! part vel new
       real, dimension(:), allocatable :: q1po,q2po,q3po ! part vel old
       real, dimension(:), allocatable :: renp       ! part Re number
       real, dimension(:,:), allocatable :: aap      ! part acc
       real, dimension(:,:), allocatable :: facc_for,drag_for,lift_for

       real, dimension(:), allocatable :: dbd12,stokes
       real, dimension(:), allocatable :: gammap
       real cpi(3), cpf(3)


 
      end module pointparticle
!====================================================
      module local_aux
       use param
       implicit none
       real, allocatable, dimension(:,:,:) :: vxc, vyc, vzc
       real, allocatable, dimension(:,:,:) :: matderxc, matderx
       real, allocatable, dimension(:,:,:) :: matderyc, matdery
       real, allocatable, dimension(:,:,:) :: matderzc, matderz
       real,allocatable,dimension(:,:,:) :: vorx, vory, vorz
       real,allocatable,dimension(:,:,:) :: vxo, vyo, vzo
      end module local_aux
