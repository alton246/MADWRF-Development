program sample_read 

! Fortran 90 version. 
!   This is a simple program to write data in the WPS intermediate 
!   format.  It is included mainly as an aid to understanding said 
!   format. 

  implicit none 

! Declarations: 

  integer, parameter :: IUNIT = 10 
  integer, parameter :: OUNIT = 11 
  integer :: ierr=0 

  integer :: IFV=5 

  character(len=24) :: HDATE = '2020-06-22_16:00:00' 
  real :: XFCST=0.000000 
  character(len=8) :: STARTLOC='SWCORNER' 
  character(len=9) :: FIELD= 'CLDMASK'
  character(len=25) :: UNITS='Unknown' 
  character(len=46) :: DESC='1-Cloud 0-No Cloud' 
  character(len=32) :: MAP_SOURCE='GOES Imagery' 
  real :: XLVL=100 
  integer :: NX=2466 
  integer :: NY=1295 
  integer :: IPROJ=0 
  real :: STARTLAT=-6.97753 
  real :: STARTLON=-99.99886 
  real :: DELTALAT= 0.01938 
  real :: DELTALON= 0.02101 
  real :: DX=2.00025 
  real :: DY=2.11138 
  real :: XLONC=-76.58453 
  real :: TRUELAT1= 12.941 
  real :: TRUELAT2= 0 
  real :: NLATS=1295
  real :: EARTH_RADIUS = 6367470. * .001 
  logical :: IS_WIND_EARTH_REL = .FALSE. 
! SLAB is an allocatable array, because we do not necessarily know in 
! advance the size of the array we need to read. 
  real, allocatable, dimension(:,:) :: SLAB 

  SLAB = reshape( (/ 