MODULE CMF_NMLIST_MOD
!==========================================================
!* PURPOSE: read CaMa-Flood Model configulations from namelist ("input_flood.nam" as default)
!
!* CONTAINS:
! -- CMF_CONFIG_NAMELIST  : read namelist for CaMa-Flood 
! -- CMF_CONFIG_CHECK     : check config conflict
!
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  29Jul 2019
!             Adapted mostly from CMF v362 CONTROL0.F90
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
!==========================================================
!* PURPOSE: Shared variables for CaMa-Flood model configulation
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
  USE PARKIND1,                ONLY: JPIM, JPRB, JPRM
  IMPLICIT NONE
  SAVE 
  ! Type definitions
  TYPE :: CaMaFiles
      LOGICAL :: LLOGOUT = .TRUE.                    !! true: log output to file
      INTEGER(KIND=JPIM) :: LOGNAM                   !! default log    file FID
      INTEGER(KIND=JPIM) :: NSETFILE                 !! input namelist file FID
      INTEGER(KIND=JPIM) :: TMPNAM                   !! temporal I/O   file FIG
      CHARACTER(LEN=256) :: CLOGOUT = './log_CaMa.txt'  !! default log file name
      CHARACTER(LEN=256) :: CSETFILE = 'input_cmf.nam'  !! input namelist file name
  END TYPE CaMaFiles
  
  TYPE :: CaMaOptions
      ! Core simulation options
      LOGICAL :: LADPSTP = .TRUE.                     !! true: use adaptive time step
      LOGICAL :: LFPLAIN = .TRUE.                     !! true: consider floodplain
      LOGICAL :: LKINE = .FALSE.                      !! true: use kinematic wave
      LOGICAL :: LFLDOUT = .TRUE.                     !! true: floodplain flow active
      LOGICAL :: LPTHOUT = .TRUE.                     !! true: activate bifurcation scheme
      LOGICAL :: LDAMOUT = .FALSE.                    !! true: activate dam operation
      LOGICAL :: LLEVEE = .FALSE.                     !! true: activate levee scheme

      ! ECMWF specific options
      LOGICAL :: LROSPLIT = .FALSE.                   !! true: split surface/sub-surface runoff
      LOGICAL :: LWEVAP = .TRUE.                      !! true: input water evaporation
      LOGICAL :: LWEVAPFIX = .TRUE.                   !! true: water balance closure
      LOGICAL :: LWEXTRACTRIV = .FALSE.               !! true: extract water from rivers
      LOGICAL :: LSLOPEMOUTH = .FALSE.                !! true: prescribe water level slope
      LOGICAL :: LGDWDLY = .FALSE.                    !! true: ground water reservoir
      LOGICAL :: LSLPMIX = .FALSE.                    !! true: mixed kinematic/inertia
  
      ! Boundary and initial conditions
      LOGICAL :: LMEANSL = .FALSE.                    !! true: mean sea level BC
      LOGICAL :: LSEALEV = .FALSE.                    !! true: variable sea level BC
      LOGICAL :: LRESTART = .FALSE.                   !! true: restart from file
      LOGICAL :: LSTOONLY = .FALSE.                   !! true: storage only restart
  
      ! Output options
      LOGICAL :: LOUTPUT = .TRUE.                     !! true: standard output
      LOGICAL :: LOUTINI = .FALSE.                    !! true: output initial storage
      LOGICAL :: LOUTINS = .FALSE.                    !! true: instantaneous discharge
  
      ! Technical options
      LOGICAL :: LGRIDMAP = .TRUE.                    !! true: XY gridded 2D map
      LOGICAL :: LLEAPYR  = .TRUE.                     !! true: neglect leap year
      LOGICAL :: LMAPEND  = .FALSE.                    !! true: map data endian conversion
      LOGICAL :: LBITSAFE = .FALSE.                    !! true: for Bit Identical (removed from v410, set in Mkinclude)
      LOGICAL :: LSTG_ES  = .FALSE.                    !! true: Vector Processor opt
  
      ! Additional schemes
      LOGICAL :: LSEDOUT     = .FALSE.                 !! true: sediment scheme
      LOGICAL :: LTRACE      = .FALSE.                 !! true: tracer scheme
      LOGICAL :: LWINFILT    = .FALSE.                 !! true: water infiltration 
      LOGICAL :: LWINFILTFIX = .FALSE.                 !! true: infiltration fix
  END TYPE CaMaOptions
  
  TYPE :: CaMaConfig
      CHARACTER(LEN=256) :: CDIMINFO="NONE"                 !! Dimension Information
      REAL(KIND=JPRB)    :: DT     = 24*60*60                       !! Time Step Length [SEC]
      INTEGER(KIND=JPIM) :: IFRQ_INP =24               !! Runoff update frequency [hour]
      REAL(KIND=JPRB)    :: DTIN   = 24*60*60               !! Input time step [SEC]
      INTEGER(KIND=JPIM) :: NX     = 1440
      INTEGER(KIND=JPIM) :: NY     = 720
      INTEGER(KIND=JPIM) :: NLFP   = 10           !! Domain dimensions
      INTEGER(KIND=JPIM) :: NXIN   = 360
      INTEGER(KIND=JPIM) :: NYIN   = 180
      INTEGER(KIND=JPIM) :: INPN   = 1      !! Input dimensions
      REAL(KIND=JPRB)    :: WEST   = -180._JPRB 
      REAL(KIND=JPRB)    :: EAST   = 180._JPRB
      REAL(KIND=JPRB)    :: NORTH  = 90._JPRB
      REAL(KIND=JPRB)    :: SOUTH  = -90._JPRB  !! Domain edges [deg]
  END TYPE CaMaConfig

  TYPE :: CaMaParams
      ! Physical parameters
      REAL(KIND=JPRB) :: PMANRIV    =    0.03_JPRB             !! Manning (river)
      REAL(KIND=JPRB) :: PMANFLD    =   0.10_JPRB                  !! Manning (floodplain)
      REAL(KIND=JPRB) :: PGRV       = 9.8_JPRB                 !! Gravity [m/s2]
      REAL(KIND=JPRB) :: PDSTMTH    =10000._JPRB                  !! Downstream distance [m]
      REAL(KIND=JPRB) :: PCADP      =0.7_JPRB                  !! CFL coefficient
      REAL(KIND=JPRB) :: PMINSLP   =1.E-5                     !! Min slope [m/m]
      
      ! Missing values
      INTEGER(KIND=JPIM) :: IMIS   =-9999_JPIM                 !! Integer undefined
      REAL(KIND=JPRM) :: RMIS       =1.E20_JPRM                !! Real undefined
      REAL(KIND=JPRB) :: DMIS    =1.E20_JPRB                   !! Double undefined
      
      ! File suffixes
      CHARACTER(LEN=256) :: CSUFBIN  ='.bin'               !! Binary suffix
      CHARACTER(LEN=256) :: CSUFVEC  ='.vec'               !! Vector suffix
      CHARACTER(LEN=256) :: CSUFPTH  ='.pth'               !! Path suffix
      CHARACTER(LEN=256) :: CSUFCDF  ='.nc'              !! NetCDF suffix
  END TYPE CaMaParams

  TYPE :: CaMaTime
      INTEGER(KIND=JPIM) :: syear = 2000    !! START YEAR
      INTEGER(KIND=JPIM) :: smon  = 1       !! START MONTH
      INTEGER(KIND=JPIM) :: sday  = 1       !! START DAY
      INTEGER(KIND=JPIM) :: shour = 0       !! START HOUR
      INTEGER(KIND=JPIM) :: eyear = 2001    !! END   YEAR
      INTEGER(KIND=JPIM) :: emon  = 1       !! END   MONTH
      INTEGER(KIND=JPIM) :: eday  = 1       !! END   DAY
      INTEGER(KIND=JPIM) :: ehour = 0       !! END   HOUR 
      INTEGER(KIND=JPIM) :: YYYY0 = 1       !! 
      INTEGER(KIND=JPIM) :: MM0   = 1       !! 
      INTEGER(KIND=JPIM) :: DD0   = 1       !! 
  END TYPE CaMaTime

  TYPE :: CaMaMaps
      ! Map file format
      LOGICAL :: LMAPCDF = .FALSE.              !! true for netCDF map input
      
      ! Binary map files
      CHARACTER(LEN=256) :: CNEXTXY = "./nextxy.bin"     !! river network nextxy
      CHARACTER(LEN=256) :: CGRAREA = "./ctmare.bin"     !! catchment area
      CHARACTER(LEN=256) :: CELEVTN = "./elevtn.bin"     !! bank top elevation
      CHARACTER(LEN=256) :: CNXTDST = "./nxtdst.bin"     !! distance to next outlet
      CHARACTER(LEN=256) :: CRIVLEN = "./rivlen.bin"     !! river channel length
      CHARACTER(LEN=256) :: CFLDHGT = "./fldhgt.bin"     !! floodplain elevation profile
      
      ! River channel parameters
      CHARACTER(LEN=256) :: CRIVWTH = "./rivwth.bin"     !! channel width
      CHARACTER(LEN=256) :: CRIVHGT = "./rivhgt.bin"     !! channel depth
      CHARACTER(LEN=256) :: CRIVMAN = "./rivman.bin"     !! river manning coefficient
      
      ! Optional maps
      CHARACTER(LEN=256) :: CPTHOUT = "./bifprm.txt"     !! bifurcation channel table
      CHARACTER(LEN=256) :: CGDWDLY = "NONE"             !! Groundwater Delay Parameter
      CHARACTER(LEN=256) :: CMEANSL = "NONE"             !! mean sea level
      CHARACTER(LEN=256) :: CMPIREG = "NONE"             !! MPI region map
      
      ! netCDF maps
      CHARACTER(LEN=256) :: CRIVCLINC = "NONE"           !! river map netcdf
      CHARACTER(LEN=256) :: CRIVPARNC = "NONE"           !! river parameter netcdf
      CHARACTER(LEN=256) :: CMEANSLNC = "NONE"           !! mean sea level netCDF  
      CHARACTER(LEN=256) :: CMPIREGNC = "NONE"           !! MPI region map in netcdf
  END TYPE CaMaMaps

  TYPE :: CaMaForcing
      ! Forcing configuration
      LOGICAL :: LINPCDF = .FALSE.              !! true: netCDF runoff forcing
      LOGICAL :: LINPEND = .FALSE.              !! true: input endian conversion
      LOGICAL :: LINTERP = .FALSE.              !! true: runoff interpolation using input matrix
      LOGICAL :: LITRPCDF = .FALSE.             !! true: netCDF input matrix file

      ! File paths and names
      CHARACTER(LEN=256) :: CINPMAT = "NONE"    !! Input matrix filename
      CHARACTER(LEN=256) :: CROFDIR = "./runoff/"!! Forcing: runoff directory
      CHARACTER(LEN=256) :: CROFPRE = "Roff____"!! Forcing: runoff prefix
      CHARACTER(LEN=256) :: CROFSUF = ".one"    !! Forcing: runoff suffix
      CHARACTER(LEN=256) :: CSUBDIR = "./runoff/"!! Forcing: sub-surface runoff directory
      CHARACTER(LEN=256) :: CSUBPRE = "Rsub____"!! Forcing: sub-surface runoff prefix
      CHARACTER(LEN=256) :: CSUBSUF = ".one"    !! Forcing: sub-surface runoff suffix
      CHARACTER(LEN=256) :: CROFCDF = "NONE"    !! Netcdf forcing file
      CHARACTER(LEN=256) :: CVNROF = "runoff"   !! NetCDF VARNAME of runoff
      CHARACTER(LEN=256) :: CVNSUB = "NONE"     !! NetCDF VARNAME of sub-surface runoff

      ! Unit conversion
      REAL(KIND=JPRB) :: DROFUNIT = 86400._JPRB*1000._JPRB  !! runoff unit conversion

      ! Time information
      INTEGER(KIND=JPIM) :: SYEARIN = 0         !! START YEAR IN NETCDF INPUT RUNOFF
      INTEGER(KIND=JPIM) :: SMONIN = 0          !! START MONTH IN NETCDF INPUT RUNOFF
      INTEGER(KIND=JPIM) :: SDAYIN = 0          !! START DAY IN NETCDF INPUT RUNOFF
      INTEGER(KIND=JPIM) :: SHOURIN = 0         !! START HOUR IN NETCDF INPUT RUNOFF
  END TYPE CaMaForcing

  TYPE :: CaMaBoundary
    ! Configuration
    LOGICAL :: LSEALEVCDF = .FALSE.              !! true: netCDF sea level boundary

    ! Plain binary data paths
    CHARACTER(LEN=256) :: CSEALEVDIR = "./sealev/"!! Sea level boundary DIRECTORY
    CHARACTER(LEN=256) :: CSEALEVPRE = "sealev"  !! Sea level boundary PREFIX
    CHARACTER(LEN=256) :: CSEALEVSUF = ".bin"    !! Sea level boundary SUFFIX

    ! netCDF data
    CHARACTER(LEN=256) :: CSEALEVCDF = "./sealev/"!! Sea level netCDF file name
    CHARACTER(LEN=256) :: CVNSEALEV = "variable" !! Sea Level netCDF variable name
    CHARACTER(LEN=256) :: CSLMAP = "./sealev/"   !! Conversion table (Sta -> XY)

    ! Time information
    INTEGER(KIND=JPIM) :: SYEARSL = 0           !! Start YEAR of netCDF sea level
    INTEGER(KIND=JPIM) :: SMONSL = 0            !! Start MONTH of netCDF sea level
    INTEGER(KIND=JPIM) :: SDAYSL = 0            !! Start DAY of netCDF sea level
    INTEGER(KIND=JPIM) :: SHOURSL = 0           !! Start HOUR of netCDF sea level

    ! For interpolation (netCDF only)
    INTEGER(KIND=JPIM) :: NLINKS = 0            !! Number of sea level station
    INTEGER(KIND=JPIM) :: NCDFSTAT = 0          !! Number of stations in netCDF

    INTEGER(KIND=JPIM) :: IFRQ_SL = 9999          !! default: dynamic sea level not used
    REAL(KIND=JPRB)                 :: DTSL                    !! SECOND IN TIME STEP [SEC]

  END TYPE CaMaBoundary
  
  TYPE :: CaMaRestart
      ! Input/Output file settings
      CHARACTER(LEN=256) :: CRESTSTO = "restart"    !! Input restart file name
      CHARACTER(LEN=256) :: CRESTDIR = "./"         !! Output restart file directory
      CHARACTER(LEN=256) :: CVNREST = "restart"     !! Output restart prefix

      ! Configuration flags
      LOGICAL :: LRESTCDF = .FALSE.                 !! true: netCDF restart file
      LOGICAL :: LRESTDBL = .TRUE.                  !! true: binary restart double precision

      ! Frequency settings
      INTEGER(KIND=JPIM) :: IFRQ_RST = 0           !! 0: only end of simulation
                                                     !! [1,2,3,6,12,24]: at selected hour
                                                     !! 30: monthly
  END TYPE CaMaRestart

  TYPE :: CaMaDam
    ! File paths and configuration
    CHARACTER(LEN=256) :: CDAMFILE = "./dam_params.csv"  !! Dam parameter file
    ! Operation flags
    LOGICAL :: LDAMTXT = .TRUE.     !! true: dam inflow-outflow txt output
    LOGICAL :: LDAMH22 = .FALSE.    !! true: Use Hanazaki 2022 scheme
    LOGICAL :: LDAMYBY = .FALSE.    !! true: Use Year-By-Year dam activation
    LOGICAL :: LiVnorm = .FALSE.    !! true: initialize dam storage with Normal Volume
  END TYPE CaMaDam

  TYPE :: CaMaLevee
      CHARACTER(LEN=256)             ::  CLEVHGT         !! LEVEE HEIGHT from RIVER
      CHARACTER(LEN=256)             ::  CLEVFRC         !! Unprotected fraction. Relative Levee distance from RIVER
  END TYPE CaMaLevee

  TYPE :: CaMaTracer
    INTEGER(KIND=JPIM)              :: NTRACE      !! number of tracer
    CHARACTER(LEN=256)              :: CTRCNAM     !! tracer name
    CHARACTER(LEN=256)              :: CTRCDIR     !! tracer input file directory
    CHARACTER(LEN=256)              :: CTRCPRE     !! tracer file prefix
    CHARACTER(LEN=256)              :: CTRCSUF     !! tracer file suffix
    
    INTEGER(KIND=JPIM)              :: IFRQ_TRIN   !! tracer input frequency (hour)
    REAL(KIND=JPRM)                 :: DTRCUNIT    !! tracer input unit conversion (DTRCUNIT=1 when tracer input file unit is [(MASS)/m2/s]. )
    LOGICAL                         :: LINPEND     !! true  for input    endian conversion
    
    LOGICAL                         :: LTRCBIF     !! true  for consider bifurcation in tracer scheme
    
    CHARACTER(LEN=256)              :: CRESTTRC               ! input restart file name
    CHARACTER(LEN=256)              :: CRESTDIR               ! output restart file directory
    CHARACTER(LEN=256)              :: CVNRSTTRC              ! output restart prefix
    LOGICAL                         :: LRESTDBL               ! true: binary restart in double precision
    INTEGER(KIND=JPIM)              :: IFRQ_RST               ! 0: only at last time, (1,2,3,...,24) hourly restart, 30: monthly restart
    
    ! output
    CHARACTER(LEN=256)              :: COUTDIR           ! OUTPUT DIRECTORY
    CHARACTER(LEN=256)              :: COUTTAG           ! Output Tag Name for each experiment
    LOGICAL                         :: LOUTVEC           ! TRUE FOR VECTORIAL OUTPUT, FALSE FOR NX,NY OUTPUT
    REAL(KIND=JPRB)                 :: DTIN_TRC
  END TYPE CaMaTracer

  TYPE :: CaMaSed
    character(len=256)              :: crocdph
    character(len=256)              :: sedD
    real(kind=JPRB)                 :: lambda           ! porosity (default:0.4)
    real(kind=JPRB)                 :: lyrdph           ! exchange layer depth
    integer(kind=JPIM)              :: nsed             ! number of sediment particle size
    real(kind=JPRB)                 :: sedDT             ! sediment timestep (s)
    real(kind=JPRB)                 :: psedD            ! density of sediment (default:2.65g/m3)
    real(kind=JPRB)                 :: DSYLUNIT         ! unit conversion factor for sediment (default:1e-6)
    real(kind=JPRB)                 :: pset             ! parameter for setting velocity
    real(kind=JPRB)                 :: pwatD            ! density of water (default:1g/m3)
    logical                         :: revEgia          ! if use Egiazoroff
    integer(kind=JPIM)              :: totlyrnum        ! number of deposition layers
    real(kind=JPRB)                 :: visKin           ! viscosity (default:1e-6)
    real(kind=JPRB)                 :: vonKar           ! von Karman coefficient (default: 0.4)
    integer(kind=JPIM)              :: psedDT           ! number of timestep within river timestep (DT/sedDT)
    logical                         :: lsedflw          ! if calculate sediment 
    integer(kind=JPIM)              :: ifrq_rst_sed
    character(len=256)              :: sedrest_infile
    character(len=256)              :: sedrest_outpre
    CHARACTER(LEN=256)              :: csedsout = ""      !! Output variable in netCDF
    character(len=256)              :: CSEDCDF="./glb_15min_sediment_00/sediment_15min/Precipitation_0p25D_2001.nc"
    character(len=256)              :: CPREPCDF="./glb_15min_sediment_00/sediment_15min/Precipitation_0p25D_2001.nc"
    character(len=256)              :: CVNSED="precipitation"
    character(len=256)              :: cslope="./glb_15min_sediment_00/sediment_15min/slope.nc4"
    character(len=256)              :: cslopevar="slope"
    character(len=256)              :: csedfrc="./glb_15min_sediment_00/sediment_15min/sedfrc.nc4"
    character(len=256)              :: csedfrcvar="sedfrc"
    character(len=256)              :: sedinput_dir="./"
    character(len=256)              :: sedinput_pre="sed"
    character(len=256)              :: sedinput_suf=".bin"
    logical                         :: LSEDCDF = .TRUE.   ! Flag to use NetCDF for sediment forcing
    integer(kind=JPIM)              :: DTIN = 86400       ! Timestep of sediment input (seconds)
    integer(kind=JPIM)              :: SYEARIN = 0        ! Start year for sediment forcing
    integer(kind=JPIM)              :: SMONIN = 0         ! Start month for sediment forcing
    integer(kind=JPIM)              :: SDAYIN = 0         ! Start day for sediment forcing
    integer(kind=JPIM)              :: SHOURIN = 0        ! Start hour for sediment forcing

  END TYPE CaMaSed

  TYPE :: CaMaOutput
      ! Output directory settings
      CHARACTER(LEN=256) :: COUTDIR = "./"     !! Output directory
      CHARACTER(LEN=256) :: CVARSOUT = ""      !! Output variable in netCDF
      CHARACTER(LEN=256) :: COUTTAG = ""       !! Output Tag Name for each experiment

      CHARACTER(LEN=256) :: CVNMOUT = ""       !! Output variable name
      CHARACTER(LEN=256) :: CSTATOUT = ""      !! Statistical variable name
      INTEGER(KIND=JPIM) :: NDLEVEL = 1        !! NETCDF DEFLATION LEVEL 

      LOGICAL :: LOUTVEC = .FALSE.    !! true: vectrized binary output
      LOGICAL :: LOUTCDF = .FALSE.    !! true: netCDF output

      ! Output format flags
      LOGICAL :: LOUTPUT = .TRUE.     !! true: use standard output
      LOGICAL :: LOUTINI = .FALSE.    !! true: output initial storage

      LOGICAL :: LMAPEND = .FALSE.    !! true: map output at last time
      LOGICAL :: LSTOONLY = .FALSE.   !! true: storage only
      LOGICAL :: LOUTPTH = .FALSE.    !! true: output bifurcation channel flow

      ! Time step settings
      INTEGER(KIND=JPIM) :: IFRQ_OUT = 24     !! output frequency (hour)
      LOGICAL :: LOUTTXT = .FALSE.             !! TRUE FOR Text output for some gauges
      CHARACTER(LEN=256) :: CGAUTXT = "NONE"   !! List of Gauges (ID, IX, IY)

  END TYPE CaMaOutput
  ! Module variables
  TYPE(CaMaFiles)    :: CMF_FILES
  TYPE(CaMaOptions)  :: CMF_OPTIONS
  TYPE(CaMaConfig)   :: CMF_CONFIG
  TYPE(CaMaParams)   :: CMF_PARAMS
  TYPE(CaMaTime)     :: CMF_TIME
  TYPE(CaMaMaps)     :: CMF_MAPS
  TYPE(CaMaForcing)  :: CMF_FORCING
  TYPE(CaMaBoundary) :: CMF_BOUNDARY
  TYPE(CaMaRestart)  :: CMF_RESTART
  TYPE(CaMaDam)      :: CMF_DAM
  TYPE(CaMaLevee)    :: CMF_LEVEE
  TYPE(CaMaTracer)   :: CMF_TRACER
  TYPE(CaMaSed)      :: CMF_SED
  TYPE(CaMaOutput)   :: CMF_OUTPUT
CONTAINS

!####################################################################
SUBROUTINE CMF_CONFIG_NMLIST
    ! Temporary variables for NRUNVER namelist
    LOGICAL :: LADPSTP, LFPLAIN, LKINE, LFLDOUT, LPTHOUT, LDAMOUT, LLEVEE
    LOGICAL :: LROSPLIT, LWEVAP, LWEVAPFIX, LWINFILT, LWINFILTFIX, LWEXTRACTRIV
    LOGICAL :: LSLOPEMOUTH, LGDWDLY, LSLPMIX, LMEANSL, LSEALEV,LOUTINS, LRESTART
    LOGICAL :: LSTOONLY, LOUTPUT, LOUTINI, LGRIDMAP, LLEAPYR, LMAPEND,LBITSAFE, LSTG_ES,LSEDOUT,LTRACE

    ! Temporary variables for NDIMTIME namelist
    CHARACTER(LEN=256) :: CDIMINFO
    REAL(KIND=JPRB) :: DT
    INTEGER(KIND=JPIM) :: IFRQ_INP

    ! Temporary variables for NPARAM namelist
    REAL(KIND=JPRB) :: PMANRIV, PMANFLD, PGRV, PDSTMTH, PCADP, PMINSLP
    INTEGER(KIND=JPIM) :: IMIS
    REAL(KIND=JPRM) :: RMIS
    REAL(KIND=JPRB) :: DMIS
    CHARACTER(LEN=256) :: CSUFBIN, CSUFVEC, CSUFPTH, CSUFCDF

    ! === Previous temporary variables remain unchanged ===

    ! Define the namelist variables
    NAMELIST /NRUNVER/    LADPSTP, LFPLAIN, LKINE, LFLDOUT, LPTHOUT, LDAMOUT, LLEVEE, &
                          LROSPLIT, LWEVAP, LWEVAPFIX, LWINFILT, LWINFILTFIX, LWEXTRACTRIV, &
                          LSLOPEMOUTH, LGDWDLY, LSLPMIX, LMEANSL, LSEALEV, LOUTINS,LRESTART, &
                          LSTOONLY, LOUTPUT, LOUTINI, LGRIDMAP, LLEAPYR, LMAPEND,LBITSAFE, LSTG_ES,LSEDOUT,LTRACE

    NAMELIST /NDIMTIME/  CDIMINFO, DT, IFRQ_INP

    NAMELIST /NPARAM/    PMANRIV, PMANFLD, PGRV, PDSTMTH, PCADP, PMINSLP, &
                        IMIS, RMIS, DMIS, CSUFBIN, CSUFVEC, CSUFPTH, CSUFCDF

    ! *** 0. SET INPUT UNIT AND OPEN FILE 
    CMF_FILES%NSETFILE = INQUIRE_FID()               !!  for namelist
    OPEN(CMF_FILES%NSETFILE, FILE=CMF_FILES%CSETFILE, STATUS="OLD", ACTION="READ")
    WRITE(CMF_FILES%LOGNAM,*) "CMF::CONFIG_NMLIST: namelist opened: ", TRIM(CMF_FILES%CSETFILE), CMF_FILES%NSETFILE 

    !*** Read NRUNVER namelist and update CMF_OPTIONS
    REWIND(CMF_FILES%NSETFILE)
    READ(CMF_FILES%NSETFILE,NML=NRUNVER)
    CMF_OPTIONS%LADPSTP = LADPSTP
    CMF_OPTIONS%LFPLAIN = LFPLAIN
    CMF_OPTIONS%LKINE = LKINE
    CMF_OPTIONS%LFLDOUT = LFLDOUT
    CMF_OPTIONS%LPTHOUT = LPTHOUT
    CMF_OPTIONS%LDAMOUT = LDAMOUT
    CMF_OPTIONS%LLEVEE = LLEVEE
    CMF_OPTIONS%LROSPLIT = LROSPLIT
    CMF_OPTIONS%LWEVAP = LWEVAP
    CMF_OPTIONS%LWEVAPFIX = LWEVAPFIX
    CMF_OPTIONS%LWINFILT = LWINFILT
    CMF_OPTIONS%LWINFILTFIX = LWINFILTFIX
    CMF_OPTIONS%LWEXTRACTRIV = LWEXTRACTRIV
    CMF_OPTIONS%LSLOPEMOUTH = LSLOPEMOUTH
    CMF_OPTIONS%LGDWDLY = LGDWDLY
    CMF_OPTIONS%LSLPMIX = LSLPMIX
    CMF_OPTIONS%LMEANSL = LMEANSL
    CMF_OPTIONS%LSEALEV = LSEALEV
    CMF_OPTIONS%LOUTINS = LOUTINS
    CMF_OPTIONS%LRESTART = LRESTART
    CMF_OPTIONS%LSTOONLY = LSTOONLY
    CMF_OPTIONS%LOUTPUT = LOUTPUT
    CMF_OPTIONS%LOUTINI = LOUTINI
    CMF_OPTIONS%LGRIDMAP = LGRIDMAP
    CMF_OPTIONS%LLEAPYR = LLEAPYR
    CMF_OPTIONS%LMAPEND = LMAPEND
    CMF_OPTIONS%LSTG_ES = LSTG_ES
    CMF_OPTIONS%LSEDOUT = LSEDOUT
    CMF_OPTIONS%LTRACE  = LTRACE




    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "=== NAMELIST, NRUNVER ==="
    WRITE(CMF_FILES%LOGNAM,*) "LADPSTP ",  CMF_OPTIONS%LADPSTP
    WRITE(CMF_FILES%LOGNAM,*) "LFPLAIN ",  CMF_OPTIONS%LFPLAIN
    WRITE(CMF_FILES%LOGNAM,*) "LKINE   ",  CMF_OPTIONS%LKINE
    WRITE(CMF_FILES%LOGNAM,*) "LFLDOUT ",  CMF_OPTIONS%LFLDOUT
    WRITE(CMF_FILES%LOGNAM,*) "LPTHOUT ",  CMF_OPTIONS%LPTHOUT
    WRITE(CMF_FILES%LOGNAM,*) "LDAMOUT ",  CMF_OPTIONS%LDAMOUT
    WRITE(CMF_FILES%LOGNAM,*) "LLEVEE  ",  CMF_OPTIONS%LLEVEE
    WRITE(CMF_FILES%LOGNAM,*) "LSEDOUT ",  CMF_OPTIONS%LSEDOUT
    WRITE(CMF_FILES%LOGNAM,*) "LTRACE  ",  CMF_OPTIONS%LTRACE
    WRITE(CMF_FILES%LOGNAM,*) "LOUTINS ",  CMF_OPTIONS%LOUTINS
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "LROSPLIT ", CMF_OPTIONS%LROSPLIT
    WRITE(CMF_FILES%LOGNAM,*) "LWEVAP   ", CMF_OPTIONS%LWEVAP
    WRITE(CMF_FILES%LOGNAM,*) "LWEVAPFIX", CMF_OPTIONS%LWEVAPFIX
    WRITE(CMF_FILES%LOGNAM,*) "LWEXTRACTRIV", CMF_OPTIONS%LWEXTRACTRIV
    WRITE(CMF_FILES%LOGNAM,*) "LGDWDLY  ", CMF_OPTIONS%LGDWDLY
    WRITE(CMF_FILES%LOGNAM,*) "LSLPMIX  ", CMF_OPTIONS%LSLPMIX
    WRITE(CMF_FILES%LOGNAM,*) "LSLOPEMOUTH ", CMF_OPTIONS%LSLOPEMOUTH
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "LMEANSL: ", CMF_OPTIONS%LMEANSL
    WRITE(CMF_FILES%LOGNAM,*) "LSEALEV: ", CMF_OPTIONS%LSEALEV
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "LRESTART ", CMF_OPTIONS%LRESTART
    WRITE(CMF_FILES%LOGNAM,*) "LSTOONLY ", CMF_OPTIONS%LSTOONLY
    WRITE(CMF_FILES%LOGNAM,*) "LOUTPUT  ", CMF_OPTIONS%LOUTPUT
    WRITE(CMF_FILES%LOGNAM,*) "LOUTINI  ", CMF_OPTIONS%LOUTINI
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "LGRIDMAP ", CMF_OPTIONS%LGRIDMAP
    WRITE(CMF_FILES%LOGNAM,*) "LLEAPYR  ", CMF_OPTIONS%LLEAPYR
    WRITE(CMF_FILES%LOGNAM,*) "LMAPEND  ", CMF_OPTIONS%LMAPEND
    WRITE(CMF_FILES%LOGNAM,*) "LBITSAFE ", CMF_OPTIONS%LBITSAFE
    WRITE(CMF_FILES%LOGNAM,*) "LSTG_ES " , CMF_OPTIONS%LSTG_ES
    WRITE(CMF_FILES%LOGNAM,*) "LSEDOUT " , CMF_OPTIONS%LSEDOUT
    WRITE(CMF_FILES%LOGNAM,*) "LTRACE " , CMF_OPTIONS%LTRACE

    !*** Read NDIMTIME namelist and update CMF_CONFIG
    REWIND(CMF_FILES%NSETFILE)
    READ(CMF_FILES%NSETFILE,NML=NDIMTIME)
    CMF_CONFIG%CDIMINFO = CDIMINFO
    CMF_CONFIG%DT = DT
    CMF_CONFIG%IFRQ_INP = IFRQ_INP
    CMF_CONFIG%DTIN  = IFRQ_INP*60*60       !! hour -> second

    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "=== NAMELIST, NCONF ==="
    WRITE(CMF_FILES%LOGNAM,*) "CDIMINFO  ", TRIM(CMF_CONFIG%CDIMINFO)
    WRITE(CMF_FILES%LOGNAM,*) "DT        ", CMF_CONFIG%DT
    WRITE(CMF_FILES%LOGNAM,*) "DTIN      ", CMF_CONFIG%DTIN
    WRITE(CMF_FILES%LOGNAM,*) "IFRQ_INP  ", CMF_CONFIG%IFRQ_INP


    !* value from CDIMINFO
    IF( CMF_CONFIG%CDIMINFO/="NONE" )THEN
      WRITE(CMF_FILES%LOGNAM,*) "CMF::CONFIG_NMLIST: read DIMINFO ", TRIM(CMF_CONFIG%CDIMINFO)

      CMF_FILES%TMPNAM=INQUIRE_FID()
      OPEN(CMF_FILES%TMPNAM,FILE=CDIMINFO,FORM='FORMATTED')
      READ(CMF_FILES%TMPNAM,*) CMF_CONFIG%NX
      READ(CMF_FILES%TMPNAM,*) CMF_CONFIG%NY
      READ(CMF_FILES%TMPNAM,*) CMF_CONFIG%NLFP
      READ(CMF_FILES%TMPNAM,*) CMF_CONFIG%NXIN
      READ(CMF_FILES%TMPNAM,*) CMF_CONFIG%NYIN
      READ(CMF_FILES%TMPNAM,*) CMF_CONFIG%INPN
      READ(CMF_FILES%TMPNAM,*) 
      IF( CMF_OPTIONS%LGRIDMAP )THEN
        READ(CMF_FILES%TMPNAM,*) CMF_CONFIG%WEST
        READ(CMF_FILES%TMPNAM,*) CMF_CONFIG%EAST
        READ(CMF_FILES%TMPNAM,*) CMF_CONFIG%NORTH
        READ(CMF_FILES%TMPNAM,*) CMF_CONFIG%SOUTH
      ENDIF
      CLOSE(CMF_FILES%TMPNAM)
    ENDIF
    !* check
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "=== DIMINFO ==="
    WRITE(CMF_FILES%LOGNAM,*) "NX,NY,NLFP     ", CMF_CONFIG%NX,  CMF_CONFIG%NY,  CMF_CONFIG%NLFP
    WRITE(CMF_FILES%LOGNAM,*) "NXIN,NYIN,INPN ", CMF_CONFIG%NXIN,CMF_CONFIG%NYIN,CMF_CONFIG%INPN
    IF( CMF_OPTIONS%LGRIDMAP ) THEN
      WRITE(CMF_FILES%LOGNAM,*) "WEST,EAST,NORTH,SOUTH ", CMF_CONFIG%WEST,CMF_CONFIG%EAST,CMF_CONFIG%NORTH,CMF_CONFIG%SOUTH
    ENDIF

    !*** Read NPARAM namelist and update CMF_PARAMS
    REWIND(CMF_FILES%NSETFILE)
    READ(CMF_FILES%NSETFILE,NML=NPARAM)
    CMF_PARAMS%PMANRIV = PMANRIV
    CMF_PARAMS%PMANFLD = PMANFLD
    CMF_PARAMS%PGRV = PGRV
    CMF_PARAMS%PDSTMTH = PDSTMTH
    CMF_PARAMS%PCADP = PCADP
    CMF_PARAMS%PMINSLP = PMINSLP
    CMF_PARAMS%IMIS = IMIS
    CMF_PARAMS%RMIS = RMIS
    CMF_PARAMS%DMIS = DMIS
    CMF_PARAMS%CSUFBIN = CSUFBIN
    CMF_PARAMS%CSUFVEC = CSUFVEC
    CMF_PARAMS%CSUFPTH = CSUFPTH
    CMF_PARAMS%CSUFCDF = CSUFCDF

    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "=== NAMELIST, NPARAM ==="
    WRITE(CMF_FILES%LOGNAM,*) "PMANRIV  ", CMF_PARAMS%PMANRIV
    WRITE(CMF_FILES%LOGNAM,*) "PMANRIV  ", CMF_PARAMS%PMANFLD
    WRITE(CMF_FILES%LOGNAM,*) "PGRV     ", CMF_PARAMS%PGRV
    WRITE(CMF_FILES%LOGNAM,*) "PDSTMTH  ", CMF_PARAMS%PDSTMTH
    WRITE(CMF_FILES%LOGNAM,*) "PCADP    ", CMF_PARAMS%PCADP
    WRITE(CMF_FILES%LOGNAM,*) "PMINSLP  ", CMF_PARAMS%PMINSLP
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "IMIS     ", CMF_PARAMS%IMIS
    WRITE(CMF_FILES%LOGNAM,*) "RMIS     ", CMF_PARAMS%RMIS
    WRITE(CMF_FILES%LOGNAM,*) "DMIS     ", CMF_PARAMS%DMIS
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "CSUFBIN  ", TRIM(CMF_PARAMS%CSUFBIN)
    WRITE(CMF_FILES%LOGNAM,*) "CSUFVEC  ", TRIM(CMF_PARAMS%CSUFVEC)
    WRITE(CMF_FILES%LOGNAM,*) "CSUFPTH  ", TRIM(CMF_PARAMS%CSUFPTH)
    WRITE(CMF_FILES%LOGNAM,*) "CSUFCDF  ", TRIM(CMF_PARAMS%CSUFCDF)


    CLOSE(CMF_FILES%NSETFILE)

END SUBROUTINE CMF_CONFIG_NMLIST

SUBROUTINE CMF_TIME_NMLIST
  ! Declare variables
  INTEGER(KIND=JPIM) :: syear, smon, sday, shour, eyear, emon, eday, ehour

  ! Define the NAMELIST
  NAMELIST /NSIMTIME/ SYEAR, SMON, SDAY, SHOUR, EYEAR, EMON, EDAY, EHOUR

  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
  
  !*** 0. Open namelist
  ! *** 0. SET INPUT UNIT AND OPEN FILE 
  CMF_FILES%NSETFILE = INQUIRE_FID()               !!  for namelist
  OPEN(CMF_FILES%NSETFILE, FILE=CMF_FILES%CSETFILE, STATUS="OLD", ACTION="READ")
  WRITE(CMF_FILES%LOGNAM,*) "CMF::CONFIG_NMLIST: namelist opened: ", TRIM(CMF_FILES%CSETFILE), CMF_FILES%NSETFILE 

  !*** Read NSIMTIME namelist and update CMF_TIME
  REWIND(CMF_FILES%NSETFILE)
  READ(CMF_FILES%NSETFILE,NML=NSIMTIME)

  CMF_TIME%SYEAR = SYEAR
  CMF_TIME%SMON = SMON
  CMF_TIME%SDAY = SDAY
  CMF_TIME%SHOUR = SHOUR
  CMF_TIME%EYEAR = EYEAR
  CMF_TIME%EMON = EMON
  CMF_TIME%EDAY = EDAY
  CMF_TIME%EHOUR = EHOUR
  
  WRITE(CMF_FILES%LOGNAM,*) "=== NAMELIST, NSIMTIME ==="
  WRITE(CMF_FILES%LOGNAM,*) "SYEAR,SMON,SDAY,SHOUR:", CMF_TIME%SYEAR,CMF_TIME%SMON,CMF_TIME%SDAY,CMF_TIME%SHOUR
  WRITE(CMF_FILES%LOGNAM,*) "EYEAR,EMON,EDAY,EHOUR:", CMF_TIME%EYEAR,CMF_TIME%EMON,CMF_TIME%EDAY,CMF_TIME%EHOUR
  
  !*** 3. close namelist
  CLOSE(CMF_FILES%NSETFILE)
  
  !*** 4. Define base date for KMIN calculation
  CMF_TIME%YYYY0=CMF_TIME%SYEAR
  CMF_TIME%MM0=1
  CMF_TIME%DD0=1
  WRITE(CMF_FILES%LOGNAM,*) "TIME_NMLIST: YYYY0 MM0 DD0 set to : ", CMF_TIME%SYEAR, CMF_TIME%SMON, CMF_TIME%SDAY
  WRITE(CMF_FILES%LOGNAM,*) "CMF::TIME_NMLIST: end: "
END SUBROUTINE CMF_TIME_NMLIST

SUBROUTINE CMF_MAPS_NMLIST
  ! Declare variables
  INTEGER :: NSETFILE
  CHARACTER(LEN=256) :: CNEXTXY, CGRAREA, CELEVTN, CNXTDST, CRIVLEN, CFLDHGT
  CHARACTER(LEN=256) :: CRIVWTH, CRIVHGT, CRIVMAN, CPTHOUT, CGDWDLY, CMEANSL
  CHARACTER(LEN=256) :: CMPIREG, CRIVCLINC, CRIVPARNC, CMEANSLNC, CMPIREGNC
  LOGICAL :: LMAPCDF

  NAMELIST/NMAP/ CNEXTXY, CGRAREA, CELEVTN, CNXTDST, CRIVLEN, CFLDHGT, &
                 CRIVWTH, CRIVHGT, CRIVMAN, CPTHOUT, CGDWDLY, CMEANSL, &
                 CMPIREG, LMAPCDF, CRIVCLINC, CRIVPARNC, CMEANSLNC, CMPIREGNC

  !================================================
  !*** 1. open namelist
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"

  NSETFILE = INQUIRE_FID()
  OPEN(NSETFILE, FILE=CMF_FILES%CSETFILE, STATUS="OLD")
  WRITE(CMF_FILES%LOGNAM,*) "CMF::MAP_NMLIST: namelist OPEN in unit: ", TRIM(CMF_FILES%CSETFILE), NSETFILE 

  !*** 2. default value
  CNEXTXY = "./nextxy.bin"
  CGRAREA = "./ctmare.bin"
  CELEVTN = "./elevtn.bin"
  CNXTDST = "./nxtdst.bin"
  CRIVLEN = "./rivlen.bin"
  CFLDHGT = "./fldhgt.bin"

  CRIVWTH = "./rivwth.bin"
  CRIVHGT = "./rivhgt.bin"
  CRIVMAN = "./rivman.bin"

  CPTHOUT = "./bifprm.txt"
  CGDWDLY = "NONE"
  CMEANSL = "NONE"

  CMPIREG = "NONE"

  LMAPCDF = .FALSE.
  CRIVCLINC = "NONE"
  CRIVPARNC = "NONE"
  CMEANSLNC = "NONE"
  CMPIREGNC = "NONE"

  !*** 3. read namelist
  REWIND(NSETFILE)
  READ(NSETFILE, NML=NMAP)
  CMF_MAPS%CNEXTXY = CNEXTXY
  CMF_MAPS%CGRAREA = CGRAREA
  CMF_MAPS%CELEVTN = CELEVTN
  CMF_MAPS%CNXTDST = CNXTDST
  CMF_MAPS%CRIVLEN = CRIVLEN
  CMF_MAPS%CFLDHGT = CFLDHGT

  CMF_MAPS%CRIVWTH = CRIVWTH
  CMF_MAPS%CRIVHGT = CRIVHGT
  CMF_MAPS%CRIVMAN = CRIVMAN
  CMF_MAPS%CPTHOUT = CPTHOUT
  CMF_MAPS%CGDWDLY = CGDWDLY
  CMF_MAPS%CMEANSL = CMEANSL
  CMF_MAPS%CMPIREG = CMPIREG
  CMF_MAPS%LMAPCDF = LMAPCDF
  CMF_MAPS%CRIVCLINC = CRIVCLINC
  CMF_MAPS%CRIVPARNC = CRIVPARNC
  CMF_MAPS%CMEANSLNC = CMEANSLNC
  CMF_MAPS%CMPIREGNC = CMPIREGNC


  WRITE(CMF_FILES%LOGNAM,*) "=== NAMELIST, NMAP ==="
  WRITE(CMF_FILES%LOGNAM,*) "LMAPCDF:   ", CMF_MAPS%LMAPCDF
  IF (CMF_MAPS%LMAPCDF) THEN
    WRITE(CMF_FILES%LOGNAM,*) "CRIVCLINC: ", TRIM(CMF_MAPS%CRIVCLINC)
    WRITE(CMF_FILES%LOGNAM,*) "CRIVPARNC: ", TRIM(CMF_MAPS%CRIVPARNC)
    IF (CMF_OPTIONS%LMEANSL) THEN
      WRITE(CMF_FILES%LOGNAM,*) "CMEANSLNC: ", TRIM(CMF_MAPS%CMEANSLNC)
    ENDIF
#ifdef UseMPI_CMF
    WRITE(CMF_FILES%LOGNAM,*) "CMPIREGNC:   ", TRIM(CMF_MAPS%CMPIREGNC)
#endif
  ELSE
    WRITE(CMF_FILES%LOGNAM,*) "CNEXTXY:   ", TRIM(CMF_MAPS%CNEXTXY)
    WRITE(CMF_FILES%LOGNAM,*) "CGRAREA:   ", TRIM(CMF_MAPS%CGRAREA)
    WRITE(CMF_FILES%LOGNAM,*) "CELEVTN:   ", TRIM(CMF_MAPS%CELEVTN)
    WRITE(CMF_FILES%LOGNAM,*) "CNXTDST:   ", TRIM(CMF_MAPS%CNXTDST)
    WRITE(CMF_FILES%LOGNAM,*) "CRIVLEN:   ", TRIM(CMF_MAPS%CRIVLEN)
    WRITE(CMF_FILES%LOGNAM,*) "CFLDHGT:   ", TRIM(CMF_MAPS%CFLDHGT)

    WRITE(CMF_FILES%LOGNAM,*) "CRIVWTH:   ", TRIM(CMF_MAPS%CRIVWTH)
    WRITE(CMF_FILES%LOGNAM,*) "CRIVHGT:   ", TRIM(CMF_MAPS%CRIVHGT)
    WRITE(CMF_FILES%LOGNAM,*) "CRIVMAN:   ", TRIM(CMF_MAPS%CRIVMAN)

    WRITE(CMF_FILES%LOGNAM,*) "CPTHOUT:   ", TRIM(CMF_MAPS%CPTHOUT)
    IF (CMF_OPTIONS%LGDWDLY) THEN
      WRITE(CMF_FILES%LOGNAM,*) "CGDWDLY:    ", TRIM(CMF_MAPS%CGDWDLY)
    ENDIF
    IF (CMF_OPTIONS%LMEANSL) THEN
      WRITE(CMF_FILES%LOGNAM,*) "CMEANSL:   ", TRIM(CMF_MAPS%CMEANSL)
    ENDIF
#ifdef UseMPI_CMF
    WRITE(CMF_FILES%LOGNAM,*) "CMPIREG:   ", TRIM(CMF_MAPS%CMPIREG)
#endif
  ENDIF

  CLOSE(NSETFILE)

  WRITE(CMF_FILES%LOGNAM,*) "CMF::MAP_NMLIST: end"
END SUBROUTINE CMF_MAPS_NMLIST

SUBROUTINE CMF_FORCING_NMLIST
    ! Temporary variables for NFORCE namelist
    LOGICAL :: LINPCDF, LINPEND, LINTERP, LITRPCDF
    CHARACTER(LEN=256) :: CINPMAT, CROFDIR, CROFPRE, CROFSUF
    CHARACTER(LEN=256) :: CSUBDIR, CSUBPRE, CSUBSUF, CROFCDF
    CHARACTER(LEN=256) :: CVNROF, CVNSUB
    REAL(KIND=JPRB) :: DROFUNIT
    INTEGER(KIND=JPIM) :: SYEARIN, SMONIN, SDAYIN, SHOURIN

    ! Define the NAMELIST
    NAMELIST /NFORCE/ LINTERP, LINPEND, LINPCDF, LITRPCDF, CINPMAT, DROFUNIT, &
                     CROFDIR, CROFPRE, CROFSUF, CSUBDIR, CSUBPRE, CSUBSUF, &
                     CROFCDF, CVNROF, CVNSUB, SYEARIN, SMONIN, SDAYIN, SHOURIN

    !*** 0. Open namelist
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"

    CMF_FILES%NSETFILE = INQUIRE_FID()
    OPEN(CMF_FILES%NSETFILE, FILE=CMF_FILES%CSETFILE, STATUS="OLD", ACTION="READ")
    WRITE(CMF_FILES%LOGNAM,*) "CMF::FORCING_NMLIST: namelist opened: ", TRIM(CMF_FILES%CSETFILE), CMF_FILES%NSETFILE

    !*** Set default values
    LINPCDF = .FALSE.
    LINPEND = .FALSE.
    LINTERP = .FALSE.
    LITRPCDF = .FALSE.
    CINPMAT = "NONE"
    DROFUNIT = 86400._JPRB * 1000._JPRB    ! defaults mm/day -> m3/m2/s

    CROFDIR = "./runoff/"
    CROFPRE = "Roff____"
    CROFSUF = ".one"
    CSUBDIR = "./runoff/"
    CSUBPRE = "Rsub____"
    CSUBSUF = ".one"
    CROFCDF = "NONE"
    CVNROF = "runoff"
    CVNSUB = "NONE"

    SYEARIN = 0
    SMONIN = 0
    SDAYIN = 0
    SHOURIN = 0

    !*** Read namelist
    REWIND(CMF_FILES%NSETFILE)
    READ(CMF_FILES%NSETFILE,NML=NFORCE)

    !*** Update CMF_FORCING settings
    CMF_FORCING%LINPCDF = LINPCDF
    CMF_FORCING%LINPEND = LINPEND
    CMF_FORCING%LINTERP = LINTERP
    CMF_FORCING%LITRPCDF = LITRPCDF
    CMF_FORCING%CINPMAT = CINPMAT
    CMF_FORCING%DROFUNIT = DROFUNIT
    CMF_FORCING%CROFDIR = CROFDIR
    CMF_FORCING%CROFPRE = CROFPRE
    CMF_FORCING%CROFSUF = CROFSUF
    CMF_FORCING%CSUBDIR = CSUBDIR
    CMF_FORCING%CSUBPRE = CSUBPRE
    CMF_FORCING%CSUBSUF = CSUBSUF
    CMF_FORCING%CROFCDF = CROFCDF
    CMF_FORCING%CVNROF = CVNROF
    CMF_FORCING%CVNSUB = CVNSUB
    CMF_FORCING%SYEARIN = SYEARIN
    CMF_FORCING%SMONIN = SMONIN
    CMF_FORCING%SDAYIN = SDAYIN
    CMF_FORCING%SHOURIN = SHOURIN
    !*** Write configuration to log
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "=== NAMELIST, NFORCE ==="
    WRITE(CMF_FILES%LOGNAM,*) "LINPCDF:   ", CMF_FORCING%LINPCDF
    WRITE(CMF_FILES%LOGNAM,*) "LINTERP:   ", CMF_FORCING%LINTERP
    WRITE(CMF_FILES%LOGNAM,*) "LITRPCDF:  ", CMF_FORCING%LITRPCDF
    WRITE(CMF_FILES%LOGNAM,*) "CINPMAT:   ", TRIM(CMF_FORCING%CINPMAT)
    WRITE(CMF_FILES%LOGNAM,*) "DROFUNIT:  ", CMF_FORCING%DROFUNIT

    IF (.NOT. CMF_FORCING%LINPCDF) THEN
        WRITE(CMF_FILES%LOGNAM,*) "CROFDIR:   ", TRIM(CMF_FORCING%CROFDIR)
        WRITE(CMF_FILES%LOGNAM,*) "CROFPRE:   ", TRIM(CMF_FORCING%CROFPRE)
        WRITE(CMF_FILES%LOGNAM,*) "CROFSUF:   ", TRIM(CMF_FORCING%CROFSUF)
    ELSE
        WRITE(CMF_FILES%LOGNAM,*) "CROFCDF:   ", TRIM(CMF_FORCING%CROFCDF)
        WRITE(CMF_FILES%LOGNAM,*) "CVNROF:    ", TRIM(CMF_FORCING%CVNROF)
        WRITE(CMF_FILES%LOGNAM,*) "CVNSUB:    ", TRIM(CMF_FORCING%CVNSUB)
        WRITE(CMF_FILES%LOGNAM,*) "SYEARIN,SMONIN,SDAYIN,SHOURIN: ", &
            CMF_FORCING%SYEARIN, CMF_FORCING%SMONIN, CMF_FORCING%SDAYIN, CMF_FORCING%SHOURIN
    END IF

    CLOSE(CMF_FILES%NSETFILE)

    WRITE(CMF_FILES%LOGNAM,*) "CMF::FORCING_NMLIST: end"

END SUBROUTINE CMF_FORCING_NMLIST

SUBROUTINE CMF_BOUNDARY_NMLIST
  ! Temporary variables for NBOUND namelist
  LOGICAL :: LSEALEVCDF
  CHARACTER(LEN=256) :: CSEALEVDIR, CSEALEVPRE, CSEALEVSUF
  CHARACTER(LEN=256) :: CSEALEVCDF, CVNSEALEV, CSLMAP
  INTEGER(KIND=JPIM) :: SYEARSL, SMONSL, SDAYSL, SHOURSL
  INTEGER(KIND=JPIM) :: NLINKS, NCDFSTAT, IFRQ_SL

  ! Define the NAMELIST
  NAMELIST /NBOUND/ LSEALEVCDF, CSEALEVDIR, CSEALEVPRE, CSEALEVSUF, &
                    CSEALEVCDF, CVNSEALEV, SYEARSL, SMONSL, SDAYSL, SHOURSL, &
                    CSLMAP, NLINKS, NCDFSTAT, IFRQ_SL

  !*** 0. Open namelist
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"

  CMF_FILES%NSETFILE = INQUIRE_FID()
  OPEN(CMF_FILES%NSETFILE, FILE=CMF_FILES%CSETFILE, STATUS="OLD", ACTION="READ")
  WRITE(CMF_FILES%LOGNAM,*) "CMF::BOUNDARY_NMLIST: namelist opened: ", TRIM(CMF_FILES%CSETFILE), CMF_FILES%NSETFILE

  !*** Set default values
  LSEALEVCDF = .FALSE.
  CSEALEVDIR = "./sealev/"
  CSEALEVPRE = "sealev"
  CSEALEVSUF = ".bin"
  CSEALEVCDF = "./sealev/"
  CVNSEALEV = "variable"
  CSLMAP = "./sealev/"
  SYEARSL = 0
  SMONSL = 0
  SDAYSL = 0
  SHOURSL = 0
  IFRQ_SL = 9999
  !*** Read namelist
  REWIND(CMF_FILES%NSETFILE)
  READ(CMF_FILES%NSETFILE,NML=NBOUND)

  !*** Update CMF_BOUNDARY settings
  CMF_BOUNDARY%LSEALEVCDF = LSEALEVCDF
  CMF_BOUNDARY%CSEALEVDIR = CSEALEVDIR
  CMF_BOUNDARY%CSEALEVPRE = CSEALEVPRE
  CMF_BOUNDARY%CSEALEVSUF = CSEALEVSUF
  CMF_BOUNDARY%CSEALEVCDF = CSEALEVCDF
  CMF_BOUNDARY%CVNSEALEV = CVNSEALEV
  CMF_BOUNDARY%CSLMAP = CSLMAP
  CMF_BOUNDARY%SYEARSL = SYEARSL
  CMF_BOUNDARY%SMONSL = SMONSL
  CMF_BOUNDARY%SDAYSL = SDAYSL
  CMF_BOUNDARY%SHOURSL = SHOURSL
  CMF_BOUNDARY%NLINKS = NLINKS
  CMF_BOUNDARY%NCDFSTAT = NCDFSTAT
  CMF_BOUNDARY%IFRQ_SL = IFRQ_SL
  CMF_BOUNDARY%DTSL  = CMF_BOUNDARY%IFRQ_SL *60          !! min  -> second

  !*** Write configuration to log
  IF (CMF_OPTIONS%LSEALEV) THEN
      WRITE(CMF_FILES%LOGNAM,*) ""
      WRITE(CMF_FILES%LOGNAM,*) "=== NAMELIST, NBOUND ==="
      WRITE(CMF_FILES%LOGNAM,*) "LSEALEVCDF: ", CMF_BOUNDARY%LSEALEVCDF
      IF (CMF_BOUNDARY%LSEALEVCDF) THEN
          WRITE(CMF_FILES%LOGNAM,*) "CSEALEVCDF: ", TRIM(CMF_BOUNDARY%CSEALEVCDF)
          WRITE(CMF_FILES%LOGNAM,*) "CVNSEALEV:  ", TRIM(CMF_BOUNDARY%CVNSEALEV)
          WRITE(CMF_FILES%LOGNAM,*) "SYEARSL:    ", CMF_BOUNDARY%SYEARSL
          WRITE(CMF_FILES%LOGNAM,*) "SMONSL:     ", CMF_BOUNDARY%SMONSL
          WRITE(CMF_FILES%LOGNAM,*) "SDAYSL:     ", CMF_BOUNDARY%SDAYSL
          WRITE(CMF_FILES%LOGNAM,*) "SHOURSL:    ", CMF_BOUNDARY%SHOURSL
          WRITE(CMF_FILES%LOGNAM,*) "CSLMAP:     ", TRIM(CMF_BOUNDARY%CSLMAP)
      ELSE
          WRITE(CMF_FILES%LOGNAM,*) "CSEALEVDIR: ", TRIM(CMF_BOUNDARY%CSEALEVDIR)
          WRITE(CMF_FILES%LOGNAM,*) "CSEALEVPRE: ", TRIM(CMF_BOUNDARY%CSEALEVPRE)
          WRITE(CMF_FILES%LOGNAM,*) "CSEALEVSUF: ", TRIM(CMF_BOUNDARY%CSEALEVSUF)
      END IF
  END IF

  CLOSE(CMF_FILES%NSETFILE)

  WRITE(CMF_FILES%LOGNAM,*) "CMF::BOUNDARY_NMLIST: end"

END SUBROUTINE CMF_BOUNDARY_NMLIST

SUBROUTINE CMF_RESTART_NMLIST
    ! Temporary variables for NRESTART namelist
    CHARACTER(LEN=256) :: CRESTSTO   ! input restart file name
    CHARACTER(LEN=256) :: CRESTDIR   ! output restart file directory
    CHARACTER(LEN=256) :: CVNREST    ! output restart prefix
    LOGICAL :: LRESTCDF             ! true: netCDF restart file
    LOGICAL :: LRESTDBL             ! true: binary restart double precision
    INTEGER(KIND=JPIM) :: IFRQ_RST  ! restart frequency
    
    NAMELIST /NRESTART/ CRESTSTO, CRESTDIR, CVNREST, LRESTCDF, LRESTDBL, IFRQ_RST

    !*** 0. Open namelist
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"

    CMF_FILES%NSETFILE = INQUIRE_FID()
    OPEN(CMF_FILES%NSETFILE, FILE=CMF_FILES%CSETFILE, STATUS="OLD", ACTION="READ")
    WRITE(CMF_FILES%LOGNAM,*) "CMF::RESTART_NMLIST: namelist opened: ", TRIM(CMF_FILES%CSETFILE), CMF_FILES%NSETFILE 

    !*** Set default values
    CRESTSTO = "restart"   ! input restart file name
    CRESTDIR = "./"       ! output restart file directory
    CVNREST = "restart"   ! output restart file prefix
    LRESTCDF = .FALSE.    ! true: netCDF restart file
    LRESTDBL = .TRUE.     ! true: binary restart double precision
    IFRQ_RST = 0         ! 0: only end of simulation
                                                     !! [1,2,3,6,12,24]: at selected hour
                                                     !! 30: monthly

    !*** Read namelist
    REWIND(CMF_FILES%NSETFILE)
    READ(CMF_FILES%NSETFILE, NML=NRESTART)

    !*** Update CMF_RESTART settings
    CMF_RESTART%CRESTSTO = CRESTSTO
    CMF_RESTART%CRESTDIR = CRESTDIR
    CMF_RESTART%CVNREST = CVNREST
    CMF_RESTART%LRESTCDF = LRESTCDF
    CMF_RESTART%LRESTDBL = LRESTDBL
    CMF_RESTART%IFRQ_RST = IFRQ_RST

    !*** Write configuration to log
    WRITE(CMF_FILES%LOGNAM,*) "=== NAMELIST, NRESTART ==="
    WRITE(CMF_FILES%LOGNAM,*) "CRESTSTO:  ", TRIM(CMF_RESTART%CRESTSTO)
    WRITE(CMF_FILES%LOGNAM,*) "CRESTDIR:  ", TRIM(CMF_RESTART%CRESTDIR)
    WRITE(CMF_FILES%LOGNAM,*) "CVNREST:   ", TRIM(CMF_RESTART%CVNREST)
    WRITE(CMF_FILES%LOGNAM,*) "LRESTCDF:  ", CMF_RESTART%LRESTCDF
    WRITE(CMF_FILES%LOGNAM,*) "LRESTDBL:  ", CMF_RESTART%LRESTDBL
    WRITE(CMF_FILES%LOGNAM,*) "IFRQ_RST:  ", CMF_RESTART%IFRQ_RST

    CLOSE(CMF_FILES%NSETFILE)

END SUBROUTINE CMF_RESTART_NMLIST

SUBROUTINE CMF_DAMOUT_NMLIST
  ! Temporary variables for NDAMOUT namelist
  CHARACTER(LEN=256) :: CDAMFILE    ! dam parameter file
  LOGICAL :: LDAMTXT                ! true: dam inflow-outflow txt output
  LOGICAL :: LDAMH22                ! true: Use Hanazaki 2022 scheme
  LOGICAL :: LDAMYBY               ! true: Use Year-By-Year dam activation
  LOGICAL :: LiVnorm               ! true: initialize dam storage with Normal Volume

  ! Define the NAMELIST
  NAMELIST /NDAMOUT/ CDAMFILE, LDAMTXT, LDAMH22, LDAMYBY, LiVnorm

  !*** 0. Open namelist
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"

  CMF_FILES%NSETFILE = INQUIRE_FID()
  OPEN(CMF_FILES%NSETFILE, FILE=CMF_FILES%CSETFILE, STATUS="OLD", ACTION="READ")
  WRITE(CMF_FILES%LOGNAM,*) "CMF::DAMOUT_NMLIST: namelist opened: ", TRIM(CMF_FILES%CSETFILE), CMF_FILES%NSETFILE

  !*** Set default values
  CDAMFILE = "./dam_params.csv"
  LDAMTXT = .TRUE.
  LDAMH22 = .FALSE.
  LDAMYBY = .FALSE.
  LiVnorm = .FALSE.

  !*** Read namelist
  REWIND(CMF_FILES%NSETFILE)
  READ(CMF_FILES%NSETFILE,NML=NDAMOUT)

  !*** Update CMF_DAM settings
  CMF_DAM%CDAMFILE = CDAMFILE
  CMF_DAM%LDAMTXT = LDAMTXT
  CMF_DAM%LDAMH22 = LDAMH22
  CMF_DAM%LDAMYBY = LDAMYBY
  CMF_DAM%LiVnorm = LiVnorm

  !*** Write configuration to log
  IF (CMF_OPTIONS%LDAMOUT) THEN
    WRITE(CMF_FILES%LOGNAM,*) "=== NAMELIST, NDAMOUT ==="
    WRITE(CMF_FILES%LOGNAM,*) "CDAMFILE: ", TRIM(CMF_DAM%CDAMFILE)
    WRITE(CMF_FILES%LOGNAM,*) "LDAMTXT:  ", CMF_DAM%LDAMTXT
    WRITE(CMF_FILES%LOGNAM,*) "LDAMH22:  ", CMF_DAM%LDAMH22
    WRITE(CMF_FILES%LOGNAM,*) "LDAMYBY:  ", CMF_DAM%LDAMYBY
    WRITE(CMF_FILES%LOGNAM,*) "LiVnorm:  ", CMF_DAM%LiVnorm
  END IF

  CLOSE(CMF_FILES%NSETFILE)

  WRITE(CMF_FILES%LOGNAM,*) "CMF::DAMOUT_NMLIST: end"

END SUBROUTINE CMF_DAMOUT_NMLIST

SUBROUTINE CMF_LEVEE_NMLIST
    ! Temporary variables for NLEVEE namelist
    CHARACTER(LEN=256) :: CLEVHGT    ! LEVEE HEIGHT from RIVER
    CHARACTER(LEN=256) :: CLEVFRC    ! Unprotected fraction. Relative Levee distance from RIVER

    ! Define the NAMELIST
    NAMELIST /NLEVEE/ CLEVHGT, CLEVFRC

    !*** 0. Open namelist
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"

    CMF_FILES%NSETFILE = INQUIRE_FID()
    OPEN(CMF_FILES%NSETFILE, FILE=CMF_FILES%CSETFILE, STATUS="OLD", ACTION="READ")
    WRITE(CMF_FILES%LOGNAM,*) "CMF::LEVEE_NMLIST: namelist opened: ", TRIM(CMF_FILES%CSETFILE), CMF_FILES%NSETFILE

    !*** Set default values
    CLEVHGT = "NONE"
    CLEVFRC = "NONE"

    !*** Read namelist
    REWIND(CMF_FILES%NSETFILE)
    READ(CMF_FILES%NSETFILE,NML=NLEVEE)

    !*** Update CMF_LEVEE settings
    CMF_LEVEE%CLEVHGT = CLEVHGT
    CMF_LEVEE%CLEVFRC = CLEVFRC

    !*** Write configuration to log
    IF (CMF_OPTIONS%LLEVEE) THEN
        WRITE(CMF_FILES%LOGNAM,*) "=== NAMELIST, NLEVEE ==="
        WRITE(CMF_FILES%LOGNAM,*) "CLEVHGT: ", TRIM(CMF_LEVEE%CLEVHGT)
        WRITE(CMF_FILES%LOGNAM,*) "CLEVFRC: ", TRIM(CMF_LEVEE%CLEVFRC)
    END IF

    CLOSE(CMF_FILES%NSETFILE)

    WRITE(CMF_FILES%LOGNAM,*) "CMF::LEVEE_NMLIST: end"

END SUBROUTINE CMF_LEVEE_NMLIST

SUBROUTINE CMF_TRACER_NMLIST
  ! Temporary variables for NTRACER namelist
  CHARACTER(LEN=256) :: CTRCNAM     ! tracer name
  CHARACTER(LEN=256) :: CTRCDIR     ! tracer input file directory
  CHARACTER(LEN=256) :: CTRCPRE     ! tracer file prefix
  CHARACTER(LEN=256) :: CTRCSUF     ! tracer file suffix
  INTEGER(KIND=JPIM) :: NTRACE      ! number of tracer
  INTEGER(KIND=JPIM) :: IFRQ_TRIN   ! tracer input frequency (hour)
  REAL(KIND=JPRM)    :: DTRCUNIT    !! tracer input unit conversion (DTRCUNIT=1 when tracer input file unit is [(MASS)/m2/s]. )
  LOGICAL            :: LINPEND     !! true for input endian conversion
  LOGICAL            :: LTRCBIF     !! true for consider bifurcation in tracer scheme
  CHARACTER(LEN=256) :: CRESTTRC    ! input restart file name
  CHARACTER(LEN=256) :: CRESTDIR    ! output restart file directory
  CHARACTER(LEN=256) :: CVNRSTTRC   ! output restart prefix
  LOGICAL            :: LRESTDBL    !! true: binary restart in double precision
  INTEGER(KIND=JPIM) :: IFRQ_RST    ! 0: only at last time, (1,2,3,...,24) hourly restart, 30: monthly restart
  CHARACTER(LEN=256) :: COUTDIR     ! output directory
  CHARACTER(LEN=256) :: COUTTAG     ! output tag name for each experiment
  LOGICAL            :: LOUTVEC     !! true for vectorial output, false for NX,NY output

  NAMELIST /NTRACER/ CTRCNAM, CTRCDIR, CTRCPRE, CTRCSUF, NTRACE, IFRQ_TRIN, DTRCUNIT, LINPEND, LTRCBIF, CRESTTRC, CRESTDIR, CVNRSTTRC, LRESTDBL, IFRQ_RST, COUTDIR, COUTTAG, LOUTVEC

  !================================================
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"

  !*** 1. open namelist
  CMF_FILES%NSETFILE = INQUIRE_FID()
  OPEN(CMF_FILES%NSETFILE, FILE=CMF_FILES%CSETFILE, STATUS="OLD", ACTION="READ")
  WRITE(CMF_FILES%LOGNAM,*) "CMF::TRACER_NMLIST: namelist opened: ", TRIM(CMF_FILES%CSETFILE), CMF_FILES%NSETFILE

  !*** 2. default values
  CTRCNAM = "NONE"
  CTRCDIR = "./"
  CTRCPRE = "tracer"
  CTRCSUF = ".dat"
  NTRACE = 1
  IFRQ_TRIN = 24
  DTRCUNIT = 1.0_JPRM
  LINPEND = .FALSE.
  LTRCBIF = .FALSE.
  CRESTTRC = "resttrc"
  CRESTDIR = "./"
  CVNRSTTRC = "resttrc"
  LRESTDBL = .TRUE.
  IFRQ_RST = 0
  COUTDIR = "./"
  COUTTAG = "_cmf"
  LOUTVEC = .FALSE.

  !*** 3. read namelist
  REWIND(CMF_FILES%NSETFILE)
  READ(CMF_FILES%NSETFILE, NML=NTRACER)

  !*** 4. Update CMF_TRACER settings
  CMF_TRACER%CTRCNAM = CTRCNAM
  CMF_TRACER%CTRCDIR = CTRCDIR
  CMF_TRACER%CTRCPRE = CTRCPRE
  CMF_TRACER%CTRCSUF = CTRCSUF
  CMF_TRACER%NTRACE = NTRACE
  CMF_TRACER%IFRQ_TRIN = IFRQ_TRIN
  CMF_TRACER%DTRCUNIT = DTRCUNIT
  CMF_TRACER%LINPEND = LINPEND
  CMF_TRACER%LTRCBIF = LTRCBIF
  CMF_TRACER%CRESTTRC = CRESTTRC
  CMF_TRACER%CRESTDIR = CRESTDIR
  CMF_TRACER%CVNRSTTRC = CVNRSTTRC
  CMF_TRACER%LRESTDBL = LRESTDBL
  CMF_TRACER%IFRQ_RST = IFRQ_RST
  CMF_TRACER%COUTDIR = COUTDIR
  CMF_TRACER%COUTTAG = COUTTAG
  CMF_TRACER%LOUTVEC = LOUTVEC
  !*** 5. Write configuration to log
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "=== NAMELIST, NTRACER ==="
  WRITE(CMF_FILES%LOGNAM,*) "CTRCNAM  ", CMF_TRACER%CTRCNAM
  WRITE(CMF_FILES%LOGNAM,*) "CTRCDIR  ", CMF_TRACER%CTRCDIR
  WRITE(CMF_FILES%LOGNAM,*) "CTRCPRE  ", CMF_TRACER%CTRCPRE
  WRITE(CMF_FILES%LOGNAM,*) "CTRCSUF  ", CMF_TRACER%CTRCSUF
  WRITE(CMF_FILES%LOGNAM,*) "NTRACE   ", CMF_TRACER%NTRACE
  WRITE(CMF_FILES%LOGNAM,*) "IFRQ_TRIN", CMF_TRACER%IFRQ_TRIN
  WRITE(CMF_FILES%LOGNAM,*) "DTRCUNIT ", CMF_TRACER%DTRCUNIT
  WRITE(CMF_FILES%LOGNAM,*) "LINPEND  ", CMF_TRACER%LINPEND
  WRITE(CMF_FILES%LOGNAM,*) "LTRCBIF  ", CMF_TRACER%LTRCBIF
  WRITE(CMF_FILES%LOGNAM,*) "CRESTTRC ", CMF_TRACER%CRESTTRC
  WRITE(CMF_FILES%LOGNAM,*) "CRESTDIR ", CMF_TRACER%CRESTDIR
  WRITE(CMF_FILES%LOGNAM,*) "CVNRSTTRC", CMF_TRACER%CVNRSTTRC
  WRITE(CMF_FILES%LOGNAM,*) "LRESTDBL ", CMF_TRACER%LRESTDBL
  WRITE(CMF_FILES%LOGNAM,*) "IFRQ_RST ", CMF_TRACER%IFRQ_RST
  WRITE(CMF_FILES%LOGNAM,*) "COUTDIR  ", CMF_TRACER%COUTDIR
  WRITE(CMF_FILES%LOGNAM,*) "COUTTAG  ", CMF_TRACER%COUTTAG
  WRITE(CMF_FILES%LOGNAM,*) "LOUTVEC  ", CMF_TRACER%LOUTVEC

  CMF_TRACER%DTIN_TRC= CMF_TRACER%IFRQ_TRIN*60.*60.  !! hour to sec
  WRITE(CMF_FILES%LOGNAM,*)   "DTIN_TRC: " , CMF_TRACER%DTIN_TRC
  CLOSE(CMF_FILES%NSETFILE)

  WRITE(CMF_FILES%LOGNAM,*) "CMF::TRACER_NMLIST: end"


END SUBROUTINE CMF_TRACER_NMLIST

SUBROUTINE CMF_SED_NMLIST
  ! Temporary variables for sediment_param namelist
  REAL(KIND=JPRB) :: lambda, lyrdph, sedDT, psedD, pset, pwatD, visKin, vonKar,DSYLUNIT
  INTEGER(KIND=JPIM) :: nsed, totlyrnum
  LOGICAL :: revEgia, LSEDCDF
  INTEGER(KIND=JPIM) :: DTIN, SYEARIN, SMONIN, SDAYIN, SHOURIN

  ! Temporary variables for sediment_map namelist
  CHARACTER(LEN=256) :: crocdph, csedfrc, sedD
  integer(kind=JPIM)              :: ifrq_rst_sed
  character(len=256)              :: sedrest_infile, sedrest_outpre
  character(len=256)              :: csedsout
  character(len=256)              :: sedinput_dir, sedinput_pre, sedinput_suf,CVNSED, CSEDCDF, cslope,CPREPCDF
  ! Define the namelist variables
  NAMELIST /NSEDIMENT/ crocdph, sedD, csedfrc, lambda, lyrdph, nsed, DSYLUNIT,sedDT, psedD, pset, pwatD, revEgia, totlyrnum, visKin, vonKar, sedrest_infile, sedrest_outpre, ifrq_rst_sed,csedsout,sedinput_dir,sedinput_pre,sedinput_suf,CPREPCDF,CSEDCDF,CVNSED,cslope,LSEDCDF,DTIN,SYEARIN,SMONIN,SDAYIN,SHOURIN

  ! *** 0. SET INPUT UNIT AND OPEN FILE 
  CMF_FILES%NSETFILE = INQUIRE_FID()               !!  for namelist
  OPEN(CMF_FILES%NSETFILE, FILE=CMF_FILES%CSETFILE, STATUS="OLD", ACTION="READ")
  WRITE(CMF_FILES%LOGNAM,*) "CMF::SED_NMLIST: namelist opened: ", TRIM(CMF_FILES%CSETFILE), CMF_FILES%NSETFILE 
  
  lambda = 0.4d0
  lyrdph = 0.00005d0
  nsed = 3
  DSYLUNIT = 1.d-6
  sedDT = 3600
  psedD = 2.65d0
  pset = 1.d0
  pwatD = 1.d0
  revEgia = .true.
  totlyrnum = 5
  visKin = 1.d-6
  vonKar = 0.4d0
  LSEDCDF = .TRUE.
  DTIN = 86400
  SYEARIN = 0
  SMONIN = 0
  SDAYIN = 0
  SHOURIN = 0
  
  !*** Read sediment namelist and update CMF_SED
  REWIND(CMF_FILES%NSETFILE)
  READ(CMF_FILES%NSETFILE, NML=NSEDIMENT)
  CMF_SED%lambda = lambda
  CMF_SED%lyrdph = lyrdph
  CMF_SED%nsed = nsed
  CMF_SED%sedDT = sedDT
  CMF_SED%psedD = psedD
  CMF_SED%DSYLUNIT = DSYLUNIT
  CMF_SED%pset = pset
  CMF_SED%pwatD = pwatD
  CMF_SED%revEgia = revEgia
  CMF_SED%totlyrnum = totlyrnum
  CMF_SED%visKin = visKin
  CMF_SED%vonKar = vonKar
  CMF_SED%crocdph = crocdph
  CMF_SED%sedD = sedD
  CMF_SED%sedrest_infile = sedrest_infile
  CMF_SED%sedrest_outpre = sedrest_outpre
  CMF_SED%ifrq_rst_sed   = ifrq_rst_sed
  CMF_SED%csedsout = csedsout
  CMF_SED%sedinput_dir = sedinput_dir
  CMF_SED%sedinput_pre = sedinput_pre
  CMF_SED%sedinput_suf = sedinput_suf
  CMF_SED%CSEDCDF = CSEDCDF
  CMF_SED%CPREPCDF = CPREPCDF
  CMF_SED%CVNSED = CVNSED
  CMF_SED%cslope = cslope
  CMF_SED%csedfrc = csedfrc
  CMF_SED%LSEDCDF = LSEDCDF
  CMF_SED%DTIN = DTIN
  CMF_SED%SYEARIN = SYEARIN
  CMF_SED%SMONIN = SMONIN
  CMF_SED%SDAYIN = SDAYIN
  CMF_SED%SHOURIN = SHOURIN
  
  !*** Write configuration to log
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "=== NAMELIST, sediment ==="
  WRITE(CMF_FILES%LOGNAM,*) "lambda    ", CMF_SED%lambda
  WRITE(CMF_FILES%LOGNAM,*) "lyrdph    ", CMF_SED%lyrdph
  WRITE(CMF_FILES%LOGNAM,*) "nsed      ", CMF_SED%nsed
  WRITE(CMF_FILES%LOGNAM,*) "sedDT     ", CMF_SED%sedDT
  WRITE(CMF_FILES%LOGNAM,*) "psedD     ", CMF_SED%psedD
  WRITE(CMF_FILES%LOGNAM,*) "pset      ", CMF_SED%pset
  WRITE(CMF_FILES%LOGNAM,*) "pwatD     ", CMF_SED%pwatD
  WRITE(CMF_FILES%LOGNAM,*) "revEgia   ", CMF_SED%revEgia
  WRITE(CMF_FILES%LOGNAM,*) "totlyrnum ", CMF_SED%totlyrnum
  WRITE(CMF_FILES%LOGNAM,*) "visKin    ", CMF_SED%visKin
  WRITE(CMF_FILES%LOGNAM,*) "vonKar    ", CMF_SED%vonKar
  WRITE(CMF_FILES%LOGNAM,*) "LSEDCDF   ", CMF_SED%LSEDCDF
  WRITE(CMF_FILES%LOGNAM,*) "DTIN      ", CMF_SED%DTIN
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "crocdph   ", TRIM(CMF_SED%crocdph)
  WRITE(CMF_FILES%LOGNAM,*) "sedD      ", CMF_SED%sedD
  WRITE(CMF_FILES%LOGNAM,*) "csedfrc   ", CMF_SED%csedfrc
  WRITE(CMF_FILES%LOGNAM,*) "sedrest_infile   ", CMF_SED%sedrest_infile
  WRITE(CMF_FILES%LOGNAM,*) "sedrest_outpre   ", CMF_SED%sedrest_outpre
  WRITE(CMF_FILES%LOGNAM,*) "ifrq_rst_sed   ", CMF_SED%ifrq_rst_sed
  WRITE(CMF_FILES%LOGNAM,*) "csedsout   ", CMF_SED%csedsout
  WRITE(CMF_FILES%LOGNAM,*) "CSEDCDF   ", TRIM(CMF_SED%CSEDCDF)
  WRITE(CMF_FILES%LOGNAM,*) "CVNSED    ", TRIM(CMF_SED%CVNSED)
  WRITE(CMF_FILES%LOGNAM,*) "CPREPCDF  ", TRIM(CMF_SED%CPREPCDF)
  CLOSE(CMF_FILES%NSETFILE)

  WRITE(CMF_FILES%LOGNAM,*) "CMF::SED_NMLIST: end"

END SUBROUTINE CMF_SED_NMLIST


SUBROUTINE CMF_OUTPUT_NMLIST
  ! Temporary variables for NOUTPUT namelist
  CHARACTER(LEN=256) :: COUTDIR           ! OUTPUT DIRECTORY
  CHARACTER(LEN=256) :: CVARSOUT          ! Comma-separated list of output variables to save 
  CHARACTER(LEN=256) :: COUTTAG           ! Output Tag Name for each experiment
  LOGICAL :: LOUTCDF           ! true for netcdf output false for binary
  INTEGER(KIND=JPIM) :: NDLEVEL           ! NETCDF DEFLATION LEVEL 
  LOGICAL :: LOUTVEC           ! TRUE FOR VECTORIAL OUTPUT, FALSE FOR NX,NY OUTPUT
  INTEGER(KIND=JPIM) :: IFRQ_OUT          ! output frequency (hour)
  LOGICAL :: LOUTTXT           ! TRUE FOR Text output for some gauges
  CHARACTER(LEN=256) :: CGAUTXT           ! List of Gauges (ID, IX, IY)

  NAMELIST/NOUTPUT/ COUTDIR,CVARSOUT,COUTTAG,LOUTCDF,NDLEVEL,LOUTVEC,IFRQ_OUT,LOUTTXT,CGAUTXT

  !================================================
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"

  !*** 1. open namelist
  CMF_FILES%NSETFILE = INQUIRE_FID()
  OPEN(CMF_FILES%NSETFILE, FILE=CMF_FILES%CSETFILE, STATUS="OLD", ACTION="READ")
  WRITE(CMF_FILES%LOGNAM,*) "CMF::OUTPUT_NMLIST: namelist opened: ", TRIM(CMF_FILES%CSETFILE), CMF_FILES%NSETFILE 

  !*** 2. default values
  COUTDIR="./"
  CVARSOUT="outflw,storge,rivdph"
  COUTTAG="_cmf"
  LOUTCDF=.FALSE.
  NDLEVEL=0
  LOUTVEC=.FALSE.
  IFRQ_OUT=24                !! daily (24h) output
  LOUTTXT=.FALSE.
  CGAUTXT="None"

  !*** 3. read namelist
  REWIND(CMF_FILES%NSETFILE)
  READ(CMF_FILES%NSETFILE,NML=NOUTPUT)

  !*** 4. Update CMF_OUTPUT settings
  CMF_OUTPUT%COUTDIR = COUTDIR
  CMF_OUTPUT%CVARSOUT = CVARSOUT  
  CMF_OUTPUT%COUTTAG = COUTTAG
  CMF_OUTPUT%LOUTCDF = LOUTCDF
  CMF_OUTPUT%NDLEVEL = NDLEVEL
  CMF_OUTPUT%LOUTVEC = LOUTVEC
  CMF_OUTPUT%IFRQ_OUT = IFRQ_OUT
  CMF_OUTPUT%LOUTTXT = LOUTTXT
  CMF_OUTPUT%CGAUTXT = CGAUTXT

  !*** 5. Write configuration to log
  WRITE(CMF_FILES%LOGNAM,*) "=== NAMELIST, NOUTPUT ==="
  WRITE(CMF_FILES%LOGNAM,*) "COUTDIR:  ", TRIM(COUTDIR)
  WRITE(CMF_FILES%LOGNAM,*) "CVARSOUT: ", TRIM(CVARSOUT)
  WRITE(CMF_FILES%LOGNAM,*) "COUTTAG:  ", TRIM(COUTTAG)
  WRITE(CMF_FILES%LOGNAM,*) "LOUTCDF:  ", LOUTCDF
  IF (LOUTCDF) THEN
    WRITE(CMF_FILES%LOGNAM,*) "NDLEVEL:  ", NDLEVEL
  ENDIF
  IF (LOUTVEC) THEN
    WRITE(CMF_FILES%LOGNAM,*) "LOUTVEC:  ", LOUTVEC
  ENDIF
  WRITE(CMF_FILES%LOGNAM,*) "IFRQ_OUT: ", IFRQ_OUT
  WRITE(CMF_FILES%LOGNAM,*) "LOUTTXT:  ", LOUTTXT
  WRITE(CMF_FILES%LOGNAM,*) "CGAUTXT:  ", TRIM(CGAUTXT)

  CLOSE(CMF_FILES%NSETFILE)

  WRITE(CMF_FILES%LOGNAM,*) "CMF::OUTPUT_NMLIST: end"

END SUBROUTINE CMF_OUTPUT_NMLIST


SUBROUTINE CMF_READ_NMLIST

  call CMF_CONFIG_NMLIST
  call CMF_TIME_NMLIST
  call CMF_MAPS_NMLIST
  call CMF_FORCING_NMLIST
  call CMF_RESTART_NMLIST

  if (CMF_OPTIONS%LSEALEV) then
      call CMF_BOUNDARY_NMLIST
  end if

  if (CMF_OPTIONS%LDAMOUT) then
      call CMF_DAMOUT_NMLIST
  end if

  if (CMF_OPTIONS%LSEDOUT) then
      call CMF_SED_NMLIST
  end if

  if (CMF_OPTIONS%LLEVEE) then
      call CMF_LEVEE_NMLIST
  end if

  call CMF_OUTPUT_NMLIST


END SUBROUTINE CMF_READ_NMLIST

!####################################################################
!####################################################################
! file I/O
!-- INQUIRE_FID : inruire unused file FID
!-- NCERROR     : netCDF I/O wrapper
!####################################################################
FUNCTION INQUIRE_FID() RESULT(FID)
  IMPLICIT NONE
  !* input/output
  INTEGER :: FID ! FILE ID
  !* local variable
  LOGICAL :: I_OPENED ! FILE ID IS ALREADY USED OR NOT?
  !================================================
  DO FID = 10, 999
    INQUIRE(FID,OPENED=I_OPENED)
    IF ( .NOT. I_OPENED ) RETURN
  ENDDO
  END FUNCTION INQUIRE_FID
  !==========================================================


END MODULE CMF_NMLIST_MOD
