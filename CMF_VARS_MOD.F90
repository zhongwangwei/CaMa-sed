module CMF_VARS_MOD
!==========================================================
!* PURPOSE: Shared time-related variables
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
  USE PARKIND1,                ONLY: JPIM, JPRB, JPRM,JPRD
  IMPLICIT NONE
  !======================================
  SAVE
  ! simulation time step
  INTEGER(KIND=JPIM)              :: KSTEP              !! time step since start
  INTEGER(KIND=JPIM)              :: NSTEPS             !! total time step (from start to end), given in CMF_TIME_INIT
  ! elapsed minute from base date (YYYY0,MM0,DD0)
  INTEGER(KIND=JPIM)              :: KMIN               !! KMIN at the start of time step
  INTEGER(KIND=JPIM)              :: KMINNEXT           !! KMIN at the end   of time step
  !
  INTEGER(KIND=JPIM)              :: KMINSTART          !! KMIN at the start of simulation
  INTEGER(KIND=JPIM)              :: KMINEND            !! KMIN at the end   of simulation
  !
  INTEGER(KIND=JPIM)              :: KMINSTAIN          !! KMIN at the start of forcing runoff  data (netCDF)
  INTEGER(KIND=JPIM)              :: KMINSTASL          !! KMIN at the start of boundary sealev data (netCDF)
  ! simulation start date:hour (KMINSTART)
  INTEGER(KIND=JPIM)              :: ISYYYYMMDD         !! date     at simulation start time
  INTEGER(KIND=JPIM)              :: ISHHMM             !! hour+min at simulation start time
  INTEGER(KIND=JPIM)              :: ISYYYY
  INTEGER(KIND=JPIM)              :: ISMM
  INTEGER(KIND=JPIM)              :: ISDD
  INTEGER(KIND=JPIM)              :: ISHOUR
  INTEGER(KIND=JPIM)              :: ISMIN
  ! simulation end   date:hour (KMINEND)
  INTEGER(KIND=JPIM)              :: IEYYYYMMDD         !! date     of simulation end time
  INTEGER(KIND=JPIM)              :: IEHHMM             !! hour+min of simulation end time
  INTEGER(KIND=JPIM)              :: IEYYYY
  INTEGER(KIND=JPIM)              :: IEMM
  INTEGER(KIND=JPIM)              :: IEDD
  INTEGER(KIND=JPIM)              :: IEHOUR
  INTEGER(KIND=JPIM)              :: IEMIN
  !*** date:hour at START of time steop (KMIN)
  INTEGER(KIND=JPIM)              :: IYYYYMMDD          !! date     at the start of time-step
  INTEGER(KIND=JPIM)              :: IYYYY              !! year     at the start of time-step
  INTEGER(KIND=JPIM)              :: IMM                !! month    at the start of time-step
  INTEGER(KIND=JPIM)              :: IDD                !! day      at the start of time-step
  INTEGER(KIND=JPIM)              :: IHHMM              !! hour+min at the start of time-step
  INTEGER(KIND=JPIM)              :: IHOUR              !! hour     at the start of time-step
  INTEGER(KIND=JPIM)              :: IMIN               !! min      at the start of time-step
  !*** date:hour at END   of time steop (KMINNEXT)
  INTEGER(KIND=JPIM)              :: JYYYYMMDD          !! date     at the end   of time-step
  INTEGER(KIND=JPIM)              :: JYYYY              !! year     at the end   of time-step
  INTEGER(KIND=JPIM)              :: JMM                !! month    at the end   of time-step
  INTEGER(KIND=JPIM)              :: JDD                !! day      at the end   of time-step
  INTEGER(KIND=JPIM)              :: JHHMM              !! hour+min at the end   of time-step
  INTEGER(KIND=JPIM)              :: JHOUR              !! hour     at the end   of time-step
  INTEGER(KIND=JPIM)              :: JMIN               !! min      at the end   of time-step
  
  !*** base time to define kmin
  INTEGER(KIND=JPIM)              :: YYYY0              !! base year
  INTEGER(KIND=JPIM)              :: MM0                !! base month
  INTEGER(KIND=JPIM)              :: DD0                !! base day


  !================================================
  !*** river network
  INTEGER(KIND=JPIM),ALLOCATABLE           ::  I2NEXTX(:,:)       !! POINT DOWNSTREAM HORIZONTAL
  INTEGER(KIND=JPIM),ALLOCATABLE           ::  I2NEXTY(:,:)       !! POINT DOWNSTREAM VERTICAL
  
  INTEGER(KIND=JPIM),ALLOCATABLE           ::  I1SEQX(:)          !! 1D SEQUENCE HORIZONTAL
  INTEGER(KIND=JPIM),ALLOCATABLE           ::  I1SEQY(:)          !! 1D SEQUENCE VERTICAL
  INTEGER(KIND=JPIM),ALLOCATABLE           ::  I1NEXT(:)          !! 1D DOWNSTREAM
  INTEGER(KIND=JPIM)                       ::  NSEQRIV            !! LENGTH OF 1D SEQUNECE FOR RIVER
  INTEGER(KIND=JPIM)                       ::  NSEQALL            !! LENGTH OF 1D SEQUNECE FOR RIVER AND MOUTH
  INTEGER(KIND=JPIM)                       ::  NSEQMAX            !! MAX OF NSEQALL (PARALLEL)
  
  INTEGER(KIND=JPIM),ALLOCATABLE           ::  I2VECTOR(:,:)      !! VECTOR INDEX
  INTEGER(KIND=JPIM),ALLOCATABLE           ::  I2REGION(:,:)      !! REGION INDEX
  INTEGER(KIND=JPIM)                       ::  REGIONALL          !! REGION TOTAL
  INTEGER(KIND=JPIM)                       ::  REGIONTHIS         !! REGION THIS CPU
  INTEGER(KIND=JPIM)                       ::  MPI_COMM_CAMA      !! MPI COMMUNICATOR
  
  !================================================
  !*** lat, lon
  REAL(KIND=JPRB),ALLOCATABLE              ::  D1LON(:)           !! longitude [degree_east]
  REAL(KIND=JPRB),ALLOCATABLE              ::  D1LAT(:)           !! latitude  [degree_north]
  
  !================================================
  !*** River + Floodplain topography (map)
  REAL(KIND=JPRB),ALLOCATABLE              ::  D2GRAREA(:,:)      !! GRID AREA [M2]
  REAL(KIND=JPRB),ALLOCATABLE              ::  D2ELEVTN(:,:)      !! ELEVATION [M]
  REAL(KIND=JPRB),ALLOCATABLE              ::  D2NXTDST(:,:)      !! DISTANCE TO THE NEXT GRID [M]
  REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVLEN(:,:)      !! RIVER LENGTH [M]
  REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVWTH(:,:)      !! RIVER WIDTH [M]
  REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVMAN(:,:)      !! RIVER MANNING COEFFICIENT
  REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVHGT(:,:)      !! RIVER HEIGHT [M]
  REAL(KIND=JPRB),ALLOCATABLE              ::  D2FLDHGT(:,:,:)    !! FLOODPLAIN HEIGHT [M]
  
  REAL(KIND=JPRB),ALLOCATABLE              ::  D2GDWDLY(:,:)      !! Ground water delay
  REAL(KIND=JPRB),ALLOCATABLE              ::  D2ELEVSLOPE(:,:)   !! River bed slope
  INTEGER(KIND=JPIM),ALLOCATABLE           ::  I2MASK(:,:)        !! Mask 
  
  !================================================
  !*** Floodplain Topography (diagnosed)
  REAL(KIND=JPRD),ALLOCATABLE              ::  P2RIVSTOMAX(:,:)   !! maximum river storage [m3]
  REAL(KIND=JPRD),ALLOCATABLE              ::  P2FLDSTOMAX(:,:,:) !! MAXIMUM FLOODPLAIN STORAGE [M3]
  
  REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVELV(:,:)      !! elevation of river bed [m3]
  REAL(KIND=JPRB),ALLOCATABLE              ::  D2FLDGRD(:,:,:)    !! FLOODPLAIN GRADIENT
  REAL(KIND=JPRB)                          ::  DFRCINC            !! FLOODPLAIN FRACTION INCREMENT [-] (1/NLFP)
  
  !================================================
  !*** Downstream boundary
  REAL(KIND=JPRB),ALLOCATABLE              ::  D2MEANSL(:,:)      !! MEAN SEA LEVEL [M]
  REAL(KIND=JPRB),ALLOCATABLE              ::  D2SEALEV(:,:)        !! sea level variation [m]
  REAL(KIND=JPRB),ALLOCATABLE              ::  D2DWNELV(:,:)        !! downstream boundary elevation [m]
  
  !*** local variables
#ifdef UseCDF_CMF
  TYPE TYPESL
  CHARACTER(LEN=256)              :: CNAME       !! Netcdf file name
  CHARACTER(LEN=256)              :: CVAR        !! Netcdf variable name 
  INTEGER(KIND=JPIM)              :: NCID        !! Netcdf file     ID
  INTEGER(KIND=JPIM)              :: NVARID      !! Netcdf variable ID
  INTEGER(KIND=JPIM)              :: NSTAID      !! Netcdf station  ID
  INTEGER(KIND=JPIM)              :: NSTART      !! start date of netcdf (KMIN)
  INTEGER(KIND=JPIM)              :: NSTEP       !! steps in netCDF
  END TYPE TYPESL
  TYPE(TYPESL)                    :: SLCDF  !!  Derived type for Sea Level boundary 

  REAL(KIND=JPRM),ALLOCATABLE     :: R1SLIN(:)  ! 1D input boundary condition (m)
  INTEGER(KIND=JPIM),ALLOCATABLE  :: I2SLMAP(:,:)

  !==========================================================
  INTEGER(KIND=JPIM)              :: NVARS              ! temporal output var number
  PARAMETER                         (NVARS=100)          ! actual   output var number
  INTEGER(KIND=JPIM)              :: NVARSOUT
  INTEGER(KIND=JPIM)              :: NVARSOUT_SED
  INTEGER(KIND=JPIM)              :: NVARSOUT_TRACER
  INTEGER(KIND=JPIM)              :: IRECOUT            ! Output file irec
  CHARACTER(LEN=256)              :: CVARSOUT          ! Comma-separated list of output variables to save 
 
  INTEGER(KIND=JPIM)              :: LECMF2LAKEC        ! Output file irec

!*** TYPE for output file    
  TYPE TVAROUT
  CHARACTER(LEN=256)              :: CVNAME             ! output variable name
  CHARACTER(LEN=256)              :: CVLNAME            ! output variable long name
  CHARACTER(LEN=256)              :: CVUNITS            ! output units
  CHARACTER(LEN=256)              :: CFILE              ! output full path file name 
  INTEGER(KIND=JPIM)              :: BINID              ! output binary output file ID
  INTEGER(KIND=JPIM)              :: NCID               ! output netCDF output file ID
  INTEGER(KIND=JPIM)              :: VARID              ! output netCDF output variable ID
  INTEGER(KIND=JPIM)              :: TIMID              ! output netCDF time   variable ID 
  INTEGER(KIND=JPIM)              :: IRECNC             ! Current time record for writting 
  END TYPE TVAROUT 
  TYPE(TVAROUT),ALLOCATABLE       :: VAROUT(:)          ! output variable TYPE set
  TYPE(TVAROUT),ALLOCATABLE       :: VAROUT_SED(:)      ! output variable TYPE set
  TYPE(TVAROUT),ALLOCATABLE       :: VAROUT_TRACER(:)   ! output variable TYPE set

#ifdef UseCDF_CMF
  TYPE TYPEROF
  CHARACTER(LEN=256)              :: CNAME       !! netCDF file name
  CHARACTER(LEN=256)              :: CVAR(3)     !! netCDF variable name
  INTEGER(KIND=JPIM)              :: NCID        !! netCDF file     ID
  INTEGER(KIND=JPIM)              :: NVARID(3)   !! netCDF variable ID
  INTEGER(KIND=JPIM)              :: NSTART      !! Start date of netNDF (in KMIN)
  END TYPE TYPEROF
  TYPE(TYPEROF)                   :: ROFCDF      !! Derived type for Runoff input

  TYPE TYPEROFSED
  CHARACTER(LEN=256)              :: CNAME       !! netCDF file name
  CHARACTER(LEN=256)              :: CVAR(1)     !! netCDF variable name
  INTEGER(KIND=JPIM)              :: NCID        !! netCDF file     ID
  INTEGER(KIND=JPIM)              :: NVARID(1)   !! netCDF variable ID
  INTEGER(KIND=JPIM)              :: NSTART      !! Start date of netNDF (in KMIN)
  END TYPE TYPEROFSED
  TYPE(TYPEROFSED)                 :: SEDCDF      !! Derived type for Sediment input
#endif

#endif

  !================================================
  !*** bifurcation channel
  INTEGER(KIND=JPIM)                       ::  NPTHOUT            !! NUMBER OF FLOODPLAIN PATH
  INTEGER(KIND=JPIM)                       ::  NPTHLEV            !! NUMBER OF FLOODPLAIN PATH LAYER
  INTEGER(KIND=JPIM),ALLOCATABLE           ::  PTH_UPST(:)        !! FLOOD PATHWAY UPSTREAM   ISEQ
  INTEGER(KIND=JPIM),ALLOCATABLE           ::  PTH_DOWN(:)        !! FLOOD PATHWAY DOWNSTREAM JSEQ
  REAL(KIND=JPRB),ALLOCATABLE              ::  PTH_DST(:)         !! FLOOD PATHWAY DISTANCE [m]
  REAL(KIND=JPRB),ALLOCATABLE              ::  PTH_ELV(:,:)         !! FLOOD PATHWAY ELEVATION [m]
  REAL(KIND=JPRB),ALLOCATABLE              ::  PTH_WTH(:,:)         !! FLOOD PATHWAY WIDTH [m]
  REAL(KIND=JPRB),ALLOCATABLE              ::  PTH_MAN(:)         !! FLOOD PATHWAY Manning
  

  !*** Levee Parameters from map
  REAL(KIND=JPRB),ALLOCATABLE    ::  D2LEVHGT(:,:)        !! LEVEE HEIGHT [M] (levee croen elevation above elevtn.bin-river elevation)
  REAL(KIND=JPRB),ALLOCATABLE    ::  D2LEVFRC(:,:)        !! Unprotected fraction = RELATIVE DISTANCE between LEVEE and RIVER [0-1].
                                                          !!  0 = just aside channel, 1 = edge of catchment

  !*** Levee stage parameter (calculated)
  REAL(KIND=JPRB),ALLOCATABLE    ::  D2BASHGT(:,:)        !! LEVEE Base height [M] (levee base elevation above elevtn.bin-river elev)
  REAL(KIND=JPRB),ALLOCATABLE    ::  D2LEVDST(:,:)        !! Absolute DISTANCE between LEVEE and RIVER [0-1]. 
                                                          !! 0 = just aside channel, 1 = edge of catchment

  REAL(KIND=JPRB),ALLOCATABLE    ::  D2LEVBASSTO(:,:)  !! MAXIMUM STORAGE under LEVEE BASE [M3]
  REAL(KIND=JPRB),ALLOCATABLE    ::  D2LEVTOPSTO(:,:)  !! MAXIMUM STORAGE at LEVEE TOP [M3] (only river side)
  REAL(KIND=JPRB),ALLOCATABLE    ::  D2LEVFILSTO(:,:) !! MAXIMUM STORAGE at LEVEE TOP [M3] (both river & protected side are filled)


  !================================================
  ! input matrix (converted from NX:NY*INPN to NSEQMAX*INPN)
  INTEGER(KIND=JPIM),ALLOCATABLE           :: INPX(:,:)        !! INPUT GRID XIN
  INTEGER(KIND=JPIM),ALLOCATABLE           :: INPY(:,:)        !! INPUT GRID YIN
  REAL(KIND=JPRB),ALLOCATABLE              :: INPA(:,:)        !! INPUT AREA
  
  ! input matrix Inverse
  INTEGER(KIND=JPIM),ALLOCATABLE           :: INPXI(:,:,:)        !! OUTPUT GRID XOUT
  INTEGER(KIND=JPIM),ALLOCATABLE           :: INPYI(:,:,:)        !! OUTPUT GRID YOUT
  REAL(KIND=JPRB),ALLOCATABLE              :: INPAI(:,:,:)        !! OUTPUT AREA
  INTEGER(KIND=JPIM)                       :: INPNI               !! MAX INPUT NUMBER for inverse interpolation
  
  DATA REGIONALL  /1/
  DATA REGIONTHIS /1/

  real(kind=JPRB),allocatable     :: d2slope(:,:)     ! floodplain slope [deg]

  !================================================
  ! Pointer was removed in v4.08 in order to keep simple codes when activating Single Precision Mode
  !*** prognostics / state variables initial conditions

  ! Dammy variable for input/output
  REAL(KIND=JPRB),ALLOCATABLE,TARGET     :: D2DAMMY(:,:)       !! Dammy Array for unused variables
  REAL(KIND=JPRB),ALLOCATABLE,TARGET     :: D2COPY(:,:)        !! Dammy Array for Float64/32 switch

  !================================================
  !*** input runoff (interporlated)
  REAL(KIND=JPRB),ALLOCATABLE,TARGET     :: D2RUNOFF(:,:)         !! input runoff             [m3/s]
  REAL(KIND=JPRB),ALLOCATABLE,TARGET     :: D2ROFSUB(:,:)         !! input sub-surface runoff [m3/s]
  REAL(KIND=JPRB),ALLOCATABLE,TARGET     :: D2WEVAP(:,:)          !! input Evaporation [m3/s]

  !================================================
  !*** river & floodpain
  ! storage variables are always in double precision
  REAL(KIND=JPRD),ALLOCATABLE,TARGET     :: P2RIVSTO(:,:)         !! river      storage [m3]
  REAL(KIND=JPRD),ALLOCATABLE,TARGET     :: P2FLDSTO(:,:)         !! floodplain storage [m3]

  REAL(KIND=JPRB),ALLOCATABLE,TARGET     :: D2RIVOUT(:,:)         !! river      outflow [m3/s]
  REAL(KIND=JPRB),ALLOCATABLE,TARGET     :: D2FLDOUT(:,:)         !! floodplain outflow [m3/s]

  !================================================
  !*** for implicit schemes of the local inertial equation
  REAL(KIND=JPRB),ALLOCATABLE,TARGET     :: D2RIVOUT_PRE(:,:)     !! river      outflow [m3/s] (prev t-step)
  REAL(KIND=JPRB),ALLOCATABLE,TARGET     :: D2RIVDPH_PRE(:,:)     !! river      depth   [m]    (prev t-step)
  REAL(KIND=JPRB),ALLOCATABLE,TARGET     :: D2FLDOUT_PRE(:,:)     !! floodplain outflow [m3/s] (prev t-step)
  REAL(KIND=JPRB),ALLOCATABLE,TARGET     :: D2FLDSTO_PRE(:,:)     !! floodplain storage [m3]   (prev t-step)

  !================================================
  !*** Groundwater Delay
  REAL(KIND=JPRD),ALLOCATABLE,TARGET     :: P2GDWSTO(:,:)         !! ground water storage  [m3]
  REAL(KIND=JPRB),ALLOCATABLE,TARGET     :: D2GDWRTN(:,:)         !! Ground water return flow [m3/s]

  !================================================
  !*** These have a different share, not part of the D2PROG array
  REAL(KIND=JPRB),ALLOCATABLE,TARGET     :: D1PTHFLW(:,:)         !! flood path outflow [m3/s]
  REAL(KIND=JPRB),ALLOCATABLE,TARGET     :: D1PTHFLW_PRE(:,:)     !! flood path outflow [m3/s] (prev t-step)

  !================================================
  !!!*** dam variables
  REAL(KIND=JPRD),ALLOCATABLE,TARGET     :: P2DAMSTO(:,:)         !! reservoir storage [m3]
  REAL(KIND=JPRD),ALLOCATABLE,TARGET     :: P2DAMINF(:,:)         !! reservoir inflow [m3/s]; discharge before operation

  !================================================
  !!!*** levee variables
  REAL(KIND=JPRD),ALLOCATABLE,TARGET     :: P2LEVSTO(:,:)         !! flood storage in protected side (storage betwen river & levee)

  !! dam variables
  INTEGER(KIND=JPIM)                     :: IDAM, NDAM  !! number of dams
  INTEGER(KIND=JPIM)                     :: NDAMX       !! exclude dams
  
  INTEGER(KIND=JPIM),ALLOCATABLE         :: DamID(:) !! Dam ID
  CHARACTER(LEN=256),ALLOCATABLE         :: DamName(:)  !! 
  INTEGER(KIND=JPIM),ALLOCATABLE         :: DamIX(:), DamIY(:)  !! IX,IY of dam grid
  REAL(KIND=JPRB),ALLOCATABLE     :: DamLon(:), DamLat(:)  !! longitude, latitude of dam body
  REAL(KIND=JPRB),ALLOCATABLE     :: upreal(:)   !! observed drainage area of reservoir
  REAL(KIND=JPRB),ALLOCATABLE     :: R_VolUpa(:) !! ratio: flood storage capacity / drainage area
  REAL(KIND=JPRB),ALLOCATABLE     :: Qf(:), Qn(:) !! Qf: flood discharge, Qn: normal discharge
  INTEGER(KIND=JPIM),ALLOCATABLE  :: DamYear(:)  !! Dam activation year
  INTEGER(KIND=JPIM),ALLOCATABLE  :: DamStat(:)  !! Dam Status, 2=old, 1=new, -1=not_yet, IMIS=out_of_domain
  
  REAL(KIND=JPRB),ALLOCATABLE     :: EmeVol(:)   !! storage volume to start emergency operation
  REAL(KIND=JPRB),ALLOCATABLE     :: FldVol(:)   !! flood control volume: exclusive for flood control
  REAL(KIND=JPRB),ALLOCATABLE     :: ConVol(:)   !! conservative volume: mainly for water supply
  REAL(KIND=JPRB),ALLOCATABLE     :: NorVol(:)   !! normal storage volume: impoundment
  
  ! internal dam param for stability
  REAL(KIND=JPRB),ALLOCATABLE     :: AdjVol(:)   !! Dam storage for stabilization 
  REAL(KIND=JPRB),ALLOCATABLE     :: Qa(:)       !! Dam outflow for stabilization
  
  
  !*** dam map
  INTEGER(KIND=JPIM),ALLOCATABLE  :: DamSeq(:)   !! coresponding ISEQ of each dam
  INTEGER(KIND=JPIM),ALLOCATABLE  :: I1DAM(:)    !! dam map: 1=dam, 10=upstream of dam, 11: dam grid & downstream is also dam, 0=other
  

  
! Pointer was removed in v4.12 in order to keep simple codes when activating Single Precision Mode
!*** prognostics / state variables initial conditions

!*** Inst. diagnostics 
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2RIVINF(:,:)           !! river      inflow   [m3/s] (from upstream)
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2RIVDPH(:,:)           !! river      depth    [m]
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2RIVVEL(:,:)           !! flow velocity       [m/s]
  
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2FLDINF(:,:)           !! floodplain inflow   [m3/s]
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2FLDDPH(:,:)           !! floodplain depth    [m]
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2FLDFRC(:,:)           !! flooded    fractipn [m2/m2]
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2FLDARE(:,:)           !! flooded    area     [m2]
  
  REAL(KIND=JPRB),ALLOCATABLE                :: D1PTHFLWSUM(:)          !! bifurcation channel flow (1D, not 2D variable), all layer sum
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2PTHOUT(:,:)           !! flood path outflow   [m3/s]
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2PTHINF(:,:)           !! flood path inflow   [m3/s]
  
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2SFCELV(:,:)           !! water surface elev  [m]    (elevtn - rivhgt + rivdph)
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2OUTFLW(:,:)           !! total outflow       [m3/s] (rivout + fldout)
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2STORGE(:,:)           !! total storage       [m3]   (rivsto + fldsto)
  
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2OUTINS(:,:)           !! instantaneous discharge [m3/s] (unrouted runoff)
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2WEVAPEX(:,:)          !! Evaporation water extracted
  
  !================================================
  !*** Average diagnostics for adaptive time step
  REAL(KIND=JPRB)                            :: NADD_adp                    !! sum DT to calculate average
  
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2RIVOUT_aAVG(:,:)       !! average river       discharge
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2OUTFLW_aAVG(:,:)       !! average total outflow       [m3/s] (rivout + fldout)  !! bugfix v362
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2FLDOUT_aAVG(:,:)       !! average floodplain  discharge
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2RIVVEL_aAVG(:,:)       !! average flow velocity
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2PTHOUT_aAVG(:,:)       !! flood pathway net outflow (2D)
  
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2GDWRTN_aAVG(:,:)       !! average ground water return flow
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2RUNOFF_aAVG(:,:)       !! average input runoff
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2ROFSUB_aAVG(:,:)       !! average input sub-surface runoff
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2WEVAPEX_aAVG(:,:)      !! average extracted water evaporation
  
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2DAMINF_aAVG(:,:)       !! average reservoir inflow [m3/s]  !!!added
  
  !*** Average diagnostics (1D) for output
  REAL(KIND=JPRB),ALLOCATABLE                :: D1PTHFLW_aAVG(:,:)       !! bifurcation channel flow (1D, not 2D variable)
  REAL(KIND=JPRB),ALLOCATABLE                :: D1PTHFLWSUM_aAVG(:)      !! bifurcation channel flow (1D, not 2D variable), all layer sum
  
  !*** Daily max diagnostics for output
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2OUTFLW_aMAX(:,:)       !! max total outflow       [m3/s] (rivout + fldout)
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2STORGE_aMAX(:,:)       !! max total outflow       [m3/s] (rivout + fldout)
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2RIVDPH_aMAX(:,:)       !! max total outflow       [m3/s] (rivout + fldout)
  
  !================================================
  !*** Average diagnostics for output
  REAL(KIND=JPRB)                            :: NADD_out                    !! sum DT to calculate average
  
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2RIVOUT_oAVG(:,:)       !! average river       discharge
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2OUTFLW_oAVG(:,:)       !! average total outflow       [m3/s] (rivout + fldout)  !! bugfix v362
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2FLDOUT_oAVG(:,:)       !! average floodplain  discharge
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2RIVVEL_oAVG(:,:)       !! average flow velocity
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2PTHOUT_oAVG(:,:)       !! flood pathway net outflow (2D)
  
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2GDWRTN_oAVG(:,:)       !! average ground water return flow
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2RUNOFF_oAVG(:,:)       !! average input runoff
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2ROFSUB_oAVG(:,:)       !! average input sub-surface runoff
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2WEVAPEX_oAVG(:,:)      !! average extracted water evaporation
  
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2DAMINF_oAVG(:,:)       !! average reservoir inflow [m3/s]  !!!added
  
  !*** Average diagnostics (1D) for output
  REAL(KIND=JPRB),ALLOCATABLE                :: D1PTHFLW_oAVG(:,:)          !! bifurcation channel flow (1D, not 2D variable)
  
  !*** Daily max diagnostics for output
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2OUTFLW_oMAX(:,:)       !! max total outflow       [m3/s] (rivout + fldout)
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2STORGE_oMAX(:,:)       !! max total outflow       [m3/s] (rivout + fldout)
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2RIVDPH_oMAX(:,:)       !! max total outflow       [m3/s] (rivout + fldout)
  
  
  !================================================
  !*** Global total
  ! discharge calculation budget
  REAL(KIND=JPRD)                 :: P0GLBSTOPRE              !! global water storage      [m3] (befre flow calculation)
  REAL(KIND=JPRD)                 :: P0GLBSTONXT              !! global water storage      [m3] (after flow calculation)
  REAL(KIND=JPRD)                 :: P0GLBSTONEW              !! global water storage      [m3] (after runoff input)
  REAL(KIND=JPRD)                 :: P0GLBRIVINF              !! global inflow             [m3] (rivinf + fldinf)
  REAL(KIND=JPRD)                 :: P0GLBRIVOUT              !! global outflow            [m3] (rivout + fldout)
  
  ! stage calculation budget
  REAL(KIND=JPRD)                 :: P0GLBSTOPRE2             !! global water storage      [m3] (befre stage calculation)
  REAL(KIND=JPRD)                 :: P0GLBSTONEW2             !! global water storage      [m3] (after stage calculation)
  REAL(KIND=JPRD)                 :: P0GLBRIVSTO              !! global river storage      [m3]
  REAL(KIND=JPRD)                 :: P0GLBFLDSTO              !! global floodplain storage [m3]
  REAL(KIND=JPRD)                 :: P0GLBLEVSTO              !! global protected-side storage [m3] (levee scheme)
  REAL(KIND=JPRD)                 :: P0GLBFLDARE              !! global flooded area       [m2]
  
  !================================================
  !!!*** levee variables
  REAL(KIND=JPRB),ALLOCATABLE,TARGET         :: D2LEVDPH(:,:)           !! flood depth in protected side (water depth betwen river & levee)


  !*** tracer
  !*** TYPE for tracer data (each tracer should have one data)
  TYPE TTRACE
  CHARACTER(LEN=256)              :: TRCNAME            ! tracer variable name
  CHARACTER(LEN=256)              :: TRCPRE             ! input  file prefix
  CHARACTER(LEN=256)              :: OUTFILE            ! output full path file name 
  INTEGER(KIND=JPIM)              :: BINID              ! tracer binary output file ID
  END TYPE TTRACE
  TYPE(TTRACE),ALLOCATABLE        :: VTRACE(:)          ! tracer variable TYPE set
  REAL(KIND=JPRB),ALLOCATABLE     :: TBUFF(:,:,:)       ! Buffer to store forcing tracer
  INTEGER(KIND=JPIM)              :: ITRACE             !! tracer id

  REAL(KIND=JPRD),ALLOCATABLE     :: P2TRCSTO(:,:)       ! Tracer Storage
  
  REAL(KIND=JPRB),ALLOCATABLE     :: D2TRCDNS(:,:)       ! Tracer Density
  REAL(KIND=JPRB),ALLOCATABLE     :: D2TRCOUT(:,:)       ! Tracer Flux  (main channel)
  REAL(KIND=JPRB),ALLOCATABLE     :: D2TRCINP(:,:)       ! Tracer input (interporlated to catchments)
  
  REAL(KIND=JPRB),ALLOCATABLE     :: D1TRCPFLW(:,:)      ! Tracer Bifurcation Path Flux     (1D: NPTHALL)
  REAL(KIND=JPRB),ALLOCATABLE     :: D2TRCPOUT(:,:)      ! Tracer Bifurcation Path Net Flux (2D: NSEQALL)
  
  !*** Average diagnostics for output
  REAL(KIND=JPRB),ALLOCATABLE     :: D2TRCOUT_oAVG(:,:)   ! Tracer Flux (main channel)
  REAL(KIND=JPRB),ALLOCATABLE     :: D2TRCDNS_oAVG(:,:)   ! Tracer Density
  REAL(KIND=JPRB),ALLOCATABLE     :: D2TRCPOUT_oAVG(:,:)  ! Tracer Bifurcation Path Net Flux (2D: NSEQALL)
  
!================
  !*** sediment
    !================================================  
  
  real(kind=JPRB),allocatable,target  :: d2sedv(:,:,:)    ! storage array for sediment variables
  real(kind=JPRB),pointer             :: d2bedout(:,:)    ! bedflow (m3/s)
  real(kind=JPRB),pointer             :: d2layer(:,:)     ! exchange layer storage (m3)
  real(kind=JPRB),pointer             :: d2netflw(:,:)    ! suspension - deposition (m3/s)
  real(kind=JPRB),pointer             :: d2sedcon(:,:)    ! suspended sediment concentration (m3/m3)
  real(kind=JPRB),pointer             :: d2sedfrc(:,:)    ! sediment distribution fraction [-]
  real(kind=JPRB),pointer             :: d2sedout(:,:)    ! suspended sediment flow (m3/s)
  real(kind=JPRB),pointer             :: d2sedinp(:,:)    ! sediment inflow (m3/s)
  
  real(kind=JPRB),allocatable,target  :: d2depv(:,:,:)    ! storage array for sediment variables
  real(kind=JPRB),pointer             :: d2seddep(:,:,:)  ! deposition storage
  
  real(kind=JPRB),allocatable,target  :: d2sedv_avg(:,:,:)    ! storage array for averaged sediment variables
  real(kind=JPRB),pointer             :: d2bedout_avg(:,:)  ! bedflow (m3/s)
  real(kind=JPRB),pointer             :: d2netflw_avg(:,:)    ! suspension - deposition (m3/s)
  real(kind=JPRB),pointer             :: d2sedout_avg(:,:)    ! suspended sediment flow (m3/s)
  real(kind=JPRB),pointer             :: d2sedinp_avg(:,:)    ! sediment inflow (m3/s)
  
  real(kind=JPRB),allocatable         :: sDiam(:)          ! sediment diameter
  real(kind=JPRB),allocatable         :: setVel(:)         ! setting velocity (m/s)
  
  integer(kind=JPIM)                  :: STEP_SED          ! number of river timesteps within sediment timestep (sedDT/DT)
  real(kind=JPRB)                     :: sadd_riv          ! sum DT to calculate river variable average for sediment
  real(kind=JPRB)                     :: sadd_out          ! sum sedDT to calculate output average
  real(kind=JPRB),allocatable         :: d2rivout_sed(:)   ! accumulate rivout at DT to average into sedDT
  real(kind=JPRB),allocatable         :: d2rivvel_sed(:)   ! accumulate rivvel at DT to average into sedDT
  real(kind=JPRB),allocatable         :: d2rivsto_pre(:)   ! save river storage from previous timestep


  real(kind=JPRB)                 :: dsylunit         ! unit conversion for sediment [m3/km2] -> [m3/m2]
  real(kind=JPRB)                 :: pyld, pyldc, pyldpc  ! parameters for sediment erosion calculation
  
  !*** ILS
  LOGICAL                     :: LLAKEIN
  REAL(KIND=JPRB),ALLOCATABLE :: D2LAKEFRC(:,:)
  REAL(KIND=JPRB),ALLOCATABLE :: D2RUNIN(:,:)
  REAL(KIND=JPRB),ALLOCATABLE :: D2RUNIN_AVG(:,:)

  
end module CMF_VARS_MOD