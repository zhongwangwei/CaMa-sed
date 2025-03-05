MODULE CMF_OUTPUT_MOD
!==========================================================
!* PURPOSE: Control CaMa-Flood standard output file (binary / netCDF)  
!
!* CONTAINS:
! -- CMF_OUTPUT_NMLIST : Read output file info from namelist
! -- CMF_OUTPUT_INIT   : Create & Open standard output files
! -- CMF_OUTPUT_WRITE  : Write output to files
! -- CMF_OUTPUT_END    : Close standard output files
!
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
! shared variables in module
USE PARKIND1,           ONLY: JPIM, JPRB, JPRM,JPRD
USE CMF_NMLIST_MOD,     ONLY: CMF_FILES, CMF_OPTIONS, CMF_CONFIG, CMF_PARAMS
USE CMF_NMLIST_MOD,     ONLY: CMF_TIME, CMF_MAPS, CMF_FORCING, CMF_BOUNDARY
USE CMF_NMLIST_MOD,     ONLY: CMF_RESTART, CMF_DAM, CMF_LEVEE, CMF_OUTPUT,CMF_TRACER,CMF_SED

IMPLICIT NONE
SAVE!
!*** local variables
INTEGER(KIND=JPIM)              :: NVARS              ! temporal output var number
PARAMETER                         (NVARS=100)          ! actual   output var number
INTEGER(KIND=JPIM)              :: NVARSOUT
INTEGER(KIND=JPIM)              :: IRECOUT            ! Output file irec

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
INTEGER(KIND=JPIM)              :: IRECNC               ! Current time record for writting 
END TYPE TVAROUT 
TYPE(TVAROUT),ALLOCATABLE       :: VAROUT(:)          ! output variable TYPE set

CONTAINS

!####################################################################
SUBROUTINE CMF_OUTPUT_WRITE
!======
USE CMF_UTILS_MOD,           ONLY: vecD2mapR
! save results to output files
! -- Called either from "MAIN/Coupler" or CMF_DRV_ADVANCE
USE CMF_VARS_MOD,       ONLY: NSEQMAX, NPTHOUT, NPTHLEV, REGIONTHIS
USE CMF_VARS_MOD,       ONLY: JYYYYMMDD, JHHMM, JHOUR, JMIN, KSTEP
USE CMF_VARS_MOD,       ONLY: P2RIVSTO,     P2FLDSTO,     P2GDWSTO, &
                            & P2DAMSTO,     P2LEVSTO,     D2COPY       !!! added
USE CMF_VARS_MOD,       ONLY: D2RIVDPH,     D2FLDDPH,     D2FLDFRC,     D2FLDARE,     D2SFCELV,     D2STORGE, &
                            & D2OUTFLW_oAVG, D2RIVOUT_oAVG, D2FLDOUT_oAVG, D2PTHOUT_oAVG, D1PTHFLW_oAVG,  &
                            & D2RIVVEL_oAVG, D2GDWRTN_oAVG, D2RUNOFF_oAVG, D2ROFSUB_oAVG, D2WEVAPEX_oAVG, &
                            & D2OUTFLW_oMAX, D2STORGE_oMAX, D2RIVDPH_oMAX, &
                            & D2DAMINF_oAVG, D2OUTINS, D2LEVDPH   !!! added
                            
#ifdef UseMPI_CMF
USE CMF_CTRL_MPI_MOD,   ONLY: CMF_MPI_AllReduce_R2MAP, CMF_MPI_AllReduce_R1PTH
#endif
IMPLICIT NONE
INTEGER(KIND=JPIM)          :: JF
REAL(KIND=JPRB),POINTER     :: D2VEC(:,:) ! point data location to output
!*** LOCAL
REAL(KIND=JPRM)             :: R2OUT(CMF_CONFIG%NX, CMF_CONFIG%NY)
REAL(KIND=JPRM)             :: R1POUT(NPTHOUT,NPTHLEV)
!================================================
WRITE(CMF_FILES%LOGNAM,*) ""
WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"

!*** 0. check date:hour with output frequency
IF ( MOD(JHOUR,CMF_OUTPUT%IFRQ_OUT)==0 .and. JMIN==0 ) THEN             ! JHOUR: end of time step , NFPPH: output frequency (hour)

  !*** 1. update IREC & calc average variable
  IRECOUT=IRECOUT+1 
  WRITE(CMF_FILES%LOGNAM,*) 'CMF::OUTPUT_WRITE: write at time: ', JYYYYMMDD, JHHMM, IRECOUT

  !*** 2. check variable name & allocate data to pointer DVEC
  DO JF=1,NVARSOUT
    SELECT CASE (VAROUT(JF)%CVNAME)
      CASE ('rivsto')
        D2COPY=P2RIVSTO  !! convert Double to Single precision when using SinglePrecisionMode 
        D2VEC => D2COPY  !!   (Storage variables are kept as Float64 in SinglePrecisionMode)
      CASE ('fldsto')
        D2COPY=P2FLDSTO
        D2VEC => D2COPY

      CASE ('rivout')
        D2VEC => D2RIVOUT_oAVG
      CASE ('rivdph')
        D2VEC => D2RIVDPH
      CASE ('rivvel')
        D2VEC => D2RIVVEL_oAVG
      CASE ('fldout')
        D2VEC => D2FLDOUT_oAVG

      CASE ('flddph')
        D2VEC => D2FLDDPH
      CASE ('fldfrc')
        D2VEC => D2FLDFRC
      CASE ('fldare')
        D2VEC => D2FLDARE
      CASE ('sfcelv')
        D2VEC => D2SFCELV

      CASE ('totout')
        D2VEC => D2OUTFLW_oAVG
      CASE ('outflw')            !!  compatibility for previous file name
        D2VEC => D2OUTFLW_oAVG
      CASE ('totsto')
        D2VEC => D2STORGE
      CASE ('storge')            !!  compatibility for previous file name
        D2VEC => D2STORGE

      CASE ('pthout')
        IF( .not. CMF_OPTIONS%LPTHOUT  ) CYCLE
        D2VEC => D2PTHOUT_oAVG
      CASE ('pthflw')
        IF( .not. CMF_OPTIONS%LPTHOUT  ) CYCLE
      CASE ('maxflw')
        D2VEC =>  D2OUTFLW_oMAX
      CASE ('maxdph')
        D2VEC =>  D2RIVDPH_oMAX
      CASE ('maxsto')
        D2VEC =>  D2STORGE_oMAX

      CASE ('outins')
        IF( .not. CMF_OPTIONS%LOUTINS ) CYCLE
        D2VEC =>  D2OUTINS

      CASE ('gwsto')
        IF( .not. CMF_OPTIONS%LGDWDLY ) CYCLE
        D2COPY=P2GDWSTO
        D2VEC =>  D2COPY
      CASE ('gdwsto')
        IF( .not. CMF_OPTIONS%LGDWDLY ) CYCLE
        D2COPY=P2GDWSTO
        D2VEC =>  D2COPY
      CASE ('gwout')
        IF( .not. CMF_OPTIONS%LGDWDLY ) CYCLE
        D2VEC =>  D2GDWRTN_oAVG
      CASE ('gdwrtn')
        IF( .not. CMF_OPTIONS%LGDWDLY ) CYCLE
        D2VEC =>  D2GDWRTN_oAVG

      CASE ('runoff')             !!  compatibility for previous file name
        D2VEC =>  D2RUNOFF_oAVG  
      CASE ('runoffsub')           !!  compatibility for previous file name
        IF( .not. CMF_OPTIONS%LROSPLIT ) CYCLE
        D2VEC =>  D2ROFSUB_oAVG  
      CASE ('rofsfc')
        D2VEC =>  D2RUNOFF_oAVG
      CASE ('rofsub')
        D2VEC =>  D2ROFSUB_oAVG
      CASE ('wevap')
        IF( .not. CMF_OPTIONS%LWEVAP ) CYCLE
        D2VEC => D2WEVAPEX_oAVG

      CASE ('damsto')   !!! added
        IF( .not. CMF_OPTIONS%LDAMOUT ) CYCLE
        D2COPY=P2DAMSTO
        D2VEC => D2COPY
      CASE ('daminf')   !!! added
        IF( .not. CMF_OPTIONS%LDAMOUT ) CYCLE
        D2VEC =>  d2daminf_oAVG

      CASE ('levsto')   !!! added
        IF( .not. CMF_OPTIONS%LLEVEE ) CYCLE
        D2COPY=P2LEVSTO
        D2VEC => D2COPY
      CASE ('levdph')   !!! added
        IF( .not. CMF_OPTIONS%LLEVEE ) CYCLE
        D2VEC =>  D2LEVDPH

      CASE DEFAULT
        WRITE(CMF_FILES%LOGNAM,*) VAROUT(JF)%CVNAME, ' Not defined in CMF_OUTPUT_MOD'

    END SELECT   !! variable name select

    IF( KSTEP==0 .and. CMF_OPTIONS%LOUTINI  )THEN  !! write storage only when LOUTINI specified
      IF ( .not. CMF_OUTPUT%LOUTCDF ) CYCLE
      IF ( VAROUT(JF)%CVNAME/='rivsto' .and. VAROUT(JF)%CVNAME/='fldsto' .and. VAROUT(JF)%CVNAME/='gwsto' ) CYCLE
    ENDIF

!! convert 1Dvector to 2Dmap
    IF( VAROUT(JF)%CVNAME/='pthflw' ) THEN  !! usual 2D map variable
      CALL vecD2mapR(D2VEC,R2OUT)             !! MPI node data is gathered by vecP2mapR
#ifdef UseMPI_CMF
      CALL CMF_MPI_AllReduce_R2MAP(R2OUT)
#endif
    ELSE
      IF( .not. CMF_OPTIONS%LPTHOUT ) CYCLE
      R1POUT(:,:)=REAL(D1PTHFLW_oAVG(:,:))
#ifdef UseMPI_CMF
      CALL CMF_MPI_AllReduce_R1PTH(R1POUT)
#endif
    ENDIF

    !*** 3. write D2VEC to output file
    IF ( CMF_OUTPUT%LOUTCDF ) THEN
      IF ( REGIONTHIS==1 ) CALL WRTE_OUTCDF  !! netCDFG
    ELSE
      IF( VAROUT(JF)%CVNAME=='pthflw' ) THEN
        IF ( REGIONTHIS==1 ) CALL WRTE_OUTPTH(VAROUT(JF)%BINID,IRECOUT,R1POUT)        !! 1D bifu channel
      ELSE
        IF( CMF_OUTPUT%LOUTVEC )THEN
          CALL WRTE_OUTVEC(VAROUT(JF)%BINID,IRECOUT,D2VEC)         !! 1D vector (optional)
        ELSE
          IF ( REGIONTHIS==1 ) CALL WRTE_OUTBIN(VAROUT(JF)%BINID,IRECOUT,R2OUT)         !! 2D map
        ENDIF
      ENDIF
    ENDIF
  END DO

  WRITE(CMF_FILES%LOGNAM,*) 'CMF::OUTPUT_WRITE: end'


ENDIF



!==========================================================
CONTAINS
!+ WRTE_OUTBIN
!+ WRTE_OUTPTH
!+ WRTE_OUTVEC
!+ WRTE_OUTCDF
!==========================================================
SUBROUTINE WRTE_OUTBIN(IFN,IREC,R2OUTDAT)
IMPLICIT NONE
!*** INPUT
INTEGER(KIND=JPIM),INTENT(IN)   :: IFN                 !! FILE NUMBER
INTEGER(KIND=JPIM),INTENT(IN)   :: IREC                !! RECORD
REAL(KIND=JPRM)                 :: R2OUTDAT(CMF_CONFIG%NX, CMF_CONFIG%NY)
!================================================
WRITE(IFN,REC=IREC) R2OUTDAT

END SUBROUTINE WRTE_OUTBIN
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE WRTE_OUTPTH(IFN,IREC,R2OUTDAT)
IMPLICIT NONE
!*** INPUT
INTEGER(KIND=JPIM),INTENT(IN)   :: IFN                 !! FILE NUMBER
INTEGER(KIND=JPIM),INTENT(IN)   :: IREC                !! RECORD
REAL(KIND=JPRM)                 :: R2OUTDAT(NPTHOUT,NPTHLEV)
!================================================
WRITE(IFN,REC=IREC) R2OUTDAT

END SUBROUTINE WRTE_OUTPTH
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE WRTE_OUTVEC(IFN,IREC,D2OUTDAT)
IMPLICIT NONE
!*** INPUT
INTEGER(KIND=JPIM),INTENT(IN)   :: IFN                 !! FILE NUMBER
INTEGER(KIND=JPIM),INTENT(IN)   :: IREC                !! RECORD
REAL(KIND=JPRB),INTENT(IN)      :: D2OUTDAT(NSEQMAX,1) !! OUTPUT DATA
!*** LOCAL
REAL(KIND=JPRM)                 :: R2OUTDAT(NSEQMAX,1)
!================================================
R2OUTDAT(:,:)=REAL(D2OUTDAT(:,:))
WRITE(IFN,REC=IREC) R2OUTDAT

END SUBROUTINE WRTE_OUTVEC
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE WRTE_OUTCDF
#ifdef UseCDF_CMF
USE NETCDF 
USE CMF_VARS_MOD,            ONLY: KMINSTART,KMINNEXT
USE CMF_UTILS_MOD,           ONLY: NCERROR
IMPLICIT NONE
REAL(KIND=JPRB)                 :: XTIME ! seconds since start of the run !

!================================================
XTIME=REAL( (KMINNEXT-KMINSTART),JPRB) *60._JPRB      !! for netCDF

! Write time variable
CALL NCERROR( NF90_PUT_VAR(VAROUT(JF)%NCID,VAROUT(JF)%TIMID,XTIME,(/VAROUT(JF)%IRECNC/)) )

! Write data variable
CALL NCERROR( NF90_PUT_VAR(VAROUT(JF)%NCID,VAROUT(JF)%VARID,R2OUT(1:CMF_CONFIG%NX,1:CMF_CONFIG%NY),&
              (/1,1,VAROUT(JF)%IRECNC/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,1/)) )

! update IREC
VAROUT(JF)%IRECNC=VAROUT(JF)%IRECNC+1

! Force sync after each write
CALL NCERROR( NF90_SYNC(VAROUT(JF)%NCID) )

! Debug message removed
#endif
END SUBROUTINE WRTE_OUTCDF 
!==========================================================

END SUBROUTINE CMF_OUTPUT_WRITE
!####################################################################





!####################################################################



!####################################################################
SUBROUTINE CMF_OUTTXT_WRTE
USE CMF_VARS_MOD,       ONLY: D2OUTFLW
USE CMF_VARS_MOD,       ONLY: IYYYYMMDD,ISYYYY
USE CMF_VARS_MOD,        ONLY: I2VECTOR
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID

! local
INTEGER(KIND=JPIM)                  :: GID, GIX, GIY, GISEQ
CHARACTER(len=256),SAVE             :: GNAME

INTEGER(KIND=JPIM),SAVE             :: IGAUGE, NGAUGE, NGAUGEX
INTEGER(KIND=JPIM),ALLOCATABLE,SAVE :: WriteID(:), WriteISEQ(:)
CHARACTER(len=9),ALLOCATABLE,SAVE   :: WriteName(:)
REAL(KIND=JPRB),ALLOCATABLE,SAVE    :: WriteOut(:)

! File IO
INTEGER(KIND=JPIM),SAVE             :: LOGOUTTXT
CHARACTER(len=4),SAVE               :: cYYYY
CHARACTER(len=256),SAVE             :: CLEN, CFMT
CHARACTER(len=256),SAVE             :: COUTTXT
LOGICAL,SAVE                        :: IsOpen
DATA IsOpen       /.FALSE./

! ======

IF( CMF_OUTPUT%LOUTTXT )THEN

  IF( .not. IsOpen)THEN
    IsOpen=.TRUE.

    NGAUGEX=0
    LOGOUTTXT=INQUIRE_FID()
    OPEN(LOGOUTTXT,FILE=CMF_OUTPUT%CGAUTXT,FORM='formatted',STATUS='old')
    READ(LOGOUTTXT,*) NGAUGE
    DO IGAUGE=1, NGAUGE
      READ(LOGOUTTXT,*) GID, GNAME, GIX, GIY
      IF( I2VECTOR(GIX,GIY)>0 )THEN
        NGAUGEX=NGAUGEX+1
      ENDIF
    END DO
    CLOSE(LOGOUTTXT)

    ALLOCATE( WriteID(NGAUGEX),WriteISEQ(NGAUGEX),WriteOut(NGAUGEX),WriteName(NGAUGEX))

    NGAUGEX=0
    OPEN(LOGOUTTXT,FILE=CMF_OUTPUT%CGAUTXT,FORM='formatted',STATUS='old')
    READ(LOGOUTTXT,*) NGAUGE
    DO IGAUGE=1, NGAUGE
      READ(LOGOUTTXT,*) GID, GNAME, GIX, GIY
      IF( I2VECTOR(GIX,GIY)>0 )THEN
        NGAUGEX=NGAUGEX+1
        WriteID(NGAUGEX)  =GID
        WriteName(NGAUGEX)=TRIM(GNAME)
        WriteISEQ(NGAUGEX)=I2VECTOR(GIX,GIY)
      ENDIF
    END DO
    CLOSE(LOGOUTTXT)

    ! ============
    WRITE(CYYYY,'(i4.4)') ISYYYY
    COUTTXT='./outtxt-'//TRIM(cYYYY)//'.txt'

    LOGOUTTXT=INQUIRE_FID()
    OPEN(LOGOUTTXT,FILE=COUTTXT,FORM='formatted')

    WRITE(CLEN,'(i0)') NGAUGE
    CFMT="(i10,"//TRIM(CLEN)//"(i10))"
    WRITE(LOGOUTTXT,CFMT) NGAUGEX, ( WriteID(IGAUGE),IGAUGE=1,NGAUGEX )

    CFMT="(i10,"//TRIM(CLEN)//"(x,a9))"
    WRITE(LOGOUTTXT,CFMT) NGAUGEX, ( WriteName(IGAUGE),IGAUGE=1,NGAUGEX )


    CFMT="(i10,"//TRIM(CLEN)//"(f10.2))"
  ENDIF

  DO IGAUGE=1, NGAUGEX
    GISEQ=WriteISEQ(IGAUGE)
    WriteOut(IGAUGE) = D2OUTFLW(GISEQ,1)
  END DO
  
  WRITE(LOGOUTTXT,CFMT) IYYYYMMDD, ( WriteOUT(IGAUGE),IGAUGE=1,NGAUGEX )

ENDIF

END SUBROUTINE CMF_OUTTXT_WRTE
!####################################################################



!####################################################################
SUBROUTINE CMF_DAMOUT_WRTE
  ! local
  USE CMF_VARS_MOD,            ONLY: IDAM, NDAM, NDAMX, DamID, DamName, DamIX, DamIY
  USE CMF_VARS_MOD,            ONLY: IYYYYMMDD,ISYYYY
  USE CMF_VARS_MOD,            ONLY: DamLon, DamLat, upreal, R_VolUpa, Qf, Qn, DamYear
  USE CMF_VARS_MOD,            ONLY: DamStat, EmeVol, FldVol, ConVol, NorVol, AdjVol, Qa, DamSeq, I1DAM
  USE CMF_UTILS_MOD,           ONLY: INQUIRE_FID
  USE CMF_VARS_MOD,            ONLY: D2RIVOUT, D2FLDOUT, P2RIVSTO, P2FLDSTO, P2DAMSTO, P2DAMINF, D2RUNOFF

  CHARACTER(len=36)          :: WriteTXT(NDAMX), WriteTXT2(NDAMX)
  REAL(KIND=JPRB)            :: DDamInf, DDamOut
  ! File IO
  INTEGER(KIND=JPIM),SAVE    :: ISEQD, JDAM
  INTEGER(KIND=JPIM),SAVE    :: LOGDAM
  CHARACTER(len=4),SAVE      :: cYYYY
  CHARACTER(len=256),SAVE    :: CLEN, CFMT
  CHARACTER(len=256),SAVE    :: DAMTXT
  LOGICAL,SAVE               :: IsOpen
  DATA IsOpen       /.FALSE./
  ! ==========================================
  
  IF( CMF_DAM%LDAMTXT ) THEN
  
    IF( .not. IsOpen )THEN
      IsOpen=.TRUE.
      WRITE(CYYYY,'(i4.4)') ISYYYY
      DAMTXT='./damtxt-'//trim(cYYYY)//'.txt'
  
      LOGDAM=INQUIRE_FID()
      OPEN(LOGDAM,FILE=DAMTXT,FORM='formatted')
  
      WRITE(CLEN,'(i0)') NDAMX
      CFMT="(i10,"//TRIM(CLEN)//"(a36))"
  
      JDAM=0
      DO IDAM=1, NDAM
        IF( DamStat(IDAM)==CMF_PARAMS%IMIS )CYCLE
        JDAM=JDAM+1
        IF( DamStat(IDAM)==-1 ) THEN   !! dam not activated yet
          WRITE(WriteTxt(JDAM), '(i12,2f12.2)') DamID(IDAM), -9., -9.
          WRITE(WriteTxt2(JDAM),'(3f12.2)') upreal(IDAM),   Qf(IDAM), Qn(IDAM)
        ELSE
          WRITE(WriteTxt(JDAM), '(i12,2f12.2)') DamID(IDAM), (FldVol(IDAM)+ConVol(IDAM))*1.E-9, ConVol(IDAM)*1.E-9
          WRITE(WriteTxt2(JDAM),'(3f12.2)') upreal(IDAM),   Qf(IDAM), Qn(IDAM)
        ENDIF
      END DO
  
      WRITE(LOGDAM,CFMT) NDAMX, (WriteTXT(JDAM) ,JDAM=1, NDAMX)
  
      IF( NDAMX > 0 )THEN
        WRITE(CLEN,'(i0)') NDAMX
        CFMT="(i10,"//TRIM(CLEN)//"(a36))"
      ELSE
        CFMT="(i10)"  ! Default format when NDAMX=0
      ENDIF
      WRITE(LOGDAM,CFMT)  "Date", (WriteTXT2(JDAM),JDAM=1, NDAMX)
    ENDIF
  
    JDAM=0
    DO IDAM=1, NDAM
      IF( DamStat(IDAM)==CMF_PARAMS%IMIS ) CYCLE
      JDAM=JDAM+1
      ISEQD=DamSeq(IDAM)
      DDamInf=P2DAMINF(ISEQD,1)
      DDamOut=D2RIVOUT(ISEQD,1) + D2FLDOUT(ISEQD,1)
      WRITE(WriteTxt(JDAM), '(3f12.2)') P2DAMSTO(ISEQD,1)*1.E-9, DDamInf, DDamOut
    END DO
  
    CFMT="(i10,"//TRIM(CLEN)//"(a36))"  
    WRITE(LOGDAM,CFMT) IYYYYMMDD, (WriteTXT(JDAM),JDAM=1, NDAMX)
  
  ENDIF
  
  END SUBROUTINE CMF_DAMOUT_WRTE
  !####################################################################
  
  


!####################################################################
SUBROUTINE CMF_RESTART_WRITE
  ! write restart files
  ! -- called CMF_from DRV_ADVANCE
  USE CMF_VARS_MOD,       ONLY: KSTEP,  NSTEPS, JYYYYMMDD, JHHMM, JDD, JHOUR, JMIN
  USE CMF_VARS_MOD,        ONLY: NPTHOUT,     NPTHLEV
  USE CMF_VARS_MOD,       ONLY: P2RIVSTO,    P2FLDSTO,    D2RIVOUT_PRE,D2FLDOUT_PRE, &
                              & D1PTHFLW_PRE,D2RIVDPH_PRE,D2FLDSTO_PRE,P2GDWSTO, &
                              & P2DAMSTO,    P2LEVSTO
  USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
  IMPLICIT NONE
  !* local variable
  INTEGER(KIND=JPIM)         :: IREST
  !================================================
  IREST=0
  
  IF ( CMF_RESTART%IFRQ_RST>=0 .and. KSTEP==NSTEPS )THEN         !! end of run
    IREST=1
  ENDIF
  
  IF ( CMF_RESTART%IFRQ_RST>=1 .and. CMF_RESTART%IFRQ_RST<=24 )THEN
    IF ( MOD(JHOUR,CMF_RESTART%IFRQ_RST)==0 .and. JMIN==0 )THEN  !! at selected hour
      IREST=1
    ENDIF
  ENDIF
  
  IF ( CMF_RESTART%IFRQ_RST==30 )THEN
    IF ( JDD==1 .and. JHOUR==0 .and. JMIN==0 )THEN  !! at start of month
      IREST=1
    ENDIF
  ENDIF
  
  
  IF( IREST==1 )THEN
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
    WRITE(CMF_FILES%LOGNAM,*) 'CMF::RESTART_WRITE: write time: ' , JYYYYMMDD, JHHMM
  
  
    IF( CMF_RESTART%LRESTCDF )THEN
      CALL WRTE_REST_CDF  !! netCDF restart write
    ELSE
      CALL WRTE_REST_BIN
    ENDIF
  END IF 
  
  CONTAINS
  !==========================================================
  !+ WRTE_REST_BIN
  !+ WRTE_REST_CDF
  !==========================================================
  SUBROUTINE WRTE_REST_BIN
  USE CMF_VARS_MOD,       ONLY: JYYYYMMDD, JHOUR
  USE CMF_VARS_MOD,        ONLY: REGIONTHIS, NSEQMAX
#ifdef UseMPI_CMF
  USE CMF_MPI_MOD,   ONLY: CMF_MPI_AllReduce_R1PTH, CMF_MPI_AllReduce_P1PTH
#endif
  IMPLICIT NONE
  ! local variable 
  INTEGER(KIND=JPIM)         :: RIREC
  CHARACTER(LEN=256)         :: CFILE,CDATE
  REAL(KIND=JPRD)            :: P2TMP(NSEQMAX,1)        !! use Real*8 for code simplicity
  REAL(KIND=JPRD)            :: P1PTH(NPTHOUT,NPTHLEV) 
  REAL(KIND=JPRM)            :: R1PTH(NPTHOUT,NPTHLEV) 
  !================================================
  !*** set file nam
  WRITE(CDATE,'(I8.8,I2.2)') JYYYYMMDD,JHOUR
  CFILE=TRIM(CMF_RESTART%CRESTDIR)//TRIM(CMF_RESTART%CVNREST)//TRIM(CDATE)//TRIM(CMF_PARAMS%CSUFBIN)
  WRITE(CMF_FILES%LOGNAM,*) 'WRTE_REST_BIN: restart file:',CFILE
  
  !*** write restart data (2D map)
  CMF_FILES%TMPNAM=INQUIRE_FID()
  
  IF( CMF_RESTART%LRESTDBL )THEN
    IF ( REGIONTHIS==1 ) OPEN(CMF_FILES%TMPNAM,FILE=CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8*CMF_CONFIG%NX*CMF_CONFIG%NY)
  ELSE
    IF ( REGIONTHIS==1 ) OPEN(CMF_FILES%TMPNAM,FILE=CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
  ENDIF
  
  RIREC=0
    P2TMP=P2RIVSTO
     CALL WRTE_BIN_MAP(P2TMP,CMF_FILES%TMPNAM,RIREC)
    P2TMP=P2FLDSTO
     CALL WRTE_BIN_MAP(P2TMP,CMF_FILES%TMPNAM,RIREC)
  !!================
  !! additional restart data for optional schemes (only write required vars)
    IF ( .not. CMF_OPTIONS%LSTOONLY )THEN           !! default restart with previous t-step outflw
      P2TMP=D2RIVOUT_PRE
       CALL WRTE_BIN_MAP(P2TMP,CMF_FILES%TMPNAM,RIREC)
      P2TMP=D2FLDOUT_PRE
       CALL WRTE_BIN_MAP(P2TMP,CMF_FILES%TMPNAM,RIREC)
      P2TMP=D2RIVDPH_PRE
       CALL WRTE_BIN_MAP(P2TMP,CMF_FILES%TMPNAM,RIREC)
      P2TMP=D2FLDSTO_PRE
       CALL WRTE_BIN_MAP(P2TMP,CMF_FILES%TMPNAM,RIREC)
    ENDIF
  
    IF ( CMF_OPTIONS%LGDWDLY ) THEN
      P2TMP=P2GDWSTO
       CALL WRTE_BIN_MAP(P2TMP,CMF_FILES%TMPNAM,RIREC)
    ENDIF
    IF ( CMF_OPTIONS%LDAMOUT ) THEN   !!! ADDED
      P2TMP=P2DAMSTO
       CALL WRTE_BIN_MAP(P2TMP,CMF_FILES%TMPNAM,RIREC)
    ENDIF
    IF ( CMF_OPTIONS%LLEVEE ) THEN   !!! ADDED
      P2TMP=P2LEVSTO
       CALL WRTE_BIN_MAP(P2TMP,CMF_FILES%TMPNAM,RIREC)
    ENDIF
  
  CLOSE(CMF_FILES%TMPNAM)
  
  !*** write restart data (1D bifucation chanenl)
  IF( CMF_OPTIONS%LPTHOUT)THEN
  
    CFILE=TRIM(CMF_RESTART%CRESTDIR)//TRIM(CMF_RESTART%CVNREST)//TRIM(CDATE)//TRIM(CMF_PARAMS%CSUFBIN)//'.pth'
    WRITE(CMF_FILES%LOGNAM,*) 'WRTE_REST: WRITE RESTART BIN:',CFILE
  
    !! Double Precision Restart
    IF( CMF_RESTART%LRESTDBL )THEN
      P1PTH(:,:)=D1PTHFLW_PRE(:,:)
#ifdef UseMPI_CMF
      CALL CMF_MPI_AllReduce_P1PTH(P1PTH)
#endif
      IF ( REGIONTHIS==1 )THEN
        OPEN(CMF_FILES%TMPNAM,FILE=CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8*NPTHOUT*NPTHLEV)
        WRITE(CMF_FILES%TMPNAM,REC=1) P1PTH
        CLOSE(CMF_FILES%TMPNAM)
      ENDIF
  
    !! Single Precision Restart
    ELSE
      R1PTH(:,:)=D1PTHFLW_PRE(:,:)
#ifdef UseMPI_CMF
      CALL CMF_MPI_AllReduce_R1PTH(R1PTH)
#endif
      IF ( REGIONTHIS==1 )THEN
        OPEN(CMF_FILES%TMPNAM,FILE=CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NPTHOUT*NPTHLEV)
        WRITE(CMF_FILES%TMPNAM,REC=1) R1PTH
        CLOSE(CMF_FILES%TMPNAM)
      ENDIF
    ENDIF
  ENDIF
  
  END SUBROUTINE WRTE_REST_BIN
  !=================
  SUBROUTINE WRTE_BIN_MAP(P2VAR,TNAM,IREC)
  USE CMF_UTILS_MOD,      ONLY: vecP2mapP,  vecP2mapR
  USE CMF_VARS_MOD,       ONLY: REGIONTHIS, NSEQMAX
#ifdef UseMPI_CMF
  USE CMF_MPI_MOD,   ONLY: CMF_MPI_AllReduce_R2MAP, CMF_MPI_AllReduce_P2MAP
#endif
  IMPLICIT NONE
  REAL(KIND=JPRD)            :: P2VAR(NSEQMAX,1)  !! use Real*8 for code simplicity
  INTEGER(KIND=JPIM)         :: TNAM,IREC
  !* local
  REAL(KIND=JPRM)            :: R2TEMP(CMF_CONFIG%NX,CMF_CONFIG%NY)
  REAL(KIND=JPRD)            :: P2TEMP(CMF_CONFIG%NX,CMF_CONFIG%NY)
  !=================
  IREC=IREC+1
  
  !! Double Precision Restart
  IF( CMF_RESTART%LRESTDBL )THEN
    CALL vecP2mapP(P2VAR,P2TEMP)  
#ifdef UseMPI_CMF
    CALL CMF_MPI_AllReduce_P2MAP(P2TEMP)
#endif
  
    IF ( REGIONTHIS==1 ) WRITE(TNAM,REC=IREC) P2TEMP
  
  !! Single Precision Restart
  ELSE
    CALL vecP2mapR(P2VAR,R2TEMP)  
#ifdef UseMPI_CMF
      CALL CMF_MPI_AllReduce_R2MAP(R2TEMP)
#endif
    IF ( REGIONTHIS==1 ) WRITE(TNAM,REC=IREC) R2TEMP
  ENDIF
  !=================
  END SUBROUTINE WRTE_BIN_MAP
  !==========================================================
  !+
  !+
  !+
  !==========================================================
  SUBROUTINE WRTE_REST_CDF
#ifdef UseCDF_CMF
  USE NETCDF
  USE CMF_VARS_MOD,       ONLY: KMINNEXT, KMINSTART, ISYYYY,ISMM,ISDD, ISHOUR, ISMIN
  USE CMF_VARS_MOD,       ONLY: JYYYYMMDD,JHOUR
  USE CMF_VARS_MOD,        ONLY: D1LON,    D1LAT,     REGIONTHIS, NSEQMAX
  USE CMF_UTILS_MOD,      ONLY: NCERROR,  vecP2mapP
#ifdef UseMPI_CMF
  USE CMF_MPI_MOD,   ONLY: CMF_MPI_AllReduce_P2MAP, CMF_MPI_AllReduce_P1PTH
#endif
  IMPLICIT NONE
  !* local variable
  CHARACTER(LEN=256)         :: CFILE, CDATE, CTIME, CVAR
  INTEGER(KIND=JPIM)         :: NCID,  VARID, LATID, LONID, TIMEID, JF, &
                                NPTHOUTID,    NPTHLEVID,    STATUS, IOUT
  REAL(KIND=JPRB)            :: XTIME ! seconds since start of the run ! 
  REAL(KIND=JPRD)            :: P2VEC(NSEQMAX,1), P2TEMP(CMF_CONFIG%NX,CMF_CONFIG%NY), P1PTH(NPTHOUT,NPTHLEV)
  !================================================
  !*** 1. set file name & tim
  XTIME=REAL( (KMINNEXT-KMINSTART),JPRB) *60._JPRB
  WRITE(CTIME,'(A14,I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)') 'seconds since ',ISYYYY,'-',ISMM,'-',ISDD,' ',ISHOUR,":",ISMIN
  
  WRITE(CDATE,'(I8.8,I2.2)') JYYYYMMDD,JHOUR
  CFILE=TRIM(CMF_RESTART%CRESTDIR)//TRIM(CMF_RESTART%CVNREST)//TRIM(CDATE)//TRIM(CMF_PARAMS%CSUFCDF)
  WRITE(CMF_FILES%LOGNAM,*) 'WRTE_REST:create RESTART NETCDF:',CFILE
  
  !============================
  !*** 2. create netCDF file
  !! Note: all restart variables are saved as Float64.
  IF( REGIONTHIS==1 )THEN   !! write restart only on master node
  
    CALL NCERROR( NF90_CREATE(CFILE,NF90_NETCDF4,NCID),'CREATING FILE:'//TRIM(CFILE) )
    
    !! dimensions 
    CALL NCERROR( NF90_DEF_DIM(NCID, 'time', NF90_UNLIMITED, TIMEID) )
    CALL NCERROR( NF90_DEF_DIM(NCID, 'lat', CMF_CONFIG%NY, LATID) )
    CALL NCERROR( NF90_DEF_DIM(NCID, 'lon', CMF_CONFIG%NX, LONID) )
    
    IF ( CMF_OPTIONS%LPTHOUT ) THEN
      CALL NCERROR( NF90_DEF_DIM(NCID, 'NPTHOUT', NPTHOUT, NPTHOUTID) )
      CALL NCERROR( NF90_DEF_DIM(NCID, 'NPTHLEV', NPTHLEV, NPTHLEVID) )
    ENDIF 
    
    !! dimentions 
    CALL NCERROR( NF90_DEF_VAR(NCID, 'lat', NF90_FLOAT, (/LATID/), VARID) )
    CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'long_name','latitude') )
    CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'units','degrees_north') )
    
    CALL NCERROR( NF90_DEF_VAR(NCID, 'lon', NF90_FLOAT, (/LONID/), VARID) )
    CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'long_name','longitude') )
    CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'units','degrees_east') )
    
    CALL NCERROR( NF90_DEF_VAR(NCID, 'time', NF90_DOUBLE, (/TIMEID/), VARID) ) 
    CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'long_name','time') )
    CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'units',CTIME) )
    
    !! variables
    CALL NCERROR( NF90_DEF_VAR(NCID, 'rivsto', NF90_DOUBLE, (/LONID,LATID,TIMEID/), &
                               VARID,DEFLATE_LEVEL=6), 'Creating Variable')
    CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'long_name',"river storage" ) )
    CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'units',"m3") )
    CALL NCERROR( NF90_PUT_ATT(NCID, VARID, '_FillValue',REAL(CMF_PARAMS%DMIS,KIND=JPRD)),'in here?' )
    
     
    CALL NCERROR( NF90_DEF_VAR(NCID, 'fldsto', NF90_DOUBLE, (/LONID,LATID,TIMEID/), &
                               VARID,DEFLATE_LEVEL=6), 'Creating Variable')  
    CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'long_name',"flood plain storage" ) )
    CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'units',"m3") )
    CALL NCERROR( NF90_PUT_ATT(NCID, VARID, '_FillValue',REAL(CMF_PARAMS%DMIS,KIND=JPRD)) )
    
    IF ( .not. CMF_OPTIONS%LSTOONLY )THEN           !! default restart with previous t-step outflw
      CALL NCERROR( NF90_DEF_VAR(NCID, 'rivout_pre', NF90_DOUBLE, (/LONID,LATID,TIMEID/),&
                                 VARID,DEFLATE_LEVEL=6), 'Creating Variable')  
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'long_name',"river outflow prev" ) )
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'units',"m3/s") )
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, '_FillValue',REAL(CMF_PARAMS%DMIS,KIND=JPRD)) )
      
      CALL NCERROR( NF90_DEF_VAR(NCID, 'fldout_pre', NF90_DOUBLE, (/LONID,LATID,TIMEID/), &
                                 VARID,DEFLATE_LEVEL=6), 'Creating Variable')  
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'long_name',"floodplain outflow prev" ) )
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'units',"m3/s") )
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, '_FillValue',REAL(CMF_PARAMS%DMIS,KIND=JPRD)) )
      
      CALL NCERROR( NF90_DEF_VAR(NCID, 'rivdph_pre', NF90_DOUBLE, (/LONID,LATID,TIMEID/), &
                                 VARID,DEFLATE_LEVEL=6), 'Creating Variable')  
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'long_name',"river depth prev" ) )
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'units',"m") )
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, '_FillValue',REAL(CMF_PARAMS%DMIS,KIND=JPRD)) )
      
      CALL NCERROR( NF90_DEF_VAR(NCID, 'fldsto_pre', NF90_DOUBLE, (/LONID,LATID,TIMEID/), &
                                 VARID,DEFLATE_LEVEL=6), 'Creating Variable')  
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'long_name',"floodplain storage prev" ) )
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'units',"m3") )
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, '_FillValue',REAL(CMF_PARAMS%DMIS,KIND=JPRD)) )
    
      !! optional variables
      IF ( CMF_OPTIONS%LPTHOUT ) THEN
        CALL NCERROR( NF90_DEF_VAR(NCID, 'pthflw_pre', NF90_DOUBLE, (/NPTHOUTID,NPTHLEVID,TIMEID/),&
                                   VARID,DEFLATE_LEVEL=6) ) 
        CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'long_name',"bifurcation outflow pre" ) )
        CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'units',"m3/s") )
      ENDIF
    ENDIF
    
    IF ( CMF_OPTIONS%LGDWDLY ) THEN
      CALL NCERROR( NF90_DEF_VAR(NCID, 'gdwsto', NF90_DOUBLE, (/LONID,LATID,TIMEID/), &
                               VARID,DEFLATE_LEVEL=6), 'Creating Variable gdwsto')  
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'long_name',"ground water storage" ) )
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'units',"m3") )
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, '_FillValue',REAL(CMF_PARAMS%DMIS,KIND=JPRD)) )
    ENDIF
    
    IF ( CMF_OPTIONS%LDAMOUT ) THEN    !!! added
      CALL NCERROR( NF90_DEF_VAR(NCID, 'damsto', NF90_DOUBLE, (/LONID,LATID,TIMEID/), &
                               VARID,DEFLATE_LEVEL=6), 'Creating Variable dasmto')  
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'long_name',"dam reservoir storage" ) )
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'units',"m3") )
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, '_FillValue',REAL(CMF_PARAMS%DMIS,KIND=JPRD)) )
    ENDIF
    
    IF ( CMF_OPTIONS%LLEVEE ) THEN    !!! added
      CALL NCERROR( NF90_DEF_VAR(NCID, 'levsto', NF90_DOUBLE, (/LONID,LATID,TIMEID/), &
                               VARID,DEFLATE_LEVEL=6), 'Creating Variable levsto')  
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'long_name',"storage exceeds levee protection" ) )
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, 'units',"m3") )
      CALL NCERROR( NF90_PUT_ATT(NCID, VARID, '_FillValue',REAL(CMF_PARAMS%DMIS,KIND=JPRD)) )
    ENDIF
  
    CALL NCERROR( NF90_ENDDEF(NCID) )
    !============================
    !*** 2. write data
    
    !! dimentions (time,lon,lat)
    CALL NCERROR( NF90_INQ_VARID(NCID,'time',VARID))
    CALL NCERROR( NF90_PUT_VAR(NCID,VARID,XTIME) )
    
    CALL NCERROR ( NF90_INQ_VARID(NCID,'lon',VARID),'getting id' )
    CALL NCERROR( NF90_PUT_VAR(NCID,VARID,D1LON))
    
    CALL NCERROR ( NF90_INQ_VARID(NCID,'lat',VARID),'getting id' )
    CALL NCERROR( NF90_PUT_VAR(NCID,VARID,D1LAT))
  
  ENDIF  !! regionthis=1: definition
  
  !! write restart variables (gather data in MPI mode)
  DO JF=1,9
    IOUT=0
    SELECT CASE(JF)
      CASE (1)
        CVAR='rivsto'
        CALL vecP2mapP(P2RIVSTO,P2TEMP)
        IOUT=1
      CASE (2)
        CVAR='fldsto'
        CALL vecP2mapP(P2FLDSTO,P2TEMP)
        IOUT=1
      CASE (3)
        CVAR='rivout_pre'
        IF( .not. CMF_OPTIONS%LSTOONLY ) THEN
          P2VEC(:,:)=D2RIVOUT_PRE(:,:)
          CALL vecP2mapP(P2VEC,P2TEMP)
          IOUT=1
        ENDIF
      CASE (4)
        CVAR='fldout_pre'
        IF( .not. CMF_OPTIONS%LSTOONLY ) THEN
          P2VEC(:,:)=D2FLDOUT_PRE(:,:)
          CALL vecP2mapP(P2VEC,P2TEMP)
          IOUT=1
        ENDIF
     CASE (5)
        CVAR='rivdph_pre'
        IF( .not. CMF_OPTIONS%LSTOONLY ) THEN
          P2VEC(:,:)=D2RIVDPH_PRE(:,:)
          CALL vecP2mapP(P2VEC,P2TEMP)
          IOUT=1
        ENDIF
      CASE (6)
        CVAR='fldsto_pre'
        IF( .not. CMF_OPTIONS%LSTOONLY ) THEN
          P2VEC(:,:)=D2FLDSTO_PRE(:,:)
          CALL vecP2mapP(P2VEC,P2TEMP)
          IOUT=1
        ENDIF
      CASE (7)
        CVAR='gdwsto'
        IF( CMF_OPTIONS%LGDWDLY ) THEN
          CALL vecP2mapP(P2GDWSTO,P2TEMP)
          IOUT=1
        ENDIF
      CASE (8)  !!! LDAMOUT
        CVAR='damsto'
        IF( CMF_OPTIONS%LDAMOUT ) THEN
          CALL vecP2mapP(P2DAMSTO,P2TEMP)  !! P2DAMSTO only allocated for LDAMOUT
          IOUT=1
        ENDIF
      CASE (9)  !!! LLEVEE
        CVAR='levsto'
        IF( CMF_OPTIONS%LLEVEE ) THEN
          CALL vecP2mapP(P2LEVSTO,P2TEMP)  !! P2DAMSTO only allocated for LDAMOUT
          IOUT=1
        ENDIF
    END SELECT
  
#ifdef UseMPI_CMF
    CALL CMF_MPI_AllReduce_P2MAP(P2TEMP)
#endif
  
    IF( IOUT==1 )THEN
      IF( REGIONTHIS==1 )THEN
        STATUS = NF90_INQ_VARID(NCID,TRIM(CVAR),VARID)  !! check VARID is defined above, write restart only when STATUS=0
        IF ( STATUS .EQ. 0 ) THEN
          CALL NCERROR( NF90_PUT_VAR(NCID,VARID,P2TEMP,(/1,1,1/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,1/)) )
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  
  IF ( CMF_OPTIONS%LPTHOUT ) THEN
    IF ( .not. CMF_OPTIONS%LSTOONLY )THEN
      P1PTH(:,:)=D1PTHFLW_PRE(:,:)  !! convert Float32 to Float64 (for Single Precison Use)
#ifdef UseMPI_CMF
      CALL CMF_MPI_AllReduce_P1PTH(P1PTH)
#endif
      IF( REGIONTHIS==1 )THEN
        CALL NCERROR( NF90_INQ_VARID(NCID,'pthflw_pre',VARID))
        CALL NCERROR( NF90_PUT_VAR(NCID,VARID,P1PTH,(/1,1,1/),(/NPTHOUT,NPTHLEV,1/)) )
      ENDIF
    ENDIF
  ENDIF
  
  IF( REGIONTHIS==1 )THEN
    CALL NCERROR( NF90_SYNC(NCID) )
    CALL NCERROR( NF90_CLOSE(NCID) )
  ENDIF
  
  WRITE(CMF_FILES%LOGNAM,*) 'WRTE_REST: WRITE RESTART NETCDF:',CFILE
  
#endif
  END SUBROUTINE WRTE_REST_CDF
  !==========================================================
  
  END SUBROUTINE CMF_RESTART_WRITE
  !####################################################################



  !+
  !####################################################################
  SUBROUTINE CMF_TRACER_OUTPUT_WRITE
  !======
  USE CMF_UTILS_MOD,      ONLY: vecD2mapR
  
  ! save results to output files
  ! -- Called either from "MAIN/Coupler" or CMF_DRV_ADVANCE
  USE CMF_VARS_MOD,        ONLY: REGIONTHIS,NSEQMAX,ITRACE,P2TRCSTO,D2TRCOUT_oAVG,D2TRCDNS_oAVG,D2TRCPOUT_oAVG
  USE CMF_VARS_MOD,       ONLY: JYYYYMMDD, JHHMM, JHOUR, JMIN
#ifdef UseMPI_CMF
  USE CMF_CTRL_MPI_MOD,   ONLY: CMF_MPI_AllReduce_R2MAP
#endif
  
  IMPLICIT NONE
  INTEGER(KIND=JPIM)          :: JF
  !*** LOCAL
  REAL(KIND=JPRM)             :: R2OUT(CMF_CONFIG%NX,CMF_CONFIG%NY)
  REAL(KIND=JPRM)             :: R2COPY(NSEQMAX,1)
  REAL(KIND=JPRB)             :: D2COPY(NSEQMAX,1)        !! Dammy Array for Float64/32 switch
  !================================================
  !*** 0. check date:hour with output frequency
  IF ( MOD(JHOUR,CMF_OUTPUT%IFRQ_OUT)==0 .and. JMIN==0 ) THEN             ! JHOUR: end of time step , NFPPH: output frequency (hour)
  
    CALL CMF_TRACER_DIAG_GETAVE
  
    !*** 1. update IREC & calc average variable
    IRECOUT=IRECOUT+1 
    WRITE(CMF_FILES%LOGNAM,*) 'CMF::TRACER_OUTPUT_WRITE: write at time: ', JYYYYMMDD, JHHMM, IRECOUT
  
    !*** 2. check variable name & allocate data to pointer DVEC
    JF=0
    DO ITRACE=1, CMF_TRACER%NTRACE
      !! storage ======  
      JF=JF+1
      D2COPY(:,1)=P2TRCSTO(:,ITRACE) !! convert Double to Single precision when using SinglePrecisionMode 
  
      IF(  CMF_TRACER%LOUTVEC )THEN
        R2COPY(:,1)=D2COPY(:,1)
        WRITE(VAROUT(JF)%BINID,REC=IRECOUT) R2COPY         !! 1D vector (optional)
      ELSE
        !! convert 1Dvector to 2Dmap
        CALL vecD2mapR(D2COPY,R2OUT)             !! MPI node data is gathered by vecP2mapR
#ifdef UseMPI_CMF
        CALL CMF_MPI_AllReduce_R2MAP(R2OUT)
#endif
        IF ( REGIONTHIS==1 ) WRITE(VAROUT(JF)%BINID,REC=IRECOUT) R2OUT         !! 2D map
      ENDIF
  
      !! flux =======
      JF=JF+1
      IF( CMF_TRACER%LOUTVEC )THEN
        R2COPY(:,1)=D2TRCOUT_oAVG(:,ITRACE)
        WRITE(VAROUT(JF)%BINID,REC=IRECOUT) R2COPY         !! 1D vector (optional)
      ELSE
        !! convert 1Dvector to 2Dmap
        CALL vecD2mapR(D2TRCOUT_oAVG(:,ITRACE),R2OUT)             !! MPI node data is gathered by vecP2mapR
#ifdef UseMPI_CMF
        CALL CMF_MPI_AllReduce_R2MAP(R2OUT)
#endif
        IF ( REGIONTHIS==1 ) WRITE(VAROUT(JF)%BINID,REC=IRECOUT) R2OUT         !! 2D map
      ENDIF
  
      !! flux =======
      JF=JF+1
      IF( CMF_TRACER%LOUTVEC )THEN
        R2COPY(:,1)=D2TRCDNS_oAVG(:,ITRACE)
        WRITE(VAROUT(JF)%BINID,REC=IRECOUT) R2COPY         !! 1D vector (optional)
      ELSE
        !! convert 1Dvector to 2Dmap
        CALL vecD2mapR(D2TRCDNS_oAVG(:,ITRACE),R2OUT)             !! MPI node data is gathered by vecP2mapR
#ifdef UseMPI_CMF
        CALL CMF_MPI_AllReduce_R2MAP(R2OUT)
#endif
        IF ( REGIONTHIS==1 ) WRITE(VAROUT(JF)%BINID,REC=IRECOUT) R2OUT         !! 2D map
      ENDIF
  
      !! bifurcation net outflow =======
      IF(  CMF_TRACER%LTRCBIF )THEN 
        JF=JF+1
        IF( CMF_TRACER%LOUTVEC )THEN
          R2COPY(:,1)=D2TRCPOUT_oAVG(:,ITRACE)
          WRITE(VAROUT(JF)%BINID,REC=IRECOUT) R2COPY         !! 1D vector (optional)
        ELSE
          !! convert 1Dvector to 2Dmap
          CALL vecD2mapR(D2TRCPOUT_oAVG(:,ITRACE),R2OUT)             !! MPI node data is gathered by vecP2mapR
#ifdef UseMPI_CMF
          CALL CMF_MPI_AllReduce_R2MAP(R2OUT)
#endif
          IF ( REGIONTHIS==1 ) WRITE(VAROUT(JF)%BINID,REC=IRECOUT) R2OUT         !! 2D map
        ENDIF
      ENDIF 
    END DO
  
    WRITE(CMF_FILES%LOGNAM,*) 'CMF::TRACER_OUTPUT_WRITE: end'
  
  ENDIF
  
  CALL CMF_TRACER_DIAG_RESET
  
  END SUBROUTINE CMF_TRACER_OUTPUT_WRITE

  !####################################################################
SUBROUTINE CMF_TRACER_DIAG_GETAVE
  USE CMF_VARS_MOD,       ONLY: NADD_out,D2TRCOUT_oAVG,D2TRCDNS_oAVG,D2TRCPOUT_oAVG
  IMPLICIT NONE
  !====================
  D2TRCOUT_oAVG(:,:)  = D2TRCOUT_oAVG(:,:)  / DBLE(NADD_out)
  D2TRCDNS_oAVG(:,:)  = D2TRCDNS_oAVG(:,:)  / DBLE(NADD_out)
  D2TRCPOUT_oAVG(:,:) = D2TRCPOUT_oAVG(:,:) / DBLE(NADD_out)
  END SUBROUTINE CMF_TRACER_DIAG_GETAVE
  !####################################################################

  
  ! @@@@@@ Tracer Diagnose (for output)
!####################################################################
SUBROUTINE CMF_TRACER_DIAG_RESET
  USE CMF_VARS_MOD,       ONLY: JYYYYMMDD, JHHMM,NADD_out,D2TRCOUT_oAVG,D2TRCDNS_oAVG,D2TRCPOUT_oAVG
  IMPLICIT NONE
  !================================================
  WRITE(CMF_FILES%LOGNAM,*) "CMF::DIAG_AVERAGE: reset", JYYYYMMDD, JHHMM
  NADD_out=0
  D2TRCOUT_oAVG(:,:) = 0._JPRB
  D2TRCDNS_oAVG(:,:) = 0._JPRB
  D2TRCPOUT_oAVG(:,:)= 0._JPRB
  END SUBROUTINE CMF_TRACER_DIAG_RESET
  !####################################################################



  !+
!####################################################################
SUBROUTINE CMF_TRACER_RESTART_WRITE
  ! write restart files
  USE CMF_VARS_MOD,       ONLY: KSTEP,  NSTEPS, JYYYYMMDD, JHHMM, JDD, JHOUR, JMIN
  USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID, vecP2mapP, vecP2mapR
  USE CMF_VARS_MOD,       ONLY: REGIONTHIS,NSEQMAX
  USE CMF_VARS_MOD,       ONLY: ITRACE,P2TRCSTO

  IMPLICIT NONE
  !* local variable
  INTEGER(KIND=JPIM)         :: IREST
  CHARACTER(LEN=256)         :: CFILE,CDATE
  REAL(KIND=JPRD)            :: P2VAR(NSEQMAX,1)  !! use Real*8 for code simplicity
  REAL(KIND=JPRM)            :: R2TEMP(CMF_CONFIG%NX,CMF_CONFIG%NY)
  REAL(KIND=JPRD)            :: P2TEMP(CMF_CONFIG%NX,CMF_CONFIG%NY)
  !================================================
  IREST=0
  
  IF ( CMF_RESTART%IFRQ_RST>=0 .and. KSTEP==NSTEPS )THEN         !! end of run
    IREST=1
  ENDIF
  
  IF ( CMF_RESTART%IFRQ_RST>=1 .and. CMF_RESTART%IFRQ_RST<=24 )THEN
    IF ( MOD(JHOUR,CMF_RESTART%IFRQ_RST)==0 .and. JMIN==0 )THEN  !! at selected hour
      IREST=1
    ENDIF
  ENDIF
  
  IF ( CMF_RESTART%IFRQ_RST==30 )THEN
    IF ( JDD==1 .and. JHOUR==0 .and. JMIN==0 )THEN  !! at start of month
      IREST=1
    ENDIF
  ENDIF
  
  
  IF( IREST==1 )THEN
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
    WRITE(CMF_FILES%LOGNAM,*) 'CMF::TRACER_RESTART_WRITE: write time: ' , JYYYYMMDD, JHHMM
  
    !*** set file nam
    WRITE(CDATE,'(I8.8,I2.2)') JYYYYMMDD,JHOUR
    CFILE=TRIM(CMF_RESTART%CRESTDIR)//TRIM(CMF_TRACER%CVNRSTTRC)//TRIM(CDATE)//TRIM(CMF_PARAMS%CSUFBIN)
    WRITE(CMF_FILES%LOGNAM,*) 'TRACER_WRTE_REST_BIN: restart file:',CFILE
  
    !*** write restart data (2D map)
    CMF_FILES%TMPNAM=INQUIRE_FID()
  
    IF( CMF_RESTART%LRESTDBL )THEN
      IF ( REGIONTHIS==1 ) OPEN(CMF_FILES%TMPNAM,FILE=CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8*CMF_CONFIG%NX*CMF_CONFIG%NY)
    ELSE
      IF ( REGIONTHIS==1 ) OPEN(CMF_FILES%TMPNAM,FILE=CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
    ENDIF
  
    DO ITRACE=1, CMF_TRACER%NTRACE
      !! Double Precision Restart
      P2VAR(:,1)=P2TRCSTO(:,ITRACE)
      IF( CMF_RESTART%LRESTDBL )THEN
        CALL vecP2mapP(P2VAR,P2TEMP)  
#ifdef UseMPI_CMF
        CALL CMF_MPI_AllReduce_P2MAP(P2TEMP)
#endif
        IF ( REGIONTHIS==1 ) WRITE(CMF_FILES%TMPNAM,REC=ITRACE) P2TEMP
  
      !! Single Precision Restart
      ELSE
        CALL vecP2mapR(P2VAR,R2TEMP)  
#ifdef UseMPI_CMF
        CALL CMF_MPI_AllReduce_R2MAP(R2TEMP)
#endif
        IF ( REGIONTHIS==1 ) WRITE(CMF_FILES%TMPNAM,REC=ITRACE) R2TEMP
      ENDIF
  
    END DO
  ENDIF
  
  END SUBROUTINE CMF_TRACER_RESTART_WRITE
  !####################################################################
  

  !==========================================================
subroutine cmf_sed_output
  use CMF_UTILS_MOD,           only: vecD2mapR
  use CMF_VARS_MOD,             only: NSEQMAX,REGIONTHIS
  use CMF_VARS_MOD,            only: JHOUR, JMIN
  use CMF_VARS_MOD,             only: d2layer, d2sedcon, d2seddep, d2bedout_avg, d2netflw_avg, &
                                     d2sedout_avg, d2sedinp_avg, d2sedv_avg, sadd_out
  USE CMF_VARS_MOD,            only: VAROUT_SED,  NVARSOUT_SED
#ifdef UseMPI_CMF
  use CMF_CTRL_MPI_MOD,        only: CMF_MPI_AllReduce_R2MAP
#endif
  
  implicit none
  save
  integer(kind=JPIM)              :: ilyr, ised
  integer(kind=JPIM)              :: jf
  real(kind=JPRB),pointer         :: d2vec(:,:) ! point data location to output
  !*** local
  real(kind=JPRM)                 :: r3out(CMF_CONFIG%NX,CMF_CONFIG%NY,CMF_SED%NSED)
  !================================================
  call sediment_restart_write
  d2sedv_avg(:,:,:) = d2sedv_avg(:,:,:) / dble(sadd_out)
  write(CMF_FILES%LOGNAM,*) 'cmf_sed_output: average ',sadd_out,' seconds'
  !*** 0. check date:hour with output frequency
  if ( mod(JHOUR,CMF_OUTPUT%IFRQ_OUT)==0 .and. JMIN==0 ) then             ! JHOUR: end of time step , nfpph: output frequency (hour)

    !*** 1. calc average variable
    write(CMF_FILES%LOGNAM,*) 'cmf::sediment_output_write: write irec: ', IRECOUT

    !*** 2. check variable name & allocate data to pointer dvec
    do jf=1,NVARSOUT_SED
      select case (VAROUT_SED(jf)%cvname)
        case ('sedout')
          d2vec => d2sedout_avg
        case ('sedcon')
          d2vec => d2sedcon
        case ('sedinp')
          d2vec => d2sedinp_avg
        case ('bedout')
          d2vec => d2bedout_avg
        case ('netflw')
          d2vec => d2netflw_avg
        case ('layer')
          d2vec => d2layer
        case default
          if ( VAROUT_SED(jf)%cvname(:6) == 'deplyr' ) then
            read(VAROUT_SED(jf)%cvname(7:8),*) ilyr
            d2vec => d2seddep(:,ilyr,:)
          else
            write(CMF_FILES%LOGNAM,*) VAROUT_SED(jf)%cvname, ' not defined in cmf_output_mod'
          endif
      end select   !! variable name select

      !! convert 1dvector to 3dmap
      r3out(:,:,:) = CMF_PARAMS%RMIS
      
      if ( .not. CMF_OUTPUT%LOUTVEC ) then
        do ised = 1, CMF_SED%NSED
          call vecD2mapR(d2vec(:,ised),r3out(:,:,ised))             !! mpi node data is gathered by vec2map
#ifdef UseMPI_CMF
          call CMF_MPI_AllReduce_R2MAP(r3out(:,:,ised))
#endif
        enddo
        
        if ( REGIONTHIS==1 ) then
          if ( CMF_OUTPUT%LOUTCDF ) then
            call wrte_outcdf
          else
            call wrte_outbin(VAROUT_SED(jf)%binid,IRECOUT,r3out)
          endif
        endif
      else 
        call wrte_outvec(VAROUT_SED(jf)%binid,IRECOUT,d2vec)
      endif
    end do

    write(CMF_FILES%LOGNAM,*) 'cmf::sediment_output_write: end'
  endif

  d2sedv_avg(:,:,:) = 0._JPRB
  sadd_out = 0._JPRB

contains
  subroutine wrte_outcdf
#ifdef UseCDF_CMF
    use NETCDF
    use CMF_VARS_MOD,            only: KMINSTART, KMINNEXT
    use CMF_UTILS_MOD,           only: NCERROR
    
    implicit none
    save
    real(kind=JPRB)                 :: xtime
    xtime = real( (KMINNEXT-KMINSTART), JPRB) *60._JPRB
    call NCERROR( nf90_put_var(VAROUT_SED(jf)%ncid,VAROUT_SED(jf)%timid,xtime,(/VAROUT_SED(jf)%irecnc/)) )

    call NCERROR( nf90_put_var(VAROUT_SED(jf)%ncid,VAROUT_SED(jf)%varid,real(r3out(1:CMF_CONFIG%NX,1:CMF_CONFIG%NY,1:CMF_SED%NSED),JPRB),&
                  (/1,1,1,VAROUT_SED(jf)%irecnc/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,CMF_SED%NSED,1/)) )
    
    ! update irec
    VAROUT_SED(jf)%irecnc=VAROUT_SED(jf)%irecnc+1
    CALL NCERROR( NF90_SYNC(VAROUT_SED(jf)%ncid) )
#endif
  end subroutine wrte_outcdf



  !==========================================================
  subroutine wrte_outbin(ifn,irec,r2outdat)
    
    implicit none
    !*** input
    save
    integer(kind=JPIM),intent(in)   :: ifn                 !! file number
    integer(kind=JPIM),intent(in)   :: irec                !! record
    real(kind=JPRM)                 :: r2outdat(CMF_CONFIG%NX,CMF_CONFIG%NY,CMF_SED%NSED)
    !================================================
    write(ifn,rec=irec) r2outdat
  end subroutine wrte_outbin
  !==========================================================
  subroutine wrte_outvec(ifn,irec,d2outdat)
    
    implicit none
    !*** input
    save
    integer(kind=JPIM),intent(in)   :: ifn                 !! file number
    integer(kind=JPIM),intent(in)   :: irec                !! record
    real(kind=JPRB),intent(in)      :: d2outdat(NSEQMAX,CMF_SED%NSED) !! output data
    !*** local
    real(kind=JPRM)                 :: r2outdat(NSEQMAX,CMF_SED%NSED)
    !================================================
    r2outdat(:,:)=real(d2outdat(:,:))
    write(ifn,rec=irec) r2outdat
  end subroutine wrte_outvec
  !==========================================================
end subroutine cmf_sed_output
!==========================================================
!+
!==========================================================
!+
!==========================================================
subroutine sediment_restart_write
  use CMF_VARS_MOD,            only: KSTEP, NSTEPS, JDD, JHHMM, JHOUR, JMIN, JYYYYMMDD,REGIONTHIS
  use CMF_VARS_MOD,            only: d2layer, d2sedcon, d2seddep
  use CMF_UTILS_MOD,           only: vecD2mapR,INQUIRE_FID
#ifdef UseMPI_CMF
  use CMF_CTRL_MPI_MOD,        only: CMF_MPI_AllReduce_R2MAP
#endif
  
  implicit none
  save
  integer(kind=JPIM)              :: irec, irest, ised, tmpnam
  real(kind=JPRM)                 :: r3final(CMF_CONFIG%NX,CMF_CONFIG%NY,CMF_SED%NSED), r2temp(CMF_CONFIG%NX,CMF_CONFIG%NY)
  character(len=256)              :: cdate, cfile

  irest = 0

  if ( CMF_RESTART%IFRQ_RST>=0 .and. KSTEP==NSTEPS ) then  !! end of run
    irest = 1
  endif

  if ( CMF_SED%ifrq_rst_sed>=1 .and. CMF_SED%ifrq_rst_sed<=24 ) then  !! at selected hour
    if ( mod(JHOUR,CMF_SED%ifrq_rst_sed)==0 .and. JMIN==0 ) then
      irest = 1
    endif
  endif

  if ( CMF_SED%ifrq_rst_sed==30 ) then  !! at end of month
    if ( JDD==1 .and. JHOUR==0 .and. JMIN==0 ) then
      irest = 1
    endif
  endif

  if ( irest==1 ) then
    write(CMF_FILES%LOGNAM,*) ""
    write(CMF_FILES%LOGNAM,*) "!---------------------!"
    write(CMF_FILES%LOGNAM,*) 'cmf::sediment_restart_write: write time: ' , JYYYYMMDD, JHHMM

    write(cdate,'(I8.8,I2.2)') JYYYYMMDD,JHOUR
    cfile=trim(CMF_RESTART%CRESTDIR)//TRIM(CMF_SED%sedrest_outpre)//TRIM(cdate)//TRIM(CMF_PARAMS%CSUFBIN)
    write(CMF_FILES%LOGNAM,*) 'wrte_rest_bin: restart file:',cfile

    !*** write restart data (2D map)
    if ( REGIONTHIS == 1 ) then
      tmpnam = INQUIRE_FID()
      open(TMPNAM,file=cfile,form='unformatted',access='direct',recl=4*CMF_CONFIG%NX*CMF_CONFIG%NY*CMF_SED%NSED)
    endif
    do irec = 1, 2
     r3final(:,:,:) = CMF_PARAMS%RMIS
     do ised = 1, CMF_SED%NSED
       select case(irec)
         case (1)
           call vecD2mapR(d2layer(:,ised),r2temp)
         case (2)
           call vecD2mapR(d2sedcon(:,ised),r2temp)
       end select
#ifdef UseMPI_CMF
       call CMF_MPI_AllReduce_R2MAP(r2temp)
#endif
       r3final(:,:,ised) = r2temp(:,:)
     enddo
     if ( REGIONTHIS == 1 ) write(tmpnam,rec=irec) r3final
    enddo

    do irec = 1, CMF_SED%totlyrnum
      r3final(:,:,:) = CMF_PARAMS%RMIS
      do ised = 1, CMF_SED%NSED
        call vecD2mapR(d2seddep(:,irec,ised),r2temp)
#ifdef UseMPI_CMF
        call CMF_MPI_AllReduce_R2MAP(r2temp)
#endif
        r3final(:,:,ised) = r2temp
      enddo
      if ( REGIONTHIS == 1 ) write(tmpnam,rec=irec+2) r3final
    enddo

    if ( REGIONTHIS == 1 ) close(tmpnam)

  endif
end subroutine sediment_restart_write
!####################################################################

END MODULE CMF_OUTPUT_MOD
