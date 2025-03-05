MODULE CMF_FORCING_MOD
!==========================================================
!* PURPOSE: Manage CaMa-Flood forcing
!
!* CONTAINS:
! -- CMF_FORCING_NMLIST : Read setting from Namelist
! -- CMF_FORCING_INIT   : Initialize forcing data file
! -- CMF_FORCING_PUT    : Put forcing data (PBUFF) to CaMa-Flood
! -- CMF_FORCING_GET    : Read forcing data from file (save as "data buffer" PBUFF)
! -- CMF_FORCING_END    : Finalize forcing data file
!
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
!
! Modifications: I. Ayan-Miguez (BSC) Apr 2023: Read inpmat.nc matrix by layers and added LECMF2LAKEC switch
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
USE PARKIND1,           ONLY: JPIM, JPRB, JPRM
USE CMF_NMLIST_MOD,     ONLY: CMF_FILES, CMF_OPTIONS, CMF_CONFIG, CMF_PARAMS
USE CMF_NMLIST_MOD,     ONLY: CMF_TIME, CMF_MAPS, CMF_FORCING, CMF_BOUNDARY
USE CMF_NMLIST_MOD,     ONLY: CMF_RESTART, CMF_DAM, CMF_LEVEE, CMF_OUTPUT,CMF_TRACER,CMF_SED
USE CMF_VARS_MOD,       ONLY: INPX, INPY, INPA, INPXI, INPYI, INPAI, INPNI
#ifdef UseCDF_CMF
USE CMF_VARS_MOD,       ONLY: ROFCDF
#endif
!============================
IMPLICIT NONE
SAVE
!*** input
REAL(KIND=JPRB),ALLOCATABLE     :: TBUFF(:,:,:)       ! Buffer to store forcing tracer

CONTAINS



!####################################################################
SUBROUTINE CMF_FORCING_GET(PBUFF)
USE CMF_UTILS_MOD,           ONLY: CMF_CheckNanB  !! check Udefined value
! read runoff from file
IMPLICIT NONE
REAL(KIND=JPRB),INTENT(INOUT)   :: PBUFF(:,:,:)

INTEGER(KIND=JPIM),SAVE         ::  IXIN, IYIN  !! FOR OUTPUT
!$OMP THREADPRIVATE                (IXIN)
!================================================
IF( CMF_FORCING%LINPCDF ) THEN
#ifdef UseCDF_CMF
  CALL CMF_FORCING_GET_CDF(PBUFF(:,:,:))
#endif
ELSE
  CALL CMF_FORCING_GET_BIN(PBUFF(:,:,:))
ENDIF 

!$OMP PARALLEL DO
DO IYIN=1,CMF_CONFIG%NYIN
  DO IXIN=1,CMF_CONFIG%NXIN
    IF( CMF_CheckNanB(PBUFF(IXIN,IYIN,1),0._JPRB) )THEN !! Check if PRUFINN(IX,IY) is NaN (Not-A-Number) ot not
      PBUFF(IXIN,IYIN,1)=CMF_PARAMS%RMIS 
    ENDIF
    PBUFF(IXIN,IYIN,1)=max(PBUFF(IXIN,IYIN,1),0._JPRB)    !! negative Runoff not assumed
  ENDDO
ENDDO
!$OMP END PARALLEL DO

IF ( CMF_OPTIONS%LROSPLIT ) THEN
!$OMP PARALLEL DO
  DO IYIN=1,CMF_CONFIG%NYIN
    DO IXIN=1,CMF_CONFIG%NXIN
      IF( CMF_CheckNanB(PBUFF(IXIN,IYIN,2),0._JPRB) )THEN !! Check if PRUFINN(IX,IY) is NaN (Not-A-Number) ot not
        PBUFF(IXIN,IYIN,2)=CMF_PARAMS%RMIS
      ENDIF
      PBUFF(IXIN,IYIN,2)=max(PBUFF(IXIN,IYIN,2),0._JPRB)    !! negative Runoff not assumed
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDIF


CONTAINS
!==========================================================
!+ CMF_FORCING_GET_BIN
!+ CMF_FORCING_GET_CDF
!==========================================================
SUBROUTINE CMF_FORCING_GET_BIN(PBUFF)
USE CMF_VARS_MOD,            ONLY: IYYYY, IMM, IDD, IHOUR, IMIN
USE CMF_UTILS_MOD,           ONLY: CONV_END,INQUIRE_FID
IMPLICIT NONE
REAL(KIND=JPRB),INTENT(OUT)     :: PBUFF(:,:,:)
!* Local variables
INTEGER(KIND=JPIM)              :: IRECINP
INTEGER(KIND=JPIM)              :: ISEC
CHARACTER(LEN=256)              :: CIFNAME             !! INPUT FILE
CHARACTER(LEN=256)              :: CDATE               !!
REAL(KIND=JPRM)                 :: R2TMP(CMF_CONFIG%NXIN,CMF_CONFIG%NYIN)
!================================================
!*** 1. calculate IREC for sub-daily runoff
ISEC    = IHOUR*60*60+IMIN*60   !! current second in a day
IRECINP = int( ISEC/CMF_CONFIG%DTIN ) +1   !! runoff irec (sub-dairy runoff)

!*** 2. set file name
WRITE(CDATE,'(I4.4,I2.2,I2.2)') IYYYY,IMM,IDD
CIFNAME=TRIM(CMF_FORCING%CROFDIR )//'/'//TRIM(CMF_FORCING%CROFPRE)//TRIM(CDATE)//TRIM(CMF_FORCING%CROFSUF)
WRITE(CMF_FILES%LOGNAM,*) "CMF::FORCING_GET_BIN:",TRIM(CIFNAME)

!*** 3. open & read runoff
CMF_FILES%TMPNAM=INQUIRE_FID()
OPEN(CMF_FILES%TMPNAM,FILE=CIFNAME,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NXIN*CMF_CONFIG%NYIN)
READ(CMF_FILES%TMPNAM,REC=IRECINP) R2TMP
CLOSE(CMF_FILES%TMPNAM)
WRITE(CMF_FILES%LOGNAM,*) "IRECINP:", IRECINP

!*** 4. copy runoff to PBUSS, endian conversion is needed
IF( CMF_FORCING%LINPEND ) CALL CONV_END(R2TMP,CMF_CONFIG%NXIN,CMF_CONFIG%NYIN)
PBUFF(:,:,1)=R2TMP(:,:)

!*** for sub-surface runoff withe LROSPLIT
PBUFF(:,:,2)=0._JPRB  !! Plain Binary subsurface runoff to be added later
IF ( CMF_OPTIONS%LROSPLIT ) THEN
  CIFNAME=TRIM(CMF_FORCING%CSUBDIR)//'/'//TRIM(CMF_FORCING%CSUBPRE)//TRIM(CDATE)//TRIM(CMF_FORCING%CSUBSUF)
  WRITE(CMF_FILES%LOGNAM,*) "CMF::FORCING_GET_BIN: (sub-surface)",TRIM(CIFNAME)

  CMF_FILES%TMPNAM=INQUIRE_FID()
  OPEN(CMF_FILES%TMPNAM,FILE=CIFNAME,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NXIN*CMF_CONFIG%NYIN)
  READ(CMF_FILES%TMPNAM,REC=IRECINP) R2TMP
  CLOSE(CMF_FILES%TMPNAM)
  WRITE(CMF_FILES%LOGNAM,*) "IRECINP:", IRECINP

  IF( CMF_FORCING%LINPEND ) CALL CONV_END(R2TMP,CMF_CONFIG%NXIN,CMF_CONFIG%NYIN)
  PBUFF(:,:,2)=R2TMP(:,:)
ENDIF

END SUBROUTINE CMF_FORCING_GET_BIN
! ================================================
!+
!+
!+
! ================================================
#ifdef UseCDF_CMF
SUBROUTINE CMF_FORCING_GET_CDF(PBUFF)
! Read forcing data from netcdf
! -- call from CMF_FORCING_GET
USE CMF_VARS_MOD,            ONLY: KMIN, IYYYYMMDD, IHHMM
USE CMF_UTILS_MOD,           ONLY: NCERROR
USE NETCDF
IMPLICIT NONE
!* Declaration of arguments 
REAL(KIND=JPRB),INTENT(OUT)     :: PBUFF(:,:,:)
!* Local variables
INTEGER(KIND=JPIM)              :: IRECINP
! ================================================
!*** 1. calculate irec
IRECINP=INT( (KMIN-ROFCDF%NSTART)*60_JPIM,JPIM ) / INT(CMF_CONFIG%DTIN,JPIM) + 1     !! (second from netcdf start time) / (input time step)
!*** 2. read runoff
CALL NCERROR( NF90_GET_VAR(ROFCDF%NCID,ROFCDF%NVARID(1),PBUFF(:,:,1),(/1,1,IRECINP/),(/CMF_CONFIG%NXIN,CMF_CONFIG%NYIN,1/)),'READING RUNOFF 1 ' )
IF ( ROFCDF%NVARID(2) .NE. -1 ) THEN
  CALL NCERROR( NF90_GET_VAR(ROFCDF%NCID,ROFCDF%NVARID(2),PBUFF(:,:,2),(/1,1,IRECINP/),(/CMF_CONFIG%NXIN,CMF_CONFIG%NYIN,1/)),'READING RUNOFF 2' )
ENDIF
WRITE(CMF_FILES%LOGNAM,*) "CMF::FORCING_GET_CDF: read runoff:",IYYYYMMDD,IHHMM,IRECINP
  
END SUBROUTINE CMF_FORCING_GET_CDF
#endif
! ================================================

END SUBROUTINE CMF_FORCING_GET
!####################################################################


!####################################################################
SUBROUTINE CMF_FORCING_COM(PBUFF)
! interporlate with inpmatI (CaMa grid -> input runoff grid), then send calling Model 
! -- called from "Main Program / Coupler" or CMF_DRV_ADVANCE
USE CMF_UTILS_MOD,           ONLY: vecD2mapD
#ifdef UseMPI_CMF
USE CMF_MPI_MOD,              ONLY: CMF_MPI_AllReduce_D2MAP
#endif
USE CMF_VARS_MOD,             ONLY: D2FLDFRC
IMPLICIT NONE 
! Declaration of arguments 
REAL(KIND=JPRB)                  :: D2MAPTMP(CMF_CONFIG%NX,CMF_CONFIG%NY)
REAL(KIND=JPRB), INTENT(OUT)     :: PBUFF(:,:,:)
!============================
CALL vecD2mapD(D2FLDFRC,D2MAPTMP)             !! MPI node data is gathered by vecP2mapR
#ifdef UseMPI_CMF
  CALL CMF_MPI_AllReduce_D2MAP(D2MAPTMP)
#endif

CALL INTERPI(D2MAPTMP,PBUFF(:,:,1))        !!  Inverse interpolation (CaMa grid -> input runoff grid)

CONTAINS
!==========================================================
!+ INTERPI
!==========================================================
SUBROUTINE INTERPI(PBUFFIN,PBUFFOUT)
! interporlate field using "input matrix inverse: from catchment to other grid"
IMPLICIT NONE
REAL(KIND=JPRB),INTENT(IN)      :: PBUFFIN(:,:)     !! CaMa-Flood variable on catchment (NX*NY)
REAL(KIND=JPRB),INTENT(OUT)     :: PBUFFOUT(:,:)    !! output on target grid = input runoff grid (NXIN * NYIN)

INTEGER(KIND=JPIM)              :: IX,IY,INP,IXIN,IYIN
! ========================================================
IF ( INPNI == -1 ) THEN
  WRITE(CMF_FILES%LOGNAM,*) "INPNI==-1, no inverse interpolation possible"
  STOP 9
ENDIF
PBUFFOUT(:,:)=1._JPRB 

DO IYIN=1, CMF_CONFIG%NYIN
  DO IXIN=1,CMF_CONFIG%NXIN
    PBUFFOUT(IXIN,IYIN)=0._JPRB
    DO INP=1,INPNI
      IX=INPXI(IXIN,IYIN,INP)
      IY=INPYI(IXIN,IYIN,INP)
      IF ( IX > 0 .AND. IY > 0 .AND. IX <= CMF_CONFIG%NX .AND. IY <= CMF_CONFIG%NY ) THEN
        PBUFFOUT(IXIN,IYIN) = PBUFFOUT(IXIN,IYIN) + PBUFFIN(IX,IY) * INPAI(IXIN,IYIN,INP)
      ENDIF
    ENDDO
  ENDDO
ENDDO



END SUBROUTINE INTERPI
!==========================================================

END SUBROUTINE CMF_FORCING_COM
!####################################################################


!####################################################################
SUBROUTINE CMF_FORCING_PUT(PBUFF)
! interporlate with inpmat, then send runoff data to CaMa-Flood 
! -- called from "Main Program / Coupler" or CMF_DRV_ADVANCE
USE CMF_VARS_MOD,            ONLY: D2RUNOFF,D2ROFSUB,D2WEVAP
IMPLICIT NONE 
! Declaration of arguments 
REAL(KIND=JPRB), INTENT(IN)     :: PBUFF(:,:,:)
!============================
! Runoff interpolation & unit conversion (mm/dt -> m3/sec)
IF (CMF_FORCING%LINTERP ) THEN ! mass conservation using "input matrix table (inpmat)"
  CALL ROFF_INTERP(PBUFF(:,:,1),D2RUNOFF)
  IF (CMF_OPTIONS%LROSPLIT) THEN
    CALL ROFF_INTERP(PBUFF(:,:,2),D2ROFSUB)
  ELSE
    D2ROFSUB(:,:) = 0._JPRB
  ENDIF
ELSE !  nearest point
  CALL CONV_RESOL(PBUFF(:,:,1),D2RUNOFF)
  IF (CMF_OPTIONS%LROSPLIT) THEN
    CALL CONV_RESOL(PBUFF(:,:,2),D2ROFSUB)
  ELSE
    D2ROFSUB(:,:) = 0._JPRB
  ENDIF
ENDIF 

IF (CMF_OPTIONS%LWEVAP) THEN
  IF ( SIZE(PBUFF,3) == 3 ) THEN
    CALL ROFF_INTERP(PBUFF(:,:,3),D2WEVAP)
  ELSE
    WRITE(CMF_FILES%LOGNAM,*)  "LWEVAP is true but evaporation not provide in input array for interpolation"
    WRITE(CMF_FILES%LOGNAM,*)  "CMF_FORCING_PUT(PBUFF), PBUFF should have 3 fields for interpolation "
    STOP 9
  ENDIF
ENDIF

CONTAINS
!==========================================================
!+ ROFF_INTERP : runoff interpolation with mass conservation using "input matrix table (inpmat)"
!+ CONV_RESOL : nearest point runoff interpolation
!==========================================================
SUBROUTINE ROFF_INTERP(PBUFFIN,PBUFFOUT)
! interporlate runoff using "input matrix"
USE CMF_VARS_MOD,             ONLY: NSEQALL,INPX,INPY,INPA
IMPLICIT NONE
REAL(KIND=JPRB),INTENT(IN)      :: PBUFFIN(:,:)     !! default [mm/dt] 
REAL(KIND=JPRB),INTENT(OUT)     :: PBUFFOUT(:,:)    !! m3/s
! SAVE for OMP
INTEGER(KIND=JPIM),SAVE  ::  ISEQ
INTEGER(KIND=JPIM),SAVE  ::  IXIN, IYIN, INPI  !! FOR OUTPUT
!$OMP THREADPRIVATE    (IXIN, IYIN, INPI)
!============================
!$OMP PARALLEL DO
DO ISEQ=1, NSEQALL
  PBUFFOUT(ISEQ,1)=0._JPRB
  DO INPI=1, CMF_CONFIG%INPN
    IXIN=INPX(ISEQ,INPI)
    IYIN=INPY(ISEQ,INPI)
    IF( IXIN>0 )THEN
      IF( IXIN > CMF_CONFIG%NXIN .OR. IYIN > CMF_CONFIG%NYIN ) THEN
        WRITE(CMF_FILES%LOGNAM,*)  "error"
        WRITE(CMF_FILES%LOGNAM,*)  'XXX',ISEQ,INPI,IXIN,IYIN
        CYCLE
      ENDIF

      IF( PBUFFIN(IXIN,IYIN).NE.CMF_PARAMS%RMIS )THEN
        PBUFFOUT(ISEQ,1) = PBUFFOUT(ISEQ,1) + PBUFFIN(IXIN,IYIN) * INPA(ISEQ,INPI) / CMF_FORCING%DROFUNIT   !! DTIN removed in v395
      ENDIF

    ENDIF
  END DO
END DO
!$OMP END PARALLEL DO
END SUBROUTINE ROFF_INTERP
!==========================================================
!+
!+
!==========================================================
SUBROUTINE CONV_RESOL(PBUFFIN,PBUFFOUT)
!! use runoff data without any interporlation. map resolution & runoff resolution should be same
USE CMF_VARS_MOD,             ONLY: NSEQALL, NSEQMAX, D2GRAREA
USE CMF_UTILS_MOD,           ONLY: mapD2vecD
IMPLICIT NONE

REAL(KIND=JPRB),INTENT(IN)      :: PBUFFIN(:,:)     !! default [mm/dt] 
REAL(KIND=JPRB),INTENT(OUT)     :: PBUFFOUT(:,:)    !! m3/s

REAL(KIND=JPRB),ALLOCATABLE     :: D2TEMP(:,:)

INTEGER(KIND=JPIM),SAVE         ::  ISEQ
! ================================================
ALLOCATE(D2TEMP(NSEQMAX,1))
CALL mapD2vecD(PBUFFIN,D2TEMP)
!$OMP PARALLEL DO
DO ISEQ=1, NSEQALL
  IF( D2TEMP(ISEQ,1).NE.CMF_PARAMS%RMIS )THEN
    PBUFFOUT(ISEQ,1) = D2TEMP(ISEQ,1) * D2GRAREA(ISEQ,1) / CMF_FORCING%DROFUNIT
    PBUFFOUT(ISEQ,1) = MAX(PBUFFOUT(ISEQ,1), 0._JPRB)
  ELSE
    PBUFFOUT(ISEQ,1)=0._JPRB
  ENDIF
END DO
!$OMP END PARALLEL DO
END SUBROUTINE CONV_RESOL
!==========================================================

END SUBROUTINE CMF_FORCING_PUT
!####################################################################


!@@@@@@ TRACER Forcing Input @@@@@@
!####################################################################
SUBROUTINE CMF_TRACER_FORC_GET
  ! tracer forcing data
  USE CMF_VARS_MOD,            ONLY: IYYYY, IMM, IDD, IHOUR, IMIN
  USE CMF_UTILS_MOD,           ONLY: CONV_END, INQUIRE_FID
  IMPLICIT NONE
  !* Local variables
  INTEGER(KIND=JPIM)              :: IRECINP
  INTEGER(KIND=JPIM)              :: ISEC
  CHARACTER(LEN=256)              :: CIFNAME             !! INPUT FILE
  CHARACTER(LEN=256)              :: CDATE               !!
  REAL(KIND=JPRM)                 :: R2TMP( CMF_CONFIG%NXIN,CMF_CONFIG%NYIN)
  INTEGER(KIND=JPIM)       :: ITRACE      !! tracer id

  !####################################################################
  !*** 1. calculate IREC for sub-daily runoff
  ISEC    = IHOUR*60*60+IMIN*60   !! current second in a day
  IRECINP = int( ISEC/CMF_TRACER%DTIN_TRC ) +1   !! runoff irec (sub-dairy runoff)
  
  DO ITRACE=1, CMF_TRACER%NTRACE
  
    !*** 2. set file name
    WRITE(CDATE,'(I4.4,I2.2,I2.2)') IYYYY,IMM,IDD
    CIFNAME=TRIM(CMF_TRACER%CTRCDIR)//'/'//TRIM(CMF_TRACER%CTRCPRE)//TRIM(CDATE)//TRIM(CMF_TRACER%CTRCSUF)
    WRITE(CMF_FILES%LOGNAM,*) "CMF::FORCING_GET_BIN:",TRIM(CIFNAME)
  
    !*** 3. open & read runoff
    CMF_FILES%TMPNAM=INQUIRE_FID()
    OPEN(CMF_FILES%TMPNAM,FILE=CIFNAME,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NXIN*CMF_CONFIG%NYIN)
    READ(CMF_FILES%TMPNAM,REC=IRECINP) R2TMP
    CLOSE(CMF_FILES%TMPNAM)
    WRITE(CMF_FILES%LOGNAM,*) "IRECINP:", IRECINP
  
    !*** 4. copy runoff to PBUSS, endian conversion is needed
    IF( CMF_TRACER%LINPEND ) CALL CONV_END(R2TMP,CMF_CONFIG%NXIN,CMF_CONFIG%NYIN)
    TBUFF(:,:,ITRACE)=R2TMP(:,:)
  
  END DO
  
  END SUBROUTINE CMF_TRACER_FORC_GET
  !####################################################################
  !+
  !+
  !+
  !####################################################################
  SUBROUTINE CMF_TRACER_FORC_INTERP
  ! interporlate with inpmat, then send runoff data to CaMa-Flood 
  ! -- called from "Main Program / Coupler" or CMF_DRV_ADVANCE
    USE CMF_VARS_MOD,            ONLY: NSEQALL, INPX, INPY, INPA,D2TRCINP
    USE CMF_UTILS_MOD,           ONLY: CMF_CheckNanB
  IMPLICIT NONE
  ! SAVE for OMP
  INTEGER(KIND=JPIM),SAVE  ::  ISEQ
  INTEGER(KIND=JPIM),SAVE  ::  IXIN, IYIN, INPI  !! FOR OUTPUT
  INTEGER(KIND=JPIM)       :: ITRACE      !! tracer id

  !$OMP THREADPRIVATE         (IXIN, IYIN, INPI)
  !============================
  
  DO ITRACE=1, CMF_TRACER%NTRACE
  
    !$OMP PARALLEL DO
      DO ISEQ=1, NSEQALL
        D2TRCINP(ISEQ,ITRACE)=0._JPRB
        DO INPI=1, CMF_CONFIG%INPN
          IXIN=INPX(ISEQ,INPI)
          IYIN=INPY(ISEQ,INPI)
          IF( IXIN>0 )THEN
            IF( IXIN > CMF_CONFIG%NXIN .OR. IYIN > CMF_CONFIG%NYIN ) THEN
              WRITE(CMF_FILES%LOGNAM,*)  "error"
              WRITE(CMF_FILES%LOGNAM,*)  'XXX',ISEQ,INPI,IXIN,IYIN
              CYCLE
            ENDIF
            IF( TBUFF(IXIN,IYIN,ITRACE).NE.CMF_PARAMS%RMIS )THEN
              D2TRCINP(ISEQ,ITRACE) = D2TRCINP(ISEQ,ITRACE) + TBUFF(IXIN,IYIN,ITRACE) * INPA(ISEQ,INPI) / CMF_TRACER%DTRCUNIT
              !! assume tracer input file unit is [ MASS/m2/s ]. If it is not "per second", DTRCUNIT should be changed. (default is 1)
            ENDIF
            IF( CMF_CheckNanB(TBUFF(IXIN,IYIN,ITRACE),0._JPRB) ) D2TRCINP(ISEQ,ITRACE)=0._JPRB  !! treat NaN runoff input 
          ENDIF
        END DO
        D2TRCINP(ISEQ,ITRACE)=MAX(D2TRCINP(ISEQ,ITRACE), 0._JPRB)
      END DO
    !$OMP END PARALLEL DO
  
  END DO
  
  END SUBROUTINE CMF_TRACER_FORC_INTERP


  !==========================================================
!+

!==========================================================
SUBROUTINE CMF_SED_FORC_GET
  ! read forcing from file
  USE CMF_VARS_MOD,            ONLY: IYYYY, IMM, IDD, IHOUR, IMIN
  USE CMF_VARS_MOD,            ONLY: NSEQMAX, NSEQALL
  USE CMF_VARS_MOD,            ONLY: D2SEDINP, D2SEDINP_AVG, D2SEDFRC
  USE CMF_VARS_MOD,            ONLY: D2FLDFRC, DSYLUNIT, PYLD, PYLDC, PYLDPC, D2SLOPE, D2GRAREA
  USE CMF_UTILS_MOD,           ONLY: CONV_END, INQUIRE_FID
  
  IMPLICIT NONE
  
  REAL(KIND=JPRB)                 :: D2TEMP(NSEQMAX,1)  ! Changed to 2D array
  !* Local variables
  INTEGER(KIND=JPIM)              :: IRECINP
  INTEGER(KIND=JPIM)              :: ISEC
  CHARACTER(LEN=256)              :: CIFNAME             !! INPUT FILE
  CHARACTER(LEN=256)              :: CDATE               !!
  REAL(KIND=JPRM)                 :: R2TMP(CMF_CONFIG%NXIN,CMF_CONFIG%NYIN)
  IF( CMF_FORCING%LINPCDF ) THEN
#ifdef UseCDF_CMF
    CALL CMF_SED_FORC_GET_CDF(D2TEMP(:,1))  ! Pass 1D slice
#endif
  ELSE
    CALL CMF_SED_FORC_GET_BIN(D2TEMP(:,1))  ! Pass 1D slice
  ENDIF

  ! Calculate sediment yield in rivers
  CALL CALC_SEDYLD(D2TEMP(:,1))  ! Pass 1D slice


CONTAINS

  SUBROUTINE CMF_SED_FORC_GET_BIN(PBUFFOUT)
    IMPLICIT NONE
    REAL(KIND=JPRB),INTENT(OUT)     :: PBUFFOUT(:)  ! Changed to 1D array
    
    !*** 1. calculate IREC for sub-daily precipitation
    ISEC    = IHOUR*60*60+IMIN*60   !! current second in a day
    IRECINP = INT( ISEC/CMF_CONFIG%DTIN) +1   !! precipitation irec (sub-daily precipitation)
  
    !*** 2. set file name
    WRITE(CDATE,'(I4.4,I2.2,I2.2)') IYYYY,IMM,IDD
    CIFNAME=TRIM(CMF_SED%SEDINPUT_DIR)//'/'//TRIM(CMF_SED%SEDINPUT_PRE)//TRIM(CDATE)//TRIM(CMF_SED%SEDINPUT_SUF)
    WRITE(CMF_FILES%LOGNAM,*) "CMF::SED_FORCING_GET_BIN:",TRIM(CIFNAME)
 
    !*** 3. open & read forcing data
    CMF_FILES%TMPNAM=INQUIRE_FID()
    
    OPEN(CMF_FILES%TMPNAM,FILE=CIFNAME,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NXIN*CMF_CONFIG%NYIN)
    READ(CMF_FILES%TMPNAM,REC=IRECINP) R2TMP
    CLOSE(CMF_FILES%TMPNAM)
  
    !*** 4. interpolate forcing grid to model grid
    CALL SEDINP_INTERP(R2TMP,PBUFFOUT)
  END SUBROUTINE CMF_SED_FORC_GET_BIN

#ifdef UseCDF_CMF
  SUBROUTINE CMF_SED_FORC_GET_CDF(PBUFFOUT)
    USE CMF_VARS_MOD,            ONLY: KMIN, SEDCDF
    USE CMF_UTILS_MOD,           ONLY: NCERROR
    USE NETCDF
    IMPLICIT NONE
    REAL(KIND=JPRB),INTENT(OUT)     :: PBUFFOUT(:)  ! Changed to 1D array
    !*** 1. calculate irec
    IRECINP=INT( (KMIN-ROFCDF%NSTART)*60_JPIM,JPIM ) / INT(CMF_CONFIG%DTIN,JPIM) + 1
    !*** 2. read sediment data
    IF (SEDCDF%NCID <= 0) THEN
      WRITE(CMF_FILES%LOGNAM,*) "ERROR: SEDCDF%NCID is not a valid NetCDF ID:", SEDCDF%NCID
      CALL NCERROR(NF90_EINVAL, 'Invalid NetCDF ID for sediment file')
    ENDIF
    
    CALL NCERROR( NF90_GET_VAR(SEDCDF%NCID,SEDCDF%NVARID(1),R2TMP,&
                  (/1,1,IRECINP/),(/CMF_CONFIG%NXIN,CMF_CONFIG%NYIN,1/)),&
                  'READING SEDIMENT' )
    WRITE(CMF_FILES%LOGNAM,*) "CMF::SED_FORCING_GET_CDF: read sediment data, IRECINP:",IRECINP
    
    !*** 3. interpolate forcing grid to model grid
    CALL SEDINP_INTERP(R2TMP,PBUFFOUT)
  END SUBROUTINE CMF_SED_FORC_GET_CDF
#endif


  subroutine SEDINP_INTERP(pbuffin,pbuffout)
    ! interporlate sediment forcing data using "input matrix"
      USE CMF_VARS_MOD,            only: NSEQALL, INPX, INPY, INPA, D2GRAREA
      implicit none
      real(kind=JPRM),intent(in)      :: pbuffin(:,:)     !! default for prcp[kg/m2/s]
      real(kind=JPRB),intent(out)     :: pbuffout(:)    !! kg/m2/s
      ! save for omp
      integer(kind=jpim),save  ::  iseq, ixin, iyin, inpi  !! for output
      !$omp threadprivate    (ixin, iyin)
      !============================
      !$omp parallel do
      do iseq=1, NSEQALL
        pbuffout(iseq)=0._JPRB
        do inpi=1, CMF_CONFIG%INPN
          ixin=INPX(iseq,inpi)
          iyin=INPY(iseq,inpi)
          if( ixin>0 )then
            if( ixin > CMF_CONFIG%NXIN .or. iyin > CMF_CONFIG%NYIN ) then
              write(CMF_FILES%LOGNAM,*)  "error"
              write(CMF_FILES%LOGNAM,*)  'xxx',iseq,inpi,ixin,iyin
              cycle
            endif
            if( pbuffin(ixin,iyin).ne.CMF_PARAMS%RMIS )then
              pbuffout(iseq) = pbuffout(iseq) + pbuffin(ixin,iyin) * INPA(iseq,inpi) / D2GRAREA(iseq,1)
            endif
          endif
        end do
        pbuffout(iseq)=max(pbuffout(iseq), 0._JPRB)
      end do
      !$omp end parallel do
    end subroutine sedinp_interp

  SUBROUTINE CALC_SEDYLD(PBUFFIN)
    IMPLICIT NONE
    REAL(KIND=JPRB), INTENT(IN)     :: PBUFFIN(:)
    REAL(KIND=JPRB)                 :: SBUFF(NSEQMAX)
    INTEGER(KIND=JPIM)              :: ISEQ
    CALL PRCP_CONVERT_SED(PBUFFIN, SBUFF)

    !$OMP PARALLEL DO 
    DO ISEQ = 1, NSEQALL
      D2SEDINP(ISEQ,:) = SBUFF(ISEQ) * D2SEDFRC(ISEQ,:)
      D2SEDINP_AVG(ISEQ,:) = D2SEDINP_AVG(ISEQ,:) + D2SEDINP(ISEQ,:) * CMF_CONFIG%DTIN
    ENDDO
    !$OMP END PARALLEL DO
  END SUBROUTINE CALC_SEDYLD

  SUBROUTINE PRCP_CONVERT_SED(PBUFFIN,PBUFFOUT)
    IMPLICIT NONE
    REAL(KIND=JPRB), INTENT(IN)     :: PBUFFIN(:)     !! kg/m2/s
    REAL(KIND=JPRB), INTENT(OUT)    :: PBUFFOUT(:)    !! m3/s
    INTEGER(KIND=JPIM)              :: I, ISEQ
    !$OMP PARALLEL DO
    DO ISEQ = 1, NSEQALL
      PBUFFOUT(ISEQ) = 0.D0
      IF ( PBUFFIN(ISEQ) * 86400.D0 <= 10.D0 ) CYCLE

      DO I = 1, CMF_CONFIG%NLFP
        IF ( D2FLDFRC(ISEQ,1) * CMF_CONFIG%NLFP > DBLE(I) ) CYCLE  ! no erosion if submerged
        PBUFFOUT(ISEQ) = PBUFFOUT(ISEQ) + PYLD * (PBUFFIN(ISEQ)*3600.D0)**PYLDPC * D2SLOPE(ISEQ,I)**PYLDC / 3600.D0 & 
          & * D2GRAREA(ISEQ,1) * MIN(DBLE(I)/DBLE(CMF_CONFIG%NLFP)-D2FLDFRC(ISEQ,1), 1.D0/DBLE(CMF_CONFIG%NLFP)) * DSYLUNIT
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO
    DO ISEQ = 1, NSEQALL
      IF ( PBUFFOUT(ISEQ) > 0.D0 ) THEN
        WRITE(CMF_FILES%LOGNAM,*) 'PBUFFOUT(',ISEQ,') = ', PBUFFOUT(ISEQ)
      ENDIF
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE PRCP_CONVERT_SED

END SUBROUTINE CMF_SED_FORC_GET


END MODULE CMF_FORCING_MOD
