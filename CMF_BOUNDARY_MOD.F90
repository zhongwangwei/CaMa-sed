MODULE CMF_BOUNDARY_MOD
!==========================================================
!* PURPOSE: Manage CaMa-Flood sea level boundary
!
!* CONTAINS:
! -- CMF_BOUNDARY_NMLIST : Read setting from namelist
! -- CMF_BOUNDARY_INIT   : Initialize boundary data file
! -- CMF_BOUNDARY_UPDATE : Update sea level boundary from file
! -- CMF_BOUNDARY_END    : Finalize   boundary data file
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
USE PARKIND1,                ONLY: JPIM, JPRB, JPRM
USE CMF_NMLIST_MOD,          ONLY: CMF_FILES, CMF_OPTIONS, CMF_CONFIG, CMF_PARAMS
USE CMF_NMLIST_MOD,          ONLY: CMF_TIME, CMF_MAPS, CMF_FORCING, CMF_BOUNDARY
USE CMF_NMLIST_MOD,          ONLY: CMF_RESTART, CMF_DAM, CMF_LEVEE, CMF_OUTPUT,CMF_TRACER,CMF_SED
!============================
IMPLICIT NONE
SAVE
  

CONTAINS


!####################################################################
SUBROUTINE CMF_BOUNDARY_UPDATE
! read runoff from file
USE CMF_VARS_MOD,            ONLY: IMIN, IYYYYMMDD, IHHMM
USE CMF_VARS_MOD,             ONLY: D2DWNELV, D2ELEVTN, D2SEALEV, D2MEANSL
IMPLICIT NONE
!* local variable
INTEGER(KIND=JPIM)              :: IUPDATE
!================================================
IUPDATE=0
IF( MOD( INT(IMIN),CMF_BOUNDARY%IFRQ_SL)==0 )THEN
  IUPDATE=1
ENDIF


IF( CMF_OPTIONS%LSEALEV .and. IUPDATE==1 ) THEN
  WRITE(CMF_FILES%LOGNAM,*) "CMF::BOUNDARY_UPDATE: update at time: ", IYYYYMMDD, IHHMM


  IF( CMF_BOUNDARY%LSEALEVCDF )THEN
#ifdef UseCDF_CMF
    CALL CMF_BOUNDARY_GET_CDF
#endif
  ELSE
    CALL CMF_BOUNDARY_GET_BIN
  ENDIF
ENDIF

IF( CMF_OPTIONS%LMEANSL ) THEN
  D2DWNELV(:,:)=D2ELEVTN(:,:) + D2MEANSL(:,:)
ELSE
  D2DWNELV(:,:)=D2ELEVTN(:,:)
ENDIF

IF( CMF_OPTIONS%LSEALEV ) THEN
  D2DWNELV(:,:)=D2DWNELV(:,:) + D2SEALEV(:,:)
ENDIF

CONTAINS
!==========================================================
!+ CMF_BOUNDARY_GET_BIN
!+ CMF_BOUNDARY_GET_CDF
!==========================================================
SUBROUTINE CMF_BOUNDARY_GET_BIN
USE CMF_VARS_MOD,            ONLY: IYYYY, IMM, IDD, IHHMM
USE CMF_VARS_MOD,             ONLY: D2SEALEV
USE CMF_UTILS_MOD,           ONLY: mapR2vecD,INQUIRE_FID
IMPLICIT NONE
CHARACTER(LEN=256)              :: CIFNAME             !! INPUT FILE
CHARACTER(LEN=256)              :: CDATE               !!
REAL(KIND=JPRM)                 :: R2TMP(CMF_CONFIG%NX,CMF_CONFIG%NY)
!====================
!*** 1. set file name
WRITE(CDATE,'(I4.4,I2.2,I2.2,I4.4)') IYYYY,IMM,IDD,IHHMM
CIFNAME = TRIM(CMF_BOUNDARY%CSEALEVDIR)//'/'//TRIM(CMF_BOUNDARY%CSEALEVPRE)//TRIM(CDATE)//TRIM(CMF_BOUNDARY%CSEALEVSUF)
WRITE(CMF_FILES%LOGNAM,*) "CMF::BOUNDARY_GET_BIN: read sealev:",TRIM(CIFNAME)

!*** open & read sea level
CMF_FILES%TMPNAM=INQUIRE_FID()
OPEN(CMF_FILES%TMPNAM,FILE=CIFNAME,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
READ(CMF_FILES%TMPNAM,REC=1) R2TMP
CLOSE(CMF_FILES%TMPNAM)
CALL mapR2vecD(R2TMP,D2SEALEV)

END SUBROUTINE CMF_BOUNDARY_GET_BIN
!==========================================================
!+
!+
!+
!==========================================================
#ifdef UseCDF_CMF
SUBROUTINE CMF_BOUNDARY_GET_CDF
USE CMF_VARS_MOD,            ONLY: KMIN
USE CMF_VARS_MOD,             ONLY: D2SEALEV, SLCDF, R1SLIN, I2SLMAP
USE CMF_UTILS_MOD,           ONLY: NCERROR, mapR2vecD
USE NETCDF
IMPLICIT NONE
!* Local variables
INTEGER(KIND=JPIM)              :: IRECSL, IX, IY, IS, ILNK
REAL(KIND=JPRM)                 :: R2TMP(CMF_CONFIG%NX,CMF_CONFIG%NY)
!===============
!*** 1. calculate irec
IRECSL = ( KMIN-SLCDF%NSTART )*60_JPIM / INT(CMF_BOUNDARY%DTSL,JPIM) + 1     !! (second from netcdf start time) / (input time step)
WRITE(CMF_FILES%LOGNAM,*) "CMF::BOUNDARY_GET_CDF:", TRIM(SLCDF%CNAME), IRECSL

!*** 2. read sea level
CALL NCERROR( NF90_GET_VAR(SLCDF%NCID,SLCDF%NVARID,R1SLIN,(/1,IRECSL/),(/CMF_BOUNDARY%NCDFSTAT,1/)),'READING SEA LEVEL' )

!*** 3. convert 1D station data -> 2D map
R2TMP(:,:)=0.E0
DO ILNK = 1, CMF_BOUNDARY%NLINKS
    IX = I2SLMAP(1,ILNK)
    IY = I2SLMAP(2,ILNK)
    IS = I2SLMAP(3,ILNK)
    R2TMP(IX,IY) = R1SLIN(IS)
END DO
CALL mapR2vecD(R2TMP,D2SEALEV)

END SUBROUTINE CMF_BOUNDARY_GET_CDF
#endif
!==========================================================

END SUBROUTINE CMF_BOUNDARY_UPDATE
!####################################################################







END MODULE CMF_BOUNDARY_MOD