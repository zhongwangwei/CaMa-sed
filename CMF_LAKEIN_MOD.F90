MODULE CMF_LAKEIN_MOD
!==========================================================
!* PURPOSE: Lake-river coupling in ILS
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
USE PARKIND1,                ONLY: JPIM, JPRM, JPRB
USE CMF_NMLIST_MOD,        ONLY: CMF_FILES, CMF_OPTIONS, CMF_CONFIG, CMF_PARAMS
USE CMF_NMLIST_MOD,        ONLY: CMF_TIME, CMF_MAPS, CMF_FORCING, CMF_BOUNDARY
USE CMF_NMLIST_MOD,        ONLY: CMF_RESTART, CMF_DAM, CMF_LEVEE, CMF_OUTPUT,CMF_TRACER,CMF_SED 
USE CMF_VARS_MOD,            ONLY: D2RUNIN_AVG, D2LAKEFRC, D2RUNIN
CONTAINS
!####################################################################
! -- CMF_CALC_LAKEIN
! -- CMF_LAKEIN_AVE
! -- CMF_LAKEIN_AVERAGE
! -- CMF_RESET_LAKEIN
!####################################################################
SUBROUTINE CMF_CALC_LAKEIN
USE CMF_VARS_MOD,         ONLY: NSEQALL, D2GRAREA
USE CMF_VARS_MOD,         ONLY: P2RIVSTO, P2FLDSTO
USE CMF_VARS_MOD,         ONLY: D2STORGE
USE CMF_VARS_MOD,          ONLY: D2LAKEFRC, D2RUNIN
IMPLICIT NONE
!*** Local
INTEGER(KIND=JPIM)         :: ISEQ
REAL(KIND=JPRB)            :: DRIVRIN, DFLDRIN
!*** Lake Parameter
REAL(KIND=JPRB)            :: RINDMP
!================================================
RINDMP = 1.D0

DO ISEQ=1, NSEQALL
  DRIVRIN = P2RIVSTO(ISEQ,1) * D2LAKEFRC(ISEQ,1) / (RINDMP * 8.64D4) * CMF_CONFIG%DT !! m3
  DFLDRIN = P2FLDSTO(ISEQ,1) * D2LAKEFRC(ISEQ,1) / (RINDMP * 8.64D4) * CMF_CONFIG%DT !! m3
  P2RIVSTO(ISEQ,1) = P2RIVSTO(ISEQ,1) - DRIVRIN
  P2FLDSTO(ISEQ,1) = P2FLDSTO(ISEQ,1) - DFLDRIN
  D2RUNIN(ISEQ,1) = (DRIVRIN + DFLDRIN) / CMF_CONFIG%DT / D2GRAREA(ISEQ,1) * 1.D3 !! kg/m2/s

  D2STORGE(ISEQ,1)=P2RIVSTO(ISEQ,1)+P2FLDSTO(ISEQ,1)
ENDDO

END SUBROUTINE CMF_CALC_LAKEIN
!####################################################################






!####################################################################
SUBROUTINE CMF_LAKEIN_AVE
USE CMF_VARS_MOD,        ONLY: NSEQALL
USE CMF_VARS_MOD,        ONLY: D2RUNIN, D2RUNIN_AVG
IMPLICIT NONE
!*** Local
INTEGER(KIND=JPIM)         :: ISEQ
!================================================
DO ISEQ=1, NSEQALL
  D2RUNIN_AVG(ISEQ,1)=D2RUNIN_AVG(ISEQ,1)+D2RUNIN(ISEQ,1)*CMF_CONFIG%DT
END DO

END SUBROUTINE CMF_LAKEIN_AVE
!####################################################################





!####################################################################
SUBROUTINE CMF_LAKEIN_AVERAGE
USE CMF_VARS_MOD,        ONLY: NADD_adp 
USE CMF_VARS_MOD,        ONLY: NSEQALL
IMPLICIT NONE
!*** Local
INTEGER(KIND=JPIM)         :: ISEQ
!================================================
DO ISEQ=1, NSEQALL
  D2RUNIN_AVG(ISEQ,1)=D2RUNIN_AVG(ISEQ,1)/DBLE(NADD_adp)
END DO

END SUBROUTINE CMF_LAKEIN_AVERAGE
!####################################################################



!####################################################################
SUBROUTINE CMF_RESET_LAKEIN
IMPLICIT NONE
!================================================
D2RUNIN_AVG(:,:)=0._JPRB

END SUBROUTINE CMF_RESET_LAKEIN
!####################################################################

END MODULE CMF_LAKEIN_MOD
