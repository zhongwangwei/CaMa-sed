MODULE CMF_TIME_MOD
!==========================================================
!* PURPOSE: Manage time-related variables in CaMa-Flood
!
!* CONTAINS:
! -- CMF_TIME_NMLIST : Read setting from namelist
! -- CMF_TIME_INIT   : Initialize    time-related variables
! -- CMF_TIME_NEXT   : Set next-step time-related variables
! -- CMF_TIME_UPDATE : Update        time-related variables
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
USE PARKIND1,                ONLY: JPIM, JPRB, JPRM
USE CMF_NMLIST_MOD,          ONLY: CMF_FILES, CMF_OPTIONS, CMF_CONFIG, CMF_PARAMS
USE CMF_NMLIST_MOD,          ONLY: CMF_TIME, CMF_MAPS, CMF_FORCING, CMF_BOUNDARY
USE CMF_NMLIST_MOD,          ONLY: CMF_RESTART, CMF_DAM, CMF_LEVEE, CMF_OUTPUT,CMF_TRACER,CMF_SED 
USE CMF_UTILS_MOD,           ONLY: MIN2DATE, DATE2MIN, SPLITDATE, SPLITHOUR
!================================================
IMPLICIT NONE


CONTAINS 

!####################################################################
SUBROUTINE CMF_TIME_NEXT
! update time-related valiable
! -- Called from CMF_DRV_ADVANCE
!================================================
USE CMF_VARS_MOD,       ONLY: KSTEP, KMIN, KMINNEXT
USE CMF_VARS_MOD,       ONLY: IYYYYMMDD, IHHMM                               !! date:hour at start of time step
USE CMF_VARS_MOD,       ONLY: JYYYYMMDD, JYYYY, JMM, JDD, JHHMM, JHOUR, JMIN !! date:hour at end   of time step
IMPLICIT NONE
!================================================
WRITE(CMF_FILES%LOGNAM,*) ""
!*** 1. Advance KMIN, KSTEP
KSTEP=KSTEP+1
KMINNEXT=KMIN+INT(CMF_CONFIG%DT/60,JPIM)

WRITE(CMF_FILES%LOGNAM,*) "CMF::TIME_NEXT: ", KSTEP, KMIN, KMINNEXT, CMF_CONFIG%DT

!*** 2. Update J-time
CALL MIN2DATE(KMINNEXT,JYYYYMMDD,JHHMM)
CALL SPLITDATE(JYYYYMMDD,JYYYY,JMM,JDD)
CALL SPLITHOUR(JHHMM,JHOUR,JMIN)

WRITE(CMF_FILES%LOGNAM,*) "Strt of Tstep: KMIN,     IYYYYMMDD, IHHMM", KMIN,     IYYYYMMDD, IHHMM
WRITE(CMF_FILES%LOGNAM,*) "End  of Tstep: KMINNEXT, JYYYYMMDD, JHHMM", KMINNEXT, JYYYYMMDD, JHHMM


END SUBROUTINE CMF_TIME_NEXT
!####################################################################





!####################################################################
SUBROUTINE CMF_TIME_UPDATE
! update time-related valiable
! -- Called from CMF_DRV_ADVANCE
!================================================
USE CMF_VARS_MOD,       ONLY: KMIN, KMINNEXT
USE CMF_VARS_MOD,       ONLY: IYYYYMMDD, IYYYY, IMM, IDD, IHHMM, IHOUR, IMIN !! date:hour at start of time step
USE CMF_VARS_MOD,       ONLY: JYYYYMMDD, JYYYY, JMM, JDD, JHHMM, JHOUR, JMIN !! date:hour at end   of time step
IMPLICIT NONE
!================================================
WRITE(CMF_FILES%LOGNAM,*) ""
WRITE(CMF_FILES%LOGNAM,*) "CMF_TIME_UPDATE:"
!*** 1. Advance KMIN, KSTEP
KMIN=KMINNEXT

!*** 2. Update I-time
IYYYYMMDD=JYYYYMMDD
IYYYY=JYYYY
IMM  =JMM
IDD  =JDD
IHHMM=JHHMM
IHOUR=JHOUR
IMIN =JMIN

WRITE(CMF_FILES%LOGNAM,*) "Current time update: KMIN, IYYYYMMDD, IHHMM", KMIN, IYYYYMMDD, IHHMM

END SUBROUTINE CMF_TIME_UPDATE
!####################################################################

END MODULE CMF_TIME_MOD
