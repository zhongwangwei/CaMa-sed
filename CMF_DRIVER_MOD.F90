MODULE CMF_DRIVER_MOD
!==========================================================
!* PURPOSE: Advance CaMa-Flood time integration  
!
!* CONTAINS:
! -- CMF_DRV_ADVANCE : Advance integration for KSPETS (given as argument)
!
!* INTERFACE:
! -- Called from "Main Program" or "Coupler"
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
USE PARKIND1,           ONLY: JPIM, JPRM, JPRB
USE CMF_NMLIST_MOD,     ONLY: CMF_FILES, CMF_OPTIONS, CMF_CONFIG, CMF_PARAMS
USE CMF_NMLIST_MOD,     ONLY: CMF_TIME, CMF_MAPS, CMF_FORCING, CMF_BOUNDARY
USE CMF_NMLIST_MOD,     ONLY: CMF_RESTART, CMF_DAM, CMF_LEVEE, CMF_OUTPUT,CMF_TRACER,CMF_SED

IMPLICIT NONE
REAL(KIND=JPRB)                 :: ZTT0, ZTT1, ZTT2   !! Time elapsed related 

CONTAINS
!####################################################################
! -- CMF_DRV_ADVANCE : Advance integration for KSPETS
!
!
!####################################################################
SUBROUTINE CMF_DRV_ADVANCE(KSTEPS)
USE CMF_VARS_MOD,            ONLY: KSTEP, JYYYYMMDD, JHHMM, JHOUR, JMIN
USE CMF_VARS_MOD,             ONLY: step_sed

USE CMF_TIME_MOD,       ONLY: CMF_TIME_NEXT, CMF_TIME_UPDATE
USE CMF_PHYSICS_MOD,    ONLY: CMF_PHYSICS_ADVANCE
USE CMF_OUTPUT_MOD,    ONLY: CMF_RESTART_WRITE
USE CMF_OUTPUT_MOD,     ONLY: CMF_OUTPUT_WRITE, CMF_OUTTXT_WRTE
USE CMF_OUTPUT_MOD,     ONLY: CMF_DAMOUT_WRTE

USE CMF_DIAG_MOD,       ONLY: CMF_DIAG_AVEMAX_OUTPUT, CMF_DIAG_GETAVE_OUTPUT, CMF_DIAG_RESET_OUTPUT
USE CMF_BOUNDARY_MOD,   ONLY: CMF_BOUNDARY_UPDATE
USE CMF_TRACER_MOD,     ONLY: CMF_TRACER_DENSITY, CMF_TRACER_FLUX
USE CMF_OUTPUT_MOD,     ONLY: CMF_TRACER_OUTPUT_WRITE, CMF_TRACER_RESTART_WRITE
USE CMF_OUTPUT_MOD,     ONLY: cmf_sed_output
USE CMF_SEDFLW_MOD,     ONLY: cmf_calc_sedflw
!$ USE OMP_LIB
IMPLICIT NONE 
SAVE
! Input argument 
INTEGER(KIND=JPIM)              :: KSTEPS             !! Number of timesteps to advance 
!* Local variables 
INTEGER(KIND=JPIM)              :: ISTEP              !! Time Step
!$ INTEGER(KIND=JPIM)           :: NTHREADS           !! OpenMP thread number
!==========================================================

!*** get OMP thread number
!$OMP PARALLEL
!$ NTHREADS=OMP_GET_MAX_THREADS()
!$OMP END PARALLEL 

!================================================
!*** START: time step loop
DO ISTEP=1,KSTEPS
  !============================
  !*** 0. get start CPU time
  CALL CPU_TIME(ZTT0)
  !$ ZTT0=OMP_GET_WTIME()

  !============================
  !*** 1. Set next time
  CALL CMF_TIME_NEXT               !! set KMINNEXT, JYYYYMMDD, JHHMM

  !*** (optional)
  IF( CMF_OPTIONS%LSEALEV )THEN
    CALL CMF_BOUNDARY_UPDATE
  ENDIF

  IF( CMF_OPTIONS%LTRACE )THEN
    CALL CMF_TRACER_DENSITY
  ENDIF

  !============================
  !*** 2. Advance model integration 
  CALL CMF_PHYSICS_ADVANCE

  !*** 2b.  Advance sediment model integration
  IF( CMF_OPTIONS%LSEDOUT .and. MOD(KSTEP,step_sed)==0 )THEN
    CALL cmf_calc_sedflw
  ENDIF
  IF( CMF_OPTIONS%LTRACE )THEN
    CALL CMF_TRACER_FLUX
  ENDIF

  CALL CMF_DIAG_AVEMAX_OUTPUT   !! average & maximum calculation for output

  CALL CPU_TIME(ZTT1)
  !$ ZTT1=OMP_GET_WTIME()

  !============================
  !*** 3. Write output file (when needed)
  IF( CMF_OPTIONS%LOUTPUT .and. MOD(JHOUR,CMF_OUTPUT%IFRQ_OUT)==0 .and. JMIN==0 )then
    !*** average variable
    CALL CMF_DIAG_GETAVE_OUTPUT !! average & maximum calculation for output: finalize

    !*** write output data
    CALL CMF_OUTPUT_WRITE

    IF ( CMF_OPTIONS%LSEDOUT ) THEN
      CALL cmf_sed_output
    ENDIF

    ! --- Optional: text file output
    IF ( CMF_OPTIONS%LDAMOUT ) THEN
    CALL CMF_OUTTXT_WRTE            !! reservoir operation
    CALL CMF_DAMOUT_WRTE            !! reservoir operation
    ENDIF
    IF( CMF_OPTIONS%LTRACE )THEN
      CALL CMF_TRACER_OUTPUT_WRITE
    ENDIF

    !*** reset variable
    CALL CMF_DIAG_RESET_OUTPUT  !! average & maximum calculation for output: reset
  ENDIF

  !============================ 
  !*** 4. Write restart file 
  CALL CMF_RESTART_WRITE
  IF( CMF_OPTIONS%LTRACE )THEN
    CALL CMF_TRACER_RESTART_WRITE
  ENDIF

  !============================ 
  !*** 5. Update current time      !! Update KMIN, IYYYYMMDD, IHHMM (to KMINNEXT, JYYYYMMDD, JHHMM)
  CALL CMF_TIME_UPDATE

  !============================
  !*** 6. Check CPU time 
  CALL CPU_TIME(ZTT2)
  !$ ZTT2=OMP_GET_WTIME()
  WRITE(CMF_FILES%LOGNAM,*) "CMF::DRV_ADVANCE END: KSTEP, time (end of Tstep):", KSTEP, JYYYYMMDD, JHHMM
  WRITE(CMF_FILES%LOGNAM,'(a,f8.1,a,f8.1,a)') "Elapsed cpu time", ZTT2-ZTT0,"Sec. // File output ", ZTT2-ZTT1, "Sec"

ENDDO
!*** END:time step loop
!================================================

END SUBROUTINE CMF_DRV_ADVANCE
!####################################################################

END MODULE CMF_DRIVER_MOD
