MODULE CMF_END_MOD
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
USE PARKIND1,                ONLY: JPIM, JPRB, JPRM
USE CMF_NMLIST_MOD,     ONLY: CMF_FILES, CMF_OPTIONS, CMF_CONFIG, CMF_PARAMS
USE CMF_NMLIST_MOD,     ONLY: CMF_TIME, CMF_MAPS, CMF_FORCING, CMF_BOUNDARY
USE CMF_NMLIST_MOD,     ONLY: CMF_RESTART, CMF_DAM, CMF_LEVEE, CMF_OUTPUT,CMF_TRACER,CMF_SED
USE CMF_VARS_MOD,             ONLY: INPX, INPY, INPA, INPXI, INPYI, INPAI, INPNI
USE CMF_DRIVER_MOD,           ONLY: ZTT0, ZTT1, ZTT2
!============================
IMPLICIT NONE
SAVE

CONTAINS

!####################################################################
SUBROUTINE CMF_END

  !$ USE OMP_LIB    
  IMPLICIT NONE 
  !==========================================================
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!******************************!"
  WRITE(CMF_FILES%LOGNAM,*) "CMF::DRV_END: finalize forcing & output modules"
  CALL CMF_FORCING_END
  IF( CMF_OPTIONS%LOUTPUT )THEN
    CALL CMF_OUTPUT_END
    IF( CMF_OPTIONS%LTRACE ) CALL CMF_TRACER_END
    IF( CMF_OPTIONS%LSEDOUT ) CALL CMF_SED_END
  ENDIF
  IF( CMF_OPTIONS%LSEALEV ) THEN
    CALL CMF_BOUNDARY_END
  ENDIF
  
  !*** get simulation end time
  CALL CPU_TIME(ZTT2)
  !$ ZTT2=OMP_GET_WTIME()
  WRITE(CMF_FILES%LOGNAM,*) "CMF::DRV_END: simulation finished in:",ZTT2-ZTT0,' Seconds'
  
  WRITE(CMF_FILES%LOGNAM,*) "CMF::DRV_END: close logfile"
  WRITE(CMF_FILES%LOGNAM,*) "CMF::===== CALCULATION END ====="
  CLOSE(CMF_FILES%LOGNAM)
  
  END SUBROUTINE CMF_END
  !####################################################################


!####################################################################
  SUBROUTINE CMF_OUTPUT_END
    ! Finalize output module (close files)
    ! -- Called from CMF_DRV_END
#ifdef UseCDF_CMF
    USE NETCDF
    USE CMF_UTILS_MOD,           ONLY: NCERROR
#endif
    USE CMF_VARS_MOD,             ONLY: REGIONTHIS,VAROUT, NVARSOUT !
    IMPLICIT NONE
    ! Local variables
    INTEGER(KIND=JPIM)              :: JF
    !================================================
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
    WRITE(CMF_FILES%LOGNAM,*) "CMF::OUTPUT_END: finalize output module"
    
    IF( REGIONTHIS==1 )THEN
      IF (CMF_OUTPUT%LOUTCDF) THEN
#ifdef UseCDF_CMF
        DO JF=1,NVARSOUT 
          CALL NCERROR( NF90_CLOSE(VAROUT(JF)%NCID))
          WRITE(CMF_FILES%LOGNAM,*) "Output netcdf output unit closed:",VAROUT(JF)%NCID
        ENDDO
#endif
      ELSE !! binary output
        DO JF=1,NVARSOUT
          CLOSE(VAROUT(JF)%BINID)
          WRITE(CMF_FILES%LOGNAM,*) "Output binary output unit closed:",VAROUT(JF)%BINID
        ENDDO
        IF( CMF_OUTPUT%LOUTVEC )THEN
          CALL WRTE_mapR2vecD  !! write map-vector conversion file
        ENDIF
      ENDIF
    ENDIF
    
    WRITE(CMF_FILES%LOGNAM,*) "CMF::OUTPUT_END: end"
    
    
    CONTAINS
    !==========================================================
    !+ WRTE_mapR2vecD
    !+
    !+
    !==========================================================
    SUBROUTINE WRTE_mapR2vecD       !! 1D sequence vector informtion required to convert MPI distributed vector output to 2D map
    USE CMF_VARS_MOD,        ONLY: I1SEQX, I1SEQY, NSEQMAX
    USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
    IMPLICIT NONE
    !* local variable
    CHARACTER(LEN=256)         :: CFILE1
    !================================================
    IF( CMF_OUTPUT%LOUTVEC )THEN
      CFILE1='./ind_xy'//TRIM(CMF_PARAMS%CSUFVEC)
    
      WRITE(CMF_FILES%LOGNAM,*) "LOUTVEC: write mapR2vecD conversion table", TRIM(CFILE1)
    
      CMF_FILES%TMPNAM=INQUIRE_FID()
      OPEN(CMF_FILES%TMPNAM,FILE=CFILE1,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NSEQMAX)
      WRITE(CMF_FILES%TMPNAM,REC=1) I1SEQX
      WRITE(CMF_FILES%TMPNAM,REC=2) I1SEQY
      CLOSE(CMF_FILES%TMPNAM)
    ENDIF
    
    END SUBROUTINE WRTE_mapR2vecD
    !================================================
    
    
    END SUBROUTINE CMF_OUTPUT_END

!###################################################################
SUBROUTINE CMF_FORCING_END
#ifdef UseCDF_CMF
USE CMF_UTILS_MOD,         ONLY: NCERROR
USE CMF_VARS_MOD,          ONLY: ROFCDF
USE NETCDF
#endif
IMPLICIT NONE
!================================================
WRITE(CMF_FILES%LOGNAM,*) ""
WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
WRITE(CMF_FILES%LOGNAM,*) "CMF::FORCING_END: Finalize forcing module"

!* Close Input netcdf 
IF( CMF_FORCING%LINPCDF ) THEN
#ifdef UseCDF_CMF
  CALL NCERROR( NF90_CLOSE(ROFCDF%NCID))
  WRITE(CMF_FILES%LOGNAM,*) "input netCDF runoff closed:",ROFCDF%NCID
#endif
ENDIF 

WRITE(CMF_FILES%LOGNAM,*) "CMF::FORCING_END: end"

END SUBROUTINE CMF_FORCING_END
!####################################################################



!####################################################################
SUBROUTINE CMF_BOUNDARY_END
#ifdef UseCDF_CMF
  USE CMF_UTILS_MOD,         ONLY: NCERROR
  USE CMF_VARS_MOD,          ONLY: SLCDF
  USE NETCDF
#endif
  IMPLICIT NONE
  !================================================
  WRITE(CMF_FILES%LOGNAM,*) "CMF::BOUNDARY_END: Finalize boundary module"
  
  IF( CMF_BOUNDARY%LSEALEVCDF )THEN
#ifdef UseCDF_CMF
    CALL NCERROR( NF90_CLOSE(SLCDF%NCID))
    WRITE(CMF_FILES%LOGNAM,*) "Input netcdf sealev closed:",SLCDF%NCID
#endif
  ENDIF 
  
  WRITE(CMF_FILES%LOGNAM,*) "CMF::BOUNDARY_END: end"
  
  END SUBROUTINE CMF_BOUNDARY_END
  !####################################################################

!####################################################################
SUBROUTINE CMF_TRACER_END
  ! Finalize output module (close files)
  ! -- Called from CMF_DRV_END
  USE CMF_VARS_MOD,          ONLY: REGIONTHIS, NVARSOUT, VAROUT
  IMPLICIT NONE
  ! Local variables
  INTEGER(KIND=JPIM)              :: JF
  !================================================
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
  WRITE(CMF_FILES%LOGNAM,*) "CMF::TRACER_OUTPUT_END: finalize output module"
  
  IF( REGIONTHIS==1 )THEN
    DO JF=1,NVARSOUT
      CLOSE(VAROUT(JF)%BINID)
      WRITE(CMF_FILES%LOGNAM,*) "Output binary output unit closed:",VAROUT(JF)%BINID
    ENDDO
  ENDIF
  
  END SUBROUTINE CMF_TRACER_END
  !####################################################################

  !==========================================================
subroutine CMF_SED_END
#ifdef UseCDF_CMF
    use NETCDF
    use CMF_UTILS_MOD,           only: NCERROR
#endif
    use CMF_VARS_MOD,             only: REGIONTHIS,NVARSOUT,VAROUT
    implicit none
    save
    integer(kind=JPIM)              :: jf
    write(CMF_FILES%LOGNAM,*) ""
    write(CMF_FILES%LOGNAM,*) "!---------------------!"
    write(CMF_FILES%LOGNAM,*) "sediment_output_end: finalize output module"
  
    if ( CMF_OUTPUT%LOUTVEC) then
      do jf = 1, nvarsout
        close(varout(jf)%binid)
      enddo
    else if ( REGIONTHIS==1 ) then
      do jf = 1, nvarsout
        if ( CMF_OUTPUT%LOUTCDF ) then
#ifdef UseCDF_CMF
          call NCERROR( nf90_close(varout(jf)%ncid) )
#endif
        else
          close(varout(jf)%binid)
        endif
      enddo
    endif 
    
    write(CMF_FILES%LOGNAM,*) 'sediment_output_end: end'
  end subroutine CMF_SED_END
    
  
END MODULE CMF_END_MOD
