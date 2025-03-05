PROGRAM MAIN_cmf
  !=======================================================================
  !> @brief CaMa-Flood default stand-alone driver
  !>
  !> Main program for CaMa-Flood model execution in stand-alone mode
  !>
  !> @author D.Yamazaki & E. Dutra (U-Tokyo/FCUL)
  !> @date Aug 2019
  !>
  !> Licensed under the Apache License, Version 2.0 (the "License");
  !> modified by Zhongwang Wei (SYSU/CHINA)
  !> @date 2025-01-05
  !=======================================================================
  
  ! Core modules
  USE PARKIND1,        ONLY: JPRB, JPRM, JPIM 
  USE CMF_NMLIST_MOD
  USE CMF_INIT_MOD
  USE CMF_UTILS_MOD
  USE CMF_VARS_MOD
  
  ! Physics and forcing modules
  USE CMF_FORCING_MOD
  USE CMF_PHYSICS_MOD
  USE CMF_DRIVER_MOD
  USE CMF_END_MOD
  IMPLICIT NONE

  ! Local variables
  INTEGER(KIND=JPIM)          :: istep              ! total time step
  INTEGER(KIND=JPIM)          :: istepadv           ! time step to be advanced
  REAL(KIND=JPRB), &
  ALLOCATABLE            :: zbuff(:,:,:)       ! Buffer for forcing runoff

  !=====================================================================
  ! 1. Initialization
  !=====================================================================
  
  ! Read namelist and initialize model
  CALL CMF_READ_NMLIST()
  CALL CMF_MAKE_INIT(0)
  ! Allocate data buffer for input forcing
  ALLOCATE(zbuff(CMF_CONFIG%NXIN, CMF_CONFIG%NYIN, 2))
  !=====================================================================
  ! 2. Main temporal loop
  !=====================================================================
  
  istepadv = INT(CMF_CONFIG%DTIN/CMF_CONFIG%DT, JPIM)

  time_loop: DO istep = 1, NSTEPS, istepadv
    ! Read and process forcing data
    CALL CMF_FORCING_GET(zbuff)
    CALL CMF_FORCING_PUT(zbuff)
    
    ! Optional tracer calculations
    IF (CMF_OPTIONS%LTRACE) THEN
      CALL CMF_TRACER_FORC_GET()
      CALL CMF_TRACER_FORC_INTERP()
    END IF

    ! Optional sediment transport
    IF (CMF_OPTIONS%LSEDOUT) THEN
      CALL CMF_SED_FORC_GET()
    END IF

    ! run the routing model
    CALL CMF_DRV_ADVANCE(istepadv)
    
  END DO time_loop

  ! Clean up
  IF (ALLOCATED(zbuff)) DEALLOCATE(zbuff)
  ! end the routing model
  CALL CMF_END

  !*** 3b. MPI specific finalization
#ifdef UseMPI_CMF
   CALL CMF_MPI_END
#endif


END PROGRAM MAIN_cmf
!####################################################################
