MODULE CMF_INIT_MOD
  USE PARKIND1,                ONLY: JPIM, JPRB, JPRM, JPRD
  USE CMF_NMLIST_MOD,     ONLY: CMF_FILES, CMF_OPTIONS, CMF_CONFIG, CMF_PARAMS
  USE CMF_NMLIST_MOD,     ONLY: CMF_TIME, CMF_MAPS, CMF_FORCING, CMF_BOUNDARY
  USE CMF_NMLIST_MOD,     ONLY: CMF_RESTART, CMF_DAM, CMF_LEVEE, CMF_OUTPUT,CMF_TRACER,CMF_SED
  USE NETCDF

  IMPLICIT NONE
  !** local variables
  SAVE
  REAL(KIND=JPRB)                 :: ZTT0, ZTT1, ZTT2   ! Time elapsed related 
  INTEGER(KIND=JPIM)              :: IRECOUT            ! Output file irec

CONTAINS
!####################################################################
SUBROUTINE CMF_TIME_INIT
  ! initialize time-related valiable
  ! -- Called from CMF_DRV_INIT
  !================================================
  USE CMF_VARS_MOD,            ONLY: KSTEP, NSTEPS, KMIN, KMINNEXT,  KMINSTART, KMINEND
  USE CMF_VARS_MOD,            ONLY: ISYYYYMMDD,ISHHMM,ISYYYY,ISMM,ISDD,ISHOUR,ISMIN     !! start date:hour
  USE CMF_VARS_MOD,            ONLY: IEYYYYMMDD,IEHHMM,IEYYYY,IEMM,IEDD,IEHOUR,IEMIN     !! start date:hour
  USE CMF_VARS_MOD,            ONLY: IYYYYMMDD, IHHMM, IYYYY, IMM, IDD, IHOUR, IMIN !! date:hour at start of time step
  USE CMF_VARS_MOD,            ONLY: JYYYYMMDD, JHHMM, JYYYY, JMM, JDD, JHOUR, JMIN !! date:hour at end   of time step
  USE CMF_NMLIST_MOD,          ONLY: CMF_TIME  ! Add this to access namelist time variables
  USE CMF_UTILS_MOD,           ONLY: MIN2DATE, DATE2MIN, SPLITDATE, SPLITHOUR

  IMPLICIT NONE
  !================================================
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
  
  WRITE(CMF_FILES%LOGNAM,*) "CMF::TIME_INIT:  initialize time variables"

  !*** 1. Start time & End Time
  ISYYYYMMDD=CMF_TIME%SYEAR*10000+CMF_TIME%SMON*100+CMF_TIME%SDAY
  ISHHMM=CMF_TIME%SHOUR*100_JPIM
  ISYYYY=CMF_TIME%SYEAR
  ISMM  =CMF_TIME%SMON  
  ISDD  =CMF_TIME%SDAY
  ISHOUR=CMF_TIME%SHOUR
  ISMIN =0_JPIM
  
  IEYYYYMMDD=CMF_TIME%EYEAR*10000+CMF_TIME%EMON*100+CMF_TIME%EDAY  !! End   time
  IEHHMM=CMF_TIME%EHOUR*100_JPIM
  IEYYYY=CMF_TIME%EYEAR
  IEMM  =CMF_TIME%EMON
  IEDD  =CMF_TIME%EDAY
  IEHOUR=CMF_TIME%EHOUR
  IEMIN =0_JPIM
  
  WRITE(CMF_FILES%LOGNAM,*) 'Start Date:',ISYYYYMMDD, ISHHMM, KMINSTART
  WRITE(CMF_FILES%LOGNAM,*) 'End   Date:',IEYYYYMMDD, IEHHMM, KMINEND
  
  !*** 2. Initialize KMIN for START & END Time
  KMINSTART=DATE2MIN(ISYYYYMMDD,ISHHMM)
  KMINEND  =DATE2MIN(IEYYYYMMDD,IEHHMM)
  
  KMIN=KMINSTART
  
  !*** 3. Calculate NSTEPS: time steps within simulation time
  KSTEP=0 
  NSTEPS=int ( ( (KMINEND-KMINSTART)*60_JPIM ) / CMF_CONFIG%DT )      !!  (End - Start) / DT
  
  WRITE(CMF_FILES%LOGNAM,*) 'NSTEPS    :',NSTEPS
  
  !*** 4. Initial time step setting
  IYYYYMMDD=ISYYYYMMDD
  CALL SPLITDATE(IYYYYMMDD,IYYYY,IMM,IDD)
  IHHMM=ISHHMM
  CALL SPLITHOUR(IHHMM,IHOUR,IMIN)
  
  ! tentatively set KMINNEXT to KMIN (just within initialization phase)
  KMINNEXT =KMIN
  JYYYYMMDD=IYYYYMMDD
  JHHMM=IHHMM
  CALL SPLITDATE(JYYYYMMDD,JYYYY,JMM,JDD)
  CALL SPLITHOUR(JHHMM,JHOUR,JMIN)
  
  WRITE(CMF_FILES%LOGNAM,*) 'Initial Time Step Date:Hour :', IYYYYMMDD,'_',IHOUR,':',IMIN
  
  !*** end 
  WRITE(CMF_FILES%LOGNAM,*) "CMF::TIME_INIT: end"
  
END SUBROUTINE CMF_TIME_INIT
!####################################################################
  
!####################################################################
SUBROUTINE CMF_RIVMAP_INIT
  ! read & set river network map 
  ! -- call from CMF_DRV_INIT
  USE CMF_VARS_MOD,        ONLY: I2NEXTX,I2NEXTY, I2REGION, REGIONALL,REGIONTHIS, &
                              & I1SEQX, I1SEQY,  I1NEXT,  I2VECTOR, D1LON,    D1LAT,      &
                              & NSEQRIV,  NSEQALL,  NSEQMAX
  USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
  USE CMF_NMLIST_MOD,ONLY: CMF_FILES, CMF_MAPS, CMF_CONFIG
  IMPLICIT NONE
  !================================================
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
  WRITE(CMF_FILES%LOGNAM,*) 'CMF::RIVMAP_INIT: river network initialization'
  
  ! *** 1. ALLOCATE ARRAYS
  ALLOCATE( I2NEXTX(CMF_CONFIG%NX,CMF_CONFIG%NY) )
  ALLOCATE( I2NEXTY(CMF_CONFIG%NX,CMF_CONFIG%NY) )
  ALLOCATE( I2REGION(CMF_CONFIG%NX,CMF_CONFIG%NY) )
  ALLOCATE( D1LON(CMF_CONFIG%NX) )
  ALLOCATE( D1LAT(CMF_CONFIG%NY) )
  
  !============================
  !*** 2a. read river network map
  WRITE(CMF_FILES%LOGNAM,*) 'CMF::RIVMAP_INIT: read nextXY & set lat lon'
  IF( CMF_MAPS%LMAPCDF )THEN
#ifdef UseCDF_CMF
    CALL READ_MAP_CDF
#endif
  ELSE
    CALL READ_MAP_BIN
  ENDIF
  
  !*** 2b. calculate river sequence & regions
  WRITE(CMF_FILES%LOGNAM,*) 'CMF::RIVMAP_INIT: calc region'
  CALL CALC_REGION
  
  !============================
  !*** 3. conversion 2D map -> 1D vector
  WRITE(CMF_FILES%LOGNAM,*) 'CMF::RIVMAP_INIT: calculate 1d river sequence'
  
  CALL CALC_1D_SEQ                                  !! 2D map to 1D vector conversion. for faster calculation
  
  WRITE(CMF_FILES%LOGNAM,*) '  NSEQRIV=',NSEQRIV
  WRITE(CMF_FILES%LOGNAM,*) '  NSEQALL=',NSEQALL
  
  !*** 3c. Write Map Data                                       !! used for combining mpi distributed output into one map
  IF( REGIONTHIS==1 )THEN
    CMF_FILES%TMPNAM=INQUIRE_FID()
    OPEN(CMF_FILES%TMPNAM,FILE='./mapdata.txt',FORM='FORMATTED')
    WRITE(CMF_FILES%TMPNAM,*) 'NX',        CMF_CONFIG%NX
    WRITE(CMF_FILES%TMPNAM,*) 'NY',        CMF_CONFIG%NY
    WRITE(CMF_FILES%TMPNAM,*) 'NLFP',      CMF_CONFIG%NLFP
    WRITE(CMF_FILES%TMPNAM,*) 'REGIONALL', REGIONALL
    WRITE(CMF_FILES%TMPNAM,*) 'NSEQMAX',   NSEQMAX
    CLOSE(CMF_FILES%TMPNAM)
  ENDIF
  
  !============================
  !*** 4.  bifurcation channel parameters
  IF(CMF_OPTIONS%LPTHOUT)THEN
    WRITE(CMF_FILES%LOGNAM,*) 'CMF::RIVMAP_INIT: read bifurcation channel setting'
    CALL READ_BIFPARAM
  ENDIF
  
  DEALLOCATE( I2NEXTX,I2NEXTY,I2REGION )
  
  WRITE(CMF_FILES%LOGNAM,*) 'CMF::RIVMAP_INIT: end'
  
  CONTAINS
  !==========================================================
  !+ READ_MAP_BIN
  !+ READ_MAP_CDF
  !+ CALC_REGION
  !+ CALC_1D_SEQ
  !+ READ_BIFPRM
  !==========================================================
  SUBROUTINE READ_MAP_BIN
  USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID, CONV_ENDI
  IMPLICIT NONE
  !* local variables
  INTEGER(KIND=JPIM),SAVE    :: IX,IY
  !==========================================================
  !*** read river map
  WRITE(CMF_FILES%LOGNAM,*)'RIVMAP_INIT: nextxy binary: ',TRIM(CMF_MAPS%CNEXTXY)
  CMF_FILES%TMPNAM=INQUIRE_FID()
  OPEN(CMF_FILES%TMPNAM,FILE=CMF_MAPS%CNEXTXY,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
  READ(CMF_FILES%TMPNAM,REC=1) I2NEXTX
  READ(CMF_FILES%TMPNAM,REC=2) I2NEXTY
  CLOSE(CMF_FILES%TMPNAM)
  
  IF ( CMF_OUTPUT%LMAPEND )THEN
    CALL CONV_ENDI(I2NEXTX,CMF_CONFIG%NX,CMF_CONFIG%NY)
    CALL CONV_ENDI(I2NEXTY,CMF_CONFIG%NX,CMF_CONFIG%NY)
  ENDIF
  
  !*** calculate lat, lon
  IF( CMF_CONFIG%WEST>=-180._JPRB .and. CMF_CONFIG%EAST<=360._JPRB .and. CMF_CONFIG%SOUTH>=-180._JPRB .and. CMF_CONFIG%NORTH<=180._JPRB )THEN  !! bugfix_v396a
!$OMP PARALLEL DO
    DO IX=1,CMF_CONFIG%NX
      D1LON(IX)=CMF_CONFIG%WEST +(DBLE(IX)-0.5D0)*(CMF_CONFIG%EAST-CMF_CONFIG%WEST)  /DBLE(CMF_CONFIG%NX)
    ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
    DO IY=1,CMF_CONFIG%NY
      D1LAT(IY)=CMF_CONFIG%NORTH-(DBLE(IY)-0.5D0)*(CMF_CONFIG%NORTH-CMF_CONFIG%SOUTH)/DBLE(CMF_CONFIG%NY)
    ENDDO
!$OMP END PARALLEL DO
  ENDIF
  
  END SUBROUTINE READ_MAP_BIN
  !==========================================================
  !+
  !+
  !+
  !==========================================================
#ifdef UseCDF_CMF
  SUBROUTINE READ_MAP_CDF
  USE CMF_UTILS_MOD  ,ONLY: NCERROR
  USE NETCDF
  IMPLICIT NONE
  !* local variables
  INTEGER(KIND=JPIM)              :: NCID,VARID
  !================================================
  WRITE(CMF_FILES%LOGNAM,*)'RIVMAP_INIT: nextxy netCDF: ', TRIM(CMF_MAPS%CRIVCLINC)
  
  CALL NCERROR (NF90_OPEN(CMF_MAPS%CRIVCLINC,NF90_NOWRITE,NCID),'opening '//TRIM(CMF_MAPS%CRIVCLINC) )
  
  !*** next xy
  CALL NCERROR ( NF90_INQ_VARID(NCID,'nextx',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,I2NEXTX),'reading data' ) 
  
  CALL NCERROR ( NF90_INQ_VARID(NCID,'nexty',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,I2NEXTY),'reading data' )
  
  !*** lat, lon
  CALL NCERROR ( NF90_INQ_VARID(NCID,'lat',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,D1LAT),'reading data' )
  
  CALL NCERROR ( NF90_INQ_VARID(NCID,'lon',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,D1LON),'reading data' )
  
  CALL NCERROR( NF90_CLOSE(NCID))
  
  END SUBROUTINE READ_MAP_CDF
#endif
  !==========================================================
  !+
  !+
  !+
  !==========================================================
  SUBROUTINE CALC_REGION    !! evenly allocate pixels to mpi nodes (updated in v4.03. MPI region given from file)
#ifdef UseCDF_CMF
  USE CMF_UTILS_MOD,           ONLY: NCERROR
  USE NETCDF
#endif
  IMPLICIT NONE
  !* local variables
  INTEGER(KIND=JPIM),ALLOCATABLE  :: REGIONGRID(:)
  !
  INTEGER(KIND=JPIM),SAVE         :: IX,IY
  INTEGER(KIND=JPIM),SAVE         :: IREGION
#ifdef UseMPI_CMF
#ifdef UseCDF_CMF
  INTEGER(KIND=JPIM)              :: NCID
  INTEGER(KIND=JPIM)              :: VARID
#endif
#endif
  !$OMP THREADPRIVATE               (IX)
  !================================================
  WRITE(CMF_FILES%LOGNAM,*) 'RIVMAP_INIT: region code'
  
  !*** read MPI region map
  REGIONALL=1
  I2REGION(:,:)=CMF_PARAMS%IMIS
  !$OMP PARALLEL DO
  DO IY=1, CMF_CONFIG%NY
    DO IX=1, CMF_CONFIG%NX
      IF( I2NEXTX(IX,IY)/=CMF_PARAMS%IMIS ) THEN
        I2REGION(IX,IY)=1
      ENDIF
    END DO
  END DO
  !$OMP END PARALLEL DO
  
  !! Use MPI: read MPI region map, allocate regions to MPI nodes
#ifdef UseMPI_CMF
    IF ( CMF_MAPS%LMAPCDF ) THEN
#ifdef UseCDF_CMF
      CALL NCERROR (NF90_OPEN(CMF_MAPS%CMPIREGNC,NF90_NOWRITE,NCID),'opening '//TRIM(CMF_MAPS%CMPIREGNC) )
      CALL NCERROR (NF90_INQ_VARID(NCID, 'mpireg',VARID),'getting id' )
      CALL NCERROR (NF90_GET_VAR(NCID,VARID,I2REGION),'reading data' )
      CALL NCERROR (NF90_CLOSE(NCID))
#endif
    ELSE
      WRITE(CMF_FILES%LOGNAM,*)'RIVMAP_INIT: read MPI region: ',TRIM(CMF_MAPS%CMPIREG)
      CMF_FILES%TMPNAM=INQUIRE_FID()
      OPEN(CMF_FILES%TMPNAM,FILE=CMF_MAPS%CMPIREG,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
      READ(CMF_FILES%TMPNAM,REC=1) I2REGION
      CLOSE(CMF_FILES%TMPNAM)
    ENDIF
  
    REGIONALL=1
!$OMP PARALLEL DO REDUCTION(max:REGIONALL)
    DO IY=1, CMF_CONFIG%NY
      DO IX=1, CMF_CONFIG%NX
        REGIONALL=MAX( REGIONALL, I2REGION(IX,IY) )
      END DO
    END DO
!$OMP END PARALLEL DO
#endif
  
  
  WRITE(CMF_FILES%LOGNAM,*)'RIVMAP_INIT: count number of grid in each region: '
  ALLOCATE(REGIONGRID(REGIONALL))
  REGIONGRID(:)=0
  !! OMP reduction operation for array might not be available in some environment
  DO IY=1, CMF_CONFIG%NY
    DO IX=1, CMF_CONFIG%NX
      IF( I2REGION(IX,IY)>0 ) THEN
        IREGION=I2REGION(IX,IY)
        REGIONGRID(IREGION)=REGIONGRID(IREGION)+1
      ENDIF
    END DO
  END DO
  
  NSEQMAX=0
  DO IREGION=1, REGIONALL
    NSEQMAX=MAX(NSEQMAX,REGIONGRID(IREGION))  !! maximum nseqall among all MPI region
  END DO
  
  WRITE(CMF_FILES%LOGNAM,*) 'CALC_REGION: REGIONALL= ', REGIONALL
  WRITE(CMF_FILES%LOGNAM,*) 'CALC_REGION: NSEQMAX='   , NSEQMAX
  WRITE(CMF_FILES%LOGNAM,*) 'CALC_REGION: NSEQALL='   , NSEQALL
  
  END SUBROUTINE CALC_REGION
  !==========================================================
  !+
  !+
  !+
  !==========================================================
  SUBROUTINE CALC_1D_SEQ
  !* local variables
  INTEGER(KIND=JPIM)              :: IX,IY,JX,JY,ISEQ,JSEQ,ISEQ1,ISEQ2,AGAIN
  
  INTEGER(KIND=JPIM),ALLOCATABLE  :: NUPST(:,:), UPNOW(:,:)
  !================================================
  WRITE(CMF_FILES%LOGNAM,*) 'RIVMAP_INIT: convert 2D map to 1D sequence'
  
  ALLOCATE( NUPST(CMF_CONFIG%NX,CMF_CONFIG%NY) )
  ALLOCATE( UPNOW(CMF_CONFIG%NX,CMF_CONFIG%NY) )
  
  ALLOCATE( I1SEQX(NSEQMAX) )
  ALLOCATE( I1SEQY(NSEQMAX) )
  ALLOCATE( I1NEXT(NSEQMAX) )
  ALLOCATE( I2VECTOR(CMF_CONFIG%NX,CMF_CONFIG%NY) )
  I1SEQX(:)=0
  I1SEQY(:)=0
  I1NEXT(:)=0
  I2VECTOR(:,:)=0
  
  ! count number of upstream 
  NUPST(:,:)=0
  UPNOW(:,:)=0
  DO IY=1, CMF_CONFIG%NY
    DO IX=1, CMF_CONFIG%NX
      IF( I2NEXTX(IX,IY).GT.0 .and. I2REGION(IX,IY)==REGIONTHIS )THEN
        JX=I2NEXTX(IX,IY)
        JY=I2NEXTY(IX,IY)
        NUPST(JX,JY)=NUPST(JX,JY)+1
      ENDIF
    END DO
  END DO
  
  ! register upmost grid in 1d sequence
  ISEQ=0
  DO IY=1, CMF_CONFIG%NY
    DO IX=1, CMF_CONFIG%NX
      IF( I2NEXTX(IX,IY).GT.0 .and. I2REGION(IX,IY)==REGIONTHIS )THEN
        IF( NUPST(IX,IY)==UPNOW(IX,IY) )THEN
          ISEQ=ISEQ+1
          I1SEQX(ISEQ)=IX
          I1SEQY(ISEQ)=IY
          I2VECTOR(IX,IY)=ISEQ
        ENDIF
      ENDIF
    END DO
  END DO
  ISEQ1=1
  ISEQ2=ISEQ
  
  AGAIN=1
  DO WHILE( AGAIN==1 )
    AGAIN=0
    JSEQ=ISEQ2
    DO ISEQ=ISEQ1, ISEQ2
      IX=I1SEQX(ISEQ)
      IY=I1SEQY(ISEQ)
      JX=I2NEXTX(IX,IY)
      JY=I2NEXTY(IX,IY)
      UPNOW(JX,JY)=UPNOW(JX,JY)+1
      IF( UPNOW(JX,JY)==NUPST(JX,JY) .and. I2NEXTX(JX,JY)>0 )THEN !! if all upstream calculated, register to 1D sequence
        JSEQ=JSEQ+1
        I1SEQX(JSEQ)=JX
        I1SEQY(JSEQ)=JY
        I2VECTOR(JX,JY)=JSEQ
        AGAIN=1
      ENDIF
    END DO
    ISEQ1=ISEQ2+1
    ISEQ2=JSEQ
  END DO
  NSEQRIV=JSEQ
  
  ISEQ=NSEQRIV
  DO IY=1, CMF_CONFIG%NY
    DO IX=1, CMF_CONFIG%NX
      IF( I2NEXTX(IX,IY).LT.0 .AND. I2NEXTX(IX,IY).NE.CMF_PARAMS%IMIS .AND. I2REGION(IX,IY)== REGIONTHIS )THEN
        ISEQ=ISEQ+1
        I1SEQX(ISEQ)=IX
        I1SEQY(ISEQ)=IY
        I2VECTOR(IX,IY)=ISEQ
      ENDIF
    END DO
  END DO
  NSEQALL=ISEQ
  
  DO ISEQ=1, NSEQALL
    IX=I1SEQX(ISEQ)
    IY=I1SEQY(ISEQ)
    IF( I2NEXTX(IX,IY)>0 )THEN
      JX=I2NEXTX(IX,IY)
      JY=I2NEXTY(IX,IY)
      I1NEXT(ISEQ)=I2VECTOR(JX,JY)
    ELSE
      I1NEXT(ISEQ)=I2NEXTX(IX,IY)
    ENDIF
  END DO
  
  DEALLOCATE(NUPST,UPNOW)
        
  END SUBROUTINE CALC_1D_SEQ
  !==========================================================
  !+
  !+
  !+
  !==========================================================
  SUBROUTINE READ_BIFPARAM    !! evenly allocate pixels to mpi nodes (not used in vcurrent version)
  USE CMF_VARS_MOD,        ONLY: NPTHOUT, NPTHLEV, PTH_UPST, PTH_DOWN,&
                              & PTH_DST, PTH_ELV, PTH_WTH,  PTH_MAN
  USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
  IMPLICIT NONE
  !* local variables
  INTEGER(KIND=JPIM)         :: IX,IY, JX,JY
  INTEGER(KIND=JPIM)         :: IPTH,  ILEV,  NPTHOUT1
  REAL(KIND=JPRB)            :: PELV,  PWTH,  PDPH
  !================================================
  WRITE(CMF_FILES%LOGNAM,*)"RIVMAP_INIT: Bifuraction channel:", TRIM(CMF_MAPS%CPTHOUT)
  
  CMF_FILES%TMPNAM=INQUIRE_FID()
  OPEN(CMF_FILES%TMPNAM,FILE=TRIM(CMF_MAPS%CPTHOUT),FORM='FORMATTED')
  READ(CMF_FILES%TMPNAM,*) NPTHOUT,NPTHLEV
  
  WRITE(CMF_FILES%LOGNAM,*) "Bifurcation channel dimantion", NPTHOUT, NPTHLEV
  
  ALLOCATE( PTH_UPST(NPTHOUT) )
  ALLOCATE( PTH_DOWN(NPTHOUT) )
  ALLOCATE( PTH_DST(NPTHOUT)  )
  ALLOCATE( PTH_ELV(NPTHOUT,NPTHLEV) )
  ALLOCATE( PTH_WTH(NPTHOUT,NPTHLEV) )
  ALLOCATE( PTH_MAN(NPTHLEV)  )
  
  NPTHOUT1=0
  DO IPTH=1, NPTHOUT
    READ(CMF_FILES%TMPNAM,*) IX, IY, JX, JY, PTH_DST(IPTH), PELV, PDPH, (PTH_WTH(IPTH,ILEV),ILEV=1,NPTHLEV)
    PTH_UPST(IPTH)=I2VECTOR(IX,IY)
    PTH_DOWN(IPTH)=I2VECTOR(JX,JY)
    IF (PTH_UPST(IPTH) > 0 .AND. PTH_DOWN(IPTH) > 0) THEN
      NPTHOUT1=NPTHOUT1+1
    ENDIF
    DO ILEV=1, NPTHLEV
      IF( ILEV==1 )THEN            !!ILEV=1: water channel bifurcation. consider bifurcation channel depth
        PWTH=PTH_WTH(IPTH,ILEV)
        IF( PWTH>0 )then
          PTH_ELV(IPTH,ILEV)=PELV - PDPH
        ELSE
          PTH_ELV(IPTH,ILEV)=1.E20
        ENDIF
      ELSE
        PWTH=PTH_WTH(IPTH,ILEV)
        IF( PWTH>0 )then
          PTH_ELV(IPTH,ILEV)=PELV + ILEV - 2.0    !! ILEV=2: bank top level 
        ELSE
          PTH_ELV(IPTH,ILEV)=1.E20
        ENDIF
      ENDIF
    END DO
  END DO
  CLOSE(CMF_FILES%TMPNAM)
  
  DO ILEV=1, NPTHLEV
    IF( ILEV==1 )THEN
      PTH_MAN(ILEV)=CMF_PARAMS%PMANRIV
    ELSE
      PTH_MAN(ILEV)=CMF_PARAMS%PMANFLD
    ENDIF
  END DO
  
  IF (NPTHOUT /= NPTHOUT1) THEN
    WRITE(CMF_FILES%LOGNAM,*)"Bifuraction channel outside of domain. Only valid:", NPTHOUT1
  ENDIF
  
  END SUBROUTINE READ_BIFPARAM
  !==========================================================
  
  END SUBROUTINE CMF_RIVMAP_INIT
  !####################################################################
  
  
  
  
  
  !####################################################################
  SUBROUTINE CMF_TOPO_INIT
  ! read & set topography map 
  ! -- call from CMF_DRV_INIT

  USE CMF_VARS_MOD,   ONLY: D2NXTDST, D2GRAREA, D2ELEVTN, D2RIVLEN, &
                          & D2RIVWTH, D2RIVHGT, D2FLDHGT, D2RIVELV, &
                          & D2FLDGRD, D2RIVMAN, P2RIVSTOMAX, P2FLDSTOMAX,  &
                          & DFRCINC,  NSEQALL,  NSEQMAX, D2MEANSL, D2DWNELV, &
                          & D2GDWDLY, I2MASK
  IMPLICIT NONE
  !================================================
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
  
  WRITE(CMF_FILES%LOGNAM,*) 'CMF::TOPO_INIT: topography map initialization'
  
  ! *** 1. ALLOCATE ARRAYS
  ALLOCATE( D2GRAREA(NSEQMAX,1) )
  ALLOCATE( D2ELEVTN(NSEQMAX,1) )
  ALLOCATE( D2NXTDST(NSEQMAX,1) )
  ALLOCATE( D2RIVLEN(NSEQMAX,1) )
  ALLOCATE( D2RIVWTH(NSEQMAX,1) )
  ALLOCATE( D2RIVHGT(NSEQMAX,1) )
  ALLOCATE( D2FLDHGT(NSEQMAX,1,CMF_CONFIG%NLFP) )
  ALLOCATE( D2RIVMAN(NSEQMAX,1) )
  ALLOCATE( D2MEANSL(NSEQMAX,1) )
  ALLOCATE( D2DWNELV(NSEQMAX,1) )
  ALLOCATE( D2GDWDLY(NSEQMAX,1) )
  ALLOCATE( I2MASK(NSEQMAX,1) )
  
  D2GRAREA(:,:)  =0._JPRB
  D2ELEVTN(:,:)  =0._JPRB
  D2NXTDST(:,:)  =0._JPRB
  D2RIVLEN(:,:)  =0._JPRB
  D2RIVWTH(:,:)  =0._JPRB
  D2RIVHGT(:,:)  =0._JPRB
  D2FLDHGT(:,:,:)=0._JPRB
  D2RIVMAN(:,:)  =0._JPRB
  D2MEANSL(:,:)  =0._JPRB
  D2DWNELV(:,:)  =0._JPRB
  D2GDWDLY(:,:)  =0._JPRB
  I2MASK(:,:)    =0._JPIM     !! mask for calculation (IFS slopemix: Kinemacti Wave for Mask=1; Reservoir: dam=2, dam upstream=1)
  
  !============================
  ! *** 2. Read topo map
  WRITE(CMF_FILES%LOGNAM,*) 'CMF::TOPO_INIT: read topography maps'
  IF ( .not. CMF_MAPS%LMAPCDF ) THEN
    CALL READ_TOPO_BIN
  ELSE
    CALL READ_TOPO_CDF
  ENDIF
  
  !============================
  ! *** 3a. Calc Channel Parameters
  WRITE(CMF_FILES%LOGNAM,*) 'TOPO_INIT: calc river channel parameters'
  
  ALLOCATE(P2RIVSTOMAX(NSEQMAX,1))
  ALLOCATE(D2RIVELV(NSEQMAX,1))
  
  IF ( CMF_OPTIONS%LFPLAIN ) THEN
    P2RIVSTOMAX(:,:) = D2RIVLEN(:,:) * D2RIVWTH(:,:) * D2RIVHGT(:,:)
  ELSE
    WRITE(CMF_FILES%LOGNAM,*) 'TOPO_INIT: no floodplain (rivstomax=1.D18)'
    P2RIVSTOMAX(:,:) = 1.E18
  ENDIF
  D2RIVELV(:,:) = D2ELEVTN(:,:) - D2RIVHGT(:,:)
  
  !*** 3b. Calc Channel Parameters
  WRITE(CMF_FILES%LOGNAM,*) 'TOPO_INIT: calc floodplain parameters'
  
  ALLOCATE(P2FLDSTOMAX(NSEQMAX,1,CMF_CONFIG%NLFP))
  ALLOCATE(D2FLDGRD(NSEQMAX,1,CMF_CONFIG%NLFP))
  CALL SET_FLDSTG
  
  !*** 3c. Calc downstream boundary
  WRITE(CMF_FILES%LOGNAM,*) 'TOPO_INIT: calc downstream boundary elevation'
  D2DWNELV(:,:)=D2ELEVTN(:,:)
  IF( CMF_OPTIONS%LMEANSL ) THEN
    D2DWNELV(:,:)=D2ELEVTN(:,:)+D2MEANSL(:,:)
  ENDIF
  
  CONTAINS
  !==========================================================
  !+ READ_TOPO_BIN
  !+ READ_TOPO_CDF
  !+ SET_FLDSTG
  !+ SET_SLOPEMIX
  !==========================================================
  SUBROUTINE READ_TOPO_BIN
  USE CMF_UTILS_MOD,       ONLY: mapR2vecD, CONV_END,  INQUIRE_FID
  IMPLICIT NONE
  !* local variables
  INTEGER(KIND=JPIM)          :: ILFP
  REAL(KIND=JPRM),ALLOCATABLE :: R2TEMP(:,:)
  REAL(KIND=JPRB),ALLOCATABLE :: D2TEMP(:,:)
  !================================================
  ALLOCATE(R2TEMP(CMF_CONFIG%NX,CMF_CONFIG%NY))
  ALLOCATE(D2TEMP(NSEQMAX,1))
  
  CMF_FILES%TMPNAM=INQUIRE_FID()
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: unit-catchment area : ',TRIM(CMF_MAPS%CGRAREA) 
  OPEN(CMF_FILES%TMPNAM,FILE=CMF_MAPS%CGRAREA,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
  READ(CMF_FILES%TMPNAM,REC=1) R2TEMP(:,:)
    IF( CMF_OPTIONS%LMAPEND ) CALL CONV_END(R2TEMP,CMF_CONFIG%NX,CMF_CONFIG%NY)
  CALL mapR2vecD(R2TEMP,D2GRAREA)
  CLOSE(CMF_FILES%TMPNAM)
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: ground elevation : ',TRIM(CMF_MAPS%CELEVTN)
  OPEN(CMF_FILES%TMPNAM,FILE=CMF_MAPS%CELEVTN,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
  READ(CMF_FILES%TMPNAM,REC=1) R2TEMP(:,:)
    IF( CMF_OPTIONS%LMAPEND ) CALL CONV_END(R2TEMP,CMF_CONFIG%NX,CMF_CONFIG%NY)
  CALL mapR2vecD(R2TEMP,D2ELEVTN)
  CLOSE(CMF_FILES%TMPNAM)
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: downstream distance : ',TRIM(CMF_MAPS%CNXTDST)
  OPEN(CMF_FILES%TMPNAM,FILE=CMF_MAPS%CNXTDST,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
  READ(CMF_FILES%TMPNAM,REC=1) R2TEMP(:,:)
    IF( CMF_OPTIONS%LMAPEND ) CALL CONV_END(R2TEMP,CMF_CONFIG%NX,CMF_CONFIG%NY)
  CALL mapR2vecD(R2TEMP,D2NXTDST)
  CLOSE(CMF_FILES%TMPNAM)
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: river channel length : ',TRIM(CMF_MAPS%CRIVLEN)
  OPEN(CMF_FILES%TMPNAM,FILE=CMF_MAPS%CRIVLEN,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
  READ(CMF_FILES%TMPNAM,REC=1) R2TEMP(:,:)
    IF( CMF_OPTIONS%LMAPEND ) CALL CONV_END(R2TEMP,CMF_CONFIG%NX,CMF_CONFIG%NY)
  CALL mapR2vecD(R2TEMP,D2RIVLEN)
  CLOSE(CMF_FILES%TMPNAM)
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: floodplain elevation profile : ',TRIM(CMF_MAPS%CFLDHGT)
  OPEN(CMF_FILES%TMPNAM,FILE=TRIM(CMF_MAPS%CFLDHGT),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
  DO ILFP=1,CMF_CONFIG%NLFP
    READ(CMF_FILES%TMPNAM,REC=ILFP) R2TEMP
    IF( CMF_OPTIONS%LMAPEND ) CALL CONV_END(R2TEMP,CMF_CONFIG%NX,CMF_CONFIG%NY)
    CALL mapR2vecD(R2TEMP,D2TEMP)
    D2FLDHGT(:,:,ILFP)= D2TEMP(:,:)
  ENDDO
  CLOSE(CMF_FILES%TMPNAM)
  
  !*** river channel / groundwater parameters)
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: river channel depth : ',TRIM(CMF_MAPS%CRIVHGT)
  OPEN(CMF_FILES%TMPNAM,FILE=CMF_MAPS%CRIVHGT,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
  READ(CMF_FILES%TMPNAM,REC=1) R2TEMP(:,:)
    IF( CMF_OPTIONS%LMAPEND ) CALL CONV_END(R2TEMP,CMF_CONFIG%NX,CMF_CONFIG%NY)
  CALL mapR2vecD(R2TEMP,D2RIVHGT)
  CLOSE(CMF_FILES%TMPNAM)
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: river channel width : ',TRIM(CMF_MAPS%CRIVWTH)
  OPEN(CMF_FILES%TMPNAM,FILE=CMF_MAPS%CRIVWTH,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
  READ(CMF_FILES%TMPNAM,REC=1) R2TEMP(:,:)
    IF( CMF_OPTIONS%LMAPEND ) CALL CONV_END(R2TEMP,CMF_CONFIG%NX,CMF_CONFIG%NY)
  CALL mapR2vecD(R2TEMP,D2RIVWTH)
  CLOSE(CMF_FILES%TMPNAM)
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: manning coefficient river: ',TRIM(CMF_MAPS%CRIVMAN)
  OPEN(CMF_FILES%TMPNAM,FILE=TRIM(CMF_MAPS%CRIVMAN),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
  READ(CMF_FILES%TMPNAM,REC=1) R2TEMP(:,:)
    IF( CMF_OPTIONS%LMAPEND ) CALL CONV_END(R2TEMP,CMF_CONFIG%NX,CMF_CONFIG%NY)
  CALL mapR2vecD(R2TEMP,D2RIVMAN)
  CLOSE(CMF_FILES%TMPNAM)
  
  IF( CMF_OPTIONS%LGDWDLY  )THEN
    WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: groundwater delay parameter: ',TRIM(CMF_MAPS%CGDWDLY)
    OPEN(CMF_FILES%TMPNAM,FILE=TRIM(CMF_MAPS%CGDWDLY),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
    READ(CMF_FILES%TMPNAM,REC=1) R2TEMP(:,:)
      IF( CMF_OPTIONS%LMAPEND ) CALL CONV_END(R2TEMP,CMF_CONFIG%NX,CMF_CONFIG%NY)
    CALL mapR2vecD(R2TEMP,D2GDWDLY)
    CLOSE(CMF_FILES%TMPNAM)
  ENDIF
  
  IF( CMF_OPTIONS%LSLPMIX )THEN
    WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: LSLPMIX only used in IFS, not availabke with binary map'
  ENDIF
  IF( CMF_OPTIONS%LSLOPEMOUTH )THEN
    WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: LSLOPEMOUTH only used in IFS, not availabke with binary map'
  ENDIF
  
  ! ==========
  
  IF( CMF_OPTIONS%LMEANSL ) THEN
    WRITE(CMF_FILES%LOGNAM, *)'TOPO_INIT: mean sea level: ', TRIM(CMF_MAPS%CMEANSL)
    OPEN(CMF_FILES%TMPNAM, FILE=CMF_MAPS%CMEANSL, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
    READ(CMF_FILES%TMPNAM, REC=1) R2TEMP(:,:)
      IF( CMF_OPTIONS%LMAPEND ) CALL CONV_END(R2TEMP,CMF_CONFIG%NX,CMF_CONFIG%NY)
    CALL mapR2vecD(R2TEMP, D2MEANSL)
    CLOSE(CMF_FILES%TMPNAM)
  ENDIF
  
  DEALLOCATE(R2TEMP)
  DEALLOCATE(D2TEMP)
  
  END SUBROUTINE READ_TOPO_BIN
  !==========================================================
  !+
  !+
  !+
  !==========================================================
  SUBROUTINE READ_TOPO_CDF
#ifdef UseCDF_CMF
  USE NETCDF 
  USE CMF_UTILS_MOD,            ONLY: NCERROR,mapR2vecD
  USE CMF_VARS_MOD,              ONLY: D2ELEVSLOPE     !! only used in ECMWF
  IMPLICIT NONE
  !* local variables
  INTEGER(KIND=JPIM)               :: NCID,VARID,STATUS
  INTEGER(KIND=JPIM)               :: ILEV
  REAL(KIND=JPRM),ALLOCATABLE      :: R2TEMP(:,:)
  REAL(KIND=JPRB),ALLOCATABLE      :: D2TEMP(:,:)
  !================================================
  ALLOCATE(R2TEMP(CMF_CONFIG%NX,CMF_CONFIG%NY))
  ALLOCATE(D2TEMP(NSEQMAX,1))
  
  !! CLIM FILE
  CALL NCERROR (NF90_OPEN(CMF_MAPS%CRIVCLINC,NF90_NOWRITE,NCID),'opening '//TRIM(CMF_MAPS%CRIVCLINC) )
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: ctmare:',TRIM(CMF_MAPS%CRIVCLINC)
  CALL NCERROR ( NF90_INQ_VARID(NCID,'ctmare',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' )
  CALL mapR2vecD(R2TEMP,D2GRAREA)
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: elevtn:',TRIM(CMF_MAPS%CRIVCLINC)
  CALL NCERROR ( NF90_INQ_VARID(NCID,'elevtn',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
  CALL mapR2vecD(R2TEMP,D2ELEVTN)
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: nxtdst:',TRIM(CMF_MAPS%CRIVCLINC)
  CALL NCERROR ( NF90_INQ_VARID(NCID,'nxtdst',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
  CALL mapR2vecD(R2TEMP,D2NXTDST)
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: rivlen:',TRIM(CMF_MAPS%CRIVCLINC)
  CALL NCERROR ( NF90_INQ_VARID(NCID,'rivlen',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
  CALL mapR2vecD(R2TEMP,D2RIVLEN)
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: fldhgt:',TRIM(CMF_MAPS%CRIVCLINC)
  CALL NCERROR ( NF90_INQ_VARID(NCID,'fldhgt',VARID),'getting id' )
  DO ILEV=1,CMF_CONFIG%NLFP
    CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP,(/1,1,ILEV/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,1/)),'reading data' ) 
    CALL mapR2vecD(R2TEMP,D2TEMP)
    D2FLDHGT(:,:,ILEV)=D2TEMP(:,:)
  ENDDO
  
  CALL NCERROR( NF90_CLOSE(NCID))
  
  IF ( CMF_OPTIONS%LSLOPEMOUTH ) THEN
    ALLOCATE( D2ELEVSLOPE(NSEQMAX,1) )
    WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: elevslope:',TRIM(CMF_MAPS%CRIVPARNC)
    STATUS = NF90_INQ_VARID(NCID,'elevslope',VARID)
    IF (STATUS /= 0 ) THEN
      WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: elevslope: not present, aborting'
      STOP 9 
    ELSE
      CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
    ENDIF 
    CALL mapR2vecD(R2TEMP,D2ELEVSLOPE)
  ENDIF
  
  !!========== 
  !! PAR FILE (river channel / groundwater parameters)
  CALL NCERROR (NF90_OPEN(CMF_MAPS%CRIVPARNC,NF90_NOWRITE,NCID),'opening '//TRIM(CMF_MAPS%CRIVPARNC) )
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: rivwth:',TRIM(CMF_MAPS%CRIVPARNC)
  CALL NCERROR ( NF90_INQ_VARID(NCID,'rivwth',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
  CALL mapR2vecD(R2TEMP,D2RIVWTH)
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: rivhgt:',TRIM(CMF_MAPS%CRIVPARNC)
  CALL NCERROR ( NF90_INQ_VARID(NCID,'rivhgt',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
  CALL mapR2vecD(R2TEMP,D2RIVHGT)
  
  WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: rivman:',TRIM(CMF_MAPS%CRIVPARNC)
  CALL NCERROR ( NF90_INQ_VARID(NCID,'rivman',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
  CALL mapR2vecD(R2TEMP,D2RIVMAN)
  
  IF ( CMF_OPTIONS%LGDWDLY ) THEN
    WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: GDWDLY:',TRIM(CMF_MAPS%CRIVPARNC)
    STATUS = NF90_INQ_VARID(NCID,'gdwdly',VARID)
    IF (STATUS /= 0 ) THEN
      WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: GDWDLY: not present, setting to zero'
      R2TEMP(:,:) = 0._JPRB
    ELSE
      CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
    ENDIF 
    CALL mapR2vecD(R2TEMP,D2GDWDLY)
  ENDIF
  
  I2MASK(:,:)=0_JPIM
  IF ( CMF_OPTIONS%LSLPMIX ) THEN
    CALL SET_SLOPEMIX
  ENDIF
  
  CALL NCERROR( NF90_CLOSE(NCID))
  
  !!========== 
  !! MEAN SEA LEVEL FILE
  IF( CMF_OPTIONS%LMEANSL ) THEN
    CALL NCERROR (NF90_OPEN(CMF_MAPS%CMEANSLNC,NF90_NOWRITE,NCID),'opening '//TRIM(CMF_MAPS%CMEANSLNC) )
    WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: rivhgt:',TRIM(CMF_MAPS%CMEANSLNC)
    CALL NCERROR ( NF90_INQ_VARID(NCID,'meansl',VARID),'getting id' )
    CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
    CALL mapR2vecD ( R2TEMP,D2MEANSL  )
    CALL NCERROR ( NF90_CLOSE(NCID) )
  ENDIF 
  
  DEALLOCATE(R2TEMP)
  DEALLOCATE(D2TEMP)
#endif
  END SUBROUTINE READ_TOPO_CDF
  !==========================================================
  !+
!####################################################################
!==========================================================
!+
!+
!+
!==========================================================
  SUBROUTINE SET_FLDSTG
    IMPLICIT NONE
    !* local variables
    INTEGER(KIND=JPIM),SAVE  ::  ISEQ, I
    REAL(KIND=JPRB),SAVE     ::  DSTONOW
    REAL(KIND=JPRB),SAVE     ::  DSTOPRE
    REAL(KIND=JPRB),SAVE     ::  DHGTPRE
    REAL(KIND=JPRB),SAVE     ::  DWTHINC
!$OMP THREADPRIVATE               (I,DSTONOW,DSTOPRE,DHGTPRE,DWTHINC)
    !================================================
    P2FLDSTOMAX(:,:,:) = 0._JPRD
    D2FLDGRD(:,:,:)    = 0._JPRB
    DFRCINC=dble(CMF_CONFIG%NLFP)**(-1.)
    !
    !$OMP PARALLEL DO
    DO ISEQ=1, NSEQALL
      DSTOPRE = P2RIVSTOMAX(ISEQ,1)
      DHGTPRE = 0._JPRB
      DWTHINC = D2GRAREA(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * DFRCINC
      DO I=1, CMF_CONFIG%NLFP
        DSTONOW = D2RIVLEN(ISEQ,1) * ( D2RIVWTH(ISEQ,1) + DWTHINC*(DBLE(I)-0.5) ) * (D2FLDHGT(ISEQ,1,I)-DHGTPRE)
        P2FLDSTOMAX(ISEQ,1,I) = DSTOPRE + DSTONOW
        D2FLDGRD(ISEQ,1,I) = (D2FLDHGT(ISEQ,1,I)-DHGTPRE) * DWTHINC**(-1.)
        DSTOPRE = P2FLDSTOMAX(ISEQ,1,I)
        DHGTPRE = D2FLDHGT(ISEQ,1,I)
      END DO
    END DO
!$OMP END PARALLEL DO
    
    !
    END SUBROUTINE SET_FLDSTG
    !==========================================================
    !+
    !+
    !+
    !==========================================================
    SUBROUTINE SET_SLOPEMIX    !! only used in IFS0
#ifdef UseCDF_CMF
    USE NETCDF 
    USE CMF_UTILS_MOD,           ONLY: NCERROR,mapI2vecI
    
    IMPLICIT NONE
    INTEGER(KIND=JPIM),ALLOCATABLE  :: I2TEMP(:,:)
    INTEGER(KIND=JPIM)              :: ISEQ, I0, I1
    INTEGER(KIND=JPIM)              :: NCID,VARID,STATUS
    
      ALLOCATE(I2TEMP(CMF_CONFIG%NX,CMF_CONFIG%NY))
      WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: mask_slope:',TRIM(CMF_MAPS%CRIVPARNC)
      STATUS =  NF90_INQ_VARID(NCID,'mask_slope',VARID)
      IF (STATUS /= 0 ) THEN
        WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: mask_slope: LSLPMIX should be set to FALSE: ABORTING!'
        STOP 9
      ENDIF 
      CALL NCERROR ( NF90_GET_VAR(NCID,VARID,I2TEMP),'reading data' ) 
      CALL mapI2vecI(I2TEMP,I2MASK)
      I0=0
      I1=0
      DO ISEQ=1,NSEQALL
        IF (I2MASK(ISEQ,1) == 1 ) THEN  !! kinematic wave applied
          I1=I1+1
        ENDIF
        IF (I2MASK(ISEQ,1) == 0 ) THEN 
          I0=I0+1
        ENDIF
      ENDDO
      WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: sum(mask==0), sum(mask==1)',I0,I1
      IF ( I0+I1 .NE. NSEQALL ) THEN 
         WRITE(CMF_FILES%LOGNAM,*)'TOPO_INIT: mask==0 + mask == 1 does not match NSEQALL.. something wrong, aborting'
         STOP 9
      ENDIF 
    
      DEALLOCATE(I2TEMP)
#endif
  END SUBROUTINE SET_SLOPEMIX
!==========================================================
END SUBROUTINE CMF_TOPO_INIT

SUBROUTINE CMF_LEVEE_INIT
  USE CMF_VARS_MOD,   ONLY: NSEQALL, NSEQMAX, D2GRAREA, D2RIVLEN, D2RIVWTH, D2FLDHGT, &
                             & D2FLDGRD, P2RIVSTOMAX, P2FLDSTOMAX, DFRCINC
  USE CMF_VARS_MOD,   ONLY: D2LEVHGT, D2LEVFRC, D2BASHGT, D2LEVDST, &
                             & D2LEVBASSTO, D2LEVTOPSTO, D2LEVFILSTO,NPTHLEV
  USE CMF_UTILS_MOD,      ONLY: mapR2vecD, INQUIRE_FID
  !
  IMPLICIT NONE
  !* local variables
  REAL(KIND=JPRM)            ::  R2TEMP(CMF_CONFIG%NX,CMF_CONFIG%NY)
  ! SAVE for OpenMP
  INTEGER(KIND=JPIM),SAVE    ::  ISEQ, I, ILEV  
  REAL(KIND=JPRB),SAVE       ::  DSTONOW,DSTOPRE,DHGTPRE,DWTHINC,DWTHPRE,DWTHNOW,DHGTNOW,DHGTDIF
!$OMP THREADPRIVATE    (I,ILEV,DSTONOW,DSTOPRE,DHGTPRE,DWTHINC,DWTHPRE,DWTHNOW,DHGTNOW,DHGTDIF)
!####################################################################
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
  WRITE(CMF_FILES%LOGNAM,*) "CMF::LEVEE_INIT: initialize levee"
  
  !********************
  ! [1] Read Levee Parameter Map
  WRITE(CMF_FILES%LOGNAM,*) "CMF::LEVEE_INIT: read levee parameter files"
  
  ALLOCATE( D2LEVHGT(NSEQMAX,1) )
  ALLOCATE( D2LEVFRC(NSEQMAX,1) )
  D2LEVHGT(:,:)   =0._JPRB
  D2LEVFRC(:,:)   =0._JPRB
  
  CMF_FILES%TMPNAM=INQUIRE_FID()
  
  WRITE(CMF_FILES%LOGNAM,*)'INIT_LEVEE: levee crown height : ',TRIM(CMF_LEVEE%CLEVHGT)
  OPEN(CMF_FILES%TMPNAM,FILE=CMF_LEVEE%CLEVHGT,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
  READ(CMF_FILES%TMPNAM,REC=1) R2TEMP(:,:)
  CALL mapR2vecD(R2TEMP,D2LEVHGT)
  CLOSE(CMF_FILES%TMPNAM)
  
  WRITE(CMF_FILES%LOGNAM,*)'INIT_LEVEE: distance from levee to river : ',TRIM(CMF_LEVEE%CLEVFRC)
  OPEN(CMF_FILES%TMPNAM,FILE=CMF_LEVEE%CLEVFRC,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
  READ(CMF_FILES%TMPNAM,REC=1) R2TEMP(:,:)
  CALL mapR2vecD(R2TEMP,D2LEVFRC)
  CLOSE(CMF_FILES%TMPNAM)
  
  !*******************************
  ! [2] Calculate Levee Stage Parameter
  WRITE(CMF_FILES%LOGNAM,*) "CMF::LEVEE_INIT: flood stage parameters considering levee"
  
  ALLOCATE( D2BASHGT(NSEQMAX,1) )
  ALLOCATE( D2LEVDST(NSEQMAX,1) )
  
  ALLOCATE( D2LEVBASSTO(NSEQMAX,1) )
  ALLOCATE( D2LEVTOPSTO(NSEQMAX,1) )
  ALLOCATE( D2LEVFILSTO(NSEQMAX,1) )
  
  P2FLDSTOMAX(:,:,:) = 0._JPRD   !! max floodplain  storage  at each layer
  D2FLDGRD(:,:,:)    = 0._JPRB   !! floodplain topo gradient of each layer
  DFRCINC=dble(CMF_CONFIG%NLFP)**(-1.)   !! fration of each layer
  
  D2LEVBASSTO(:,:)= 0._JPRB      !! storage at levee base     (levee protection start)
  D2LEVTOPSTO(:,:)= 0._JPRB      !! storage at levee top      (levee protection end)
  D2LEVFILSTO(:,:)= 0._JPRB      !! storage when levee filled (protected-side depth reach levee top)
  
  !$OMP PARALLEL DO
  DO ISEQ=1, NSEQALL
    IF( D2LEVHGT(ISEQ,1)<=0._JPRB )THEN
      D2LEVHGT(ISEQ,1)=0._JPRB
      D2LEVFRC(ISEQ,1)=1._JPRB   !! If no levee, all area is unprotected/
    ENDIF
    D2LEVFRC(ISEQ,1)=MAX(0._JPRB,MIN(1._JPRB,D2LEVFRC(ISEQ,1)))
  END DO
!$OMP END PARALLEL DO
  
!$OMP PARALLEL DO
  DO ISEQ=1, NSEQALL
  ! calculate floodplain parameters (without levee, same as SET_FLDSTG)
    DSTOPRE = P2RIVSTOMAX(ISEQ,1)
    DHGTPRE = 0._JPRB
    DWTHINC = D2GRAREA(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * DFRCINC  !! width increlment for each layer
    DO I=1, CMF_CONFIG%NLFP
      DSTONOW = D2RIVLEN(ISEQ,1) * ( D2RIVWTH(ISEQ,1) + DWTHINC*(DBLE(I)-0.5) ) * (D2FLDHGT(ISEQ,1,I)-DHGTPRE)  !! storage increment
      P2FLDSTOMAX(ISEQ,1,I) = DSTOPRE + DSTONOW
      D2FLDGRD(ISEQ,1,I) = (D2FLDHGT(ISEQ,1,I)-DHGTPRE) * DWTHINC**(-1.)
      DSTOPRE = P2FLDSTOMAX(ISEQ,1,I)
      DHGTPRE = D2FLDHGT(ISEQ,1,I)
    END DO
  
  ! Levee parameters calculation
    IF( D2LEVHGT(ISEQ,1) == 0._JPRB )THEN ! Grid without levee, treat everything as unprotected
      D2BASHGT(ISEQ,1) = 1.E18
      D2LEVDST(ISEQ,1) = 1.E18
      D2LEVBASSTO(ISEQ,1) = 1.E18
      D2LEVTOPSTO(ISEQ,1) = 1.E18
      D2LEVFILSTO(ISEQ,1) = 1.E18
    ELSE  !! levee exist
      !!*********
      !! [1] levee base storage & levee top storage (water only in river side)
  
      DSTOPRE = P2RIVSTOMAX(ISEQ,1)
      DHGTPRE = 0._JPRB
      DWTHPRE = 0._JPRB
      D2LEVDST(ISEQ,1) = D2LEVFRC(ISEQ,1) * DWTHINC*CMF_CONFIG%NLFP !! distance from channel to levee [m]
  
      ILEV=INT( D2LEVFRC(ISEQ,1)*CMF_CONFIG%NLFP )+1 !! which layer levee exist
      IF( ILEV>=2 )THEN
        DSTOPRE = P2FLDSTOMAX(ISEQ,1,ILEV-1)
        DHGTPRE = D2FLDHGT(ISEQ,1,ILEV-1)
        DWTHPRE = DWTHINC * (ILEV-1)
      ENDIF
  
      IF( ILEV<=CMF_CONFIG%NLFP )THEN
        !! levee in floodplain layer ILEV
        DWTHNOW = D2LEVDST(ISEQ,1) - DWTHPRE
        DHGTNOW = DWTHNOW * D2FLDGRD(ISEQ,1,ILEV) !! levee height above lower floodplain profile point
        D2BASHGT(ISEQ,1) = DHGTNOW + DHGTPRE
        D2LEVHGT(ISEQ,1) = max( D2LEVHGT(ISEQ,1), D2BASHGT(ISEQ,1) ) !! levee height >= base height
  
        DSTONOW = ( DWTHNOW*0.5 + DWTHPRE + D2RIVWTH(ISEQ,1) ) * DHGTNOW * D2RIVLEN(ISEQ,1) 
        D2LEVBASSTO(ISEQ,1) = DSTOPRE + DSTONOW
  
        DHGTDIF = D2LEVHGT(ISEQ,1) - D2BASHGT(ISEQ,1)
        D2LEVTOPSTO(ISEQ,1) = D2LEVBASSTO(ISEQ,1) + ( D2LEVDST(ISEQ,1)+D2RIVWTH(ISEQ,1) ) * DHGTDIF * D2RIVLEN(ISEQ,1)
      ELSE
        !! levee on the floodplain edge (ILEV=NLEV+1)
        D2BASHGT(ISEQ,1) = DHGTPRE
        D2LEVHGT(ISEQ,1) = max( D2LEVHGT(ISEQ,1), D2BASHGT(ISEQ,1) ) !! levee height >= base height
  
        D2LEVBASSTO(ISEQ,1) = DSTOPRE
  
        DHGTDIF = D2LEVHGT(ISEQ,1) - D2BASHGT(ISEQ,1)
        D2LEVTOPSTO(ISEQ,1) = D2LEVBASSTO(ISEQ,1) + ( D2LEVDST(ISEQ,1)+D2RIVWTH(ISEQ,1) ) * DHGTDIF * D2RIVLEN(ISEQ,1)
      ENDIF
  
      !!*********
      !! [2] levee fill storage (water in both river side & protected side)
      I=1
      DSTOPRE = P2RIVSTOMAX(ISEQ,1)
      DWTHPRE = D2RIVWTH(ISEQ,1)
      DHGTPRE = 0._JPRB
  
      !! check which layer levee top belongs
      DO WHILE( D2LEVHGT(ISEQ,1) > D2FLDHGT(ISEQ,1,I) .AND. I<=CMF_CONFIG%NLFP )
        DSTOPRE = P2FLDSTOMAX(ISEQ,1,I)
        DWTHPRE = DWTHPRE + DWTHINC
        DHGTPRE = D2FLDHGT(ISEQ,1,I)
        I=I+1
        IF( I>CMF_CONFIG%NLFP ) EXIT
      END DO
  
      !! calculate levee fill volume
      IF( I<=CMF_CONFIG%NLFP )THEN 
        !! levee top height collesponds to layer I
        DHGTNOW = D2LEVHGT(ISEQ,1) - DHGTPRE
        DWTHNOW = DHGTNOW * D2FLDGRD(ISEQ,1,I)**(-1.)
  
        DSTONOW = ( DWTHNOW*0.5 + DWTHPRE ) * DHGTNOW * D2RIVLEN(ISEQ,1)
        D2LEVFILSTO(ISEQ,1) = DSTOPRE + DSTONOW
      ELSE
        !! levee higher than catchment boundary height
        DHGTNOW = D2LEVHGT(ISEQ,1) - DHGTPRE
        DSTONOW = DWTHPRE * DHGTNOW * D2RIVLEN(ISEQ,1)
        D2LEVFILSTO(ISEQ,1) = DSTOPRE + DSTONOW
      ENDIF
    ENDIF
  
  END DO
!$OMP END PARALLEL DO
  
END SUBROUTINE CMF_LEVEE_INIT

SUBROUTINE CMF_OUTPUT_INIT
! Initialize output module (create/open files)
! -- Called from CMF_DRV_INIT
USE CMF_VARS_MOD,            ONLY: ISYYYY, ISMM,   ISDD,   ISHOUR, ISMIN
USE CMF_VARS_MOD,            ONLY: NSEQMAX,NPTHOUT,NPTHLEV,REGIONTHIS,NVARS
USE CMF_UTILS_MOD,           ONLY: INQUIRE_FID
USE CMF_OUTPUT_MOD,          ONLY: NVARSOUT,VAROUT
IMPLICIT NONE
!* Local variables 

CHARACTER(LEN=256)              :: CTIME, CTMP
INTEGER(KIND=JPIM)              :: JF,J,J0
CHARACTER(LEN=256)              :: CVNAMES(NVARS)
!================================================


WRITE(CMF_FILES%LOGNAM,*) ""
WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"

WRITE(CMF_FILES%LOGNAM,*) "CMF::OUTPUT_INIT: check output variables"
!! Start by finding out # of output variables 
NVARSOUT=0
J0=1

DO J=1,LEN(TRIM(CMF_OUTPUT%CVARSOUT))
  IF( (J>J0) .AND. (CMF_OUTPUT%CVARSOUT(J:J) .EQ. ',') ) THEN
    CTMP=TRIM(ADJUSTL(CMF_OUTPUT%CVARSOUT(J0:J-1)))
    IF (LEN(CTMP) > 0 ) THEN
      NVARSOUT=NVARSOUT+1
      CVNAMES(NVARSOUT)=CTMP
    ENDIF
    J0=J+1
  ENDIF
ENDDO
! Last one 
IF ( J0 <= LEN(TRIM(CMF_OUTPUT%CVARSOUT)) ) THEN
  J=LEN(TRIM(CMF_OUTPUT%CVARSOUT))
  CTMP=TRIM(ADJUSTL(CMF_OUTPUT%CVARSOUT(J0:J)))
  IF (LEN(CTMP) > 0 ) THEN
     NVARSOUT=NVARSOUT+1
     CVNAMES(NVARSOUT)=CTMP
  ENDIF
ENDIF 

IF ( NVARSOUT == 0 ) THEN
  WRITE(CMF_FILES%LOGNAM,*) "CMF::OUTPUT_INIT: No output files will be produced!"
  RETURN
ENDIF 

ALLOCATE(VAROUT(NVARSOUT))
WRITE(CTIME,'(A14,I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)') 'seconds since ',ISYYYY,'-',ISMM,'-',ISDD,' ',ISHOUR,":",ISMIN
!* Loop on variables and create files 
DO JF=1,NVARSOUT
  WRITE(CMF_FILES%LOGNAM,*) "Creating output for variable:", TRIM( CVNAMES(JF) )
  SELECT CASE (CVNAMES(JF))
  CASE ('rivout')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='river discharge'
    VAROUT(JF)%CVUNITS='m3/s'
  CASE ('rivsto')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='river storage'
    VAROUT(JF)%CVUNITS='m3'
  CASE ('rivdph')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='river depth'
    VAROUT(JF)%CVUNITS='m'
  CASE ('rivvel')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='river velocity'
    VAROUT(JF)%CVUNITS='m/s'

  CASE ('fldout')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='floodplain discharge'
    VAROUT(JF)%CVUNITS='m3/s'
  CASE ('fldsto')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='floodplain storage'
    VAROUT(JF)%CVUNITS='m3'
  CASE ('flddph')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='floodplain depth'
    VAROUT(JF)%CVUNITS='m'  
  CASE ('fldfrc')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='flooded fraction'
    VAROUT(JF)%CVUNITS='0-1'  
  CASE ('fldare')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='flooded area'
    VAROUT(JF)%CVUNITS='m2'

  CASE ('sfcelv')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='water surface elevation'
    VAROUT(JF)%CVUNITS='m'
  CASE ('totout')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='discharge (river+floodplain)'
    VAROUT(JF)%CVUNITS='m3/s'
  CASE ('outflw')                   !! comparability for previous output name
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='discharge (river+floodplain)'
    VAROUT(JF)%CVUNITS='m3/s'
  CASE ('totsto')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='total storage (river+floodplain)'
    VAROUT(JF)%CVUNITS='m3'
  CASE ('storge')                   !! comparability for previous output name
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='total storage (river+floodplain)'
    VAROUT(JF)%CVUNITS='m3'

  CASE ('pthflw')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='bifurcation channel discharge'
    VAROUT(JF)%CVUNITS='m3/s'
  CASE ('pthout')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='net bifurcation discharge'
    VAROUT(JF)%CVUNITS='m3/s'

  CASE ('maxsto')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='daily maximum storage'
    VAROUT(JF)%CVUNITS='m3'  
  CASE ('maxflw')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='daily maximum discharge'
    VAROUT(JF)%CVUNITS='m3/s' 
  CASE ('maxdph')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='daily maximum river depth'
    VAROUT(JF)%CVUNITS='m' 

  CASE ('runoff')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='Surface runoff'
    VAROUT(JF)%CVUNITS='m3/s' 
  CASE ('runoffsub')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='sub-surface runoff'
    VAROUT(JF)%CVUNITS='m3/s' 

  CASE ('damsto')   !!! added
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='reservoir storage'
    VAROUT(JF)%CVUNITS='m3' 
  CASE ('daminf')   !!! added
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='reservoir inflow'
    VAROUT(JF)%CVUNITS='m3/s' 

  CASE ('levsto')   !!! added
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='protected area storage'
    VAROUT(JF)%CVUNITS='m3' 
  CASE ('levdph')   !!! added
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='protected area depth'
    VAROUT(JF)%CVUNITS='m' 

  CASE ('gdwsto')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='ground water storage'
    VAROUT(JF)%CVUNITS='m3'
  CASE ('gwsto')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='ground water storage'
    VAROUT(JF)%CVUNITS='m3'
  CASE ('gwout')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='ground water discharge'
    VAROUT(JF)%CVUNITS='m3/s'  

  CASE ('wevap')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='water evaporation'
    VAROUT(JF)%CVUNITS='m3/s'
  CASE ('outins')
    VAROUT(JF)%CVNAME=CVNAMES(JF)
    VAROUT(JF)%CVLNAME='instantaneous discharge'
    VAROUT(JF)%CVUNITS='m3/s' 

  CASE DEFAULT
    WRITE(CMF_FILES%LOGNAM,*) trim(CVNAMES(JF)), ' Not defined in CMF_CREATE_OUTCDF_MOD'
  END SELECT

  VAROUT(JF)%BINID=INQUIRE_FID()

  IF( CMF_OUTPUT%LOUTCDF )THEN
    IF( REGIONTHIS==1 )THEN
      CALL CREATE_OUTCDF
    ENDIF
  ELSE
    CALL CREATE_OUTBIN
  ENDIF
END DO

IRECOUT=0  ! Initialize Output record to 1 (shared in netcdf & binary)


WRITE(CMF_FILES%LOGNAM,*) "CMF::OUTPUT_INIT: end"


CONTAINS
!==========================================================
!+ CREATE_OUTBIN
!+ CREATE_OUTCDF
!==========================================================
SUBROUTINE CREATE_OUTBIN
IMPLICIT NONE
!================================================
IF( TRIM(VAROUT(JF)%CVNAME)=='pthflw' ) THEN   !! bifurcation channel
  IF( REGIONTHIS==1 )THEN
    VAROUT(JF)%CFILE=TRIM(CMF_OUTPUT%COUTDIR)//TRIM(VAROUT(JF)%CVNAME)//TRIM(CMF_OUTPUT%COUTTAG)//TRIM(CMF_PARAMS%CSUFPTH)
    OPEN(VAROUT(JF)%BINID,FILE=VAROUT(JF)%CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NPTHOUT*NPTHLEV)
    WRITE(CMF_FILES%LOGNAM,*) "output file opened in unit: ", TRIM(VAROUT(JF)%CFILE), VAROUT(JF)%BINID
  ENDIF
ELSEIF( CMF_OUTPUT%LOUTVEC )THEN   !!  1D land only output
  VAROUT(JF)%CFILE=TRIM(CMF_OUTPUT%COUTDIR)//TRIM(VAROUT(JF)%CVNAME)//TRIM(CMF_OUTPUT%COUTTAG)//TRIM(CMF_PARAMS%CSUFVEC)
  OPEN(VAROUT(JF)%BINID,FILE=VAROUT(JF)%CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NSEQMAX)
  WRITE(CMF_FILES%LOGNAM,*) "output file opened in unit: ", TRIM(VAROUT(JF)%CFILE), VAROUT(JF)%BINID
ELSE                   !!  2D default map output
  IF( REGIONTHIS==1 )THEN
    VAROUT(JF)%CFILE=TRIM(CMF_OUTPUT%COUTDIR)//TRIM(VAROUT(JF)%CVNAME)//TRIM(CMF_OUTPUT%COUTTAG)//TRIM(CMF_PARAMS%CSUFBIN)
    WRITE(CMF_FILES%LOGNAM,*) "  -- ", TRIM(VAROUT(JF)%CFILE)
    OPEN(VAROUT(JF)%BINID,FILE=VAROUT(JF)%CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
    WRITE(CMF_FILES%LOGNAM,*) "output file opened in unit: ", TRIM(VAROUT(JF)%CFILE), VAROUT(JF)%BINID
  ENDIF
ENDIF
END SUBROUTINE CREATE_OUTBIN
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE CREATE_OUTCDF
#ifdef UseCDF_CMF
USE CMF_VARS_MOD,             ONLY: D1LON, D1LAT
USE CMF_UTILS_MOD,           ONLY: NCERROR
USE NETCDF
IMPLICIT NONE
INTEGER(KIND=JPIM)  :: TIMEID,VARID,LATID,LONID
!============
VAROUT(JF)%IRECNC=1 ! initialize record current writting record to 1 

!============
VAROUT(JF)%CFILE=TRIM(CMF_OUTPUT%COUTDIR)//'o_'//TRIM(VAROUT(JF)%CVNAME)//TRIM(CMF_OUTPUT%COUTTAG)//TRIM(CMF_PARAMS%CSUFCDF)
! Create file 
CALL NCERROR( NF90_CREATE(VAROUT(JF)%CFILE,NF90_NETCDF4,VAROUT(JF)%NCID),&
              'CREATING FILE:'//TRIM(VAROUT(JF)%CFILE) )
!=== set dimension ===
CALL NCERROR( NF90_DEF_DIM(VAROUT(JF)%NCID, 'time', NF90_UNLIMITED, TIMEID) )
CALL NCERROR( NF90_DEF_DIM(VAROUT(JF)%NCID, 'lat', CMF_CONFIG%NY, LATID) )
CALL NCERROR( NF90_DEF_DIM(VAROUT(JF)%NCID, 'lon', CMF_CONFIG%NX, LONID) )

!=== define variables ===
CALL NCERROR( NF90_DEF_VAR(VAROUT(JF)%NCID, 'lat', NF90_FLOAT, (/LATID/), VARID) )
CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VARID, 'long_name','latitude') )
CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VARID, 'units','degrees_north') )

CALL NCERROR( NF90_DEF_VAR(VAROUT(JF)%NCID, 'lon', NF90_FLOAT, (/LONID/), VARID) )
CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VARID, 'long_name','longitude') )
CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VARID, 'units','degrees_east') )

CALL NCERROR( NF90_DEF_VAR(VAROUT(JF)%NCID, 'time', NF90_DOUBLE, (/TIMEID/), VAROUT(JF)%TIMID) ) 
CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VAROUT(JF)%TIMID, 'long_name','time') )
CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VAROUT(JF)%TIMID, 'units',CTIME) )

!===
CALL NCERROR( NF90_DEF_VAR(VAROUT(JF)%NCID, VAROUT(JF)%CVNAME, NF90_FLOAT, &
              (/LONID,LATID,TIMEID/), VAROUT(JF)%VARID,DEFLATE_LEVEL=CMF_OUTPUT%NDLEVEL ),     &
              'Creating Variable')

CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VAROUT(JF)%VARID, 'long_name', TRIM(VAROUT(JF)%CVLNAME)) )
CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VAROUT(JF)%VARID, 'units',     TRIM(VAROUT(JF)%CVUNITS)) )
CALL NCERROR( NF90_PUT_ATT(VAROUT(JF)%NCID, VAROUT(JF)%VARID, '_FillValue',CMF_PARAMS%RMIS) )

CALL NCERROR( NF90_ENDDEF(VAROUT(JF)%NCID) )

!=== put lon lat info ===
CALL NCERROR ( NF90_INQ_VARID(VAROUT(JF)%NCID,'lon',VARID),'getting id' )
CALL NCERROR( NF90_PUT_VAR(VAROUT(JF)%NCID,VARID,D1LON))

CALL NCERROR ( NF90_INQ_VARID(VAROUT(JF)%NCID,'lat',VARID),'getting id' )
CALL NCERROR( NF90_PUT_VAR(VAROUT(JF)%NCID,VARID,D1LAT))

WRITE(CMF_FILES%LOGNAM,*) 'CFILE: ',TRIM(VAROUT(JF)%CFILE),' CVAR:',TRIM(VAROUT(JF)%CVNAME),&
                ' CLNAME: ',TRIM(VAROUT(JF)%CVLNAME),' CUNITS: ',TRIM(VAROUT(JF)%CVUNITS)
WRITE(CMF_FILES%LOGNAM,*) 'OPEN IN UNIT: ',VAROUT(JF)%NCID
#endif
END SUBROUTINE CREATE_OUTCDF
!==========================================================

END SUBROUTINE CMF_OUTPUT_INIT


!####################################################################
SUBROUTINE CMF_FORCING_INIT(LECMF2LAKEC)
  ! Initialize/open netcdf input 
  ! -- called from "Main Program / Coupler"
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: LECMF2LAKEC   !! Lake coupling: this is currently only used in ECMWF
  
  !================================================
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
  
  WRITE(CMF_FILES%LOGNAM,*) "CMF::FORCING_INIT: Initialize runoff/sediment forcing file (only for netCDF)" 
  IF(  CMF_FORCING%LINPCDF ) THEN
#ifdef UseCDF_CMF
    CALL CMF_FORCING_INIT_CDF
#endif
  ENDIF
  IF( CMF_FORCING%LINTERP ) THEN
    IF( CMF_FORCING%LITRPCDF )THEN
#ifdef UseCDF_CMF
      IF(PRESENT(LECMF2LAKEC)) THEN
        CALL CMF_INPMAT_INIT_CDF(LECMF2LAKEC)
      ELSE
        CALL CMF_INPMAT_INIT_CDF
      ENDIF
#endif
    ELSE
      CALL CMF_INPMAT_INIT_BIN
    ENDIF
  ENDIF 
  
  WRITE(CMF_FILES%LOGNAM,*) "CMF::FORCING_INIT: end" 
  
  CONTAINS
  !==========================================================
  !+ CMF_FORCING_INIT_CDF : open netCDF forcing
  !+ CMF_INPMAT_INIT_CDF      :  open runcrunoff interporlation matrix (inpmat)
  !+ CMF_INPMAT_INIT_BIN      :  open runoff interporlation matrix (inpmat)
  !==========================================================
#ifdef UseCDF_CMF
  SUBROUTINE CMF_FORCING_INIT_CDF
  USE CMF_VARS_MOD,            ONLY: KMINSTAIN, KMINSTART, KMINEND
  USE CMF_UTILS_MOD,           ONLY: NCERROR,   DATE2MIN
  USE CMF_VARS_MOD,            ONLY: ROFCDF, SEDCDF
  USE NETCDF
  IMPLICIT NONE

  !* Local Variables 
  INTEGER(KIND=JPIM)              :: NTIMEID,NCDFSTP
  INTEGER(KIND=JPIM)              :: KMINENDIN
  !================================================
  !*** 1. calculate KMINSTAINP (start KMIN for forcing)
  KMINSTAIN=DATE2MIN(CMF_FORCING%SYEARIN*10000+CMF_FORCING%SMONIN*100+CMF_FORCING%SDAYIN,CMF_FORCING%SHOURIN*100)
  
  !*** 2. Initialize Type for Runoff CDF:
  ROFCDF%CNAME=TRIM(CMF_FORCING%CROFCDF)
  ROFCDF%CVAR(1)=TRIM(CMF_FORCING%CVNROF)
  ROFCDF%CVAR(2)=TRIM(CMF_FORCING%CVNSUB)
  IF ( .not. CMF_OPTIONS%LROSPLIT ) THEN
    ROFCDF%CVAR(2)="NONE"
    ROFCDF%NVARID(2)=-1
  ENDIF
  IF ( .NOT.  CMF_OPTIONS%LWEVAP  ) THEN
    ROFCDF%CVAR(3)="NONE"
    ROFCDF%NVARID(3)=-1
  ENDIF 
  
  ROFCDF%NSTART=KMINSTAIN
  WRITE(CMF_FILES%LOGNAM,*) "CMF::FORCING_INIT_CDF:", TRIM(ROFCDF%CNAME), TRIM(ROFCDF%CVAR(1))
  
  !*** 3. Open netCDF ruoff file
  CALL NCERROR( NF90_OPEN(TRIM(ROFCDF%CNAME),NF90_NOWRITE,ROFCDF%NCID),'OPENING :'//ROFCDF%CNAME )
  CALL NCERROR( NF90_INQ_VARID(ROFCDF%NCID,TRIM(ROFCDF%CVAR(1)),ROFCDF%NVARID(1)) )
  
  IF ( CMF_OPTIONS%LROSPLIT ) THEN
    CALL NCERROR( NF90_INQ_VARID(ROFCDF%NCID,ROFCDF%CVAR(2),ROFCDF%NVARID(2)) )
  ENDIF 
  IF ( CMF_OPTIONS%LWEVAP ) THEN
    CALL NCERROR( NF90_INQ_VARID(ROFCDF%NCID,ROFCDF%CVAR(3),ROFCDF%NVARID(3)) )
  ENDIF
  CALL NCERROR( NF90_INQ_DIMID(ROFCDF%NCID,'time',NTIMEID),'GETTING TIME ID FORCING RUNOFF')
  CALL NCERROR( NF90_INQUIRE_DIMENSION(NCID=ROFCDF%NCID,DIMID=NTIMEID,LEN=NCDFSTP),'GETTING TIME LENGTH')
  
  WRITE(CMF_FILES%LOGNAM,*) "CMF::FORCING_INIT_CDF: CNAME,NCID,VARID", TRIM(ROFCDF%CNAME),ROFCDF%NCID,ROFCDF%NVARID(1)
  
  !*** 4. check runoff forcing time 
  IF ( KMINSTART .LT. KMINSTAIN ) THEN 
    WRITE(CMF_FILES%LOGNAM,*) "Run start earlier than forcing data", TRIM(ROFCDF%CNAME), KMINSTART, KMINSTAIN
    STOP 9
  ENDIF
  
  KMINENDIN=KMINSTAIN + NCDFSTP*INT(CMF_CONFIG%DTIN /60,JPIM)
  IF ( KMINEND .GT. KMINENDIN  ) THEN 
    WRITE(CMF_FILES%LOGNAM,*) "Run end later than forcing data", TRIM(ROFCDF%CNAME), KMINEND, KMINENDIN
    STOP 9
  ENDIF
  
  END SUBROUTINE CMF_FORCING_INIT_CDF
#endif
  !==========================================================
  !+
  !+
  !+
  !==========================================================
#ifdef UseCDF_CMF
  SUBROUTINE CMF_INPMAT_INIT_CDF(LECMF2LAKEC)
  USE CMF_VARS_MOD,            ONLY: INPX,INPA,INPY,INPNI,INPXI,INPYI,INPAI
  USE CMF_VARS_MOD,            ONLY: NSEQMAX
  USE CMF_UTILS_MOD,           ONLY: INQUIRE_FID, NCERROR, mapD2vecD, mapI2vecI
  USE NETCDF
  IMPLICIT NONE
  
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: LECMF2LAKEC  !! for lake coupling: currently only used in ECMWF
  
  INTEGER(KIND=JPIM),ALLOCATABLE  :: I2TMP(:,:,:)
  REAL(KIND=JPRB),ALLOCATABLE     :: D2TMP(:,:,:)
  
  INTEGER(KIND=JPIM)              :: INPI
  INTEGER(KIND=JPIM)              :: VARID
  INTEGER(KIND=JPIM)              :: ISTATUS,VDIMIDS(1)
  INTEGER(KIND=JPIM)              :: NCID               ! output netCDF output file ID

  ! SAVE for OpenMP
  INTEGER(KIND=JPIM),SAVE         :: IX,IY,ILEV
  REAL(KIND=JPRB),SAVE            :: ZTMP
!$OMP THREADPRIVATE               (IY,ILEV,ZTMP)
  !================================================
  !*** 1. allocate input matrix variables
  WRITE(CMF_FILES%LOGNAM,*) 'NX, NY, INPN =', CMF_CONFIG%NX, CMF_CONFIG%NY, CMF_CONFIG%INPN
  ALLOCATE( INPX(NSEQMAX,CMF_CONFIG%INPN),INPY(NSEQMAX,CMF_CONFIG%INPN),INPA(NSEQMAX,CMF_CONFIG%INPN) )
  
  !*** 2. Read Input Matrix
  WRITE(CMF_FILES%LOGNAM,*) 'INPUT MATRIX netCDF', CMF_FORCING%CINPMAT 
  
  CALL NCERROR (NF90_OPEN(CMF_FORCING%CINPMAT ,NF90_NOWRITE,NCID),'opening '//TRIM(CMF_FORCING%CINPMAT ) )
  !** input matrix area
  ALLOCATE( D2TMP(CMF_CONFIG%NX,CMF_CONFIG%NY,CMF_CONFIG%INPN) )
  WRITE(CMF_FILES%LOGNAM,*)'INIT_MAP: inpa:',TRIM(CMF_FORCING%CINPMAT )
  CALL NCERROR ( NF90_INQ_VARID(NCID,'inpa',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,D2TMP,(/1,1,1/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,CMF_CONFIG%INPN/)),'reading data' ) 
  DO INPI=1, CMF_CONFIG%INPN
    CALL mapD2vecD(D2TMP(:,:,INPI:INPI),INPA(:,INPI:INPI))
  END DO
  DEALLOCATE( D2TMP )
  
  !** input matrix IXIN
  ALLOCATE( I2TMP(CMF_CONFIG%NX,CMF_CONFIG%NY,CMF_CONFIG%INPN) )
  
  WRITE(CMF_FILES%LOGNAM,*)'INIT_MAP: inpx:',TRIM(CMF_FORCING%CINPMAT )
  CALL NCERROR ( NF90_INQ_VARID(NCID,'inpx',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,I2TMP,(/1,1,1/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,CMF_CONFIG%INPN/)),'reading data' ) 
  DO INPI=1, CMF_CONFIG%INPN
    CALL mapI2vecI(I2TMP(:,:,INPI:INPI),INPX(:,INPI:INPI))
  END DO
  
  !** input matrix IYIN
  WRITE(CMF_FILES%LOGNAM,*)'INIT_MAP: inpy:',TRIM(CMF_FORCING%CINPMAT )
  CALL NCERROR ( NF90_INQ_VARID(NCID,'inpy',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,I2TMP,(/1,1,1/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,CMF_CONFIG%INPN/)),'reading data' ) 
  DO INPI=1, CMF_CONFIG%INPN
    CALL mapI2vecI(I2TMP(:,:,INPI:INPI),INPY(:,INPI:INPI))
  END DO
  
  DEALLOCATE( I2TMP )
  
  !================================================
  !*** Check if inverse information is available  (only used in ECMWF/IFS v4.07)
  IF(PRESENT(LECMF2LAKEC) .AND. (LECMF2LAKEC .NE. 0)) THEN
  
    ISTATUS = NF90_INQ_VARID(NCID, 'levI', VARID)
    IF ( ISTATUS /= 0 ) THEN
      WRITE(CMF_FILES%LOGNAM,*) "Could not find levI variable in inpmat.nc: inverse interpolation not available"
      INPNI=-1  ! Not available 
    ELSE
      !* Find levels dimension
      CALL NCERROR( NF90_INQUIRE_VARIABLE(NCID,VARID,dimids=VDIMIDS),'getting levI dimensions ')
      CALL NCERROR( NF90_INQUIRE_DIMENSION(NCID,VDIMIDS(1),len=INPNI),'getting time len ')
      WRITE(CMF_FILES%LOGNAM,*) 'Allocating INP*I: NXIN, NYIN, INPNI =', CMF_CONFIG%NXIN, CMF_CONFIG%NYIN, INPNI
      ALLOCATE( INPXI(CMF_CONFIG%NXIN,CMF_CONFIG%NYIN,INPNI),INPYI(CMF_CONFIG%NXIN,CMF_CONFIG%NYIN,INPNI),INPAI(CMF_CONFIG%NXIN,CMF_CONFIG%NYIN,INPNI) )
    
      WRITE(CMF_FILES%LOGNAM,*)'INIT_MAP: inpaI:',TRIM(CMF_FORCING%CINPMAT )
      CALL NCERROR ( NF90_INQ_VARID(NCID,'inpaI',VARID),'getting id' )
      CALL NCERROR ( NF90_GET_VAR(NCID,VARID,INPAI,(/1,1,1/),(/CMF_CONFIG%NXIN,CMF_CONFIG%NYIN,INPNI/)),'reading data' ) 
  
      WRITE(CMF_FILES%LOGNAM,*)'INIT_MAP: inpx:',TRIM(CMF_FORCING%CINPMAT )
      CALL NCERROR ( NF90_INQ_VARID(NCID,'inpxI',VARID),'getting id' )
      CALL NCERROR ( NF90_GET_VAR(NCID,VARID,INPXI,(/1,1,1/),(/CMF_CONFIG%NXIN,CMF_CONFIG%NYIN,INPNI/)),'reading data' ) 
  
      WRITE(CMF_FILES%LOGNAM,*)'INIT_MAP: inpy:',TRIM(CMF_FORCING%CINPMAT )
      CALL NCERROR ( NF90_INQ_VARID(NCID,'inpyI',VARID),'getting id' )
      CALL NCERROR ( NF90_GET_VAR(NCID,VARID,INPYI,(/1,1,1/),(/CMF_CONFIG%NXIN,CMF_CONFIG%NYIN,INPNI/)),'reading data' )
    
      !! We normalize INPAI here as it is used to interpolate flood fraction (Input Area Inversed)
      WRITE(CMF_FILES%LOGNAM,*) 'INPAI normalization'
  !$OMP PARALLEL DO
      DO IX=1,CMF_CONFIG%NXIN
        DO IY=1,CMF_CONFIG%NYIN
          ZTMP=0._JPRB
          DO ILEV=1,INPNI
            ZTMP=ZTMP+INPAI(IX,IY,ILEV)
          ENDDO
          IF (ZTMP > 0._JPRB) THEN
            DO ILEV=1,INPNI
              INPAI(IX,IY,ILEV) = INPAI(IX,IY,ILEV) / ZTMP
            ENDDO
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    ENDIF
  ENDIF
  
  END SUBROUTINE CMF_INPMAT_INIT_CDF
#endif
  !==========================================================
  !+
  !+
  !+
  !==========================================================
  SUBROUTINE CMF_INPMAT_INIT_BIN
    USE CMF_UTILS_MOD,           ONLY: INQUIRE_FID, CONV_END,  CONV_ENDI, mapR2vecD, mapI2vecI
    USE CMF_VARS_MOD,            ONLY: NSEQMAX,INPX,INPY,INPA
  IMPLICIT NONE
  INTEGER(KIND=JPIM)              :: INPI
  INTEGER(KIND=JPIM),ALLOCATABLE  :: I2TMP(:,:)
  REAL(KIND=JPRM),ALLOCATABLE     :: R2TMP(:,:)
  !================================================
  !*** 1. allocate input matrix variables
  WRITE(CMF_FILES%LOGNAM,*) 'NX, NY, INPN =', CMF_CONFIG%NX, CMF_CONFIG%NY, CMF_CONFIG%INPN
  ALLOCATE( INPX(NSEQMAX,CMF_CONFIG%INPN),INPY(NSEQMAX,CMF_CONFIG%INPN),INPA(NSEQMAX,CMF_CONFIG%INPN) )
  
  !*** 2. Read Input Matrix
  WRITE(CMF_FILES%LOGNAM,*) 'INPUT MATRIX binary', CMF_FORCING%CINPMAT
  
  ALLOCATE( I2TMP(CMF_CONFIG%NX,CMF_CONFIG%NY) )
  ALLOCATE( R2TMP(CMF_CONFIG%NX,CMF_CONFIG%NY) )
  
  CMF_FILES%TMPNAM=INQUIRE_FID()
!OPEN(TMPNAM,FILE=CINPMAT,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY*INPN)
!READ(TMPNAM,REC=1) INPX
!READ(TMPNAM,REC=2) INPY
!READ(TMPNAM,REC=3) R2TMP
  
  OPEN(CMF_FILES%TMPNAM,FILE=CMF_FORCING%CINPMAT,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
  DO INPI=1, CMF_CONFIG%INPN
    READ(CMF_FILES%TMPNAM,REC=       INPI) I2TMP
     CALL mapI2vecI(I2TMP,INPX(:,INPI:INPI))
    READ(CMF_FILES%TMPNAM,REC=  CMF_CONFIG%INPN+INPI) I2TMP
     CALL mapI2vecI(I2TMP,INPY(:,INPI:INPI))
    READ(CMF_FILES%TMPNAM,REC=2*CMF_CONFIG%INPN+INPI) R2TMP
     CALL mapR2vecD( R2TMP,INPA(:,INPI:INPI))
  END DO
  
  CLOSE(CMF_FILES%TMPNAM)
  DEALLOCATE(I2TMP,R2TMP)
  
  END SUBROUTINE CMF_INPMAT_INIT_BIN
  !==========================================================
  
  END SUBROUTINE CMF_FORCING_INIT
  !####################################################################
SUBROUTINE CMF_BOUNDARY_INIT
USE CMF_VARS_MOD,             ONLY: NSEQMAX, D2SEALEV

IMPLICIT NONE
!####################################################################
WRITE(CMF_FILES%LOGNAM,*) ""
WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
WRITE(CMF_FILES%LOGNAM,*) "CMF::BOUNDARY_INIT: initialize boundary" 

ALLOCATE( D2SEALEV(NSEQMAX,1) )

IF( CMF_BOUNDARY%LSEALEVCDF )THEN
#ifdef UseCDF_CMF
  CALL CMF_BOUNDARY_INIT_CDF    !! initialize sea level boundary (netCDF only)
#endif
ENDIF

WRITE(CMF_FILES%LOGNAM,*) "CMF::BOUNDARY_INIT: end" 

#ifdef UseCDF_CMF
CONTAINS
!==========================================================
!+ CMF_BOUNDARY_INIT_CDF
!==========================================================
SUBROUTINE CMF_BOUNDARY_INIT_CDF
USE CMF_VARS_MOD,           ONLY: I1NEXT, I2VECTOR
USE CMF_VARS_MOD,           ONLY: KMINSTASL, KMINSTART, KMINEND
USE CMF_UTILS_MOD,           ONLY: NCERROR, INQUIRE_FID, DATE2MIN
USE NETCDF
IMPLICIT NONE
!* Local Variables 
INTEGER(KIND=JPIM)              :: NTIMEID,NCDFSTP
INTEGER(KIND=JPIM)              :: KMINENDSL
INTEGER(KIND=JPIM)              :: IX, IY, IS, ILNK, ISEQ
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
#endif
! ===============================================
!*** 1. calculate KMINSTASL (START KMIN for boundary)
KMINSTASL = DATE2MIN(CMF_BOUNDARY%SYEARSL*10000+CMF_BOUNDARY%SMONSL *100+CMF_BOUNDARY%SDAYSL,CMF_BOUNDARY%SHOURSL*100)

!*** 2. Initialize Type for sea level CDF
SLCDF%CNAME=TRIM(CMF_BOUNDARY%CSEALEVCDF)
SLCDF%CVAR=TRIM(CMF_BOUNDARY%CVNSEALEV)
SLCDF%NSTART=KMINSTASL
WRITE(CMF_FILES%LOGNAM,*) "CMF::BOUNRARY_INIT_CDF:", SLCDF%CNAME, SLCDF%NSTART

!*** Open netCDF sea level File 
CALL NCERROR( NF90_OPEN(SLCDF%CNAME,NF90_NOWRITE,SLCDF%NCID),'OPENING :'//SLCDF%CNAME )
CALL NCERROR( NF90_INQ_VARID(SLCDF%NCID,SLCDF%CVAR,SLCDF%NVARID) )
CALL NCERROR( NF90_INQ_DIMID(SLCDF%NCID,'time',NTIMEID),'GETTING TIME ID Sea Level Boundary')
CALL NCERROR( NF90_INQUIRE_DIMENSION(NCID=SLCDF%NCID,DIMID=NTIMEID,LEN=NCDFSTP),'GETTING TIME LENGTH')
CALL NCERROR( NF90_INQ_DIMID(SLCDF%NCID, 'stations', SLCDF%NSTAID ), 'GETTING STATION ID' ) 
CALL NCERROR( NF90_INQUIRE_DIMENSION(SLCDF%NCID, DIMID=SLCDF%NSTAID, LEN=CMF_BOUNDARY%NCDFSTAT ), 'GETTING STATION NUMBER' )
ALLOCATE( R1SLIN(CMF_BOUNDARY%NCDFSTAT)    ) ! 1D input boundary condition (m)

WRITE(CMF_FILES%LOGNAM,*) "CMF::BOUNDARY_INIT_CDF: CNAME,NCID,VARID", TRIM(SLCDF%CNAME),SLCDF%NCID,SLCDF%NVARID

!*** 4. check sealev forcing time 
IF ( KMINSTART .LT. KMINSTASL ) THEN 
    WRITE(CMF_FILES%LOGNAM,*) "Run start earlier than boundary data", TRIM(SLCDF%CNAME), KMINSTART, KMINSTASL
    STOP 9
ENDIF

KMINENDSL=KMINSTASL + NCDFSTP*INT(CMF_BOUNDARY%DTSL/60,JPIM)
IF ( KMINEND .GT. KMINENDSL  ) THEN 
    WRITE(CMF_FILES%LOGNAM,*) "Run end later than sealev data", TRIM(SLCDF%CNAME), KMINEND, KMINENDSL
    STOP 9
ENDIF

!*** 4. conversion table

!! suggested new mapping format with X Y STATION columns
!! read formated  mapping file and check if at river outlet and in NETCDF
CMF_FILES%TMPNAM=INQUIRE_FID()
OPEN(CMF_FILES%TMPNAM,FILE=CMF_BOUNDARY%CSLMAP,FORM='FORMATTED')
READ(CMF_FILES%TMPNAM,*) CMF_BOUNDARY%NLINKS

WRITE(CMF_FILES%LOGNAM,*) "Dynamic sea level links", CMF_BOUNDARY%NLINKS

ALLOCATE( I2SLMAP(3,CMF_BOUNDARY%NLINKS) ) ! conversion matrix (X Y STATION )
DO ILNK=1, CMF_BOUNDARY%NLINKS
  READ(CMF_FILES%TMPNAM,*) IX, IY, IS
  ! check if links with river outlet cells
  ISEQ=I2VECTOR(IX,IY)
  IF( ISEQ>0 )THEN
    IF( I1NEXT(ISEQ) .NE. -9 ) THEN
        WRITE(CMF_FILES%LOGNAM,*) "Sealev link not at river outlet cell", IX, IY
        STOP 9
    ! check if station index in netcdf
    ELSEIF (IS .LT. 1 .or. IS .GT. CMF_BOUNDARY%NCDFSTAT) THEN
        WRITE(CMF_FILES%LOGNAM,*) "Sealev link outside netcdf index", IS
        STOP 9
    ENDIF
  ELSE
    WRITE(CMF_FILES%LOGNAM,*) "Sealev link outside land grids", IX,IY
    STOP 9
  ENDIF

  I2SLMAP(1,ILNK) = IX
  I2SLMAP(2,ILNK) = IY
  I2SLMAP(3,ILNK) = IS

END DO
CLOSE(CMF_FILES%TMPNAM)

END SUBROUTINE CMF_BOUNDARY_INIT_CDF
#endif

END SUBROUTINE CMF_BOUNDARY_INIT


SUBROUTINE CMF_PROG_INIT
  USE CMF_VARS_MOD,             ONLY: NSEQMAX, NPTHOUT, NPTHLEV
  USE CMF_VARS_MOD,             ONLY: D2RUNOFF,     D2ROFSUB,     &
                                   & P2RIVSTO,     P2FLDSTO,     D2RIVOUT,     D2FLDOUT,     &
                                   & D2RIVOUT_PRE, D2FLDOUT_PRE, D2RIVDPH_PRE, D2FLDSTO_PRE, &
                                   & D1PTHFLW,     D1PTHFLW_PRE, P2GDWSTO,     D2GDWRTN,     &
                                   & P2DAMSTO,     P2DAMINF,     P2LEVSTO ,    D2WEVAP,      &   !! optional
                                   & D2DAMMY,      D2COPY
  IMPLICIT NONE
  !================================================
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
  
  WRITE(CMF_FILES%LOGNAM,*) "CMF::PROG_INIT: prognostic variable initialization"
  
  !*** 1. ALLOCATE 
  ! runoff input
  ALLOCATE( D2RUNOFF(NSEQMAX,1)     )
  ALLOCATE( D2ROFSUB(NSEQMAX,1)     )
  
  ! river+floodplain storage
  ALLOCATE( P2RIVSTO(NSEQMAX,1)     )
  ALLOCATE( P2FLDSTO(NSEQMAX,1)     )
  
  ! discharge calculation
  ALLOCATE( D2RIVOUT(NSEQMAX,1)     )
  ALLOCATE( D2FLDOUT(NSEQMAX,1)     )
  ALLOCATE( D2RIVOUT_PRE(NSEQMAX,1)     )
  ALLOCATE( D2FLDOUT_PRE(NSEQMAX,1)     )
  ALLOCATE( D2RIVDPH_PRE(NSEQMAX,1)     )
  ALLOCATE( D2FLDSTO_PRE(NSEQMAX,1)     )
  
  D2RUNOFF(:,:)=0._JPRB
  D2ROFSUB(:,:)=0._JPRB
  
  P2RIVSTO(:,:)=0._JPRD
  P2FLDSTO(:,:)=0._JPRD
  
  D2RIVOUT(:,:)=0._JPRB
  D2FLDOUT(:,:)=0._JPRB
  D2RIVOUT_PRE(:,:)=0._JPRB
  D2FLDOUT_PRE(:,:)=0._JPRB
  D2RIVDPH_PRE(:,:)=0._JPRB
  D2FLDSTO_PRE(:,:)=0._JPRB
  
  IF( CMF_OPTIONS%LPTHOUT ) THEN  !! additional prognostics for bifurcation scheme
    ALLOCATE( D1PTHFLW(NPTHOUT,NPTHLEV)     )
    ALLOCATE( D1PTHFLW_PRE(NPTHOUT,NPTHLEV) )
    D1PTHFLW(:,:)=0._JPRB
    D1PTHFLW_PRE(:,:)=0._JPRB
  ENDIF
  IF( CMF_OPTIONS%LDAMOUT ) THEN  !! additional prognostics for reservoir operation
    ALLOCATE( P2DAMSTO(NSEQMAX,1)     )
    ALLOCATE( P2DAMINF(NSEQMAX,1)     )
    P2DAMSTO(:,:)=0._JPRD
    P2DAMINF(:,:)=0._JPRD
  ENDIF
  IF( CMF_OPTIONS%LLEVEE ) THEN  !! additional prognostics for LLEVEE
    ALLOCATE( P2LEVSTO(NSEQMAX,1)     )
    P2LEVSTO(:,:)=0._JPRD
  ENDIF
  
  !! Used in ECMWF
  IF( CMF_OPTIONS%LWEVAP ) THEN  !! additional prognostics for LLEVEE
    ALLOCATE( D2WEVAP(NSEQMAX,1)     )
    D2WEVAP(:,:)=0._JPRB
  ENDIF
  
  !! keep these variables even when LGDWDLY is not used.
  ALLOCATE( P2GDWSTO(NSEQMAX,1)     )
  ALLOCATE( D2GDWRTN(NSEQMAX,1)     )
  P2GDWSTO(:,:)=0._JPRD
  D2GDWRTN(:,:)=0._JPRB
  
  !! dammy variable for data handling
  ALLOCATE( D2DAMMY(NSEQMAX,1)) !! Float64/32 switch (Dammy for unused var)
  ALLOCATE( D2COPY(NSEQMAX,1))  !! Float64/32 switch (Dammy for output)
  D2DAMMY(:,:)=0._JPRB
  D2COPY(:,:) =0._JPRB
  
  !============================
  !***  2. set initial water surface elevation to sea surface level
  WRITE(CMF_FILES%LOGNAM,*) 'PROG_INIT: fill channels below downstream boundary'
  CALL STORAGE_SEA_SURFACE
  
  
  WRITE(CMF_FILES%LOGNAM,*) "CMF::PROG_INIT: end"
  
  CONTAINS
  !==========================================================
  !+ STORAGE_SEA_SURFACE: set initial storage, assuming water surface not lower than downstream sea surface elevation
  !+
  !+
  ! ==================================================
  SUBROUTINE STORAGE_SEA_SURFACE
  ! set initial storage, assuming water surface not lower than downstream sea surface elevation
  USE CMF_VARS_MOD,  ONLY: NSEQRIV,  NSEQALL,  I1NEXT
  USE CMF_VARS_MOD,  ONLY: D2DWNELV, D2RIVELV,D2RIVHGT,D2RIVWTH,D2RIVLEN,P2RIVSTOMAX
  IMPLICIT NONE
  ! local variables
  INTEGER(KIND=JPIM)   :: ISEQ, JSEQ
  !
  REAL(KIND=JPRB),SAVE :: DSEAELV, DDPH
  !$OMP THREADPRIVATE    (DSEAELV, DDPH)
  !!=================
  ! For River Mouth Grid
  !$OMP PARALLEL DO
  DO ISEQ=NSEQRIV+1,NSEQALL
    DSEAELV=D2DWNELV(ISEQ,1) !! downstream boundary elevation
  
    !! set initial water level to sea level if river bed is lower than sea level
    DDPH=MAX( DSEAELV-D2RIVELV(ISEQ,1),0._JPRB )
    DDPH=MIN( DDPH,D2RIVHGT(ISEQ,1) )
    P2RIVSTO(ISEQ,1)=DDPH*D2RIVLEN(ISEQ,1)*D2RIVWTH(ISEQ,1)
    P2RIVSTO(ISEQ,1)=MIN( P2RIVSTO(ISEQ,1),P2RIVSTOMAX(ISEQ,1) )
    D2RIVDPH_PRE(ISEQ,1)=DDPH
  END DO
  !$OMP END PARALLEL DO
  
  !! For Usual River Grid (from downstream to upstream). OMP cannot be applied
  DO ISEQ=NSEQRIV,1, -1
    JSEQ=I1NEXT(ISEQ)
    DSEAELV=D2RIVELV(JSEQ,1)+D2RIVDPH_PRE(JSEQ,1)
  
    !! set initial water level to sea level if river bed is lower than sea level
    DDPH=MAX( DSEAELV-D2RIVELV(ISEQ,1),0._JPRB )
    DDPH=MIN( DDPH,D2RIVHGT(ISEQ,1) )
  
    P2RIVSTO(ISEQ,1)=DDPH*D2RIVLEN(ISEQ,1)*D2RIVWTH(ISEQ,1)
    P2RIVSTO(ISEQ,1)=MIN(  P2RIVSTO(ISEQ,1),P2RIVSTOMAX(ISEQ,1) )
    D2RIVDPH_PRE(ISEQ,1)=DDPH
  END DO
  
    
  ! old version before v4.02 (too slow)
  !DO ISEQ=1, NSEQALL
  !  JSEQ=ISEQ
  !  DO WHILE( I1NEXT(JSEQ)>0 )
  !    KSEQ=JSEQ
  !    JSEQ=I1NEXT(KSEQ)
  !  END DO
  !
  !  DSEAELV=D2DWNELV(JSEQ,1) !! downstream boundary elevation
  !  !! set initial water level to sea level if river bed is lower than sea level
  !  DDPH=MAX( DSEAELV-D2RIVELV(ISEQ,1),0._JPRB )
  !  DDPH=MIN( DDPH,D2RIVHGT(ISEQ,1) )
  !  P2RIVSTO(ISEQ,1)=DDPH*D2RIVLEN(ISEQ,1)*D2RIVWTH(ISEQ,1)
  !END DO
      
  END SUBROUTINE STORAGE_SEA_SURFACE
  ! ==================================================
  
  END SUBROUTINE CMF_PROG_INIT
  !####################################################################
  
  
  
  
  
  
  !####################################################################
  SUBROUTINE CMF_DIAG_INIT
  
  USE CMF_VARS_MOD,        ONLY: NSEQMAX,NPTHOUT,NPTHLEV
  USE CMF_VARS_MOD,       ONLY: D2RIVINF, D2RIVDPH, D2RIVVEL, D2FLDINF, D2FLDDPH, D2FLDFRC, D2FLDARE, &
                              & D2PTHOUT, D2PTHINF, D2SFCELV, D2OUTFLW, D2STORGE, D2OUTINS, D2LEVDPH, &
                              & D1PTHFLWSUM,   D2WEVAPEX
  USE CMF_VARS_MOD,       ONLY: D2RIVOUT_oAVG, D2FLDOUT_oAVG, D2OUTFLW_oAVG, D2RIVVEL_oAVG, D2PTHOUT_oAVG, &
                              & D2GDWRTN_oAVG, D2RUNOFF_oAVG, D2ROFSUB_oAVG, D1PTHFLW_oAVG, D2WEVAPEX_oAVG,&
                              & D2DAMINF_oAVG, D2STORGE_oMAX, D2OUTFLW_oMAX, D2RIVDPH_oMAX, NADD_out
  USE CMF_VARS_MOD,       ONLY: D2RIVOUT_aAVG, D2FLDOUT_aAVG, D2OUTFLW_aAVG, D2RIVVEL_aAVG, D2PTHOUT_aAVG, &
                              & D2GDWRTN_aAVG, D2RUNOFF_aAVG, D2ROFSUB_aAVG, D1PTHFLW_aAVG, D2WEVAPEX_aAVG,&
                              & D2DAMINF_aAVG, D2STORGE_aMAX, D2OUTFLW_aMAX, D2RIVDPH_aMAX, NADD_adp,      &
                              & D1PTHFLWSUM_aAVG
  IMPLICIT NONE
  !================================================
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
  
  WRITE(CMF_FILES%LOGNAM,*) "CMF::DIAG_INIT: initialize diagnostic variables"
  
  !*** 1. snapshot 2D diagnostics
  ALLOCATE(D2RIVINF(NSEQMAX,1))
  ALLOCATE(D2RIVDPH(NSEQMAX,1))
  ALLOCATE(D2RIVVEL(NSEQMAX,1))
  ALLOCATE(D2FLDINF(NSEQMAX,1))
  ALLOCATE(D2FLDDPH(NSEQMAX,1))
  ALLOCATE(D2FLDFRC(NSEQMAX,1))
  ALLOCATE(D2FLDARE(NSEQMAX,1))
  ALLOCATE(D2PTHOUT(NSEQMAX,1))
  ALLOCATE(D2PTHINF(NSEQMAX,1))
  ALLOCATE(D2SFCELV(NSEQMAX,1))
  ALLOCATE(D2OUTFLW(NSEQMAX,1))
  ALLOCATE(D2STORGE(NSEQMAX,1))
  D2RIVINF(:,:)=0._JPRB
  D2RIVDPH(:,:)=0._JPRB
  D2RIVVEL(:,:)=0._JPRB
  D2FLDINF(:,:)=0._JPRB
  D2FLDDPH(:,:)=0._JPRB
  D2FLDFRC(:,:)=0._JPRB
  D2FLDARE(:,:)=0._JPRB
  D2PTHOUT(:,:)=0._JPRB
  D2PTHINF(:,:)=0._JPRB
  D2SFCELV(:,:)=0._JPRB
  D2OUTFLW(:,:)=0._JPRB
  D2STORGE(:,:)=0._JPRB
  
  ALLOCATE(D1PTHFLWSUM(NPTHOUT))
  D1PTHFLWSUM(:)=0._JPRB
  
  IF ( CMF_OPTIONS%LLEVEE  )THEN
    ALLOCATE(D2LEVDPH(NSEQMAX,1))
    D2LEVDPH(:,:)=0._JPRB
  ENDIF
  IF ( CMF_OPTIONS%LWEVAP  )THEN
    ALLOCATE(D2WEVAPEX(NSEQMAX,1))
    D2WEVAPEX(:,:)=0._JPRB
  ENDIF
  IF ( CMF_OPTIONS%LOUTINS  )THEN
    ALLOCATE(D2OUTINS (NSEQMAX,1))
    D2OUTINS (:,:)=0._JPRB
  ENDIF
  
  !============================
  !*** 2a. time-average 2D diagnostics for adaptive time step
  NADD_adp=0
  
  ALLOCATE(D2RIVOUT_aAVG(NSEQMAX,1))
  ALLOCATE(D2FLDOUT_aAVG(NSEQMAX,1))
  ALLOCATE(D2OUTFLW_aAVG(NSEQMAX,1))
  ALLOCATE(D2RIVVEL_aAVG(NSEQMAX,1))
  ALLOCATE(D2PTHOUT_aAVG(NSEQMAX,1))
  ALLOCATE(D2GDWRTN_aAVG(NSEQMAX,1))
  ALLOCATE(D2RUNOFF_aAVG(NSEQMAX,1))
  ALLOCATE(D2ROFSUB_aAVG(NSEQMAX,1))
  D2RIVOUT_aAVG(:,:)=0._JPRB
  D2FLDOUT_aAVG(:,:)=0._JPRB
  D2OUTFLW_aAVG(:,:)=0._JPRB
  D2RIVVEL_aAVG(:,:)=0._JPRB
  D2PTHOUT_aAVG(:,:)=0._JPRB
  D2GDWRTN_aAVG(:,:)=0._JPRB
  D2RUNOFF_aAVG(:,:)=0._JPRB
  D2ROFSUB_aAVG(:,:)=0._JPRB
  
  IF ( CMF_OPTIONS%LDAMOUT ) THEN
    ALLOCATE(D2DAMINF_aAVG(NSEQMAX,1))
    D2DAMINF_aAVG(:,:)=0._JPRB
  ENDIF
  IF ( CMF_OPTIONS%LWEVAP ) THEN
    ALLOCATE(D2WEVAPEX_aAVG(NSEQMAX,1))
    D2WEVAPEX_aAVG(:,:)=0._JPRB
  ENDIF
  
  !*** 2b time-average 1D Diagnostics (bifurcation channel) for adaptive time step
  ALLOCATE(D1PTHFLW_aAVG(NPTHOUT,NPTHLEV))
  ALLOCATE(D1PTHFLWSUM_aAVG(NPTHOUT))
  D1PTHFLW_aAVG(:,:)  = 0._JPRB 
  D1PTHFLWSUM_aAVG(:) = 0._JPRB 
  
  !*** 2c. Maximum 2D Diagnostics 
  
  ALLOCATE(D2STORGE_aMAX(NSEQMAX,1))
  ALLOCATE(D2OUTFLW_aMAX(NSEQMAX,1))
  ALLOCATE(D2RIVDPH_aMAX(NSEQMAX,1))
  D2STORGE_aMAX(:,:)=0._JPRB
  D2OUTFLW_aMAX(:,:)=0._JPRB
  D2RIVDPH_aMAX(:,:)=0._JPRB
  
  !============
  !*** 3a. time-average 2D diagnostics for output
  NADD_out=0
  
  ALLOCATE(D2RIVOUT_oAVG(NSEQMAX,1))
  ALLOCATE(D2FLDOUT_oAVG(NSEQMAX,1))
  ALLOCATE(D2OUTFLW_oAVG(NSEQMAX,1))
  ALLOCATE(D2RIVVEL_oAVG(NSEQMAX,1))
  ALLOCATE(D2PTHOUT_oAVG(NSEQMAX,1))
  ALLOCATE(D2GDWRTN_oAVG(NSEQMAX,1))
  ALLOCATE(D2RUNOFF_oAVG(NSEQMAX,1))
  ALLOCATE(D2ROFSUB_oAVG(NSEQMAX,1))
  D2RIVOUT_oAVG(:,:)=0._JPRB
  D2FLDOUT_oAVG(:,:)=0._JPRB
  D2OUTFLW_oAVG(:,:)=0._JPRB
  D2RIVVEL_oAVG(:,:)=0._JPRB
  D2PTHOUT_oAVG(:,:)=0._JPRB
  D2GDWRTN_oAVG(:,:)=0._JPRB
  D2RUNOFF_oAVG(:,:)=0._JPRB
  D2ROFSUB_oAVG(:,:)=0._JPRB
  
  IF ( CMF_OPTIONS%LDAMOUT ) THEN
    ALLOCATE(D2DAMINF_oAVG(NSEQMAX,1))
    D2DAMINF_oAVG(:,:)=0._JPRB
  ENDIF
  IF ( CMF_OPTIONS%LWEVAP ) THEN
    ALLOCATE(D2WEVAPEX_oAVG(NSEQMAX,1))
    D2WEVAPEX_oAVG(:,:)=0._JPRB
  ENDIF
  
  !*** ab time-average 1D Diagnostics (bifurcation channel)
  ALLOCATE(D1PTHFLW_oAVG(NPTHOUT,NPTHLEV))
  D1PTHFLW_oAVG(:,:) = 0._JPRB 
  
  !*** 3c. Maximum 2D Diagnostics 
  
  ALLOCATE(D2STORGE_oMAX(NSEQMAX,1))
  ALLOCATE(D2OUTFLW_oMAX(NSEQMAX,1))
  ALLOCATE(D2RIVDPH_oMAX(NSEQMAX,1))
  D2STORGE_oMAX(:,:)=0._JPRB
  D2OUTFLW_oMAX(:,:)=0._JPRB
  D2RIVDPH_oMAX(:,:)=0._JPRB
  
  WRITE(CMF_FILES%LOGNAM,*) "CMF::DIAG_INIT: end"
  
  END SUBROUTINE CMF_DIAG_INIT
!####################################################################

!####################################################################
  SUBROUTINE CMF_RESTART_INIT
    ! read restart file
    ! -- call from CMF_DRV_INIT
    USE CMF_VARS_MOD,   ONLY: P2RIVSTO,    P2FLDSTO,    D2RIVOUT,    D2FLDOUT,    P2GDWSTO, &
                            & D2RIVOUT_PRE,D2FLDOUT_PRE,D2RIVDPH_PRE,D2FLDSTO_PRE,&
                            & D1PTHFLW,    D1PTHFLW_PRE, &
                            & P2DAMSTO,    P2LEVSTO      !!! added
    IMPLICIT NONE
    ! ===========
    P2RIVSTO(:,:)=0._JPRD
    P2FLDSTO(:,:)=0._JPRD
    D2RIVOUT(:,:)=0._JPRB
    D2FLDOUT(:,:)=0._JPRB
    
    D2RIVOUT_PRE(:,:)=0._JPRB
    D2FLDOUT_PRE(:,:)=0._JPRB
    D2RIVDPH_PRE(:,:)=0._JPRB
    D2FLDSTO_PRE(:,:)=0._JPRB
    
    IF( CMF_OPTIONS%LPTHOUT )THEN
      D1PTHFLW(:,:)=0._JPRB
      D1PTHFLW_PRE(:,:)=0._JPRB
    ENDIF
    IF( CMF_OPTIONS%LDAMOUT ) then
      P2DAMSTO(:,:)=0._JPRD   !!! added LDAMOUT
    ENDIF
    IF( CMF_OPTIONS%LLEVEE ) then
      P2LEVSTO(:,:)=0._JPRD   !!! added LLEVEE
    ENDIF
    IF( CMF_OPTIONS%LGDWDLY ) then
      P2GDWSTO(:,:)=0._JPRD
    ENDIF
    
    IF ( CMF_RESTART%LRESTCDF ) THEN
      CALL READ_REST_CDF
    ELSE
      CALL READ_REST_BIN
    ENDIF
    
    IF( CMF_OPTIONS%LSTOONLY )THEN          !!  storage only restart
      D2FLDSTO_PRE(:,:)=P2FLDSTO(:,:)
    ENDIF
    
    CONTAINS
    !==========================================================
    !+ READ_REST_BIN
    !+ READ_REST_CDF
    !+
    !==========================================================
    SUBROUTINE READ_REST_BIN
    USE CMF_VARS_MOD,             ONLY: NSEQMAX,NPTHOUT, NPTHLEV
    USE CMF_UTILS_MOD,           ONLY: INQUIRE_FID, mapR2vecD
    IMPLICIT NONE
    !*** LOCAL
    INTEGER(KIND=JPIM)              :: RIREC
    REAL(KIND=JPRD)                 :: P2VEC(NSEQMAX,1)
    REAL(KIND=JPRM)                 :: R1PTH(NPTHOUT,NPTHLEV)
    REAL(KIND=JPRD)                 :: P1PTH(NPTHOUT,NPTHLEV)
    CHARACTER(LEN=256)              :: CFILE
    !================================================
    CFILE=TRIM(CMF_RESTART%CRESTSTO)
    WRITE(CMF_FILES%LOGNAM,*)'READ_REST: read restart binary: ', TRIM(CFILE)
    
    CMF_FILES%TMPNAM=INQUIRE_FID()
    
    IF( CMF_RESTART%LRESTDBL )THEN
      OPEN(CMF_FILES%TMPNAM,FILE=CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8*CMF_CONFIG%NX*CMF_CONFIG%NY)
    ELSE
      OPEN(CMF_FILES%TMPNAM,FILE=CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
    ENDIF
    
    RIREC=0
    CALL READ_BIN_MAP(P2VEC,CMF_FILES%TMPNAM,RIREC)
     P2RIVSTO=P2VEC
    CALL READ_BIN_MAP(P2VEC,CMF_FILES%TMPNAM,RIREC)
     P2FLDSTO=P2VEC
    
    !! additional restart data for optional schemes
    IF ( .not. CMF_OPTIONS%LSTOONLY )THEN           !! default restart with previous t-step outflw
      CALL READ_BIN_MAP(P2VEC,CMF_FILES%TMPNAM,RIREC)
       D2RIVOUT_PRE=P2VEC
      CALL READ_BIN_MAP(P2VEC,CMF_FILES%TMPNAM,RIREC)
       D2FLDOUT_PRE=P2VEC
      CALL READ_BIN_MAP(P2VEC,CMF_FILES%TMPNAM,RIREC)
       D2RIVDPH_PRE=P2VEC
      CALL READ_BIN_MAP(P2VEC,CMF_FILES%TMPNAM,RIREC)
       D2FLDSTO_PRE=P2VEC
    ENDIF
    IF ( CMF_OPTIONS%LGDWDLY ) THEN
      CALL READ_BIN_MAP(P2VEC,CMF_FILES%TMPNAM,RIREC)
       P2GDWSTO=P2VEC
    ENDIF
    IF ( CMF_OPTIONS%LDAMOUT ) THEN      !!! added LDAMOUT
      CALL READ_BIN_MAP(P2VEC,CMF_FILES%TMPNAM,RIREC)
       P2DAMSTO=P2VEC
    ENDIF
    IF ( CMF_OPTIONS%LLEVEE ) THEN      !!! added LLEVEE
      CALL READ_BIN_MAP(P2VEC,CMF_FILES%TMPNAM,RIREC)
       P2LEVSTO=P2VEC
    ENDIF
    CLOSE(CMF_FILES%TMPNAM)
    
    IF( CMF_OPTIONS%LPTHOUT )THEN
      IF( .not. CMF_OPTIONS%LSTOONLY )THEN
        CFILE=TRIM(CMF_RESTART%CRESTSTO)//'.pth'
        WRITE(CMF_FILES%LOGNAM,*)'READ_REST: read restart binary: ', TRIM(CFILE)
    
        IF( CMF_RESTART%LRESTDBL )THEN
          OPEN(CMF_FILES%TMPNAM,FILE=CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8*NPTHOUT*NPTHLEV)
          READ(CMF_FILES%TMPNAM,REC=1) P1PTH
          CLOSE(CMF_FILES%TMPNAM)
          D1PTHFLW_PRE(:,:)=P1PTH(:,:)
        ELSE
          OPEN(CMF_FILES%TMPNAM,FILE=CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NPTHOUT*NPTHLEV)
          READ(CMF_FILES%TMPNAM,REC=1) R1PTH
          CLOSE(CMF_FILES%TMPNAM)
          D1PTHFLW_PRE(:,:)=R1PTH(:,:)
        ENDIF
      ELSE
        D1PTHFLW_PRE(:,:)=0._JPRB
      ENDIF
    ENDIF
    
    END SUBROUTINE READ_REST_BIN
    
    SUBROUTINE READ_BIN_MAP(P2VAR,TNAM,IREC)
    USE CMF_UTILS_MOD,      ONLY: mapP2vecP
    USE CMF_VARS_MOD,       ONLY: NSEQMAX
    IMPLICIT NONE
    REAL(KIND=JPRD)            :: P2VAR(NSEQMAX,1)
    INTEGER(KIND=JPIM)         :: TNAM, IREC
    !* local
    REAL(KIND=JPRM)            :: R2TEMP(CMF_CONFIG%NX,CMF_CONFIG%NY)
    REAL(KIND=JPRD)            :: P2TEMP(CMF_CONFIG%NX,CMF_CONFIG%NY)
    !=================
    IREC=IREC+1
    
    !=== Double Precision Restart ===
    IF( CMF_RESTART%LRESTDBL )THEN
      READ(TNAM,REC=IREC) P2TEMP
    !=== Single Precision Restart (convert to double precision once) ===
    ELSE
      READ(TNAM,REC=IREC) R2TEMP
      P2TEMP=R2TEMP
    ENDIF
    CALL mapP2vecP(P2TEMP,P2VAR)
    
    !=================
    END SUBROUTINE READ_BIN_MAP
    !======
    !+
    !+
    !+
    !==========================================================
    SUBROUTINE READ_REST_CDF
#ifdef UseCDF_CMF
    USE NETCDF
    USE CMF_VARS_MOD,      ONLY: NPTHOUT, NPTHLEV, PTH_UPST, PTH_DOWN
    USE CMF_UTILS_MOD,    ONLY: NCERROR, mapP2vecP, mapP2vecD
    IMPLICIT NONE
    ! local variables
    INTEGER(KIND=JPIM)    ::  NCID,VARID
    INTEGER(KIND=JPIM)    ::  IPTH
    CHARACTER(LEN=256)    ::  CFILE
    REAL(KIND=JPRD)       ::  P2TEMP(CMF_CONFIG%NX,CMF_CONFIG%NY), P1PTH(NPTHOUT,NPTHLEV)  !! NetCDF restart is in Double Precision
    !================================================
    CFILE=TRIM(CMF_RESTART%CRESTSTO)
    WRITE(CMF_FILES%LOGNAM,*)'READ_REST: read restart netcdf: ', TRIM(CFILE)
    
    CALL NCERROR( NF90_OPEN(CFILE,NF90_NOWRITE,NCID), 'OPENING '//CFILE)
    
    CALL NCERROR( NF90_INQ_VARID(NCID,'rivsto',VARID))
    CALL NCERROR( NF90_GET_VAR(NCID,VARID,P2TEMP,(/1,1,1/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,1/) ) )
    CALL mapP2vecP(P2TEMP,P2RIVSTO)
    
    CALL NCERROR( NF90_INQ_VARID(NCID,'fldsto',VARID))
    CALL NCERROR( NF90_GET_VAR(NCID,VARID,P2TEMP,(/1,1,1/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,1/) ) )
    CALL mapP2vecP(P2TEMP,P2FLDSTO)
    
    
    IF( .NOT. CMF_OPTIONS%LSTOONLY )THEN
      CALL NCERROR( NF90_INQ_VARID(NCID,'rivout_pre',VARID))
      CALL NCERROR( NF90_GET_VAR(NCID,VARID,P2TEMP,(/1,1,1/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,1/) ) )
      CALL mapP2vecD(P2TEMP,D2RIVOUT_PRE)
      D2RIVOUT=D2RIVOUT_PRE
    
      CALL NCERROR( NF90_INQ_VARID(NCID,'fldout_pre',VARID))
      CALL NCERROR( NF90_GET_VAR(NCID,VARID,P2TEMP,(/1,1,1/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,1/) ) )
      CALL mapP2vecD(P2TEMP,D2FLDOUT_PRE)
      D2FLDOUT=D2FLDOUT_PRE
    
      CALL NCERROR( NF90_INQ_VARID(NCID,'rivdph_pre',VARID))
      CALL NCERROR( NF90_GET_VAR(NCID,VARID,P2TEMP,(/1,1,1/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,1/) ) )
      CALL mapP2vecD(P2TEMP,D2RIVDPH_PRE)
    
      CALL NCERROR( NF90_INQ_VARID(NCID,'fldsto_pre',VARID))
      CALL NCERROR( NF90_GET_VAR(NCID,VARID,P2TEMP,(/1,1,1/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,1/) ) )
      CALL mapP2vecD(P2TEMP,D2FLDSTO_PRE)
    ENDIF
    
    IF ( CMF_OPTIONS%LGDWDLY ) THEN
      CALL NCERROR( NF90_INQ_VARID(NCID,'gdwsto',VARID))
      CALL NCERROR( NF90_GET_VAR(NCID,VARID,P2TEMP,(/1,1,1/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,1/) ) )
      CALL mapP2vecP(P2TEMP,P2GDWSTO)
    ENDIF
    
    IF ( CMF_OPTIONS%LDAMOUT ) THEN    !!! added
      CALL NCERROR( NF90_INQ_VARID(NCID,'damsto',VARID))
      CALL NCERROR( NF90_GET_VAR(NCID,VARID,P2TEMP,(/1,1,1/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,1/) ) )
      CALL mapP2vecP(P2TEMP,P2DAMSTO)
    ENDIF
    
    IF ( CMF_OPTIONS%LLEVEE ) THEN    !!! added
      CALL NCERROR( NF90_INQ_VARID(NCID,'levsto',VARID))
      CALL NCERROR( NF90_GET_VAR(NCID,VARID,P2TEMP,(/1,1,1/),(/CMF_CONFIG%NX,CMF_CONFIG%NY,1/) ) )
      CALL mapP2vecP(P2TEMP,P2LEVSTO)
    ENDIF
    
    IF ( CMF_OPTIONS%LPTHOUT .AND. .NOT. CMF_OPTIONS%LSTOONLY ) THEN
      CALL NCERROR( NF90_INQ_VARID(NCID,'pthflw_pre',VARID))
      CALL NCERROR( NF90_GET_VAR(NCID,VARID,P1PTH,(/1,1,1/),(/NPTHOUT,NPTHLEV,1/) ) )
      DO IPTH=1,NPTHOUT
        IF (PTH_UPST(IPTH)>0 .AND. PTH_DOWN(IPTH)>0 ) THEN
          D1PTHFLW_PRE(IPTH,:)=P1PTH(IPTH,:)
        ELSE
          D1PTHFLW_PRE(IPTH,:)=0._JPRB
        ENDIF
      END DO
    ENDIF
    
    CALL NCERROR( NF90_CLOSE(NCID) )
    
#endif
    END SUBROUTINE READ_REST_CDF
    !==========================================================
    
    END SUBROUTINE CMF_RESTART_INIT
    !####################################################################
    


!####################################################################
  SUBROUTINE CMF_DAMOUT_INIT
    ! reed setting from CMF_DAM%CDAMFILE
    USE CMF_UTILS_MOD,          ONLY: NCERROR, INQUIRE_FID, DATE2MIN
    USE CMF_VARS_MOD,           ONLY: DamID, DamName, DamIX, DamIY, DamLon, DamLat, upreal, Qf, Qn
    USE CMF_VARS_MOD,           ONLY: NDAM,ISYYYY,DamYear, DamStat, DamSeq, FldVol, ConVol, EmeVol, NorVol, AdjVol, Qa, R_VolUpa, I1DAM,DamYear
    USE CMF_VARS_MOD,           ONLY: NSEQMAX,NDAMX,IDAM
    USE CMF_VARS_MOD,           ONLY: I2VECTOR,I1NEXT,I2MASK,NSEQALL,P2RIVSTO,P2FLDSTO,P2DAMSTO,P2DAMINF,NPTHOUT,PTH_DOWN,PTH_UPST,NPTHLEV,PTH_ELV

    IMPLICIT NONE
    INTEGER(KIND=JPIM)         :: NDAMFILE
    INTEGER(KIND=JPIM)         :: ISEQ, JSEQ
    INTEGER(KIND=JPIM)         :: IX, IY
    REAL(KIND=JPRB)            :: FldVol_mcm, ConVol_mcm, TotVol_mcm !! from file in Million Cubic Metter
    REAL(KIND=JPRB)            :: Qsto, Vyr
    
    INTEGER(KIND=JPIM)         :: IPTH, ILEV, ISEQP, JSEQP
    !####################################################################
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
    WRITE(CMF_FILES%LOGNAM,*) "CMF::DAMOUT_INIT: initialize dam", trim(CMF_DAM%CDAMFILE) 
    
    !==========
    NDAMFILE=INQUIRE_FID()
    OPEN(NDAMFILE,FILE=CMF_DAM%CDAMFILE,STATUS="OLD")
    READ(NDAMFILE,*) NDAM
    READ(NDAMFILE,*)        !! skip header
    
    WRITE(CMF_FILES%LOGNAM,*) "CMF::DAMOUT_INIT: number of dams", NDAM
    
    !! === ALLOCATE ===
    !! from CMF_DAM%CDAMFILE
    ALLOCATE(DamID(NDAM),DamName(NDAM))
    ALLOCATE(DamIX(NDAM),DamIY(NDAM),DamLon(NDAM),DamLat(NDAM))
    ALLOCATE(upreal(NDAM))
    ALLOCATE(Qf(NDAM),Qn(NDAM))
    ALLOCATE(DamYear(NDAM),DamStat(NDAM))
    
    !! calculate from CMF_DAM%CDAMFILE
    ALLOCATE(DamSeq(NDAM))
    ALLOCATE(FldVol(NDAM),ConVol(NDAM),EmeVol(NDAM),NorVol(NDAM))
    
    !! for outflw stability
    ALLOCATE(AdjVol(NDAM),Qa(NDAM))
    
    !! H22scheme parameter (FldVol/Upreal)
    ALLOCATE(R_VolUpa(NDAM))
    
    !! dam map, dam variable
    ALLOCATE(I1DAM(NSEQMAX))
    !! =================
    DamSeq(:) = CMF_PARAMS%IMIS
    DamStat(:)= CMF_PARAMS%IMIS
    I1DAM(:)=0
    NDAMX=0
    !! read dam parameters
    DO IDAM = 1, NDAM
      IF(  CMF_DAM%LDAMYBY) THEN
        READ(NDAMFILE,*) DamID(IDAM), DamName(IDAM), DamLat(IDAM), DamLon(IDAM), upreal(IDAM), &
         DamIX(IDAM), DamIY(IDAM), FldVol_mcm, ConVol_mcm, TotVol_mcm, Qn(IDAM), Qf(IDAM), DamYear(IDAM)
      ELSE
        READ(NDAMFILE,*) DamID(IDAM), DamName(IDAM), DamLat(IDAM), DamLon(IDAM), upreal(IDAM), &
         DamIX(IDAM), DamIY(IDAM), FldVol_mcm, ConVol_mcm, TotVol_mcm, Qn(IDAM), Qf(IDAM)
      ENDIF
    
      !! storage parameter --- from Million Cubic Meter to m3
      FldVol(IDAM) = FldVol_mcm * 1.E6                  ! Flood control storage capacity: exclusive for flood control
      ConVol(IDAM) = ConVol_mcm * 1.E6
    
      EmeVol(IDAM) = ConVol(IDAM) + FldVol(IDAM) * 0.95     ! storage to start emergency operation
    
      IX=DamIX(IDAM)
      IY=DamIY(IDAM)
      IF (IX<=0 .or. IX > CMF_CONFIG%NX .or. IY<=0 .or. IY > CMF_CONFIG%NY ) cycle
    
      ISEQ=I2VECTOR(IX,IY)
      IF( I1NEXT(ISEQ)==-9999 .or. ISEQ<=0 ) cycle
      NDAMX=NDAMX+1
    
      DamSeq(IDAM) =ISEQ
      DamStat(IDAM)=2
    
      I1DAM(ISEQ)=1
      I2MASK(ISEQ,1)=2   !! reservoir grid. skipped for adaptive time step
    
      IF( CMF_DAM%LDAMH22 )THEN    !! Hanazaki 2022 scheme 
        NorVol(IDAM)   = ConVol(IDAM) * 0.5    ! normal storage
        R_VolUpa(NDAM) = FldVol(IDAM) * 1.E-6 / upreal(IDAM)
    
      ELSE  !! Yamazaki&Funato scheme (paper in prep)
        Vyr =Qn(IDAM)*(365.*24.*60.*60.)                    !! Annual inflow -> assume dry period inflow is 1/8 of annual flow 
        Qsto=(ConVol(IDAM)*0.7+Vyr/4.)/(180.*24.*60.*60.)   !! possible mean outflow in dry period (6month, ConVol*0.7 + Inflow)
        Qn(IDAM)=min(Qn(IDAM),Qsto)*1.5                     !! Outflow at normal volume (*1.5 is parameter to decide outflw balance)
    
        AdjVol(IDAM)=ConVol(IDAM)  + FldVol(IDAM)*0.1       !! AdjVol is for outflow stability (result is not so sensitive)
        Qa(IDAM)=( Qn(IDAM)+Qf(IDAM) )*0.5                  !! Qa is also for stability
      ENDIF
    
      !! Year-by-Year scheme. If dam is not yet constructed
      IF( CMF_DAM%LDAMYBY )THEN
        IF( ISYYYY==DamYear(IDAM) )THEN
          DamStat(IDAM)=1   !! new this year
        ELSEIF( ISYYYY<DamYear(IDAM) .and. DamYear(IDAM)>0 )THEN
          DamStat(IDAM)=-1  !! not yet activated
          I1DAM(ISEQ)=-1
          FldVol(IDAM)=0._JPRB
          ConVol(IDAM)=0._JPRB
        ENDIF
      ENDIF
    
    END DO
    CLOSE(NDAMFILE)
    
    WRITE(CMF_FILES%LOGNAM,*) "CMF::DAMOUT_INIT: allocated dams:", NDAMX 
    !==========
    
    !! mark upstream of dam grid, for applying kinematic wave routine to suppress storage buffer effect.
    DO ISEQ=1, NSEQALL
      IF( I1DAM(ISEQ)<=0 .and. I1NEXT(ISEQ)>0 )THEN !! if target is non-dam grid
        JSEQ=I1NEXT(ISEQ)
        IF( I1DAM(JSEQ)==1 .or. I1DAM(JSEQ)==11 )THEN !! if downstream is dam
          I1DAM(ISEQ)=10            !! mark upstream of dam grid by "10"
          I2MASK(ISEQ,1)=1   !! reservoir upstream grid. skipped for adaptive time step
        ENDIF
      ENDIF
    
      IF( I1DAM(ISEQ)==1 .and. I1NEXT(ISEQ)>0 )THEN !! if target is dam grid
        JSEQ=I1NEXT(ISEQ)
        IF( I1DAM(JSEQ)==1 .or. I1DAM(JSEQ)==11 )THEN !! if downstream is dam
          I1DAM(ISEQ)=11            !! mark upstream of dam grid by "11"
          I2MASK(ISEQ,1)=2   !! reservoir grid (cascading). skipped for adaptive time step
        ENDIF
      ENDIF
    END DO
    
    !! Initialize dam storage
    IF( .not.  CMF_OPTIONS%LRESTART  )THEN  !! Initialize without restart data
      P2DAMSTO(:,1)=0._JPRD
      DO IDAM=1, NDAM
        IF( DamStat(IDAM)==CMF_PARAMS%IMIS )CYCLE !
        ISEQ=DamSeq(IDAM)
        IF( DamStat(IDAM)==-1 )THEN !! Dam not yet constructed
          P2DAMSTO(ISEQ,1)= P2RIVSTO(ISEQ,1)+P2FLDSTO(ISEQ,1)
        ELSE
          P2DAMSTO(ISEQ,1)=P2RIVSTO(ISEQ,1)+P2FLDSTO(ISEQ,1)
          IF( P2DAMSTO(ISEQ,1)<ConVol(IDAM) )THEN  !! If Normal Volume > initial storage, replace
            P2DAMSTO(ISEQ,1)=ConVol(IDAM)
            P2RIVSTO(ISEQ,1)=ConVol(IDAM)  
            P2FLDSTO(ISEQ,1)=0._JPRD
          ENDIF
        ENDIF
      END DO
    ELSE       !! if from restart file
      IF( CMF_DAM%LDAMYBY )THEN   !! for restart with year-by-year option, set damsto for newly constructed dam 
        DO IDAM=1, NDAM
          IF( DamStat(IDAM)==1 )THEN !! Dam newly activated from this year
            ISEQ=DamSeq(IDAM)
            P2DAMSTO(ISEQ,1)=P2RIVSTO(ISEQ,1)+P2FLDSTO(ISEQ,1)
            IF( CMF_DAM%LiVnorm .and. P2DAMSTO(ISEQ,1)<ConVol(IDAM) )THEN  !! If Initialize Vnormal option & Vnor>Riv+Fld sto, replace
              P2DAMSTO(ISEQ,1)=ConVol(IDAM)
              P2RIVSTO(ISEQ,1)=ConVol(IDAM)  
              P2FLDSTO(ISEQ,1)=0._JPRD
            ENDIF
          ENDIF
        END DO
      ENDIF
    ENDIF
    
    !! Initialize dam inflow
    DO ISEQ=1, NSEQALL
      P2DAMINF(ISEQ,1)=0._JPRD
    END DO
    
    !! Stop bifurcation at dam & dam-upstream grids
    IF( CMF_OPTIONS%LPTHOUT )THEN
      DO IPTH=1, NPTHOUT
        ISEQP=PTH_UPST(IPTH)
        JSEQP=PTH_DOWN(IPTH)
        IF( ISEQP<=0 .or. JSEQP<=0) CYCLE
        IF( I1DAM(ISEQP)>0 .or. I1DAM(JSEQP)>0 )THEN
          DO ILEV=1, NPTHLEV
            PTH_ELV(IPTH,ILEV)=1.E20  !! no bifurcation
          END DO
        ENDIF
      END DO
    ENDIF
    
    END SUBROUTINE CMF_DAMOUT_INIT
    !####################################################################
    
  !####################################################################
SUBROUTINE CMF_TRACER_INIT
  ! tracer variable initialization
  USE CMF_VARS_MOD, ONLY: VTRACE,TBUFF,ITRACE,P2TRCSTO,D2TRCDNS,D2TRCOUT,D2TRCINP,D1TRCPFLW,D2TRCPOUT,D2TRCPOUT_oAVG,D2TRCDNS_oAVG,D2TRCOUT_oAVG
  USE CMF_VARS_MOD, ONLY: Nadd_out
  USE CMF_VARS_MOD, ONLY: NSEQMAX,NPTHOUT,NSEQALL,NSEQALL,NSEQMAX
  IMPLICIT NONE
  
  CHARACTER(LEN=256)              :: CTMP
  INTEGER(KIND=JPIM)              :: J,J0
  !####################################################################
  WRITE(CMF_FILES%LOGNAM,*) ""
  WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
  WRITE(CMF_FILES%LOGNAM,*) "CMF::TRACER_INIT: initialize tracer" 
  
  ALLOCATE(VTRACE(CMF_TRACER%NTRACE ))  !! tracer data variables
  ALLOCATE(TBUFF(CMF_CONFIG%NXIN,CMF_CONFIG%NYIN,CMF_TRACER%NTRACE ))  !! tracer input buffer file
  
  !==========
  WRITE(CMF_FILES%LOGNAM,*) "  Check Tracer Names"
  ITRACE=0
  J0=1
  DO J=1,LEN(TRIM( CMF_TRACER%CTRCNAM))
    IF( (J>J0) .AND. ( CMF_TRACER%CTRCNAM(J:J) .EQ. ',') ) THEN
      CTMP=TRIM(ADJUSTL( CMF_TRACER%CTRCNAM(J0:J-1)))
      IF (LEN(CTMP) > 0 ) THEN
        ITRACE=ITRACE+1
        VTRACE(ITRACE)%TRCNAME=CTMP
        WRITE(CMF_FILES%LOGNAM,*) ITRACE, trim(VTRACE(ITRACE)%TRCNAME)
      ENDIF
      J0=J+1
    ENDIF
  ENDDO
  ! Last one 
  IF ( J0 <= LEN(TRIM( CMF_TRACER%CTRCNAM)) ) THEN
    J=LEN(TRIM( CMF_TRACER%CTRCNAM))
    CTMP=TRIM(ADJUSTL( CMF_TRACER%CTRCNAM(J0:J)))
    IF (LEN(CTMP) > 0 ) THEN
        ITRACE=ITRACE+1
        VTRACE(ITRACE)%TRCNAME=CTMP
        WRITE(CMF_FILES%LOGNAM,*) ITRACE, trim(VTRACE(ITRACE)%TRCNAME)
    ENDIF
  ENDIF 
  IF( ITRACE/=CMF_TRACER%NTRACE )THEN
    WRITE(CMF_FILES%LOGNAM,*) 'ERROR: Tracer Name number do not match with NTRACE', trim( CMF_TRACER%CTRCNAM), CMF_TRACER%NTRACE
    stop
  ENDIF
  
  !==========
  WRITE(CMF_FILES%LOGNAM,*) "  Check Tracer Input File Prefix"
  ITRACE=0
  J0=1
  DO J=1,LEN(TRIM(CMF_TRACER%CTRCPRE))
    IF( (J>J0) .AND. (CMF_TRACER%CTRCPRE(J:J) .EQ. ',') ) THEN
      CTMP=TRIM(ADJUSTL(CMF_TRACER%CTRCPRE(J0:J-1)))
      IF (LEN(CTMP) > 0 ) THEN
        ITRACE=ITRACE+1
        VTRACE(ITRACE)%TRCPRE=CTMP
        WRITE(CMF_FILES%LOGNAM,*) ITRACE, trim(VTRACE(ITRACE)%TRCPRE)
      ENDIF
      J0=J+1
    ENDIF
  ENDDO
  ! Last one 
  IF ( J0 <= LEN(TRIM(CMF_TRACER%CTRCPRE)) ) THEN
    J=LEN(TRIM(CMF_TRACER%CTRCPRE))
    CTMP=TRIM(ADJUSTL(CMF_TRACER%CTRCPRE(J0:J)))
    IF (LEN(CTMP) > 0 ) THEN
        ITRACE=ITRACE+1
        VTRACE(ITRACE)%TRCPRE=CTMP
        WRITE(CMF_FILES%LOGNAM,*) ITRACE, trim(VTRACE(ITRACE)%TRCPRE)
    ENDIF
  ENDIF 
  IF( ITRACE/=CMF_TRACER%NTRACE )THEN
    WRITE(CMF_FILES%LOGNAM,*) 'ERROR: Tracer Prefix number do not match with NTRACE', trim( CMF_TRACER%CTRCNAM), CMF_TRACER%NTRACE
    stop
  ENDIF
  
  !==========
  WRITE(CMF_FILES%LOGNAM,*) "  Allocate Tracer Variables"
  ALLOCATE( P2TRCSTO(NSEQMAX,CMF_TRACER%NTRACE) )
  ALLOCATE( D2TRCDNS(NSEQMAX,CMF_TRACER%NTRACE) )
  ALLOCATE( D2TRCOUT(NSEQMAX,CMF_TRACER%NTRACE) )
  ALLOCATE( D2TRCINP(NSEQMAX,CMF_TRACER%NTRACE) )
  P2TRCSTO(:,:)=0._JPRB
  D2TRCDNS(:,:)=0._JPRB
  D2TRCOUT(:,:)=0._JPRB
  D2TRCINP(:,:)=0._JPRB
  
  ALLOCATE( D1TRCPFLW(NPTHOUT,CMF_TRACER%NTRACE) )
  ALLOCATE( D2TRCPOUT(NSEQMAX,CMF_TRACER%NTRACE) )
  D1TRCPFLW(:,:)=0._JPRB
  D2TRCPOUT(:,:)=0._JPRB
  
  Nadd_out=0._JPRB
  ALLOCATE( D2TRCDNS_oAVG(NSEQMAX,CMF_TRACER%NTRACE) )
  ALLOCATE( D2TRCOUT_oAVG(NSEQMAX,CMF_TRACER%NTRACE) )
  ALLOCATE( D2TRCPOUT_oAVG(NSEQMAX,CMF_TRACER%NTRACE) )
  D2TRCDNS_oAVG(:,:)=0._JPRB
  D2TRCOUT_oAVG(:,:)=0._JPRB
  D2TRCPOUT_oAVG(:,:)=0._JPRB
  END SUBROUTINE CMF_TRACER_INIT
  !####################################################################
  !+  
!@@@@@@ TRACER Output @@@@@@
!####################################################################
  SUBROUTINE CMF_TRACER_OUTPUT_INIT
    ! Initialize tracer output module (create/open files)
    ! -- Called from CMF_DRV_INIT
    USE CMF_VARS_MOD,            ONLY: ISYYYY, ISMM,   ISDD,   ISHOUR, ISMIN,NSEQMAX
    USE CMF_VARS_MOD,             ONLY: REGIONTHIS,VTRACE
    USE CMF_UTILS_MOD,           ONLY: INQUIRE_FID
    USE CMF_VARS_MOD,             ONLY: NVARSOUT_TRACER,VAROUT_TRACER
    IMPLICIT NONE
    !* Local variables 
    CHARACTER(LEN=256)              :: CTIME
    INTEGER(KIND=JPIM)              :: JF,ITRACE
    CHARACTER(LEN=256)              :: CVNAMES(NVARSOUT_TRACER)
    !================================================
    WRITE(CMF_FILES%LOGNAM,*) ""
    WRITE(CMF_FILES%LOGNAM,*) "!---------------------!"
    
    WRITE(CMF_FILES%LOGNAM,*) "CMF::TRACER_OUTPUT_INIT: set variable names"
    !! 
    NVARSOUT_TRACER=0
    DO ITRACE=1, CMF_TRACER%NTRACE
      NVARSOUT_TRACER=NVARSOUT_TRACER+1
      CVNAMES(NVARSOUT_TRACER)=TRIM(VTRACE(ITRACE)%TRCNAME)//'_sto'    !! storage
    
      NVARSOUT_TRACER=NVARSOUT_TRACER+1
      CVNAMES(NVARSOUT_TRACER)=TRIM(VTRACE(ITRACE)%TRCNAME)//'_out'    !! flux (main channel)
    
      NVARSOUT_TRACER=NVARSOUT_TRACER+1
      CVNAMES(NVARSOUT_TRACER)=TRIM(VTRACE(ITRACE)%TRCNAME)//'_dns'    !  density
    
      IF( CMF_TRACER%LTRCBIF )THEN
        NVARSOUT_TRACER=NVARSOUT_TRACER+1
        CVNAMES(NVARSOUT_TRACER)=TRIM(VTRACE(ITRACE)%TRCNAME)//'_bifout'  !! net flux (bifurcation channel out-in)
      ENDIF
    
    END DO
    
    ALLOCATE(VAROUT_TRACER(NVARSOUT_TRACER))
    WRITE(CTIME,'(A14,I4.4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)') 'seconds since ',ISYYYY,'-',ISMM,'-',ISDD,' ',ISHOUR,":",ISMIN
    
    !* Loop on variables and create files 
    DO JF=1,NVARSOUT_TRACER
      WRITE(CMF_FILES%LOGNAM,*) "Creating output for variable:", TRIM( CVNAMES(JF) )
      VAROUT_TRACER(JF)%CVNAME=CVNAMES(JF)
      VAROUT_TRACER(JF)%BINID=INQUIRE_FID()
    
      IF( CMF_TRACER%LOUTVEC )THEN   !!  1D land only output
        VAROUT_TRACER(JF)%CFILE=TRIM(CMF_TRACER%COUTDIR)//TRIM(VAROUT_TRACER(JF)%CVNAME)//TRIM(CMF_TRACER%COUTTAG)//TRIM(CMF_PARAMS%CSUFVEC)
        OPEN(VAROUT_TRACER(JF)%BINID,FILE=VAROUT_TRACER(JF)%CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NSEQMAX)
        WRITE(CMF_FILES%LOGNAM,*) "output file opened in unit: ", TRIM(VAROUT_TRACER(JF)%CFILE), VAROUT_TRACER(JF)%BINID
      ELSE                   !!  2D default map output
        IF( REGIONTHIS==1 )THEN
          VAROUT_TRACER(JF)%CFILE=TRIM(CMF_TRACER%COUTDIR)//TRIM(VAROUT_TRACER(JF)%CVNAME)//TRIM(CMF_TRACER%COUTTAG)//TRIM(CMF_PARAMS%CSUFBIN)
          WRITE(CMF_FILES%LOGNAM,*) "  -- ", TRIM(VAROUT_TRACER(JF)%CFILE)
          OPEN(VAROUT_TRACER(JF)%BINID,FILE=VAROUT_TRACER(JF)%CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
          WRITE(CMF_FILES%LOGNAM,*) "output file opened in unit: ", TRIM(VAROUT_TRACER(JF)%CFILE), VAROUT_TRACER(JF)%BINID
        ENDIF
      ENDIF
    END DO
    
    IRECOUT=0  ! Initialize Output record to 1 (shared in netcdf & binary)
    
    END SUBROUTINE CMF_TRACER_OUTPUT_INIT
    !####################################################################
!####################################################################
!+
!+
!+
!####################################################################
    SUBROUTINE CMF_TRACER_RESTART_INIT
      ! read restart file
      USE CMF_UTILS_MOD,           ONLY: INQUIRE_FID, mapP2vecP
      USE CMF_VARS_MOD,            ONLY: NSEQMAX,P2TRCSTO
      IMPLICIT NONE
      !*** LOCAL
      REAL(KIND=JPRM)                 :: R2TEMP(CMF_CONFIG%NX,CMF_CONFIG%NY)
      REAL(KIND=JPRD)                 :: P2TEMP(CMF_CONFIG%NX,CMF_CONFIG%NY)
      REAL(KIND=JPRD)                 :: P2VEC(NSEQMAX,1)
      CHARACTER(LEN=256)              :: CFILE
      INTEGER(KIND=JPIM)              :: ITRACE
      !================================================
      CFILE=TRIM(CMF_TRACER%CRESTTRC )
      WRITE(CMF_FILES%LOGNAM,*)'CMF::TRACER_RESTART_INIT: read restart binary: ', TRIM(CFILE)
      
      CMF_FILES%TMPNAM=INQUIRE_FID()
      
      IF( CMF_TRACER%LRESTDBL )THEN
        OPEN(CMF_FILES%TMPNAM,FILE=CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8*CMF_CONFIG%NX*CMF_CONFIG%NY)
        DO ITRACE=1, CMF_TRACER%NTRACE
          READ(CMF_FILES%TMPNAM,REC=ITRACE) P2TEMP
          CALL mapP2vecP(P2TEMP,P2VEC)
          P2TRCSTO(:,ITRACE)=P2VEC(:,1)
        END DO
      ELSE
        OPEN(CMF_FILES%TMPNAM,FILE=CFILE,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
        DO ITRACE=1, CMF_TRACER%NTRACE
          READ(CMF_FILES%TMPNAM,REC=ITRACE) R2TEMP
          P2TEMP(:,:)=R2TEMP(:,:)
          CALL mapP2vecP(P2TEMP,P2VEC)
          P2TRCSTO(:,ITRACE)=P2VEC(:,1)
        END DO
      ENDIF
      
      END SUBROUTINE CMF_TRACER_RESTART_INIT
      !####################################################################
!+
!==========================================================
subroutine cmf_sed_init
  implicit none
  call sediment_vars_init
  call sediment_map_init
  call sediment_input_init
  if ( CMF_OPTIONS%LOUTPUT  ) then
    call sediment_output_init
  endif
      
  call sediment_restart_init
            
  END subroutine cmf_sed_init
  subroutine sediment_vars_init
    use CMF_VARS_MOD,             only: d2bedout, d2netflw, d2seddep, &
    d2bedout_avg, d2netflw_avg,   &
    d2sedout, d2sedcon, d2sedinp, &
    d2sedout_avg, d2sedinp_avg, d2layer, &
    d2sedv, d2sedv_avg, d2depv,   &
    step_sed, SEDCDF
    use CMF_VARS_MOD,             only: NSEQMAX
      implicit none
      INTEGER(KIND=JPIM)    ::  NCID,VARID

      ! Allocate and initialize the SEDCDF structure
      SEDCDF%NCID = 0
      SEDCDF%NVARID(:) = 0
      SEDCDF%CVAR(:) = ""
      SEDCDF%CNAME = ""
      SEDCDF%NSTART = 0

      if ( mod(CMF_SED%sedDT,CMF_CONFIG%DT ) /= 0 ) then
        write(CMF_FILES%LOGNAM,*) 'sedDT ',CMF_SED%sedDT,'is not a multiple of DT',CMF_CONFIG%DT 
        stop
      endif
      print *, "CMF::SEDIMENT_VARS_INIT: sedDT,DT", CMF_SED%sedDT,CMF_CONFIG%DT
      step_sed = int(CMF_SED%sedDT/CMF_CONFIG%DT )
          
      allocate(d2sedv(NSEQMAX,CMF_SED%nsed ,6))
      d2sedv(:,:,:) = 0._JPRB
      d2sedout => d2sedv(:,:,1)
      d2sedcon => d2sedv(:,:,2)
      d2sedinp => d2sedv(:,:,3)
      d2bedout => d2sedv(:,:,4)
      d2netflw => d2sedv(:,:,5)
      d2layer  => d2sedv(:,:,6)
          
      allocate(d2depv(NSEQMAX,CMF_SED%totlyrnum,CMF_SED%nsed ))
      d2depv(:,:,:) = 0._JPRB
      d2seddep => d2depv
          
      allocate(d2sedv_avg(NSEQMAX,CMF_SED%nsed ,4))
      d2sedv_avg(:,:,:) = 0._JPRB
      d2sedout_avg => d2sedv_avg(:,:,1)
      d2sedinp_avg => d2sedv_avg(:,:,2)
      d2bedout_avg => d2sedv_avg(:,:,3)
      d2netflw_avg => d2sedv_avg(:,:,4)
    end subroutine sediment_vars_init
            !==================================
    subroutine sediment_map_init
      use CMF_UTILS_MOD,           only: mapR2vecD, splitchar,INQUIRE_FID
      use CMF_VARS_MOD,            only: d2sedfrc, setVel,sDiam,REGIONTHIS,MPI_COMM_CAMA,NSEQMAX,NSEQALL
      use CMF_SEDPAR_MOD,          only: calc_settingVelocity
      USE CMF_UTILS_MOD  ,ONLY: NCERROR
      USE NETCDF
      implicit none
      integer(kind=JPIM)              :: i, ierr, ised, iseq, tmpnam
      real(kind=JPRM)                 :: r2temp(CMF_CONFIG%NX,CMF_CONFIG%NY), sTmp1(CMF_SED%nsed)
      character(len=256)              :: ctmp(20)
      INTEGER(KIND=JPIM)    ::  NCID,VARID
       
      !------------------------!
      ! get sediment diameters !
      !------------------------!
      ctmp(:) = '-999'
      call splitchar( CMF_SED%sedD,ctmp)
      ised = 0
      allocate(sDiam(CMF_SED%nsed))
      do i = 1, CMF_SED%nsed
        if ( ctmp(i) /= '-999' ) then
          ised = ised + 1
          read(ctmp(i),*) sDiam(ised)
        endif
      enddo
      if ( ised /= CMF_SED%nsed ) then
        write(CMF_FILES%LOGNAM,*) 'nsed and sedD do not match',ised,CMF_SED%nsed
        stop
      endif
      write(CMF_FILES%LOGNAM,*) ised,' grain sizes: ',sDiam(:)
      
      !----------------------------!
      ! calculate setting velocity !
      !----------------------------!
      allocate(setVel(CMF_SED%nsed))
      setVel(:) = calc_settingVelocity()
      !-----------------------------!
      ! read sediment fraction file !
      !-----------------------------!
      allocate(d2sedfrc(NSEQMAX,CMF_SED%nsed))
      if ( REGIONTHIS == 1 ) then
        CALL NCERROR (NF90_OPEN(CMF_SED%csedfrc,NF90_NOWRITE,NCID),'opening '//TRIM(CMF_SED%csedfrc) )
        print*, 'opening ', TRIM(CMF_SED%csedfrc)
        CALL NCERROR (NF90_INQ_VARID(NCID, 'sedfrc',VARID),'getting id' )
        CALL NCERROR (NF90_GET_VAR(NCID,VARID,r2temp),'reading data' )
        CALL NCERROR (NF90_CLOSE(NCID),'closing '//TRIM(CMF_SED%csedfrc))
      endif
      do ised = 1, CMF_SED%nsed
        if ( REGIONTHIS == 1 )  then
#ifdef UseMPI_CMF
          call MPI_Bcast(r2temp(1,1),CMF_CONFIG%NX*CMF_CONFIG%NY,mpi_real4,0,MPI_COMM_CAMA,ierr)
#endif
          call mapR2vecD(r2temp,d2sedfrc(:,ised))
        endif 
      enddo
      ! adjust if any fractions are negative or if sum is not equal to 1
      if ( CMF_SED%nsed == 1 ) then
        d2sedfrc(:,:) = 1.d0
      else
!$omp parallel do
        do iseq = 1, NSEQALL
          if ( minval(d2sedfrc(iseq,:)) < 0.d0 .or. sum(d2sedfrc(iseq,:)) == 0.d0 ) then
              d2sedfrc(iseq,:) = 1.d0 / dble(CMF_SED%nsed)
          else if ( sum(d2sedfrc(iseq,:)) /= 1.d0 ) then
              d2sedfrc(iseq,:) = d2sedfrc(iseq,:) / sum(d2sedfrc(iseq,:))
          endif
        enddo
!$omp end parallel do
      endif
    end subroutine sediment_map_init
    !==================================
      
    subroutine sediment_input_init
#ifdef UseMPI_CMF
      use CMF_VARS_MOD,             only: MPI_COMM_CAMA
#endif
      use CMF_UTILS_MOD,           only: INQUIRE_FID, mapR2vecD
      use CMF_VARS_MOD,            only: d2slope,NSEQMAX,REGIONTHIS,SEDCDF,KMINSTAIN
      USE CMF_UTILS_MOD,           ONLY: NCERROR, DATE2MIN
      USE NETCDF
      implicit none
      INTEGER(KIND=JPIM)    ::  NCID,VARID
 
      save
      integer                       :: ierr, tmpnam, i
      real(kind=jprm)               :: r2temp(CMF_CONFIG%NX,CMF_CONFIG%NY)
      allocate(d2slope(NSEQMAX,CMF_CONFIG%NLFP))
      if ( REGIONTHIS == 1 ) then
        CALL NCERROR (NF90_OPEN(CMF_SED%cslope,NF90_NOWRITE,NCID),'opening '//TRIM(CMF_SED%cslope) )
        CALL NCERROR (NF90_INQ_VARID(NCID, 'slope',VARID),'getting id' )
        CALL NCERROR (NF90_GET_VAR(NCID,VARID,r2temp),'reading data' )
        CALL NCERROR (NF90_CLOSE(NCID),'closing '//TRIM(CMF_SED%cslope))
      endif
      do i = 1, CMF_CONFIG%NLFP
        if ( REGIONTHIS == 1 ) then
#ifdef UseMPI_CMF
        call MPI_Bcast(r2temp(1,1), CMF_CONFIG%NX*CMF_CONFIG%NY,mpi_real4,0,MPI_COMM_CAMA,ierr)
#endif
        call mapR2vecD(r2temp,d2slope(:,i))
        endif
      enddo
      if ( REGIONTHIS == 1 ) close(CMF_FILES%TMPNAM)
      
      ! Initialize sediment forcing NetCDF file
      if (CMF_SED%LSEDCDF) then
        call sediment_forcing_init
      endif
    end subroutine sediment_input_init
    
    !====================================================
    ! Initialize sediment forcing NetCDF file
    !====================================================
    subroutine sediment_forcing_init
#ifdef UseCDF_CMF
      use CMF_VARS_MOD,            only: KMINSTAIN, KMINSTART, KMINEND, SEDCDF
      use CMF_UTILS_MOD,           only: NCERROR, DATE2MIN
      use NETCDF
      implicit none

      ! Local Variables 
      INTEGER(KIND=JPIM)              :: NTIMEID, NCDFSTP
      INTEGER(KIND=JPIM)              :: KMINENDIN
      CHARACTER(LEN=256)              :: CMSG
      !================================================
      
      ! Calculate start time for sediment forcing
      IF (CMF_SED%SYEARIN > 0) THEN
        KMINSTAIN=DATE2MIN(CMF_SED%SYEARIN*10000+CMF_SED%SMONIN*100+CMF_SED%SDAYIN,CMF_SED%SHOURIN*100)
      ELSE
        ! If not specified, use the same as runoff forcing
        KMINSTAIN=DATE2MIN(CMF_FORCING%SYEARIN*10000+CMF_FORCING%SMONIN*100+CMF_FORCING%SDAYIN,CMF_FORCING%SHOURIN*100)
      ENDIF
      
      ! Initialize Type for Sediment CDF
      SEDCDF%CNAME=TRIM(CMF_SED%CPREPCDF)
      SEDCDF%CVAR(1)=TRIM(CMF_SED%CVNSED)
      SEDCDF%NSTART=KMINSTAIN
      
      WRITE(CMF_FILES%LOGNAM,*) "CMF::SEDIMENT_FORCING_INIT:", TRIM(SEDCDF%CNAME), TRIM(SEDCDF%CVAR(1))
      print*, "CMF::SEDIMENT_FORCING_INIT:", TRIM(SEDCDF%CNAME), TRIM(SEDCDF%CVAR(1))
      ! Open NetCDF sediment file
      CALL NCERROR(NF90_OPEN(TRIM(SEDCDF%CNAME),NF90_NOWRITE,SEDCDF%NCID), &
                  'OPENING SEDIMENT FORCING:'//TRIM(SEDCDF%CNAME))
      
      CALL NCERROR(NF90_INQ_VARID(SEDCDF%NCID,TRIM(SEDCDF%CVAR(1)),SEDCDF%NVARID(1)), &
                  'GETTING VARID FOR '//TRIM(SEDCDF%CVAR(1)))
      
      CALL NCERROR(NF90_INQ_DIMID(SEDCDF%NCID,'time',NTIMEID),'GETTING TIME ID FORCING SEDIMENT')
      CALL NCERROR(NF90_INQUIRE_DIMENSION(NCID=SEDCDF%NCID,DIMID=NTIMEID,LEN=NCDFSTP),'GETTING TIME LENGTH')
      
      WRITE(CMF_FILES%LOGNAM,*) "CMF::SEDIMENT_FORCING_INIT: CNAME,NCID,VARID", &
                               TRIM(SEDCDF%CNAME),SEDCDF%NCID,SEDCDF%NVARID(1)

      ! Check sediment forcing time
      IF (KMINSTART < KMINSTAIN) THEN 
        WRITE(CMF_FILES%LOGNAM,*) "Run start earlier than sediment forcing data", &
                                 TRIM(SEDCDF%CNAME), KMINSTART, KMINSTAIN
        STOP 9
      ENDIF
      
      KMINENDIN=KMINSTAIN + NCDFSTP*INT(CMF_SED%DTIN/60,JPIM)
      IF (KMINEND > KMINENDIN) THEN 
        WRITE(CMF_FILES%LOGNAM,*) "Run end later than sediment forcing data", &
                                 TRIM(SEDCDF%CNAME), KMINEND, KMINENDIN
        STOP 9
      ENDIF
#endif
    end subroutine sediment_forcing_init
                  !==========================================================
                
    subroutine sediment_output_init
            use CMF_UTILS_MOD,           only: INQUIRE_FID,splitchar      
            use CMF_VARS_MOD,            only: REGIONTHIS
            use CMF_VARS_MOD,            only: NVARSOUT_SED,VAROUT_SED
            implicit none
            save
            integer(kind=JPIM)              :: jf, j
            integer(kind=JPIM)              :: nvars, nsetfile
            parameter                         (nvars=30)
            character(len=256)              :: cvnames(nvars), fName 
      
              
            !---------------------------!
            ! get output variable names !
            !---------------------------!
            cvnames(:) = 'none'
            nvarsout_sed = 0
            call splitchar(CMF_SED%csedsout,cvnames)
            do j = 1, nvars
              if ( cvnames(j) /= 'none' ) then
                nvarsout_sed = nvarsout_sed + 1
              endif
            enddo
                
            if ( nvarsout_sed == 0 ) then
              write(CMF_FILES%LOGNAM,*) "cmf::sed_output_init: no output files will be produced!"
              return
            endif
          
            allocate(VAROUT_SED(nvarsout_sed))
                
            !* loop on variables and create files
            do jf=1,nvarsout_sed
              write(CMF_FILES%LOGNAM,*) "creating output for variable:", trim( cvnames(jf) )
              select case (cvnames(jf))
                case ('sedout')
                  VAROUT_SED(jf)%cvname=cvnames(jf)
                  VAROUT_SED(jf)%cvlname='suspended sediment flow'
                  VAROUT_SED(jf)%cvunits='m3/s'
                case ('sedcon')
                  VAROUT_SED(jf)%cvname=cvnames(jf)
                  VAROUT_SED(jf)%cvlname='suspended sediment concentration'
                  VAROUT_SED(jf)%cvunits='m3/m3'
                case ('sedinp')
                  VAROUT_SED(jf)%cvname=cvnames(jf)
                  VAROUT_SED(jf)%cvlname='sediment inflow from land'
                  VAROUT_SED(jf)%cvunits='m3/s'
                case ('bedout')
                  VAROUT_SED(jf)%cvname=cvnames(jf)
                  VAROUT_SED(jf)%cvlname='bedload'
                  VAROUT_SED(jf)%cvunits='m3/s'
                case ('netflw')
                  VAROUT_SED(jf)%cvname=cvnames(jf)
                  VAROUT_SED(jf)%cvlname='net entrainment flow'
                  VAROUT_SED(jf)%cvunits='m3/s'
                case ('layer')
                  VAROUT_SED(jf)%cvname=cvnames(jf)
                  VAROUT_SED(jf)%cvlname='exchange layer volume'
                  VAROUT_SED(jf)%cvunits='m3'
                case default  ! should only be seddep
                  if ( cvnames(jf)(:6) == 'deplyr' ) then
                    VAROUT_SED(jf)%cvname=cvnames(jf)
                    VAROUT_SED(jf)%cvlname='river bed volume (vertical layer)'
                    VAROUT_SED(jf)%cvunits='m3'
                  else
                    write(CMF_FILES%LOGNAM,*) trim(cvnames(jf)), 'Not defined in sediment output init'
                  endif
              end select
              VAROUT_SED(jf)%binid=INQUIRE_FID()
                
              if ( trim(VAROUT_SED(jf)%cvname(:6)) == 'deplyr' ) then
                fName = trim(VAROUT_SED(jf)%cvname)//'_'//trim(CMF_OUTPUT%COUTTAG)
              else
                fName = trim(VAROUT_SED(jf)%cvname)//trim(CMF_OUTPUT%COUTTAG)
              endif
                
              if ( CMF_OUTPUT%LOUTCDF ) then
                  if ( REGIONTHIS==1 ) then
                    call create_outcdf_sed
                  else
                  call create_outbin_sed
                  endif
              endif 
            enddo
            contains  
            subroutine create_outcdf_sed
#ifdef UseCDF_CMF
                          use CMF_VARS_MOD,            only: ISYYYY, ISMM, ISDD, ISHOUR, ISMIN
                          use CMF_VARS_MOD,            only: D1LON, D1LAT, REGIONTHIS
                          use CMF_UTILS_MOD,           only: NCERROR
                          use CMF_VARS_MOD,            only: sDiam
                          USE NETCDF
                          IMPLICIT NONE
                          INTEGER(KIND=JPIM)  :: TIMEID,VARID,LATID,LONID,sedid,ierr
                          character(len=256)         :: ctime    ! Time units string
                                  
                              
                          VAROUT_SED(jf)%irecnc = 1
                          
                          VAROUT_SED(jf)%cfile = trim(CMF_OUTPUT%COUTDIR)//trim(fName)//trim(CMF_PARAMS%CSUFCDF)
                          
                          ! Create netCDF file
                          ierr = nf90_create(VAROUT_SED(jf)%cfile, nf90_netcdf4, VAROUT_SED(jf)%ncid)
                          call NCERROR(ierr, 'creating file:'//trim(VAROUT_SED(jf)%cfile))
                          
                          ! Define dimensions
                          ierr = nf90_def_dim(VAROUT_SED(jf)%ncid, 'time', nf90_unlimited, TIMEID)
                          call NCERROR(ierr, 'defining time dimension')
                          
                          ierr = nf90_def_dim(VAROUT_SED(jf)%ncid, 'lat', CMF_CONFIG%NY, latid)
                          call NCERROR(ierr, 'defining lat dimension')
                          
                          ierr = nf90_def_dim(VAROUT_SED(jf)%ncid, 'lon', CMF_CONFIG%NX, lonid)
                          call NCERROR(ierr, 'defining lon dimension')
                          
                          ierr = nf90_def_dim(VAROUT_SED(jf)%ncid, 'sedD', CMF_SED%nsed, sedid)
                          call NCERROR(ierr, 'defining sedD dimension')
                          
                          ! Define variables
                          ierr = nf90_def_var(VAROUT_SED(jf)%ncid, 'sedD', nf90_double, (/sedid/), VARID)
                          call NCERROR(ierr, 'defining sedD variable')
                          
                          ierr = nf90_put_att(VAROUT_SED(jf)%ncid, varid, 'long_name', 'sediment grain size')
                          call NCERROR(ierr, 'adding sedD long_name attribute')
                          
                          ierr = nf90_put_att(VAROUT_SED(jf)%ncid, varid, 'units', 'meters')
                          call NCERROR(ierr, 'adding sedD units attribute')
                          
                          ! Define lat/lon variables
                          ierr = nf90_def_var(VAROUT_SED(jf)%ncid, 'lat', nf90_float, (/latid/), VARID)
                          call NCERROR(ierr, 'defining lat variable')
                          
                          ierr = nf90_put_att(VAROUT_SED(jf)%ncid, varid, 'long_name', 'latitude')
                          call NCERROR(ierr, 'adding lat long_name attribute')
                          
                          ierr = nf90_put_att(VAROUT_SED(jf)%ncid, varid, 'units', 'degrees_north')
                          call NCERROR(ierr, 'adding lat units attribute')
                          
                          ierr = nf90_def_var(VAROUT_SED(jf)%ncid, 'lon', nf90_float, (/lonid/), VARID)
                          call NCERROR(ierr, 'defining lon variable')
                          
                          ierr = nf90_put_att(VAROUT_SED(jf)%ncid, varid, 'long_name', 'longitude')
                          call NCERROR(ierr, 'adding lon long_name attribute')
                          
                          ierr = nf90_put_att(VAROUT_SED(jf)%ncid, varid, 'units', 'degrees_east')
                          call NCERROR(ierr, 'adding lon units attribute')
                          
                          ! Define time variable
                          write(ctime,'(a14,i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') &
                                'seconds since ',ISYYYY,'-',ISMM,'-',ISDD,' ',ISHOUR,":",ISMIN
                          
                          ierr = nf90_def_var(VAROUT_SED(jf)%ncid, 'time', nf90_double, (/timeid/), VAROUT_SED(jf)%timid)
                          call NCERROR(ierr, 'defining time variable')
                          
                          ierr = nf90_put_att(VAROUT_SED(jf)%ncid, VAROUT_SED(jf)%timid, 'long_name', 'time')
                          call NCERROR(ierr, 'adding time long_name attribute')
                          
                          ierr = nf90_put_att(VAROUT_SED(jf)%ncid, VAROUT_SED(jf)%timid, 'units', trim(ctime))
                          call NCERROR(ierr, 'adding time units attribute')
                          
                          ! Define main variable
                          ierr = nf90_def_var(VAROUT_SED(jf)%ncid, VAROUT_SED(jf)%cvname, nf90_float, &
                                              (/lonid,latid,sedid,timeid/), VAROUT_SED(jf)%varid, &
                                              deflate_level=CMF_OUTPUT%NDLEVEL)
                          call NCERROR(ierr, 'creating variable')
                          
                          ierr = nf90_put_att(VAROUT_SED(jf)%ncid, VAROUT_SED(jf)%varid, 'long_name', trim(VAROUT_SED(jf)%cvlname))
                          call NCERROR(ierr, 'adding variable long_name attribute')
                          
                          ierr = nf90_put_att(VAROUT_SED(jf)%ncid, VAROUT_SED(jf)%varid, 'units', trim(VAROUT_SED(jf)%cvunits))
                          call NCERROR(ierr, 'adding variable units attribute')
                          
                          ierr = nf90_put_att(VAROUT_SED(jf)%ncid, VAROUT_SED(jf)%varid, '_FillValue', CMF_PARAMS%RMIS)

                          call NCERROR(ierr, 'adding variable _FillValue attribute')
                          
                          ! End define mode
                          ierr = nf90_enddef(VAROUT_SED(jf)%ncid)
                          call NCERROR(ierr, 'ending define mode')
                          
                          ! Put dimension data
                          !ierr = nf90_put_var(VAROUT_SED(jf)%ncid, VAROUT_SED(jf)%varid, sDiam)
                          ierr = nf90_inq_varid(VAROUT_SED(jf)%ncid, 'sedD', varid)
                          ierr = nf90_put_var(VAROUT_SED(jf)%ncid, varid, sDiam)    
                          call NCERROR(ierr, 'writing sedD data')
                          
                          ierr = nf90_inq_varid(VAROUT_SED(jf)%ncid, 'lon', varid)
                          ierr = nf90_put_var(VAROUT_SED(jf)%ncid, varid, D1LON)
                          call NCERROR(ierr, 'writing lon data')
                          
                          ierr = nf90_inq_varid(VAROUT_SED(jf)%ncid, 'lat', varid)
                          ierr = nf90_put_var(VAROUT_SED(jf)%ncid, varid, D1LAT)
                          call NCERROR(ierr, 'writing lat data')
                          
                          write(CMF_FILES%LOGNAM,*) 'cfile: ',trim(VAROUT_SED(jf)%cfile),' cvar:',trim(VAROUT_SED(jf)%cvname),&
                                                    ' clname: ',trim(VAROUT_SED(jf)%cvlname),' cunits: ',trim(VAROUT_SED(jf)%cvunits)
                          write(CMF_FILES%LOGNAM,*) 'open in unit: ',VAROUT_SED(jf)%ncid
#endif
                        end subroutine create_outcdf_sed
                                    
                        subroutine create_outbin_sed
                          use CMF_VARS_MOD,             only: NSEQMAX, REGIONALL
                          implicit none
                              
                            if ( CMF_OUTPUT%LOUTVEC ) then
                                  VAROUT_SED(jf)%cfile=trim(CMF_OUTPUT%COUTDIR)//trim(fName)//trim(CMF_PARAMS%CSUFVEC)
                                  open(VAROUT_SED(jf)%binid,file=VAROUT_SED(jf)%cfile,form='unformatted',access='direct',recl=4*NSEQMAX*CMF_SED%nsed )
                                else
                                  if ( REGIONTHIS==1 ) then
                                    VAROUT_SED(jf)%cfile=trim(CMF_OUTPUT%COUTDIR)//trim(fName)//trim(CMF_PARAMS%CSUFBIN)
                                    open(VAROUT_SED(jf)%binid,file=VAROUT_SED(jf)%cfile,form='unformatted',access='direct',recl=4*CMF_CONFIG%NX*CMF_CONFIG%NY* CMF_SED%nsed )
                                  endif
                            endif
                            write(CMF_FILES%LOGNAM,*) "output file opened in unit: ", TRIM(VAROUT_SED(JF)%CFILE), VAROUT_SED(JF)%BINID
                          end subroutine create_outbin_sed
  end subroutine sediment_output_init
        

                  !==========================================================
            !####################################################################
  subroutine sediment_restart_init
    use CMF_VARS_MOD,             only: D2RIVLEN, D2RIVWTH,NSEQALL,NSEQMAX,REGIONTHIS
    use CMF_VARS_MOD,            only: P2RIVSTO,d2layer,d2seddep,d2sedcon
    use CMF_VARS_MOD,             only: d2sedfrc, d2rivsto_pre, &
                                      d2rivout_sed, d2rivvel_sed, sadd_riv, sadd_out
    use CMF_UTILS_MOD,           only: mapR2vecD, inquire_fid

    implicit none
    save
#ifdef UseMPI_CMF
              integer(kind=JPIM)             :: ierr
#endif
              integer(kind=JPIM)             :: ilyr, irec, ised, iseq, tmpnam, nsetfile
              real(kind=JPRM)                :: r2temp(CMF_CONFIG%NX,CMF_CONFIG%NY)
      
      
              if ( CMF_SED%sedrest_infile == "" ) then  ! set layer/bedload if no restart file
      !$omp parallel do
                do iseq = 1, NSEQALL
                  d2layer(iseq,:) = CMF_SED%lyrdph * D2RIVWTH(iseq,1) * D2RIVLEN(iseq,1) * d2sedfrc(iseq,:)
                  do ilyr = 1, CMF_SED%totlyrnum-1
                    d2seddep(iseq,ilyr,:) = d2layer(iseq,:)
                  enddo
                  d2seddep(iseq,CMF_SED%totlyrnum,:) = ( max(10.d0-CMF_SED%lyrdph*CMF_SED%totlyrnum,0.d0) ) * D2RIVWTH(iseq,1) * D2RIVLEN(iseq,1) * d2sedfrc(iseq,:)
                enddo
      !$omp end parallel do
      
              else
                if ( REGIONTHIS == 1 ) then
                  tmpnam = INQUIRE_FID()
                  open(tmpnam,file=CMF_SED%sedrest_infile,form='unformatted',access='direct',recl=4*CMF_CONFIG%NX*CMF_CONFIG%NY)
                endif
                do irec = 1, 2
                  do ised = 1, CMF_SED%nsed
                    if ( REGIONTHIS == 1 ) read(tmpnam,rec=(irec-1)*CMF_SED%nsed+ised) r2temp
#ifdef UseMPI_CMF
                    call MPI_Bcast(r2temp(1,1),CMF_CONFIG%NX*CMF_CONFIG%NY,mpi_real4,0,MPI_COMM_CAMA,ierr)
#endif
                    select case(irec)
                      case (1)
                        call mapR2vecD(r2temp,d2layer(:,ised))
                      case (2)
                        call mapR2vecD(r2temp,d2sedcon(:,ised))
                    end select
                  enddo
                enddo
      
                do irec = 1, CMF_SED%totlyrnum
                  do ised = 1, CMF_SED%nsed
                    if ( REGIONTHIS == 1 ) read(tmpnam,rec=(irec+1)*CMF_SED%nsed+ised) r2temp
#ifdef UseMPI_CMF
                    call MPI_Bcast(r2temp(1,1),CMF_CONFIG%NX*CMF_CONFIG%NY,mpi_real4,0,MPI_COMM_CAMA,ierr)
#endif
                    call mapR2vecD(r2temp,d2seddep(:,irec,ised))
                  enddo
                enddo
                if ( REGIONTHIS == 1 ) close(tmpnam)
                write(CMF_FILES%LOGNAM,*) 'read restart sediment',maxval(d2seddep(:,CMF_SED%totlyrnum,:))
              endif
      
              allocate(d2rivsto_pre(NSEQMAX), d2rivout_sed(NSEQMAX), d2rivvel_sed(NSEQMAX))
              sadd_riv = 0.d0
              sadd_out = 0.d0
              d2rivsto_pre(:) = P2RIVSTO(:,1)
              d2rivout_sed(:) = 0.d0
              d2rivvel_sed(:) = 0.d0
            end subroutine sediment_restart_init

      
      

      !####################################################################


SUBROUTINE CMF_MAKE_INIT(LECMF2LAKEC)
  USE CMF_FLDSTG_MOD, ONLY: CMF_PHYSICS_FLDSTG
  USE CMF_OUTFLW_MOD, ONLY: CMF_CALC_OUTPRE
  USE CMF_OUTPUT_MOD, ONLY: CMF_OUTPUT_WRITE
  IMPLICIT NONE
  INTEGER(KIND=JPIM), OPTIONAL :: LECMF2LAKEC !! for lake coupling, currently only used in ECMWF

  CALL CMF_TIME_INIT
  CALL CMF_RIVMAP_INIT
  CALL CMF_TOPO_INIT
  IF( CMF_OPTIONS%LLEVEE) THEN
    CALL CMF_LEVEE_INIT
  ENDIF
  IF( CMF_OPTIONS%LOUTPUT )THEN
    CALL CMF_OUTPUT_INIT
  ENDIF

  !*** 3b. Initialize forcing data
  IF(PRESENT(LECMF2LAKEC)) THEN
    CALL CMF_FORCING_INIT(LECMF2LAKEC)
  ELSE
    CALL CMF_FORCING_INIT()
  ENDIF

  IF( CMF_OPTIONS%LSEALEV )THEN
    CALL CMF_BOUNDARY_INIT
  ENDIF

  CALL CMF_PROG_INIT
  CALL CMF_DIAG_INIT

  CALL CMF_PHYSICS_FLDSTG

  IF(CMF_OPTIONS%LRESTART)THEN
    CALL CMF_RESTART_INIT
  ENDIF

  IF( CMF_OPTIONS%LDAMOUT )THEN
    CALL CMF_DAMOUT_INIT
  ENDIF

  IF( CMF_OPTIONS%LTRACE )THEN
    CALL CMF_TRACER_INIT
    IF( CMF_OPTIONS%LOUTPUT ) CALL CMF_TRACER_OUTPUT_INIT
    IF( CMF_OPTIONS%LRESTART) CALL CMF_TRACER_RESTART_INIT
  ENDIF

  IF( CMF_OPTIONS%LSEDOUT )THEN
    CALL cmf_sed_init
  ENDIF

  IF( CMF_OPTIONS%LRESTART .AND. CMF_OPTIONS%LSTOONLY )THEN
    WRITE(CMF_FILES%LOGNAM,*) "CMF::DRV_INIT: (5a) set flood stage at initial condition"
    CALL CMF_PHYSICS_FLDSTG
    CALL CMF_CALC_OUTPRE  
  ENDIF
  
  IF ( CMF_OPTIONS%LOUTINI .AND. CMF_OPTIONS%LOUTPUT ) THEN
    WRITE(CMF_FILES%LOGNAM,*) "CMF::DRV_INIT: (5b) write initial condition"
    CALL CMF_OUTPUT_WRITE
  ENDIF
  END SUBROUTINE CMF_MAKE_INIT
  
END MODULE CMF_INIT_MOD
