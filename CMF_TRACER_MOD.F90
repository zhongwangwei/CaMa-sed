MODULE CMF_TRACER_MOD

USE PARKIND1,              ONLY: JPIM,   JPRB,    JPRD,      JPRM
USE CMF_NMLIST_MOD,     ONLY: CMF_FILES, CMF_OPTIONS, CMF_CONFIG, CMF_PARAMS
USE CMF_NMLIST_MOD,     ONLY: CMF_TIME, CMF_MAPS, CMF_FORCING, CMF_BOUNDARY
USE CMF_NMLIST_MOD,     ONLY: CMF_RESTART, CMF_DAM, CMF_LEVEE, CMF_OUTPUT,CMF_TRACER,CMF_SED
USE CMF_VARS_MOD,          ONLY: NSEQALL, NSEQMAX ,P2TRCSTO,P2RIVSTO,P2FLDSTO,D2TRCDNS,D2TRCPOUT
USE CMF_VARS_MOD,          ONLY: NADD_out,ITRACE,D2TRCINP
CONTAINS


!@@@@@@ TRACER PHYSICS @@@@@@
!####################################################################
SUBROUTINE CMF_TRACER_DENSITY
IMPLICIT NONE
! SAVE for OpenMP
INTEGER(KIND=JPIM),SAVE    :: ISEQ
!============================

DO ITRACE=1, CMF_TRACER%NTRACE
!$OMP PARALLEL DO
  DO ISEQ=1, NSEQALL
    D2TRCDNS(ISEQ,ITRACE)=P2TRCSTO(ISEQ,ITRACE) / max( (P2RIVSTO(ISEQ,1)+P2FLDSTO(ISEQ,1)), 1.D-6 )
  END DO
!$OMP END PARALLEL DO
END DO

END SUBROUTINE CMF_TRACER_DENSITY
!####################################################################
!+
!+
!+
!####################################################################
SUBROUTINE CMF_TRACER_FLUX
! Calculate Tracer Physics
USE CMF_VARS_MOD,       ONLY: D2OUTFLW_aAVG, D1PTHFLWSUM_aAVG,NSEQMAX
USE CMF_VARS_MOD,       ONLY: I1NEXT,D2TRCOUT,NPTHOUT,PTH_UPST,PTH_DOWN,D1TRCPFLW,D1PTHFLWSUM_aAVG
USE CMF_VARS_MOD,       ONLY: NSEQALL,ITRACE,P2TRCSTO,NSEQRIV

IMPLICIT NONE
REAL(KIND=JPRD)            :: P2STOOUT(NSEQMAX)                      !! total outflow from a grid     [m3]
REAL(KIND=JPRD)            :: P2TRCINF(NSEQMAX)                      !! 
REAL(KIND=JPRB)            :: D2RATE(NSEQMAX)                        !! outflow correction

! SAVE for OpenMP
INTEGER(KIND=JPIM),SAVE    :: ISEQ, JSEQ, IPTH
REAL(KIND=JPRB),SAVE       :: OUT_R1, OUT_R2, DIUP, DIDW, ISEQP, JSEQP
!$OMP THREADPRIVATE     (JSEQ,OUT_R1, OUT_R2, DIUP, DIDW, ISEQP, JSEQP)
!============================

! ****** 1. calculate flux
D2TRCPOUT(:,:)=0


DO ITRACE=1, CMF_TRACER%NTRACE
  P2TRCINF(:) = 0._JPRD
  P2STOOUT(:) = 0._JPRD
  D2RATE(:)   = 1._JPRB

!$OMP PARALLEL DO
  DO ISEQ=1, NSEQALL
    JSEQ=I1NEXT(ISEQ)
    IF( D2OUTFLW_aAVG(ISEQ,1)>= 0 )THEN
      D2TRCOUT(ISEQ,ITRACE)=D2TRCDNS(ISEQ,ITRACE) * D2OUTFLW_aAVG(ISEQ,1)
    ELSE !! reverse flow
      IF( JSEQ>0 )THEN
        D2TRCOUT(ISEQ,ITRACE)=D2TRCDNS(JSEQ,ITRACE) * D2OUTFLW_aAVG(ISEQ,1)  !! use downstream density
      ELSE
        D2TRCOUT(ISEQ,ITRACE)=0._JPRB  !! from ocean (assume no flux)
      ENDIF
    ENDIF
  END DO
!$OMP END PARALLEL DO

  IF(  CMF_TRACER%LTRCBIF ) THEN
!$OMP PARALLEL DO
    DO IPTH=1, NPTHOUT  
      ISEQP=PTH_UPST(IPTH)
      JSEQP=PTH_DOWN(IPTH)
      IF (ISEQP<=0 .OR. JSEQP<=0 ) CYCLE  !! Avoid calculation outside of domain
    
      IF( D1PTHFLWSUM_aAVG(IPTH)>=0. )THEN
        D1TRCPFLW(IPTH,ITRACE)=D2TRCDNS(ISEQP,ITRACE) * D1PTHFLWSUM_aAVG(IPTH)
      ELSE !! reverse flow
        D1TRCPFLW(IPTH,ITRACE)=D2TRCDNS(JSEQP,ITRACE) * D1PTHFLWSUM_aAVG(IPTH)
      ENDIF
    END DO
!$OMP END PARALLEL DO
  ENDIF

! ****** 2. calculate total outflow from each catchment

!! flux adjustment for mass balance
  !! for normal cells ---------
#ifndef NoAtom_CMF
!$OMP PARALLEL DO
#endif
  DO ISEQ=1, NSEQRIV                                                    !! for normalcells
    JSEQ=I1NEXT(ISEQ) ! next cell's pixel
    OUT_R1 = max(  D2TRCOUT(ISEQ,ITRACE),0._JPRB )
    OUT_R2 = max( -D2TRCOUT(ISEQ,ITRACE),0._JPRB )
    DIUP=OUT_R1*CMF_CONFIG%DT
    DIDW=OUT_R2*CMF_CONFIG%DT
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
    P2STOOUT(ISEQ) = P2STOOUT(ISEQ) + DIUP 
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
    P2STOOUT(JSEQ) = P2STOOUT(JSEQ) + DIDW 
  END DO
#ifndef NoAtom_CMF
!$OMP END PARALLEL DO
#endif
  
  !! for river mouth grids ------------
!$OMP PARALLEL DO
  DO ISEQ=NSEQRIV+1, NSEQALL
    OUT_R1 = max( D2TRCOUT(ISEQ,ITRACE), 0._JPRB )
    P2STOOUT(ISEQ) = P2STOOUT(ISEQ) + OUT_R1*CMF_CONFIG%DT
  END DO
!$OMP END PARALLEL DO

  IF(  CMF_TRACER%LTRCBIF ) THEN
#ifndef NoAtom_CMF
!$OMP PARALLEL DO
#endif
    DO IPTH=1, NPTHOUT  
      ISEQP=PTH_UPST(IPTH)
      JSEQP=PTH_DOWN(IPTH)
      IF (ISEQP<=0 .OR. JSEQP<=0 ) CYCLE  !! Avoid calculation outside of domain
      OUT_R1 = max(  D1TRCPFLW(IPTH,ITRACE),0._JPRB )
      OUT_R2 = max( -D1TRCPFLW(IPTH,ITRACE),0._JPRB )
      DIUP=OUT_R1*CMF_CONFIG%DT
      DIDW=OUT_R2*CMF_CONFIG%DT
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
      P2STOOUT(ISEQP) = P2STOOUT(ISEQP) + DIUP 
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
      P2STOOUT(JSEQP) = P2STOOUT(JSEQP) + DIDW 
    END DO
#ifndef NoAtom_CMF
!$OMP END PARALLEL DO
#endif
  ENDIF

  !! calculate modification rate
!$OMP PARALLEL DO
  DO ISEQ=1, NSEQALL
    IF ( P2STOOUT(ISEQ) > 1.E-8 ) THEN
      D2RATE(ISEQ) = min( P2TRCSTO(ISEQ,ITRACE) * P2STOOUT(ISEQ)**(-1.), 1._JPRD )
    ENDIF
  END DO
!$OMP END PARALLEL DO

  !============================
  !****** 3. modify outflow

  !! normal pixels------
#ifndef NoAtom_CMF
  !$OMP PARALLEL DO  !! No OMP Atomic for bit-identical simulation (set in Mkinclude)
#endif
  DO ISEQ=1, NSEQRIV ! for normal pixels
    JSEQ=I1NEXT(ISEQ)
    IF( D2TRCOUT(ISEQ,ITRACE) >= 0._JPRB )THEN
      D2TRCOUT(ISEQ,ITRACE) = D2TRCOUT(ISEQ,ITRACE)*D2RATE(ISEQ)
    ELSE
      D2TRCOUT(ISEQ,ITRACE) = D2TRCOUT(ISEQ,ITRACE)*D2RATE(JSEQ)
    ENDIF
#ifndef NoAtom_CMF
  !$OMP ATOMIC
#endif
    P2TRCINF(JSEQ) = P2TRCINF(JSEQ) + D2TRCOUT(ISEQ,ITRACE)             !! total inflow to a grid (from upstream)
  END DO
#ifndef NoAtom_CMF
  !$OMP END PARALLEL DO
#endif
  
  !! river mouth-----------------
  !$OMP PARALLEL DO
  DO ISEQ=NSEQRIV+1, NSEQALL
    D2TRCOUT(ISEQ,ITRACE) = D2TRCOUT(ISEQ,ITRACE)*D2RATE(ISEQ)
  END DO
  !$OMP END PARALLEL DO


  !! bifurcation channel
  IF(  CMF_TRACER%LTRCBIF ) THEN
#ifndef NoAtom_CMF
    !$OMP PARALLEL DO
#endif
    DO IPTH=1, NPTHOUT  
      ISEQP=PTH_UPST(IPTH)
      JSEQP=PTH_DOWN(IPTH)
      IF (ISEQP<=0 .OR. JSEQP<=0 ) CYCLE  !! Avoid calculation outside of domain

      IF( D1TRCPFLW(IPTH,ITRACE) >= 0._JPRB )THEN
        D1TRCPFLW(IPTH,ITRACE)  = D1TRCPFLW(IPTH,ITRACE)*D2RATE(ISEQP)
      ELSE
        D1TRCPFLW(IPTH,ITRACE)  = D1TRCPFLW(IPTH,ITRACE)*D2RATE(JSEQP)  !! reverse flow
      ENDIF

#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
      D2TRCPOUT(ISEQP,ITRACE) = D2TRCPOUT(ISEQP,ITRACE) + D1TRCPFLW(IPTH,ITRACE)
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
      D2TRCPOUT(JSEQP,ITRACE) = D2TRCPOUT(JSEQP,ITRACE) - D1TRCPFLW(IPTH,ITRACE)
    END DO
#ifndef NoAtom_CMF
    !$OMP END PARALLEL DO
#endif
  ENDIF

  !============================
  !*** 4. calculate next step storage

!$OMP PARALLEL DO
  DO ISEQ=1, NSEQALL
    P2TRCSTO(ISEQ,ITRACE) = P2TRCSTO(ISEQ,ITRACE) + P2TRCINF(ISEQ)*CMF_CONFIG%DT - D2TRCOUT(ISEQ,ITRACE)*CMF_CONFIG%DT
    P2TRCSTO(ISEQ,ITRACE) = P2TRCSTO(ISEQ,ITRACE) - D2TRCPOUT(ISEQ,ITRACE)*CMF_CONFIG%DT
    P2TRCSTO(ISEQ,ITRACE) = P2TRCSTO(ISEQ,ITRACE) + D2TRCINP(ISEQ,ITRACE) *CMF_CONFIG%DT
    P2TRCSTO(ISEQ,ITRACE) = max ( P2TRCSTO(ISEQ,ITRACE), 0._JPRB )
  END DO
!$OMP END PARALLEL DO

END DO

CALL CMF_TRACER_DIAG_AVEADD

END SUBROUTINE CMF_TRACER_FLUX
!####################################################################


!





!####################################################################
SUBROUTINE CMF_TRACER_DIAG_AVEADD
USE CMF_VARS_MOD,       ONLY: NADD_out,D2TRCOUT_oAVG,D2TRCDNS_oAVG,D2TRCPOUT_oAVG
USE CMF_VARS_MOD,       ONLY: D2TRCOUT,D2TRCDNS,D2TRCPOUT
IMPLICIT NONE
!====================
NADD_out=NADD_out+CMF_CONFIG%DT 
D2TRCOUT_oAVG(:,:)  = D2TRCOUT_oAVG(:,:)  + D2TRCOUT(:,:) * CMF_CONFIG%DT 
D2TRCDNS_oAVG(:,:)  = D2TRCDNS_oAVG(:,:)  + D2TRCDNS(:,:) * CMF_CONFIG%DT 
D2TRCPOUT_oAVG(:,:) = D2TRCPOUT_oAVG(:,:) + D2TRCPOUT(:,:)* CMF_CONFIG%DT 
END SUBROUTINE CMF_TRACER_DIAG_AVEADD
!####################################################################
!

END MODULE CMF_TRACER_MOD
