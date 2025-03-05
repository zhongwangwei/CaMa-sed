MODULE CMF_FLDSTG_MOD
!==========================================================
!* PURPOSE: call CaMa-Flood physics
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
USE CMF_NMLIST_MOD
USE PARKIND1,           ONLY: JPIM, JPRB, JPRD

CONTAINS 
SUBROUTINE CMF_CALC_FLDSTG_DEF
  USE CMF_VARS_MOD,        ONLY: NSEQALL, D2GRAREA, D2RIVLEN, D2RIVWTH, D2RIVELV
  USE CMF_VARS_MOD,        ONLY: P2RIVSTOMAX, P2FLDSTOMAX, D2FLDGRD, DFRCINC
  USE CMF_VARS_MOD,        ONLY: P2RIVSTO, P2FLDSTO
  USE CMF_VARS_MOD,        ONLY: D2RIVDPH, D2FLDDPH, D2FLDFRC, D2FLDARE, D2SFCELV
  USE CMF_VARS_MOD,        ONLY: P0GLBSTOPRE2, P0GLBSTONEW2, P0GLBRIVSTO, P0GLBFLDSTO, P0GLBFLDARE
  IMPLICIT NONE
  
  !*** LOCAL
  INTEGER(KIND=JPIM),SAVE    :: ISEQ, I
  REAL(KIND=JPRD),SAVE       :: DSTOALL, DSTONOW, DSTOPRE, DWTHNOW, DWTHPRE, DDPHPRE, DWTHINC
  !$OMP THREADPRIVATE        (I,DSTOALL, DSTONOW, DSTOPRE, DWTHNOW, DWTHPRE, DDPHPRE, DWTHINC)
  !================================================
  P0GLBSTOPRE2=0._JPRD
  P0GLBSTONEW2=0._JPRD
  P0GLBRIVSTO =0._JPRD
  P0GLBFLDSTO =0._JPRD
  P0GLBFLDARE =0._JPRD
  
  ! Estimate water depth and flood extent from water storage
  !   Solution for Equations (1) and (2) in [Yamazaki et al. 2011 WRR].
  
  !$OMP PARALLEL DO REDUCTION(+:P0GLBSTOPRE2,P0GLBSTONEW2,P0GLBRIVSTO,P0GLBFLDSTO,P0GLBFLDARE)
  DO ISEQ=1, NSEQALL
  !
    DSTOALL = P2RIVSTO(ISEQ,1) + P2FLDSTO(ISEQ,1)
  
    IF( DSTOALL > P2RIVSTOMAX(ISEQ,1) )THEN
      I=1
      DSTOPRE = P2RIVSTOMAX(ISEQ,1)
      DWTHPRE = D2RIVWTH(ISEQ,1)
      DDPHPRE = 0._JPRB
      DWTHINC = D2GRAREA(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * DFRCINC
      DO WHILE( DSTOALL > P2FLDSTOMAX(ISEQ,1,I) .AND. I<=CMF_CONFIG%NLFP)
        DSTOPRE = P2FLDSTOMAX(ISEQ,1,I)
        DWTHPRE = DWTHPRE + DWTHINC
        DDPHPRE = DDPHPRE + D2FLDGRD(ISEQ,1,I) * DWTHINC
        I=I+1
        IF( I>CMF_CONFIG%NLFP ) EXIT
      END DO
      IF( I>CMF_CONFIG%NLFP )THEN
        DSTONOW = DSTOALL - DSTOPRE
        DWTHNOW = 0._JPRB
        D2FLDDPH(ISEQ,1) = DDPHPRE + DSTONOW * DWTHPRE**(-1.) * D2RIVLEN(ISEQ,1)**(-1.)
      ELSE
        DSTONOW =  DSTOALL - DSTOPRE
        DWTHNOW = -DWTHPRE + &
  &      ( DWTHPRE**2. + 2._JPRB * DSTONOW * D2RIVLEN(ISEQ,1)**(-1.) * D2FLDGRD(ISEQ,1,I)**(-1.) )**0.5
        D2FLDDPH(ISEQ,1) = DDPHPRE + D2FLDGRD(ISEQ,1,I) * DWTHNOW
      ENDIF
      P2RIVSTO(ISEQ,1) = P2RIVSTOMAX(ISEQ,1) + D2RIVLEN(ISEQ,1) * D2RIVWTH(ISEQ,1) * D2FLDDPH(ISEQ,1)
      P2RIVSTO(ISEQ,1) = MIN(P2RIVSTO(ISEQ,1),DSTOALL)
  
      D2RIVDPH(ISEQ,1) = P2RIVSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * D2RIVWTH(ISEQ,1)**(-1.)
  !
      P2FLDSTO(ISEQ,1) = DSTOALL - P2RIVSTO(ISEQ,1)
      P2FLDSTO(ISEQ,1) = MAX( P2FLDSTO(ISEQ,1), 0._JPRD )
      D2FLDFRC(ISEQ,1) = (-D2RIVWTH(ISEQ,1) + DWTHPRE + DWTHNOW ) * (DWTHINC*CMF_CONFIG%NLFP)**(-1.)  !! bugfix 191113, (10._JPRB -> CMF_CONFIG%NLFP)
      D2FLDFRC(ISEQ,1) = MAX( D2FLDFRC(ISEQ,1),0._JPRB)
      D2FLDFRC(ISEQ,1) = MIN( D2FLDFRC(ISEQ,1),1._JPRB)
      D2FLDARE(ISEQ,1) = D2GRAREA(ISEQ,1)*D2FLDFRC(ISEQ,1)
    ELSE
      P2RIVSTO(ISEQ,1) = DSTOALL
      D2RIVDPH(ISEQ,1) = DSTOALL * D2RIVLEN(ISEQ,1)**(-1.) * D2RIVWTH(ISEQ,1)**(-1.)
      D2RIVDPH(ISEQ,1) = MAX( D2RIVDPH(ISEQ,1), 0._JPRB )
      P2FLDSTO(ISEQ,1) = 0._JPRD
      D2FLDDPH(ISEQ,1) = 0._JPRB
      D2FLDFRC(ISEQ,1) = 0._JPRB
      D2FLDARE(ISEQ,1) = 0._JPRB
    ENDIF
    D2SFCELV(ISEQ,1)     = D2RIVELV(ISEQ,1) + D2RIVDPH(ISEQ,1)
  
    P0GLBSTOPRE2     = P0GLBSTOPRE2 + DSTOALL
    P0GLBSTONEW2     = P0GLBSTONEW2 + P2RIVSTO(ISEQ,1) + P2FLDSTO(ISEQ,1)
    P0GLBRIVSTO      = P0GLBRIVSTO  + P2RIVSTO(ISEQ,1)
    P0GLBFLDSTO      = P0GLBFLDSTO  + P2FLDSTO(ISEQ,1)
    P0GLBFLDARE      = P0GLBFLDARE  + D2FLDARE(ISEQ,1)
  
  END DO
  !$OMP END PARALLEL DO
  
  END SUBROUTINE CMF_CALC_FLDSTG_DEF
  !####################################################################
  !
  !
  !
  !
  !####################################################################
  SUBROUTINE CMF_OPT_FLDSTG_ES
  ! ==========
  ! Optional code for Earth Simulator (Vector Processor)
  ! Specify option: LSTG_ES=.TRUE.
  ! Faster computation on vector prosessor by avoiding IF-THEN function. (note this code will be slow on Scaler Processor)
  ! ==========
  USE CMF_VARS_MOD,        ONLY: NSEQALL, D2GRAREA, D2RIVLEN, D2RIVWTH, D2RIVELV, D2RIVHGT
  USE CMF_VARS_MOD,        ONLY: P2RIVSTOMAX, P2FLDSTOMAX, D2FLDGRD, DFRCINC
  USE CMF_VARS_MOD,        ONLY: P2RIVSTO, P2FLDSTO
  USE CMF_VARS_MOD,        ONLY: D2RIVDPH, D2FLDDPH, D2FLDFRC, D2FLDARE, D2SFCELV
  USE CMF_VARS_MOD,        ONLY: P0GLBSTOPRE2, P0GLBSTONEW2, P0GLBRIVSTO, P0GLBFLDSTO, P0GLBFLDARE
  IMPLICIT NONE
  
  !*** LOCAL
  REAL(KIND=JPRD)            :: D2STODWN(NSEQALL,1)
  REAL(KIND=JPRD)            :: D2WTHPRE(NSEQALL,1), D2WTHINC(NSEQALL,1)
  
  ! SAVE for OpenMP
  INTEGER(KIND=JPIM),SAVE    :: ISEQ, I
  REAL(KIND=JPRD),SAVE       :: DSTOALL, DSTONOW, DWTHNOW
  !$OMP THREADPRIVATE          (DSTOALL, DSTONOW, DWTHNOW)
  !================================================
  P0GLBRIVSTO=0._JPRD
  P0GLBFLDSTO=0._JPRD
  P0GLBFLDARE=0._JPRD
  P0GLBSTOPRE2=0._JPRD
  P0GLBSTONEW2=0._JPRD
  
  ! [1] Assume all waters in river channel
  !$OMP PARALLEL DO REDUCTION(+:P0GLBSTOPRE2)
  DO ISEQ=1, NSEQALL
    DSTOALL = P2RIVSTO(ISEQ,1) + P2FLDSTO(ISEQ,1)
  
    P2RIVSTO(ISEQ,1) = DSTOALL
    D2RIVDPH(ISEQ,1) = DSTOALL * D2RIVLEN(ISEQ,1)**(-1.) * D2RIVWTH(ISEQ,1)**(-1.)
    D2RIVDPH(ISEQ,1) = MAX( D2RIVDPH(ISEQ,1), 0._JPRB )
    P2FLDSTO(ISEQ,1) = 0._JPRD
    D2FLDDPH(ISEQ,1) = 0._JPRB
    D2FLDFRC(ISEQ,1) = 0._JPRB
    D2FLDARE(ISEQ,1) = 0._JPRB
  
    D2STODWN(ISEQ,1) = P2RIVSTOMAX(ISEQ,1)
    D2WTHPRE(ISEQ,1) = D2RIVWTH(ISEQ,1)
    D2WTHINC(ISEQ,1) = D2GRAREA(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * DFRCINC
  
    P0GLBSTOPRE2     = P0GLBSTOPRE2 + DSTOALL
  END DO
  !$OMP END PARALLEL DO
  
  ! [2] Check floodplain level from I=1 to CMF_CONFIG%NLFP. Make I-CMF_CONFIG%NLFP loop outside for parallel computing (SIMD/Vector)
  DO I=1, CMF_CONFIG%NLFP
  
  !$OMP PARALLEL DO
    DO ISEQ=1, NSEQALL
      DSTOALL = P2RIVSTO(ISEQ,1) + P2FLDSTO(ISEQ,1)
  
      DSTONOW = DSTOALL - D2STODWN(ISEQ,1)
      DSTONOW = MAX( DSTONOW, 0._JPRD )
      DWTHNOW = -D2WTHPRE(ISEQ,1) + &
  &    ( D2WTHPRE(ISEQ,1)**2._JPRB + 2._JPRB * DSTONOW * D2RIVLEN(ISEQ,1)**(-1.) * D2FLDGRD(ISEQ,1,I)**(-1.) )**0.5
      DWTHNOW = MIN( DWTHNOW, D2WTHINC(ISEQ,1) )
      DWTHNOW = MAX( DWTHNOW, 0.D0 )              !! modify v4.04
  
      D2FLDDPH(ISEQ,1) = D2FLDDPH(ISEQ,1) + D2FLDGRD(ISEQ,1,I) * DWTHNOW
      D2FLDFRC(ISEQ,1) = D2FLDFRC(ISEQ,1) + DWTHNOW/D2WTHINC(ISEQ,1) * CMF_CONFIG%NLFP**(-1.)
  
      !! Update downside floodplain step storage/depth/width 
      D2STODWN(ISEQ,1) = P2FLDSTOMAX(ISEQ,1,I)
      D2WTHPRE(ISEQ,1) = D2WTHPRE(ISEQ,1) + D2WTHINC(ISEQ,1)
    END DO
    !$OMP END PARALLEL DO
  END DO
  
  !! [3] flood extent saturated case
  !$OMP PARALLEL DO
  DO ISEQ=1, NSEQALL
    DSTOALL = P2RIVSTO(ISEQ,1) + P2FLDSTO(ISEQ,1)
    DSTONOW = DSTOALL - D2STODWN(ISEQ,1)
    DSTONOW = MAX( DSTONOW, 0._JPRD )
    D2FLDDPH(ISEQ,1) = D2FLDDPH(ISEQ,1) + DSTONOW * D2WTHPRE(ISEQ,1)**(-1.) * D2RIVLEN(ISEQ,1)**(-1.)
  END DO
  !$OMP END PARALLEL DO
  
  !! [4] Floodplain stage diagnose
  !$OMP PARALLEL DO
  DO ISEQ=1, NSEQALL
  !  IF( D2FLDDPH(ISEQ,1)>0 )THEN !! bugfix v4.04
    IF( D2FLDDPH(ISEQ,1)>1.D-5 )THEN !! bugfix v4.04, to avoid false positive FLDDPH due to rounding error.
      DSTOALL = P2RIVSTO(ISEQ,1) + P2FLDSTO(ISEQ,1)
  
      D2RIVDPH(ISEQ,1) = D2RIVHGT(ISEQ,1) + D2FLDDPH(ISEQ,1)
      P2RIVSTO(ISEQ,1) = D2RIVLEN(ISEQ,1) * D2RIVWTH(ISEQ,1) * D2RIVDPH(ISEQ,1)
      P2RIVSTO(ISEQ,1) = MIN( P2RIVSTO(ISEQ,1), DSTOALL )           !! modify v4.04
  !
      P2FLDSTO(ISEQ,1) = DSTOALL - P2RIVSTO(ISEQ,1)
      P2FLDSTO(ISEQ,1) = MAX( P2FLDSTO(ISEQ,1), 0._JPRD )
  
      D2FLDFRC(ISEQ,1) = MAX( D2FLDFRC(ISEQ,1),0._JPRB)
      D2FLDFRC(ISEQ,1) = MIN( D2FLDFRC(ISEQ,1),1._JPRB)
      D2FLDARE(ISEQ,1) = D2GRAREA(ISEQ,1)*D2FLDFRC(ISEQ,1)
    ENDIF
  END DO
  !$OMP END PARALLEL DO
  
  !$OMP PARALLEL DO REDUCTION(+:P0GLBSTONEW2,P0GLBRIVSTO,P0GLBFLDSTO,P0GLBFLDARE)
  DO ISEQ=1, NSEQALL
    D2SFCELV(ISEQ,1) = D2RIVELV(ISEQ,1) + D2RIVDPH(ISEQ,1)
    P0GLBSTONEW2      = P0GLBSTONEW2+ P2RIVSTO(ISEQ,1) + P2FLDSTO(ISEQ,1)
    P0GLBRIVSTO       = P0GLBRIVSTO + P2RIVSTO(ISEQ,1)
    P0GLBFLDSTO       = P0GLBFLDSTO + P2FLDSTO(ISEQ,1)
    P0GLBFLDARE       = P0GLBFLDARE + D2FLDARE(ISEQ,1)
  END DO
  !$OMP END PARALLEL DO
  
  END SUBROUTINE CMF_OPT_FLDSTG_ES
  !####################################################################
  SUBROUTINE CMF_LEVEE_FLDSTG
    ! ================================================
    ! calculate river and floodplain staging considering levee
    ! ================================================
    USE CMF_VARS_MOD,        ONLY: NSEQALL
    USE CMF_VARS_MOD,        ONLY: D2GRAREA, D2RIVLEN, D2RIVWTH, D2RIVELV, P2RIVSTOMAX, P2FLDSTOMAX, D2FLDGRD, DFRCINC, D2FLDHGT
    USE CMF_VARS_MOD,        ONLY: P2RIVSTO, P2FLDSTO
    USE CMF_VARS_MOD,        ONLY: D2RIVDPH, D2FLDDPH, D2FLDFRC, D2FLDARE, D2SFCELV
    USE CMF_VARS_MOD,        ONLY: P0GLBSTOPRE2, P0GLBSTONEW2, P0GLBRIVSTO, P0GLBFLDSTO, P0GLBLEVSTO, P0GLBFLDARE
    USE CMF_VARS_MOD,        ONLY: D2LEVBASSTO, D2LEVTOPSTO, D2LEVFILSTO, D2LEVFRC, D2LEVDST, D2BASHGT
    USE CMF_VARS_MOD,        ONLY: D2LEVHGT
    
    !! levee specific data
    USE CMF_VARS_MOD   ,ONLY: P2LEVSTO  !! flood storage in protected side (P2FLDSTO for storage betwen river & levee)
    USE CMF_VARS_MOD   ,ONLY: D2LEVDPH  !! flood depth in protected side   (D2FLDDPH for water depth betwen river & levee)
    IMPLICIT NONE
    
    !*** LOCAL
    ! Save for OpenMP
    INTEGER(KIND=JPIM),SAVE ::  ISEQ, I, ILEV
    REAL(KIND=JPRD),SAVE    ::  DSTOALL, DSTONOW, DSTOPRE, DWTHNOW, DWTHPRE, DDPHPRE, DDPHNOW, DWTHINC, DSTOADD
!$OMP THREADPRIVATE (I,ILEV,DSTOALL, DSTONOW, DSTOPRE, DWTHNOW, DWTHPRE, DDPHPRE, DDPHNOW, DWTHINC, DSTOADD)
    !!==============================
    P0GLBSTOPRE2=0._JPRD
    P0GLBSTONEW2=0._JPRD
    P0GLBRIVSTO =0._JPRD
    P0GLBFLDSTO =0._JPRD
    P0GLBLEVSTO =0._JPRD
    P0GLBFLDARE =0._JPRD
    
!$OMP PARALLEL DO REDUCTION(+:P0GLBSTOPRE2,P0GLBSTONEW2,P0GLBRIVSTO,P0GLBFLDSTO,P0GLBLEVSTO,P0GLBFLDARE)
    DO ISEQ=1, NSEQALL
    !
      DSTOALL = P2RIVSTO(ISEQ,1) + P2FLDSTO(ISEQ,1) + P2LEVSTO(ISEQ,1)
      DWTHINC = D2GRAREA(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * DFRCINC    !! width of each layer [m]
      IF( DSTOALL > P2RIVSTOMAX(ISEQ,1) )THEN
        !**********
        ! [Case-1] Water surface is under levee base (all water is between river-levee)
        IF( DSTOALL < D2LEVBASSTO(ISEQ,1) )THEN 
          I=1
          DSTOPRE = P2RIVSTOMAX(ISEQ,1)
          DWTHPRE = D2RIVWTH(ISEQ,1)
          DDPHPRE = 0._JPRB
    
          ! which layer current water level is
          DO WHILE( DSTOALL > P2FLDSTOMAX(ISEQ,1,I) .AND. I<=CMF_CONFIG%NLFP )
            DSTOPRE = P2FLDSTOMAX(ISEQ,1,I)
            DWTHPRE = DWTHPRE + DWTHINC
            DDPHPRE = DDPHPRE + D2FLDGRD(ISEQ,1,I) * DWTHINC
            I=I+1
            IF( I>CMF_CONFIG%NLFP ) EXIT
          END DO
    
          ! water depth at unprotected area
          IF( I<=CMF_CONFIG%NLFP )THEN
            DSTONOW =  DSTOALL - DSTOPRE
            DWTHNOW = -DWTHPRE + ( DWTHPRE**2. + 2. * DSTONOW * D2RIVLEN(ISEQ,1)**(-1.) * D2FLDGRD(ISEQ,1,I)**(-1.) )**0.5
            D2FLDDPH(ISEQ,1) = DDPHPRE + D2FLDGRD(ISEQ,1,I) * DWTHNOW
          ELSE
            DSTONOW = DSTOALL - DSTOPRE
            DWTHNOW = 0._JPRB
            D2FLDDPH(ISEQ,1) = DDPHPRE + DSTONOW * DWTHPRE**(-1.) * D2RIVLEN(ISEQ,1)**(-1.)
          ENDIF
    
          P2RIVSTO(ISEQ,1) = P2RIVSTOMAX(ISEQ,1) + D2RIVLEN(ISEQ,1) * D2RIVWTH(ISEQ,1) * D2FLDDPH(ISEQ,1)
          D2RIVDPH(ISEQ,1) = P2RIVSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * D2RIVWTH(ISEQ,1)**(-1.)
    !
          P2FLDSTO(ISEQ,1) = DSTOALL - P2RIVSTO(ISEQ,1)
          P2FLDSTO(ISEQ,1) = MAX( P2FLDSTO(ISEQ,1), 0._JPRD )
          D2FLDFRC(ISEQ,1) = (-D2RIVWTH(ISEQ,1) + DWTHPRE + DWTHNOW ) * (DWTHINC*CMF_CONFIG%NLFP)**(-1.)
          D2FLDFRC(ISEQ,1) = MAX( D2FLDFRC(ISEQ,1),0._JPRB )
          D2FLDFRC(ISEQ,1) = MIN( D2FLDFRC(ISEQ,1),1._JPRB )
          D2FLDARE(ISEQ,1) = D2GRAREA(ISEQ,1)*D2FLDFRC(ISEQ,1)
    !
          P2LEVSTO(ISEQ,1) = 0._JPRD  !! no flooding in protected area
          D2LEVDPH(ISEQ,1) = 0._JPRB
    
        !**********
        ! [Case-2]  River-side water surface is under levee crown (water only in river side)
        ELSEIF( DSTOALL < D2LEVTOPSTO(ISEQ,1) )THEN 
    
          DSTONOW = DSTOALL - D2LEVBASSTO(ISEQ,1)
          DWTHNOW = D2LEVDST(ISEQ,1) + D2RIVWTH(ISEQ,1)
          D2FLDDPH(ISEQ,1) = D2BASHGT(ISEQ,1) + DSTONOW * DWTHNOW**(-1.) * D2RIVLEN(ISEQ,1)**(-1.)
    
          P2RIVSTO(ISEQ,1) = P2RIVSTOMAX(ISEQ,1) + D2RIVLEN(ISEQ,1) * D2RIVWTH(ISEQ,1) * D2FLDDPH(ISEQ,1)
          D2RIVDPH(ISEQ,1) = P2RIVSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * D2RIVWTH(ISEQ,1)**(-1.)
      !
          P2FLDSTO(ISEQ,1) = DSTOALL - P2RIVSTO(ISEQ,1)
          P2FLDSTO(ISEQ,1) = MAX( P2FLDSTO(ISEQ,1), 0._JPRD )
          D2FLDFRC(ISEQ,1) = D2LEVFRC(ISEQ,1)
          D2FLDARE(ISEQ,1) = D2GRAREA(ISEQ,1)*D2FLDFRC(ISEQ,1)
      ! 
          P2LEVSTO(ISEQ,1) = 0._JPRD  !! no flooding in protected area
          D2LEVDPH(ISEQ,1) = 0._JPRB
    
        !**********
        ! [Case-3] River side is full, protected side is under levee crown height (Water both in river side & protected side)
        ELSEIF( DSTOALL < D2LEVFILSTO(ISEQ,1) )THEN 
          ! river side stage = levee height
          D2FLDDPH(ISEQ,1) = D2LEVHGT(ISEQ,1)
          P2RIVSTO(ISEQ,1) = P2RIVSTOMAX(ISEQ,1) + D2RIVLEN(ISEQ,1) * D2RIVWTH(ISEQ,1) * D2FLDDPH(ISEQ,1)
          D2RIVDPH(ISEQ,1) = P2RIVSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * D2RIVWTH(ISEQ,1)**(-1.)
    
          P2FLDSTO(ISEQ,1) = D2LEVTOPSTO(ISEQ,1) - P2RIVSTO(ISEQ,1)
          P2FLDSTO(ISEQ,1) = MAX( P2FLDSTO(ISEQ,1), 0._JPRD )
    
          !! protected side storate calculation
          P2LEVSTO(ISEQ,1) = DSTOALL - P2RIVSTO(ISEQ,1) - P2FLDSTO(ISEQ,1)
          P2LEVSTO(ISEQ,1) = MAX( P2LEVSTO(ISEQ,1), 0._JPRD )
    
          !!****
          !! protected side stage calculation
          ILEV=INT( D2LEVFRC(ISEQ,1)*CMF_CONFIG%NLFP )+1 !! levee relative distance -> floodplain layer with levee
          DSTOPRE = D2LEVTOPSTO(ISEQ,1)
          DWTHPRE = 0._JPRB
          DDPHPRE = 0._JPRB
          !! which layer current water level is
          I=ILEV
          DO WHILE( I<=CMF_CONFIG%NLFP )
            DSTOADD = ( D2LEVDST(ISEQ,1)+D2RIVWTH(ISEQ,1) ) * ( D2LEVHGT(ISEQ,1)-D2FLDHGT(ISEQ,1,I) ) * D2RIVLEN(ISEQ,1) 
            IF( DSTOALL < P2FLDSTOMAX(ISEQ,1,I) + DSTOADD ) EXIT
            DSTOPRE = P2FLDSTOMAX(ISEQ,1,I) + DSTOADD
            DWTHPRE = DWTHINC*I - D2LEVDST(ISEQ,1)
            DDPHPRE = D2FLDHGT(ISEQ,1,I) - D2BASHGT(ISEQ,1)
            I=I+1
            IF( I>CMF_CONFIG%NLFP ) EXIT
          END DO
    
          IF( I<=CMF_CONFIG%NLFP )THEN
            DSTONOW = DSTOALL - DSTOPRE
            DWTHNOW = -DWTHPRE + ( DWTHPRE**2. + 2. * DSTONOW*D2RIVLEN(ISEQ,1)**(-1.) * D2FLDGRD(ISEQ,1,I)**(-1.) )**0.5
            DDPHNOW = DWTHNOW * D2FLDGRD(ISEQ,1,I)
            D2LEVDPH(ISEQ,1) = D2BASHGT(ISEQ,1) + DDPHPRE + DDPHNOW
    
            D2FLDFRC(ISEQ,1) = ( DWTHPRE + D2LEVDST(ISEQ,1) ) * (DWTHINC*CMF_CONFIG%NLFP)**(-1.)
            D2FLDFRC(ISEQ,1) = MAX( D2FLDFRC(ISEQ,1),0._JPRB)
            D2FLDFRC(ISEQ,1) = MIN( D2FLDFRC(ISEQ,1),1._JPRB)
            D2FLDARE(ISEQ,1) = D2GRAREA(ISEQ,1)*D2FLDFRC(ISEQ,1)
          ELSE
            DSTONOW = DSTOALL - DSTOPRE
            DDPHNOW = DSTONOW * DWTHPRE**(-1.) * D2RIVLEN(ISEQ,1)**(-1.)
            D2LEVDPH(ISEQ,1) = D2BASHGT(ISEQ,1) + DDPHPRE + DDPHNOW
    
            D2FLDFRC(ISEQ,1) = 1._JPRB
            D2FLDARE(ISEQ,1) = D2GRAREA(ISEQ,1)*D2FLDFRC(ISEQ,1)
          ENDIF
    
        !**********
        ! [Case-4] Water level above levee crown (Both river side and protected side exceed levee crown height)
        ELSE 
          I=1
          DSTOPRE = P2RIVSTOMAX(ISEQ,1)
          DWTHPRE = D2RIVWTH(ISEQ,1)
          DDPHPRE = 0._JPRB
          DO WHILE( DSTOALL > P2FLDSTOMAX(ISEQ,1,I) .AND. I<=CMF_CONFIG%NLFP)
            DSTOPRE = P2FLDSTOMAX(ISEQ,1,I)
            DWTHPRE = DWTHPRE + DWTHINC
            DDPHPRE = DDPHPRE + D2FLDGRD(ISEQ,1,I) * DWTHINC
            I=I+1
            IF( I>CMF_CONFIG%NLFP ) EXIT
          END DO
    
          IF( I<=CMF_CONFIG%NLFP )THEN
            DSTONOW =  DSTOALL - DSTOPRE
            DWTHNOW = -DWTHPRE + ( DWTHPRE**2. + 2. * DSTONOW * D2RIVLEN(ISEQ,1)**(-1.) * D2FLDGRD(ISEQ,1,I)**(-1.) )**0.5
            D2FLDDPH(ISEQ,1) = DDPHPRE + D2FLDGRD(ISEQ,1,I) * DWTHNOW
          ELSE
            DSTONOW = DSTOALL - DSTOPRE
            DWTHNOW = 0._JPRB
            D2FLDDPH(ISEQ,1) = DDPHPRE + DSTONOW * DWTHPRE**(-1.) * D2RIVLEN(ISEQ,1)**(-1.)
          ENDIF
    
          D2FLDFRC(ISEQ,1) = (-D2RIVWTH(ISEQ,1) + DWTHPRE + DWTHNOW ) * (DWTHINC*CMF_CONFIG%NLFP)**(-1.)
          D2FLDARE(ISEQ,1) = D2GRAREA(ISEQ,1)*D2FLDFRC(ISEQ,1)
    
          !! river channel storage
          P2RIVSTO(ISEQ,1) = P2RIVSTOMAX(ISEQ,1) + D2RIVLEN(ISEQ,1) * D2RIVWTH(ISEQ,1) * D2FLDDPH(ISEQ,1)
          D2RIVDPH(ISEQ,1) = P2RIVSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * D2RIVWTH(ISEQ,1)**(-1.)
    !
          DSTOADD = ( D2FLDDPH(ISEQ,1)-D2LEVHGT(ISEQ,1) ) * (D2LEVDST(ISEQ,1)+D2RIVWTH(ISEQ,1)) * D2RIVLEN(ISEQ,1)
          P2FLDSTO(ISEQ,1) = D2LEVTOPSTO(ISEQ,1) + DSTOADD - P2RIVSTO(ISEQ,1)
          P2FLDSTO(ISEQ,1) = MAX( P2FLDSTO(ISEQ,1), 0._JPRD )
    
          P2LEVSTO(ISEQ,1) = DSTOALL - P2RIVSTO(ISEQ,1) - P2FLDSTO(ISEQ,1)
          P2LEVSTO(ISEQ,1) = MAX( P2LEVSTO(ISEQ,1), 0._JPRD )
          D2LEVDPH(ISEQ,1) = D2FLDDPH(ISEQ,1)
        ENDIF
    
      ! [Case-0] Water only in river channel
      ELSE
        P2RIVSTO(ISEQ,1) = DSTOALL
        D2RIVDPH(ISEQ,1) = DSTOALL * D2RIVLEN(ISEQ,1)**(-1.) * D2RIVWTH(ISEQ,1)**(-1.)
        D2RIVDPH(ISEQ,1) = MAX( D2RIVDPH(ISEQ,1), 0._JPRB )
        P2FLDSTO(ISEQ,1) = 0._JPRD
        D2FLDDPH(ISEQ,1) = 0._JPRB
        D2FLDFRC(ISEQ,1) = 0._JPRB
        D2FLDARE(ISEQ,1) = 0._JPRB
        P2LEVSTO(ISEQ,1) = 0._JPRD
        D2LEVDPH(ISEQ,1) = 0._JPRB
      ENDIF
      D2SFCELV(ISEQ,1)     = D2RIVELV(ISEQ,1) + D2RIVDPH(ISEQ,1)
    
      P0GLBSTOPRE2     = P0GLBSTOPRE2 + DSTOALL
      P0GLBSTONEW2     = P0GLBSTONEW2 + P2RIVSTO(ISEQ,1) + P2FLDSTO(ISEQ,1) + P2LEVSTO(ISEQ,1)
      P0GLBRIVSTO      = P0GLBRIVSTO  + P2RIVSTO(ISEQ,1)
      P0GLBFLDSTO      = P0GLBFLDSTO  + P2FLDSTO(ISEQ,1)
      P0GLBLEVSTO      = P0GLBLEVSTO  + P2LEVSTO(ISEQ,1)
      P0GLBFLDARE      = P0GLBFLDARE  + D2FLDARE(ISEQ,1)
    END DO
    !$OMP END PARALLEL DO
    
    END SUBROUTINE CMF_LEVEE_FLDSTG
    !####################################################################


!###############################################################
SUBROUTINE CMF_PHYSICS_FLDSTG
! flood stage scheme selecter
IMPLICIT NONE

IF( CMF_OPTIONS%LLEVEE )THEN
  CALL CMF_LEVEE_FLDSTG  !! levee floodstage (Vector åprocessor option not available)
ELSE
  IF( CMF_OPTIONS%LSTG_ES )THEN
    CALL CMF_OPT_FLDSTG_ES  !! Alternative subroutine optimized for vector processor
  ELSE 
    CALL CMF_CALC_FLDSTG_DEF     !! Default
  ENDIF
ENDIF

END SUBROUTINE CMF_PHYSICS_FLDSTG
!###############################################################

END MODULE CMF_FLDSTG_MOD
