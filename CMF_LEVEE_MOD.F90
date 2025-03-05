MODULE CMF_LEVEE_MOD
!==========================================================
!* PURPOSE: CaMa-Flood flood stage with levee scheme
!
! (C) Y. Tanaka & D.Yamazaki (U-Tokyo)  Mar 2022
!
!* CONTAINS:
! -- CMF_LEV_NMLIST      : Read setting from namelist
! -- CMF_LEV_INIT        : Initialize levee scheme
! -- CMF_CALC_FLDSTG_LEV : Calculate flood stage considering levee
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
USE PARKIND1,                ONLY: JPIM, JPRB, JPRM, JPRD
USE CMF_NMLIST_MOD,        ONLY: CMF_FILES, CMF_OPTIONS, CMF_CONFIG, CMF_PARAMS
USE CMF_NMLIST_MOD,        ONLY: CMF_TIME, CMF_MAPS, CMF_FORCING, CMF_BOUNDARY
USE CMF_NMLIST_MOD,        ONLY: CMF_RESTART, CMF_DAM, CMF_LEVEE, CMF_OUTPUT,CMF_TRACER,CMF_SED 
USE CMF_VARS_MOD,          ONLY: D2LEVHGT, D2LEVFRC, D2BASHGT, D2LEVDST, D2LEVBASSTO, D2LEVTOPSTO, D2LEVFILSTO
!============================
IMPLICIT NONE
SAVE


CONTAINS


!####################################################################
SUBROUTINE CMF_LEVEE_FLDSTG
! ================================================
! calculate river and floodplain staging considering levee
! ================================================
USE CMF_VARS_MOD,          ONLY: NSEQALL
USE CMF_VARS_MOD,          ONLY: D2GRAREA, D2RIVLEN, D2RIVWTH, D2RIVELV, P2RIVSTOMAX, P2FLDSTOMAX, D2FLDGRD, DFRCINC, D2FLDHGT
USE CMF_VARS_MOD,          ONLY: P2RIVSTO, P2FLDSTO
USE CMF_VARS_MOD,          ONLY: D2RIVDPH, D2FLDDPH, D2FLDFRC, D2FLDARE, D2SFCELV
USE CMF_VARS_MOD,          ONLY: P0GLBSTOPRE2, P0GLBSTONEW2, P0GLBRIVSTO, P0GLBFLDSTO, P0GLBLEVSTO, P0GLBFLDARE

!! levee specific data
USE CMF_VARS_MOD,          ONLY: P2LEVSTO  !! flood storage in protected side (P2FLDSTO for storage betwen river & levee)
USE CMF_VARS_MOD,          ONLY: D2LEVDPH  !! flood depth in protected side   (D2FLDDPH for water depth betwen river & levee)
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




!####################################################################
SUBROUTINE CMF_LEVEE_OPT_PTHOUT
! realistic bifurcation considering levee
USE PARKIND1,           ONLY: JPIM, JPRB, JPRD
USE CMF_VARS_MOD,        ONLY: NSEQALL, NSEQMAX, NPTHOUT, NPTHLEV, PTH_UPST, PTH_DOWN, PTH_DST, &
                            & PTH_ELV, PTH_WTH, PTH_MAN, I2MASK
USE CMF_VARS_MOD,        ONLY: D2ELEVTN, D2RIVELV
USE CMF_VARS_MOD,       ONLY: D1PTHFLW, D1PTHFLW_PRE, D2RIVDPH_PRE
USE CMF_VARS_MOD,       ONLY: D2LEVDPH, D2SFCELV
IMPLICIT NONE
!*** Local
REAL(KIND=JPRB)    ::  D2SFCELV_LEV(NSEQMAX,1)                  !! water surface elev protected [m]
REAL(KIND=JPRB)    ::  D2SFCELV_PRE(NSEQMAX,1)                  !! water surface elev (t-1) [m] (for stable calculation)

! SAVE for OpenMP
INTEGER(KIND=JPIM),SAVE ::  IPTH, ILEV, ISEQ, ISEQP, JSEQP
REAL(KIND=JPRB),SAVE    ::  DSLOPE, DFLW, DOUT_PRE, DFLW_PRE, DFLW_IMP
!$OMP THREADPRIVATE        (DSLOPE, DFLW, DOUT_PRE, DFLW_PRE, DFLW_IMP, ILEV, ISEQP, JSEQP)
!================================================
!$OMP PARALLEL DO
DO ISEQ=1, NSEQALL
  IF( D2LEVFRC(ISEQ,1)<1.0 )THEN
    D2SFCELV_LEV(ISEQ,1) = D2ELEVTN(ISEQ,1)+D2LEVDPH(ISEQ,1)  !! levee exist, calculate pthout based on levee protected depth
  ELSE
    D2SFCELV_LEV(ISEQ,1) = D2SFCELV(ISEQ,1)
  ENDIF

  D2SFCELV_PRE(ISEQ,1) = D2RIVELV(ISEQ,1)+D2RIVDPH_PRE(ISEQ,1)
END DO
!$OMP END PARALLEL DO

D1PTHFLW(:,:) = CMF_PARAMS%PGRV
!$OMP PARALLEL DO
DO IPTH=1, NPTHOUT  
  ISEQP=PTH_UPST(IPTH)
  JSEQP=PTH_DOWN(IPTH)
  !! Avoid calculation outside of domain
  IF (ISEQP<=0 .OR. JSEQP<=0 ) CYCLE
  IF (I2MASK(ISEQP,1) == 1 .OR. I2MASK(JSEQP,1) == 1 ) CYCLE  !! I2MASK is for kinematic-inertial mixed flow scheme. 

!! [1] for channel bifurcation, use river surface elevation  
  DSLOPE  = (D2SFCELV(ISEQP,1)-D2SFCELV(JSEQP,1)) * PTH_DST(IPTH)**(-1.)
  DSLOPE = max(-0.005_JPRB,min(0.005_JPRB,DSLOPE))                                    !! v390 stabilization

  ILEV=1 !! for river channek
    DFLW = MAX(D2SFCELV(ISEQP,1),D2SFCELV(JSEQP,1)) - PTH_ELV(IPTH,ILEV) 
    DFLW = MAX(DFLW,0._JPRB)

    DFLW_PRE = MAX(D2SFCELV_PRE(ISEQP,1),D2SFCELV_PRE(JSEQP,1)) - PTH_ELV(IPTH,ILEV)
    DFLW_PRE = MAX(DFLW_PRE,0._JPRB)

    DFLW_IMP = (DFLW*DFLW_PRE)**0.5                                       !! semi implicit flow depth
    IF( DFLW_IMP<=0._JPRB ) DFLW_IMP=DFLW

    IF( DFLW_IMP>1.E-5 )THEN                         !! local inertial equation, see [Bates et al., 2010, J.Hydrol.]
      DOUT_PRE = D1PTHFLW_PRE(IPTH,ILEV) * PTH_WTH(IPTH,ILEV)**(-1.)           !! outflow (t-1) [m2/s] (unit width)
      D1PTHFLW(IPTH,ILEV) = PTH_WTH(IPTH,ILEV) * ( DOUT_PRE + CMF_PARAMS%PGRV*CMF_CONFIG%DT*DFLW_IMP*DSLOPE ) &
                         * ( 1. + CMF_PARAMS%PGRV*CMF_CONFIG%DT*PTH_MAN(ILEV)**2. * abs(DOUT_PRE)*DFLW_IMP**(-7./3.) )**(-1.)
    ELSE
      D1PTHFLW(IPTH,ILEV) = 0._JPRB
    ENDIF

!! [1] for overland bifurcation, use levee protected surface elevation
  IF( NPTHLEV<=1 ) CYCLE

  DSLOPE  = (D2SFCELV_LEV(ISEQP,1)-D2SFCELV_LEV(JSEQP,1)) * PTH_DST(IPTH)**(-1.)
  DSLOPE = max(-0.005_JPRB,min(0.005_JPRB,DSLOPE))      

  DO ILEV=2, NPTHLEV
    DFLW = MAX(D2SFCELV_LEV(ISEQP,1),D2SFCELV_LEV(JSEQP,1)) - PTH_ELV(IPTH,ILEV) 
    DFLW = MAX(DFLW,0._JPRB)

    DFLW_IMP=DFLW  !! do not consider implicit flow depth for overland bifurcation
    IF( DFLW_IMP>1.E-5 )THEN                         !! local inertial equation, see [Bates et al., 2010, J.Hydrol.]
      DOUT_PRE = D1PTHFLW_PRE(IPTH,ILEV) * PTH_WTH(IPTH,ILEV)**(-1.)            !! outflow (t-1) [m2/s] (unit width)
      D1PTHFLW(IPTH,ILEV) = PTH_WTH(IPTH,ILEV) * ( DOUT_PRE + CMF_PARAMS%PGRV*CMF_CONFIG%DT*DFLW_IMP*DSLOPE ) &
                         * ( 1. + CMF_PARAMS%PGRV*CMF_CONFIG%DT*PTH_MAN(ILEV)**2. * abs(DOUT_PRE)*DFLW_IMP**(-7./3.) )**(-1.)
    ELSE
      D1PTHFLW(IPTH,ILEV) = 0._JPRB
    ENDIF
  END DO
END DO
!$OMP END PARALLEL DO

END SUBROUTINE CMF_LEVEE_OPT_PTHOUT
!################################################################

END MODULE CMF_LEVEE_MOD
