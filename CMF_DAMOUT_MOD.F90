MODULE CMF_DAMOUT_MOD
!==========================================================
!* PURPOSE: CaMa-Flood reservoir operation scheme (under development)
!
! (C) R. Hanazaki & D.Yamazaki (U-Tokyo)  Feb 2020
!
!* CONTAINS:
! -- CMF_DEM_NMLIST  : Read setting from namelist
! -- CMF_DAM_INIT    : Initialize dam data
! -- CMF_CALC_DAMOUT : Calculate inflow and outflow at dam
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
USE CMF_UTILS_MOD,           ONLY: INQUIRE_FID
USE CMF_VARS_MOD,            ONLY: I2VECTOR, I1NEXT, NSEQALL, NSEQRIV, NSEQMAX
USE CMF_VARS_MOD,            ONLY: NPTHOUT,  NPTHLEV, PTH_UPST, PTH_DOWN, PTH_ELV, I2MASK!! bifurcation pass
USE CMF_VARS_MOD,            ONLY: D2RIVOUT, D2FLDOUT, P2RIVSTO, P2FLDSTO, P2DAMSTO, P2DAMINF, D2RUNOFF
USE CMF_VARS_MOD,            ONLY: D2RIVINF, D2FLDINF
USE CMF_VARS_MOD,            ONLY: IYYYYMMDD,ISYYYY
USE CMF_VARS_MOD,            ONLY: IDAM, NDAM, NDAMX, DamID, DamName, DamIX, DamIY
USE CMF_VARS_MOD,            ONLY: DamLon, DamLat, upreal, R_VolUpa, Qf, Qn, DamYear
USE CMF_VARS_MOD,            ONLY: DamStat, EmeVol, FldVol, ConVol, NorVol, AdjVol, Qa, DamSeq, I1DAM
!============================
IMPLICIT NONE
SAVE


CONTAINS
!####################################################################
!* CONTAINS:
! -- CMF_DEMOUT_NMLIST  : Read setting from namelist
! -- CMF_DAMOUT_INIT    : Initialize dam data
! -- CMF_DAMOUT_CALC    : Calculate inflow and outflow at dam
! -- CMF_DAMOUT_WATBAL  : Calculate water balance at dam
! -- CMF_DAMOUT_WRITE   : Write dam-related variables in text file
!####################################################################


!####################################################################
SUBROUTINE CMF_DAMOUT_CALC
! local
IMPLICIT NONE
! SAVE for OMP
INTEGER(KIND=JPIM),SAVE    :: ISEQD
!** dam variables
REAL(KIND=JPRB),SAVE       :: DamVol
REAL(KIND=JPRB),SAVE       :: DamInflow
REAL(KIND=JPRB),SAVE       :: DamOutflw           !! Total outflw 
REAL(KIND=JPRB),SAVE       :: DamOutTmp           !! Total outflw 
!$OMP THREADPRIVATE    (ISEQD,DamVol,DamInflow,DamOutflw,DamOutTmp)
!====================
!CONTAINS
!+ UPDATE_INFLOW: replace dam upstream with kinamatic wave, calculate inflow to dam
!+ MODIFY_OUTFLW: modify outflw to avoid negative storage
!+
!==========================================================

!* (1) Replace discharge in upstream grids with kinematic outflow
!     to avoid storage buffer effect (Shin et al., 2019, WRR)
! ------  rivout at upstream grids of dam, rivinf to dam grids are updated.
CALL UPDATE_INFLOW


!* (2) Reservoir Operation
!====================================
!     -- compare DamVol against storage level (NorVol, ConVol, EmeVol) & DamInflow against Qf
!$OMP PARALLEL DO
DO IDAM=1, NDAM
  IF( DamStat(IDAM)<=0 ) CYCLE  !! no calculation for dams not activated

  !! *** 2a update dam volume and inflow -----------------------------------
  ISEQD=DamSeq(IDAM)
  DamVol    = P2DAMSTO(ISEQD,1)    
  DamInflow = P2DAMINF(ISEQD,1)

  !! *** 2b Reservoir Operation          ------------------------------
  !===========================
  IF( CMF_DAM%LDAMH22 )THEN !! Hanazaki 2022 scheme
    !! case1: Water
    IF( DamVol <= NorVol(IDAM) )THEN
      DamOutflw = Qn(IDAM) * (DamVol / ConVol(IDAM) )
    !! case2: water supply
    ELSEIF( NorVol(IDAM)<DamVol .and. DamVol<=ConVol(IDAM) )THEN
      IF( Qf(IDAM)<=DamInflow )THEN
        DamOutflw = Qn(IDAM)*0.5 +   (DamVol-NorVol(IDAM))/( ConVol(IDAM)-NorVol(IDAM))      * (Qf(IDAM) - Qn(IDAM))
      ELSE
        DamOutflw = Qn(IDAM)*0.5 + (((DamVol-NorVol(IDAM))/( EmeVol(IDAM)-NorVol(IDAM)))**2) * (Qf(IDAM) - Qn(IDAM))
      ENDIF  
    !! case3: flood control
    ELSEIF( ConVol(IDAM)<DamVol .and. DamVol<EmeVol(IDAM) ) THEN
      IF( Qf(IDAM) <= DamInflow ) THEN
        DamOutflw = Qf(IDAM) + max((1. - R_VolUpa(IDAM)/0.2),0._JPRB) &
          * (DamVol-ConVol(IDAM))/(EmeVol(IDAM)-ConVol(IDAM)) * (DamInflow-Qf(IDAM))
      !! pre- and after flood control
      ELSE
        DamOutflw = Qn(IDAM)*0.5 + (((DamVol-NorVol(IDAM))/(EmeVol(IDAM)-NorVol(IDAM)))**2)* (Qf(IDAM) - Qn(IDAM))
      ENDIF
    !! case4: emergency operation
    ELSE
      DamOutflw = max(DamInflow, Qf(IDAM))
    ENDIF

  !===========================
  ELSE  !! (not LDAMH22) improved reservoir operation 'Yamazaki & Funato'

    !! Case 1: water use
    IF( DamVol<=ConVol(IDAM) )THEN
      DamOutflw = Qn(IDAM) * (DamVol/ConVol(IDAM))**0.5
    !! case 2: water excess (just avobe ConVol, for outflow stability)
    ELSEIF( DamVol>ConVol(IDAM) .and. DamVol<=AdjVol(IDAM) ) THEN
      DamOutflw = Qn(IDAM) + ( (DamVol-ConVol(IDAM)) / (AdjVol(IDAM)-ConVol(IDAM)) )**3.0 * (Qa(IDAM) - Qn(IDAM))
    !! case 3: water excess
    ELSEIF( DamVol>AdjVol(IDAM) .and. DamVol<=EmeVol(IDAM) ) THEN
      !! (flood period)
      IF( DamInflow >= Qf(IDAM) ) THEN
        !!figure left side No.2
        DamOutflw = Qn(IDAM) + ( (DamVol-ConVol(IDAM)) / (EmeVol(IDAM)-ConVol(IDAM)) )**1.0 * (DamInflow- Qn(IDAM))
        DamOutTmp = Qa(IDAM) + ( (DamVol-AdjVol(IDAM)) / (EmeVol(IDAM)-AdjVol(IDAM)) )**0.1 * (Qf(IDAM) - Qa(IDAM))
        DamOutflw = max( DamOutflw,DamOutTmp )
      !! (non-flood period)
      ELSE
        DamOutflw = Qa(IDAM) + ( (DamVol-AdjVol(IDAM)) / (EmeVol(IDAM)-AdjVol(IDAM)) )**0.1 * (Qf(IDAM) - Qa(IDAM))
      ENDIF
    !! case 4: emergency operation(no.1)
    ELSEIF( DamVol>EmeVol(IDAM) )THEN
      !! (flood period)
      IF( DamInflow >= Qf(IDAM) ) THEN
        DamOutflw = DamInflow
      !! (non-flood period)
      ELSE
        DamOutflw = Qf(IDAM)
      ENDIF
    ENDIF
  ENDIF

  !! *** 2c flow limitter
  DamOutflw = min( DamOutflw, DamVol/CMF_CONFIG%DT, real(P2RIVSTO(ISEQD,1)+P2FLDSTO(ISEQD,1),JPRB)/CMF_CONFIG%DT )
  DamOutflw = max( DamOutflw, 0._JPRB )

  !! update CaMa variables  (treat all outflow as RIVOUT in dam grid, no fldout)
  D2RIVOUT(ISEQD,1) = DamOutflw
  D2FLDOUT(ISEQD,1) = 0._JPRB
END DO
!$OMP END PARALLEL DO
!====================================

CONTAINS
!==========================================================
!+ UPDATE_INFLOW: replace dam upstream with kinamatic wave, calculate inflow to dam
!+ MODIFY_OUTFLW: modify outflw to avoid negative storage
!+
!==========================================================
SUBROUTINE UPDATE_INFLOW
USE CMF_VARS_MOD,        ONLY: D2RIVLEN, D2RIVMAN, D2ELEVTN, D2NXTDST, D2RIVWTH
USE CMF_VARS_MOD,        ONLY: D2RIVOUT_PRE, D2FLDOUT_PRE
USE CMF_VARS_MOD,       ONLY: D2RIVDPH, D2RIVVEL, D2FLDDPH
IMPLICIT NONE
! SAVE for OpenMP
INTEGER(KIND=JPIM),SAVE    :: ISEQ, JSEQ
REAL(KIND=JPRB),SAVE       :: DSLOPE,DAREA,DVEL,DSLOPE_F,DARE_F,DVEL_F
!$OMP THREADPRIVATE     (JSEQ,DSLOPE,DAREA,DVEL,DSLOPE_F,DARE_F,DVEL_F)
!============================

!*** 1a. reset outflw & dam inflow
!$OMP PARALLEL DO
DO ISEQ=1, NSEQALL
  IF( I1DAM(ISEQ)>0 )THEN  !! if dam grid or upstream of dam, reset variables
    D2RIVOUT(ISEQ,1) = 0._JPRB
    D2FLDOUT(ISEQ,1) = 0._JPRB
    P2DAMINF(ISEQ,1) = 0._JPRD
  ENDIF
END DO
!$OMP END PARALLEL DO

!*** 1b. calculate dam inflow, using previous tstep discharge
#ifndef NoAtom_CMF
!$OMP PARALLEL DO  !! No OMP Atomic for bit-identical simulation (set in Mkinclude)
#endif
DO ISEQ=1, NSEQALL
  IF( I1DAM(ISEQ)==10 .or. I1DAM(ISEQ)==11 )THEN  !! if dam grid or upstream of dam
    JSEQ=I1NEXT(ISEQ)
#ifndef NoAtom_CMF
!$OMP ATOMIC
#endif
    P2DAMINF(JSEQ,1) = P2DAMINF(JSEQ,1) + D2RIVOUT_PRE(ISEQ,1) + D2FLDOUT_PRE(ISEQ,1) 
  ENDIF
END DO
#ifndef NoAtom_CMF
!$OMP END PARALLEL DO
#endif

!*** 1c. discharge for upstream grids of dams
!$OMP PARALLEL DO  !! No OMP Atomic for bit-identical simulation (set in Mkinclude)
DO ISEQ=1, NSEQRIV
  IF( I1DAM(ISEQ)==10 )THEN  !! if downstream is DAM
    JSEQ   = I1NEXT(ISEQ)
    ! === river flow
    DSLOPE = (D2ELEVTN(ISEQ,1)-D2ELEVTN(JSEQ,1)) * D2NXTDST(ISEQ,1)**(-1.)
    DSLOPE = max(DSLOPE,CMF_PARAMS%PMINSLP)

    DVEL   = D2RIVMAN(ISEQ,1)**(-1.) * DSLOPE**0.5 * D2RIVDPH(ISEQ,1)**(2./3.)
    DAREA  = D2RIVWTH(ISEQ,1) * D2RIVDPH(ISEQ,1)

    D2RIVVEL(ISEQ,1) = DVEL
    D2RIVOUT(ISEQ,1) = DAREA * DVEL
    D2RIVOUT(ISEQ,1) = MIN( D2RIVOUT(ISEQ,1), real(P2RIVSTO(ISEQ,1),JPRB)/CMF_CONFIG%DT )
    !=== floodplain flow
    DSLOPE_F = min( 0.005_JPRB,DSLOPE )    !! set min [instead of using weirequation for efficiency]
    DVEL_F   = CMF_PARAMS%PMANFLD**(-1.) * DSLOPE_F**0.5 * D2FLDDPH(ISEQ,1)**(2./3.)
    DARE_F   = P2FLDSTO(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.)
    DARE_F   = MAX( DARE_F - D2FLDDPH(ISEQ,1)*D2RIVWTH(ISEQ,1), 0._JPRB )   !!remove above river channel     area

    D2FLDOUT(ISEQ,1) = DARE_F * DVEL_F
    D2FLDOUT(ISEQ,1) = MIN(  D2FLDOUT(ISEQ,1)*1._JPRD, P2FLDSTO(ISEQ,1)/CMF_CONFIG%DT )
  ENDIF
END DO
!$OMP END PARALLEL DO

END SUBROUTINE UPDATE_INFLOW
!==========================================================
END SUBROUTINE CMF_DAMOUT_CALC
!####################################################################
!
!
!
!####################################################################
SUBROUTINE CMF_DAMOUT_WATBAL
IMPLICIT NONE
! SAVE for OMP
INTEGER(KIND=JPIM),SAVE    :: ISEQD
!*** water balance
REAL(KIND=JPRB),SAVE       :: DamInflow
REAL(KIND=JPRB),SAVE       :: DamOutflw           !! Total outflw 
REAL(KIND=JPRD),SAVE       :: GlbDAMSTO, GlbDAMSTONXT, GlbDAMINF, GlbDAMOUT, DamMiss
!$OMP THREADPRIVATE    (ISEQD,DamInflow,DamOutflw)
! ==========================================
!* 4) update reservoir storage and check water DamMiss --------------------------
GlbDAMSTO    = 0._JPRB
GlbDAMSTONXT = 0._JPRB
GlbDAMINF    = 0._JPRB
GlbDAMOUT    = 0._JPRB

!$OMP PARALLEL DO REDUCTION(+:GlbDAMSTO, GlbDAMSTONXT, GlbDAMINF, GlbDAMOUT)
DO IDAM=1, NDAM
  IF( DamStat(IDAM)==CMF_PARAMS%IMIS ) CYCLE  !! do not calculate for dams outside of calculation domain
  ISEQD = DamSeq(IDAM)

  DamInflow = D2RIVINF(ISEQD,1) + D2FLDINF(ISEQD,1) + D2RUNOFF(ISEQD,1)
  DamOutflw = D2RIVOUT(ISEQD,1) + D2FLDOUT(ISEQD,1)
!!P2DAMINF(ISEQD,1)=DamInflow   !! if water balance needs to be checked in the output file, P2DAMINF should be updated.

  GlbDAMSTO = GlbDAMSTO + P2DAMSTO(ISEQD,1)
  GlbDAMINF = GlbDAMINF + DamInflow*CMF_CONFIG%DT
  GlbDAMOUT = GlbDAMOUT + DamOutflw*CMF_CONFIG%DT

  P2DAMSTO(ISEQD,1) = P2DAMSTO(ISEQD,1) + DamInflow * CMF_CONFIG%DT - DamOutflw * CMF_CONFIG%DT

  GlbDAMSTONXT = GlbDAMSTONXT + P2DAMSTO(ISEQD,1)
END DO
!$OMP END PARALLEL DO

DamMiss = GlbDAMSTO-GlbDAMSTONXT+GlbDAMINF-GlbDAMOUT
WRITE(CMF_FILES%LOGNAM,*) "CMF::DAM_CALC: DamMiss at all dams:", DamMiss*1.D-9

END SUBROUTINE CMF_DAMOUT_WATBAL
!####################################################################
!
!
!

END MODULE CMF_DAMOUT_MOD
