module CMF_SEDPAR_MOD
!==========================================================
!* PURPOSE: parameters for sediment transport
! (C) M.Hatono  (Hiroshima-U)  May 2021
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
  use PARKIND1,                only: JPIM, JPRB
  use CMF_NMLIST_MOD,          only: CMF_PARAMS,CMF_SED
  use CMF_VARS_MOD,            only: d2sedfrc,sDiam

contains
!####################################################################
! -- csVel
! -- sVel
! -- susVel
!####################################################################
function calc_settingVelocity() result(setVel)
  implicit none
  save
  real(kind=JPRB)                 ::  setVel( CMF_SED%nsed) ! setting velocity [m/s]
  real(kind=JPRB)                 ::  sTmp( CMF_SED%nsed)

  sTmp(:) = 6.d0 * CMF_SED%visKin / sDiam(:)
  setVel(:) = CMF_SED%pset * ( sqrt( 2.d0/3.d0*(CMF_SED%psedD-CMF_SED%pwatD)/CMF_SED%pwatD*CMF_PARAMS%PGRV*sDiam(:) &
                             + sTmp(:)*sTmp(:) ) - sTmp(:) )
end function calc_settingVelocity
!=====================================================

function calc_criticalShearVelocity(diam) result(csVel)
  implicit none
  save
  real(kind=JPRB)                 ::  csVel ! critical shear velocity[(cm/s)^2]
  real(kind=JPRB), intent(in)     ::  diam ![m]
  real(kind=JPRB)                 ::  cA, cB
  !========
  cB = 1.d0
  if ( diam >= 0.00303d0 ) then
    cA = 80.9d0
  else if ( diam >= 0.00118d0 ) then
    cA = 134.6d0
    cB = 31.d0 / 32.d0
  else if ( diam >= 0.000565d0 ) then
    cA = 55.d0
  else if ( diam >= 0.000065d0 ) then
    cA = 8.41d0
    cB = 11.d0 / 32.d0
  else
    cA = 226.d0
  endif
  
  csVel = cA * ( diam*100.d0 ) ** cB
  return
end function calc_criticalShearVelocity
!=====================================================

function calc_shearVelocity(rivvel,rivdph) result(sVel)
  implicit none
  save
  real(kind=JPRB)                 ::  sVel ! shear velocity[m/s]
  real(kind=JPRB), intent(in)     ::  rivvel, rivdph 
  !========

  sVel = sqrt ( CMF_PARAMS%PGRV * CMF_PARAMS%PMANRIV**2.d0 * rivvel**2.d0 * rivdph**(-1.d0/3.d0) )  !bug fix 2022/11/22
  return
end function calc_shearVelocity
!=====================================================

function calc_suspendVelocity(csVel,sVel,setVel) result(susVel) ! Uchida and Fukuoka (2019) Eq.44
  implicit none
  save
  real(kind=JPRB)                 ::  susVel( CMF_SED%nsed)   ! suspend velocity [m/s]
!!  real(kind=JPRB), intent(in)     ::  csVel(nsed), sVel, setVel(nsed)    ! crit shear, shear velocity, setting velocity [m/s]
  real(kind=JPRB), intent(in)     ::  csVel(:), sVel, setVel(:)    ! crit shear, shear velocity, setting velocity [m/s]
  integer(kind=JPIM)              ::  ised
  real(kind=JPRB)                 ::  alpha, a, cB, sTmp
  !========
  
  !--------------!
  ! set constant !
  !--------------!
  alpha = CMF_SED%vonKar / 6.d0
  a = 0.08d0
  cB = 1.d0 - CMF_SED%lambda
  
  !-----------------------!
  ! calc suspend velocity !
  !-----------------------!
  susVel(:) = 0.d0
  do ised = 1, CMF_SED%nsed
    if ( csVel(ised) > sVel ) cycle
    sTmp = setVel(ised) / alpha / sVel
    susVel(ised) = max( setVel(ised) * cB / (1.d0+sTmp) * (1.d0-a*sTmp) / (1.d0+(1.d0-a)*sTmp), 0.d0 )
  enddo
end function calc_suspendVelocity

!####################################################################

end module CMF_SEDPAR_MOD
