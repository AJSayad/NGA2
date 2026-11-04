!> Gas-isentropic mechanical relaxation
module relax_isentropic_class
   use precision,           only: WP
   use thermorelax_class,   only: thermorelax,RELAX_OK,RELAX_FAILED,RELAX_BAD_LIQUID,RELAX_BAD_GAS,RELAX_DEGENERATE
   use material_class,      only: material
   use ideal_gas_class,     only: ideal_gas
   use stiffened_gas_class, only: stiffened_gas
   implicit none
   private

   public :: relax_isentropic

   type, extends(thermorelax) :: relax_isentropic
      class(ideal_gas), pointer :: gas=>null()   !< Gas (ideal or stiffened)
      class(material),  pointer :: liq=>null()   !< Liquid (any equation of state)
      real(WP) :: pinfG=0.0_WP                   !< Gas stiffening pressure (zero for an ideal gas)
      integer  :: nbracket=60                    !< Maximum halvings of the expanding phase fraction to bracket the root
      integer  :: nbisect=200                    !< Maximum bisection steps (stops when the bracket cannot shrink)
   contains
      procedure :: initialize
      procedure :: apply
   end type relax_isentropic

contains

   !> Initialize: store the equations of state, and the gas stiffening pressure if the gas is a stiffened gas
   subroutine initialize(this,gas,liq)
      implicit none
      class(relax_isentropic),  intent(inout) :: this
      class(ideal_gas), target, intent(in)    :: gas
      class(material),  target, intent(in)    :: liq
      this%gas=>gas; this%liq=>liq
      this%pinfG=0.0_WP
      select type (gas)
      class is (stiffened_gas)
         this%pinfG=gas%pinf
      end select
   end subroutine initialize

   !> Relax a mixed cell of liquid fraction VF and phasic state Q(1:4)=[mL,mG,mL*eL,mG*eG] to pL-pG=Pjump
   subroutine apply(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_isentropic), intent(inout) :: this
      real(WP),                intent(in)    :: dt
      real(WP),                intent(inout) :: VF
      real(WP), dimension(:),  intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      integer,  optional,      intent(out)   :: ierr
      real(WP) :: A,lo,hi,mid,flo,fhi,fmid,YG,YL
      integer :: n
      ! Mixed cells with positive phasic masses and internal energies only
      if (VF.le.0.0_WP.or.VF.ge.1.0_WP.or.any(Q(1:4).le.0.0_WP)) then; if (present(ierr)) ierr=RELAX_DEGENERATE; return; end if
      ! Gas isentrope through the current state, (pG+pinfG)*(1-VF)**gamma=A
      A=(this%gas%get_p_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=[1.0_WP])+this%pinfG)*(1.0_WP-VF)**this%gas%gamma
      if (A.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_GAS; return; end if
      ! Bracket the root: halve the fraction of the phase at the higher pressure, which expands, until the sign changes
      lo=VF; hi=VF; flo=dp(VF); fhi=flo
      if (flo.gt.0.0_WP) then
         do n=1,this%nbracket
            hi=1.0_WP-0.5_WP*(1.0_WP-hi); fhi=dp(hi)
            if (fhi.le.0.0_WP) exit
            lo=hi; flo=fhi
         end do
      else if (flo.lt.0.0_WP) then
         do n=1,this%nbracket
            lo=0.5_WP*lo; flo=dp(lo)
            if (flo.ge.0.0_WP) exit
            hi=lo; fhi=flo
         end do
      end if
      if (flo.lt.0.0_WP.or.fhi.gt.0.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      ! Bisection to round-off
      do n=1,this%nbisect
         mid=0.5_WP*(lo+hi)
         if (mid.le.lo.or.mid.ge.hi) exit
         fmid=dp(mid)
         if (fmid.gt.0.0_WP) then
            lo=mid
         else
            hi=mid
         end if
      end do
      mid=0.5_WP*(lo+hi)
      ! Relaxed state: the gas on its isentrope, the liquid with the rest less the surface energy
      YG=gas_content(mid)
      YL=Q(3)+Q(4)-YG-Pjump*(mid-VF)
      if (YL.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      if (this%liq%get_c_from_rho_e(rho=Q(1)/mid,e=YL/Q(1),y=[1.0_WP]).le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      Q(3)=YL; Q(4)=YG; VF=mid
      if (present(ierr)) ierr=RELAX_OK
   contains
      !> Gas internal energy on its isentrope at liquid fraction aL: aG*(pG+gamma*pinfG)/(gamma-1)+mG*qG with aG=1-aL
      real(WP) function gas_content(aL) result(Y)
         implicit none
         real(WP), intent(in) :: aL
         real(WP) :: aG
         aG=1.0_WP-aL
         Y=(A*aG**(1.0_WP-this%gas%gamma)+(this%gas%gamma-1.0_WP)*this%pinfG*aG)/(this%gas%gamma-1.0_WP)+Q(2)*this%gas%q
      end function gas_content
      !> Pressure difference pL-pG-Pjump at liquid fraction aL, the gas on its isentrope and the liquid with the rest
      real(WP) function dp(aL)
         implicit none
         real(WP), intent(in) :: aL
         real(WP) :: Y1
         Y1=Q(3)+Q(4)-gas_content(aL)-Pjump*(aL-VF)
         dp=this%liq%get_p_from_rho_e(rho=Q(1)/aL,e=Y1/Q(1),y=[1.0_WP])-(A*(1.0_WP-aL)**(-this%gas%gamma)-this%pinfG)-Pjump
      end function dp
   end subroutine apply

end module relax_isentropic_class
