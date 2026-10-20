!!****if* source/Simulation/SimulationMain/StirFromFile/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!! SYNOPSIS
!!  Simulation_init()
!!
!! ARGUMENTS
!!
!! DESCRIPTION
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!
!!***

subroutine Simulation_init()

  use Driver_interface, ONLY: Driver_getMype
  use Simulation_data
  use RadTrans_data, ONLY: rt_speedlt, rt_radconst, rt_boundary_T_face
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, RuntimeParameters_set
  use Eos_data, ONLY : eos_gasConstant, eos_gamma

  implicit none

#include "Flash.h"
#include "constants.h"

  integer :: dr_myPE
  ! Dimensional constants
  real :: rho0,T0,cs0,L0,sigma_a

  call RuntimeParameters_get('sim_rholeft', sim_rholeft)
  call RuntimeParameters_get('sim_rhoright', sim_rhoright)
  call RuntimeParameters_get('sim_Tleft', sim_Tleft)
  call RuntimeParameters_get('sim_Tright', sim_Tright)
  call RuntimeParameters_get('sim_eddvalue',sim_eddvalue)
  call RuntimeParameters_get('sim_vleft',sim_vleft)
  call RuntimeParameters_get('sim_vright',sim_vright)
  call RuntimeParameters_get('sim_P0',sim_P0)
  call RuntimeParameters_get('sim_kappa',sim_kappa)
  call RuntimeParameters_get('sim_sigma',sim_sigma)




  !! Dimensionalise inputs
  cs0 = SQRT(rt_speedlt**2/(3.*sim_sigma*sim_kappa))
  T0 = cs0**2/(eos_gamma*eos_gasConstant)
  rho0 = rt_radconst*T0**4/(sim_P0*cs0**2)
  L0 = 1.0

  ! rho0 = 1.0
  ! T0 = ((sim_P0*rho0*eos_gamma*eos_gasConstant)/rt_radconst)**(1./3.)
  ! cs0 = SQRT(eos_gamma*eos_gasConstant*T0)
  ! L0 = 1.0 

  sim_rholeft = sim_rholeft*rho0
  sim_rhoright = sim_rhoright*rho0
  sim_Tleft = sim_Tleft*T0
  sim_Tright = sim_Tright*T0
  sim_vleft = sim_vleft*cs0
  sim_vright = sim_vright*cs0

  ! Dimensional sigma_a = sigma_a0*cs0/(L0*c)
  sigma_a =  sim_sigma*cs0/(L0*rt_speedlt)
  !The absorption coefficient and the two boundary temperatures are derived
  !from the dimensionless Lowrie & Edwards parameters above, so the simulation
  !hands them to the solver rather than duplicating the arithmetic in the
  !parfile. Previously this was done by a private copy of cool_vet.F90.
  !
  !dust_iropac_value goes through the runtime parameter because Opacity_init
  !runs after Simulation_init and would overwrite a direct assignment;
  !dust_opacity_ir="Constant" makes it a total extinction coefficient.
  call RuntimeParameters_set("dust_iropac_value", sigma_a)

  !RadTrans_init has already run, so these can be set directly.
  rt_boundary_T_face(LOW, IAXIS)  = sim_Tleft
  rt_boundary_T_face(HIGH,IAXIS)  = sim_Tright

  call Driver_getMype(GLOBAL_COMM, dr_myPE)

end subroutine Simulation_init
