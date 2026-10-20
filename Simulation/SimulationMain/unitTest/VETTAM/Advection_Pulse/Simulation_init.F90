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
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "Flash.h"
#include "constants.h"

  integer :: dr_myPE

  call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
  call RuntimeParameters_get('sim_T0', sim_T0)
  call RuntimeParameters_get('sim_T1', sim_T1)
  call RuntimeParameters_get('sim_nu', sim_nu)
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('sim_eddvalue',sim_eddvalue)
  call RuntimeParameters_get('sim_vx',sim_vx)
  call Driver_getMype(GLOBAL_COMM, dr_myPE)

end subroutine Simulation_init
