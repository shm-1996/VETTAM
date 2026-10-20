!!****if* source/Simulation/SimulationMain/StirFromFile/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!   use Simulation_data
!!
!! DESCRIPTION
!!  Stores the local data for Simulation setup: StirFromFile
!!
!!***


module Simulation_data

  implicit none
#include "Flash.h"
#include "constants.h"
  !! *** Runtime Parameters *** !!
  real, save    :: sim_rhoAmbient, sim_gamma, sim_erAmbient, &
                   sim_eddvalue, sim_nu, sim_mu, sim_vx
  character(len=MAX_STRING_LENGTH), save :: sim_type

end module Simulation_data


