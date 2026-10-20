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

  !! *** Runtime Parameters *** !!
  real, save    :: sim_rhoAmbient, sim_gamma,  &
                   sim_eddvalue, sim_nu,sim_T0,sim_T1, sim_vx

end module Simulation_data


