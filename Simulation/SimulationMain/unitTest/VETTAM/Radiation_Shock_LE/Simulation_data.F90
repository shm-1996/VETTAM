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
  real, save    :: sim_rholeft,sim_rhoright,sim_Tleft,sim_Tright,sim_vleft,sim_vright,&
                   sim_eddvalue,sim_P0,sim_sigma,sim_kappa

end module Simulation_data


