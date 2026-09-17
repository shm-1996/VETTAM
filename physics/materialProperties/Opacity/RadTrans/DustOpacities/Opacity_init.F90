!!****f* source/physics/materialProperties/Opacity/RadTrans/Semenov
!!
!! NAME
!!
!!  Opacity_init
!!
!!
!! SYNOPSIS
!!
!!  call Opacity_init()
!!
!! DESCRIPTION
!!
!! Initialiazed data for the Semonov opacity model using run time
!! parameters.
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!
!!***
subroutine Opacity_init()
  use SemenovOpacities, ONLY: init_opacities
  use Opacity_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

  ! Opacity stuff
  !Opacity types for each band
  call RuntimeParameters_get("dusttoGasRatio", dusttoGasRatio)
  call RuntimeParameters_get("dust_opacity_ir",dust_opacity_ir)
  call RuntimeParameters_get("dust_opacity_fuv",dust_opacity_fuv)
  call RuntimeParameters_get("dust_opacity_euv",dust_opacity_euv)
  !Constant opacity values if fixed opacity used
  call RuntimeParameters_get("dust_iropac_value",dust_iropac_value)
  call RuntimeParameters_get("dust_fuvopac_value",dust_fuvopac_value)
  call RuntimeParameters_get("dust_euvopac_value",dust_euvopac_value)
  call RuntimeParameters_get("dust_lwopac_value",dust_lwopac_value)
  call RuntimeParameters_get("dust_peopac_value",dust_peopac_value)

  !Initialize the Semenov IR opacity routines for usage
  call init_opacities()

end subroutine Opacity_init
