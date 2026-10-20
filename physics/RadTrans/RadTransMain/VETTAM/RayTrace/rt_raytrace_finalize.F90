!!****if* source/physics/RadTrans/RadTransMain/VETTAM/RayTrace/rt_raytrace_finalize
!!
!! NAME
!!
!!  rt_raytrace_finalize
!!
!! SYNOPSIS
!!
!!  call rt_raytrace_finalize ()
!!
!! DESCRIPTION
!!
!!  Cleans up the rt unit.
!!
!! ARGUMENTS
!!
!!***
subroutine rt_raytrace_finalize ()

  implicit none
#include "Flash.h"
#include "constants.h"
#ifdef RAYTRACE_3DRT
  call rt_finalize_raytrace_3drt
#endif

  return
end subroutine rt_raytrace_finalize
