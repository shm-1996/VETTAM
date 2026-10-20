!!****if* source/physics/RadTrans/RadTransMain/VETTAM/RadTrans_init
!!
!! NAME
!!
!!  RadTrans_init
!!
!! SYNOPSIS
!!
!!  RadTrans_init()
!!
!! DESCRIPTION
!! 
!! Initialise the VETTAM module.
!!
!! ARGUMENTS
!!
!!   
!!
!! AUTOGENROBODOC
!!
!!
!!***

subroutine RadTrans_init()

  use RadTrans_data
  use Driver_interface, ONLY            : Driver_abortFlash, Driver_getMype, Driver_getComm
  use Driver_data, ONLY                 : dr_globalMe
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

#include "constants.h"
#include "Flash.h"
#include "petsc/finclude/petscksp.h"

  implicit none
  PetscErrorCode :: ierr
  INTEGER        :: bc
!==============================================================================

  call RuntimeParameters_get ("useRadTrans", rt_useRadTrans)
  call RuntimeParameters_get("rt_debug",rt_debug)

  
  call RuntimeParameters_get("rt_T_min",rt_T_min)

  !Flag for HLLE Correction (see Jiang et al 2013)
  call RuntimeParameters_get("rt_hlle_correction",use_hlle_correction)

  ! these are not used in the iterative implementation (no timestep based on cooling)
  call RuntimeParameters_get("rt_compute_Dt",rt_compute_Dt)
  call RuntimeParameters_get("rt_Dt_type",rt_Dt_type)
  call RuntimeParameters_get("rt_dtFactor",rt_dtFactor)

  ! Limits to erad/flux
  call RuntimeParameters_get("rt_smalle",rt_smalle)
  call RuntimeParameters_get("rt_smallf",rt_smallf)

  ! PETSc Stuff
  call RuntimeParameters_get("rt_rtol", rt_rtol)
  call RuntimeParameters_get("rt_abstol", rt_abstol)
  call RuntimeParameters_get("rt_dtol", rt_dtol)
  call RuntimeParameters_get("rt_maxits", rt_maxits)
  call RuntimeParameters_get("rt_minits", rt_minits)
  call RuntimeParameters_get("rt_preconditioner", rt_preconditioner_str)
  call RuntimeParameters_get("rt_linearsolver", rt_linearsolver_str)

  !Multiband stuff
  call RuntimeParameters_get("rt_freqbands",rt_freqbands)
  call RuntimeParameters_get("rt_band1",rt_band1)
  call RuntimeParameters_get("rt_band2",rt_band2)
  call RuntimeParameters_get("rt_band3",rt_band3)
  call RuntimeParameters_get("rt_band4",rt_band4)
  call RuntimeParameters_get("rt_band5",rt_band5)

  !Eddington Tensor computation flag
  call RuntimeParameters_get("rt_etens",rt_etens)

  ! Abort early if M1 closure is requested with multiple bands but per-band storage variables
  ! (BnER, BnFX, BnFY, BnFZ) are not compiled in. These variables are required because
  ! ERAD_VAR/MOH*_VAR are overwritten each band, so the M1 tensor must read from the
  ! per-band persistent storage instead. Define them in your simulation Config.
  if ((trim(rt_etens) .eq. 'm1' .or. trim(rt_etens) .eq. 'mom1') .and. rt_freqbands .gt. 1) then
#ifndef B1ER_VAR
    call Driver_abortFlash("[VETTAM] rt_etens=m1 with rt_freqbands>1 requires B1ER, B1FX[,B1FY,B1FZ] in simulation Config.")
#endif
    if (rt_freqbands .ge. 2) then
#ifndef B2ER_VAR
      call Driver_abortFlash("[VETTAM] rt_etens=m1 with rt_freqbands>=2 requires B2ER, B2FX[,B2FY,B2FZ] in simulation Config.")
#endif
    endif
    if (rt_freqbands .ge. 3) then
#ifndef B3ER_VAR
      call Driver_abortFlash("[VETTAM] rt_etens=m1 with rt_freqbands>=3 requires B3ER, B3FX[,B3FY,B3FZ] in simulation Config.")
#endif
    endif
    if (rt_freqbands .ge. 4) then
#ifndef B4ER_VAR
      call Driver_abortFlash("[VETTAM] rt_etens=m1 with rt_freqbands>=4 requires B4ER, B4FX[,B4FY,B4FZ] in simulation Config.")
#endif
    endif
    if (rt_freqbands .ge. 5) then
#ifndef B5ER_VAR
      call Driver_abortFlash("[VETTAM] rt_etens=m1 with rt_freqbands>=5 requires B5ER, B5FX[,B5FY,B5FZ] in simulation Config.")
#endif
    endif
  endif

  !O(v/c) terms correction flag
  call RuntimeParameters_get("rt_ovcterms",rt_ovcterms)

  !Hydro
  call RuntimeParameters_get("rt_update_hydro",rt_update_hydro)  

  !Flags for sink heating/momentum contributions
  call RuntimeParameters_get("rt_sinkheat",rt_sinkheat)
  call RuntimeParameters_get("rt_sinkmom",rt_sinkmom)

  !Implicit sink heating or not
  call RuntimeParameters_get("rt_sink_implicit",rt_sink_implicit)
  
  !Picard temperature correction
  call RuntimeParameters_get("rt_picard_correction",rt_picard_correction)

  !Flag for whether to use VETTAM if there are no sinks
  call RuntimeParameters_get("noSink_VETTAM",noSink_VETTAM)

  !Parameters to control frequency of VETTAM (only use rt every dt that captures the signal speed rt_dt_maxvel; default 0)
  call RuntimeParameters_get("rt_dt_maxvel",rt_dt_maxvel)

  if(rt_update_hydro) then
    !Picard iteration stuff
    call RuntimeParameters_get("rt_hydro_type",rt_hydro_type)  
    call RuntimeParameters_get("rt_temp_type",rt_temp_type)  
    if(rt_hydro_type .gt. 1) then
      call RuntimeParameters_get("rt_picard_rtol",rt_picard_rtol)
      call RuntimeParameters_get("rt_picard_abstol",rt_picard_abstol)
      call RuntimeParameters_get("rt_picard_theta",rt_picard_theta)  
      call RuntimeParameters_get("rt_picard_maxits",rt_picard_maxits) 
      if(rt_temp_type .eq. 1) then
        call RuntimeParameters_get("rt_nr_rtol",rt_nr_rtol)
        call RuntimeParameters_get("rt_nr_dtol",rt_nr_dtol)
        call RuntimeParameters_get("rt_nr_maxits",rt_nr_maxits)
      endif
    endif
  endif

  

  ! Boundary conditions
  call RuntimeParameters_get("xl_rad_BC", xl_rad_BC)
  rad_BC(LOW,IAXIS) = xl_rad_BC
  call RuntimeParameters_get("xr_rad_BC", xr_rad_BC)
  rad_BC(HIGH,IAXIS) = xr_rad_BC
!#if NDIM> CONSTANT_ONE
#if K2D
  call RuntimeParameters_get("yl_rad_BC", yl_rad_BC)
  rad_BC(LOW,JAXIS) = yl_rad_BC
  call RuntimeParameters_get("yr_rad_BC", yr_rad_BC)
  rad_BC(HIGH,JAXIS) = yr_rad_BC
#if NDIM>2  
  call RuntimeParameters_get("zl_rad_BC", zl_rad_BC)
  rad_BC(LOW,KAXIS) = zl_rad_BC
  call RuntimeParameters_get("zr_rad_BC", zr_rad_BC)
  rad_BC(HIGH,KAXIS) = zr_rad_BC
#endif
#endif 
  call RuntimeParameters_get("rt_boundary_T", rt_boundary_T)

  call RuntimeParameters_get("fradx_boundary_value", fradx_boundary_value)
#if NDIM>1  
  call RuntimeParameters_get("frady_boundary_value", frady_boundary_value)
#if NDIM>2
  call RuntimeParameters_get("fradz_boundary_value", fradz_boundary_value)
#endif
#endif  

  ! Obtain physical constants
  call PhysicalConstants_get("speed of light",rt_speedlt)
  call PhysicalConstants_get("Stefan-Boltzmann",rt_boltz)
  ! Radiation constant (often: aR)
  rt_radconst = 4. * rt_boltz/rt_speedlt 

  if (.not. rt_useRadTrans .and. dr_globalMe .eq. MASTER_PE) then
     write(6,*) 'WARNING:  You have included the RadTrans unit but have set '
     write(6,*) '   the runtime parameter useRadTrans to FALSE'
     write(6,*) '   No cooling will occur but RadTrans_init will continue.'
  end if
  
  PETSC_COMM_WORLD=MPI_COMM_WORLD
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  IF(ierr.NE.0) call Driver_abortFlash("Error while initializing Petsc.")

  blockCount = 0
  alloced = .FALSE.

  !Legacy variables
  call Driver_getMype(MESH_COMM,rt_meshMe)
  call RuntimeParameters_get("meshCopyCount", rt_meshCopyCount)  
  call Driver_getMype(MESH_ACROSS_COMM, rt_acrossMe)
  call Driver_getComm(MESH_ACROSS_COMM, rt_acrossComm)
  call Driver_getComm(GLOBAL_COMM,rt_globalComm)

#ifdef RAYTRACE_3DRT
  !Initialise variables for ray-tracer
  call rt_raytrace_init()
#endif

  return

end subroutine RadTrans_init
