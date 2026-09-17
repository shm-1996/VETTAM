# VETTAM Radiation-Hydrodynamics Module

VETTAM couples multi-band radiation hydrodynamics (RHD) to FLASH by solving the two-moment equations of RHD closed with a Variable Eddington Tensor (VET) obtained from a ray-trace solution. The temporal update is fully implicit and occurs at the cadence of the FLASH hydro (or other unit) timescale, and utilizes the PETSc sparse matrix solvers.

**Reference:** Menon et al. 2022 (MNRAS 512, 401) for VETTAM details and validation.

## Contents

- [Overview](#overview)
- [Architecture](#architecture)
- [User Features](#user-features)
- [Runtime Parameters](#runtime-parameters)
- [Build & Setup](#build--setup)
- [Outputs & Logging](#outputs--logging)
- [Contact & Citation](#contact--citation)

## Architecture

**Core workflow**

- `RadTrans_init.F90`: Read runtime parameters, set constants, initialise PETSc, configure radiation boundary conditions.
- `RadTrans.F90`: Main driver unit, looping over active bands; calls the internal routines for each component of the unit.
- `RadTrans_finalize.F90`: Deallocate memory for PETSc allocations and finalize.

**Solver path**

- `rt_petsc.F90`: Allocate memory for PETSc vectors/matrices (`Petsc_alloc`), configure KSP/PC, reallocate on AMR mesh changes (`grid_changed`), and subcycle solves on non-convergence (`Radiation_step`).
- `rt_fillMatrix.F90`: Assemble the linear system using current opacities, Eddington tensors, and wave speeds; applies boundary conditions and guard-cell data.

**Physics support**

- Opacities: `rt_setOpacity.F90` (dust/gas, band-dependent; includes multispecies/ionization cross sections).
- Eddington tensors: `rt_setEddTensor.F90` for closure and wave speeds.
- Dust coupling: `rt_dustTemperature.F90`, `rt_dustTerms.F90` for dust temperature evolution and momentum/energy exchange.
- Sinks/sources: Include contribution from sink particles with `rt_sinkInject.F90` (Gaussian source), `rt_sinkHydro.F90` (gas coupling), optional ray-trace paths for sinks.
- Optional submodules: `PhotoIonization/rt_ionise.F90`, `RayTrace/` (3D VET or direct transport), `LWBand/` (LW dissociation).
- Radiation variables (`erad`, `moh*`, `tau*`, `xxed`, etc.) live in the main solution arrays and are reused per band; define `BnER`, `BnFX`, … in your simulation `Config` if you need persistent per-band storage.

## User Features

- **Bands:** Up to 5 (`rt_band1`…`rt_band5`). LW can run via ray tracing or moments depending on `lwdiss_type` and `rt_etens`/ray-trace flags.
- **Ray tracing:** Enable `RAYTRACE_3DRT` for VET computation or direct-source transport. `rt_etens='vet'` requires ray tracing; otherwise execution aborts.
- **Sinks:** With `SINK_PART_TYPE`, choose Gaussian sources or ray tracing (`sink_useRayTrace`). Control heating/momentum via `rt_sinkheat`, `rt_sinkmom`, and implicit treatment via `rt_sink_implicit`.
- **Hydro coupling:** Toggle with `rt_update_hydro`. Temperature updates via Picard or Newton-Raphson (`rt_hydro_type`, `rt_temp_type`). O(v/c) terms via `rt_ovcterms`.
- **Call throttling:** `rt_dt_maxvel` skips radiation updates if `tsinceLastCall < hy_cfl * dx / rt_dt_maxvel` (default resolves light speed).
- **Boundaries:** Face-specific BCs (`xl_rad_BC`, …) with fixed-T/flux options (`rt_boundary_T`, `frad*_boundary_value`).

## Runtime Parameters

Key parameters from `Config` (defaults in parentheses):

| Category | Parameters |
| --- | --- |
| **Control** | `rt_compute_Dt` (FALSE), `rt_debug` (FALSE), `rt_etens` ("vet"), `rt_dt_maxvel` (3e10), `rt_T_min` (0.0) |
| **Timestep** | `rt_Dt_type` (2), `rt_dtFactor` (1.0), `rt_hlle_correction` (TRUE) |
| **Solver** | `rt_rtol` (1e-8), `rt_abstol` (1e-20), `rt_dtol` (1e5), `rt_maxits` (1000), `rt_minits` (1), `rt_preconditioner` ("asm"), `rt_linearsolver` ("gmres") |
| **Floors** | `rt_smalle` (0.0), `rt_smallf` (1e-20) |
| **Hydro coupling** | `rt_update_hydro` (TRUE), `rt_hydro_type` (2), `rt_temp_type` (1), `rt_picard_rtol` (1e-3), `rt_picard_abstol` (1e-10), `rt_picard_maxits` (1000), `rt_picard_theta` (1.0), `rt_nr_rtol` (1e-4), `rt_nr_dtol` (1e5), `rt_nr_maxits` (10000), `rt_picard_correction` (TRUE), `rt_ovcterms` (TRUE) |
| **Sinks** | `rt_sink_implicit` (FALSE), `rt_sinkheat` (TRUE), `rt_sinkmom` (TRUE), `sink_useRayTrace` (FALSE), `sigma_star` (2.0), `noSink_VETTAM` (TRUE) |
| **Boundaries** | `xl_rad_BC`, `xr_rad_BC`, `yl_rad_BC`, `yr_rad_BC`, `zl_rad_BC`, `zr_rad_BC` (all 2 = Outflow), `rt_boundary_T` (10.0), `frad*_boundary_value` (0.0) |
| **Bands** | `rt_freqbands` (3), `rt_band1`–`rt_band5` (IR, FUV, EUV, PE, LYMAN_WERNER) |

## Build & Setup

- **Dependencies:** Requires PETSc and FLASH Driver; dust opacity tables via `physics/materialProperties/Opacity/RadTrans/DustOpacities`. Conflicts with legacy MGD radiation module.
- **Enabling VETTAM:** Include `physics/RadTrans/RadTransMain/VETTAM` in your simulation setup and provide PETSc paths. Define extra band-specific output variables in your simulation `Config` if needed.
- **Ray tracing flags:** Build with `RAYTRACE_3DRT` for VETs or direct transport. LW dissociation (`FSHM_VAR`) and photoionization (`UEUV_VAR`) may also use ray tracing; ensure matching compile-time flags (`SINK_RT` for sink-raytrace).
- **Runtime tuning:** Adjust `rt_rtol`/`rt_abstol` for accuracy/speed; inspect per-band convergence via `Print_PETSCConvergence`. Use `rt_dt_maxvel` to lower call frequency when light-speed resolution is unnecessary.

## Outputs & Logging

- Per-band KSP convergence (cumulative/min/max iterations, reasons) written to the FLASH logfile. Additional PETSc output can be enabled through `-ksp_*` runtime options or `rt_debug`.
- Guard-cell masks for Eddington tensors/opacities are logged once (`Logfile_stampVarMask`); hydro guard fills log when enabled.

## Contact & Citation

- Shyam Harimohan Menon — smenon@flatironinstitute.org
- Please cite **Menon et al. 2022, MNRAS 512, 401** when using VETTAM.
