!! Photochemical rates from the LW band: Photodissociation of H2 and CO, and the Photoionization of C
!! AUTHOR
!!  Shyam Harimohan Menon (2024)
!!***
!! ARGUMENTS
!!
!!  blockCount : The number of blocks in the list
!!  blockList(:) : The list of blocks on which to apply the cooling operator
!!  dt : the current timestep
!! 
!!
!!***

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
SUBROUTINE rt_lw(blockCount_,blockList_,dt)
  use Grid_interface
  use Eos_interface
  use Grid_data, ONLY: gr_smallx
  use Hydro_data, ONLY: hy_useHydro
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use rt_ionisedata, ONLY: h2A, hA
  use rt_lwdata
#if defined C_SPEC || defined CO_SPEC
  use Multispecies_interface, ONLY : Multispecies_getProperty
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
#endif

  implicit none

  integer, intent(IN)                        :: blockCount_
  integer, dimension(blockCount_), intent(IN) :: blockList_
  real, intent(IN)                           :: dt
  integer                                    :: b,blockID, nAllParticles, myPE, i, j, k
  integer,dimension(2,MDIM)                  :: blkLimits, blkLimitsGC
  real,dimension(:,:,:,:),pointer            :: solnData
  logical                                    :: gcmask(NUNK_VARS), gcMaskLogged = .false.
  logical, save                              :: first_call = .true.
  real, dimension(2,MDIM)                    :: boundBox
#if defined C_SPEC || defined CO_SPEC
  real                                       :: mH
#endif

  if(first_call) then
#ifndef H2_SPEC
  call Driver_abortFlash("[rt_lw]: The Multispecies variable H2_SPEC needs to be defined for this implementation.")
#endif
    !Read runtime params
    call RuntimeParameters_get("useH2Dissociate",useH2Dissociate)
    call RuntimeParameters_get("useHIshield",useHIshield)
    call RuntimeParameters_get("lwdiss_type",lwdiss_type)
    call RuntimeParameters_get("bfive",bfive)
    call RuntimeParameters_get("fpump",fpump)
    call RuntimeParameters_get("energyPerDissociation",energyPerDissociation)
    call RuntimeParameters_get("avgEnergyLW",avgEnergyLW)
    call RuntimeParameters_get("lw_lshield_max",lw_lshield_max)

    if(lwdiss_type .eq. "raytrace" .or. lwdiss_type .eq. "hybrid") then
#ifndef RAYTRACE_3DRT
      call Driver_abortFlash("Ray-tracer needs to be compiled for the LW dissociation to work.")
#endif
    else if(lwdiss_type .eq. "noraytrace") then
      !Default cap on the local shielding length: the smallest domain extent
      if(lw_lshield_max .le. 0.0) then
        call Grid_getDomainBoundBox(boundBox)
        lw_lshield_max = minval(boundBox(HIGH,1:NDIM) - boundBox(LOW,1:NDIM))
      end if
    else
      call Driver_abortFlash("[rt_lw]: lwdiss_type is unrecognisable")
    end if

    !cA and coA are used to convert mass density to number density (n = rho*X/A) for the
    !C and CO column densities in the ray-tracers, so store them as masses in grams (as
    !rt_ionise does for hA and h2A), not as mass numbers in amu
#if defined C_SPEC || defined CO_SPEC
    call PhysicalConstants_get("proton mass",mH)
#endif
#ifdef C_SPEC
    call Multispecies_getProperty(C_SPEC,A,cA)
    cA = cA * mH
#endif

#ifdef CO_SPEC
    call Multispecies_getProperty(CO_SPEC,A,coA)
    coA = coA * mH
#endif
    first_call = .false.
    !Return to RadTrans.F90 (where the first call will be made from; so the block of code above is basically all for initialising)
    return
  end if

  call Timers_start("LWRates")

  SELECT CASE(lwdiss_type)
    CASE("raytrace")
      !Get the LW flux and shielding factor with the RT
      call RayTrace_LW()

    CASE("hybrid") ! This uses the solution from the moment eqns for the flux, and the shielding factor calculated with the sink RT
      !Get the shielding factors with the RT
      call RayTrace_LW()      
      !Get the total dust-attenuated LW flux which we have from the Petsc solve
      call Compute_LWFL()
    CASE("noraytrace") ! Dust-attenuated flux from the moment eqns, and the shielding factor from local column densities
      !Get the shielding factors from local column densities N = n * rho/|grad rho|
      call Compute_LocalShield()
      !Get the total dust-attenuated LW flux which we have from the Petsc solve
      call Compute_LWFL()
    CASE DEFAULT
      call Driver_abortFlash("[rt_lw]: lwdiss_type is unrecognisable")

  END SELECT

  !Now combine fshield from RT and net (dust-attenuated) flux from the moment solve to get shielded LW flux
  call Compute_Rates()

  call Timers_stop("LWRates")

CONTAINS 

SUBROUTINE Compute_LWFL()
  use rt_lwmodule
  implicit none
  real :: fshield, totalflux

  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)
    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

          !Get the total dust-attenuated LW flux which we have from the Petsc solve
          totalflux = solnData(MOHX_VAR,i,j,k)**2
#if NDIM>1
          totalflux = totalflux + solnData(MOHY_VAR,i,j,k)**2
#if NDIM>2
          totalflux = totalflux + solnData(MOHZ_VAR,i,j,k)**2
#endif
#endif
          totalflux = SQRT(totalflux)

          !Now set the self-shielded flux = totalflux * fshield
          solnData(LWFL_VAR,i,j,k) = totalflux
          
        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID,solnData)
  end do



END SUBROUTINE Compute_LWFL


SUBROUTINE Compute_LocalShield()
  !Local approximation to the line shielding: column densities N_X = n_X * L with the Sobolev-like
  !length L = rho/|grad rho| (Gnedin et al. 2009). The shielding factors are stored in the same containers
  !as the ray-traced ones (FSHM, FSHA, IC, DICO), so Compute_Rates treats both identically.
  use rt_lwmodule
  implicit none
  real :: del(MDIM), grad_rho, lshield, dx_min
  real :: nH2, nH, nC, nCO, fsh_H2, fsh_H, fsh_C, fsh_CO

  !Centred differences need one layer of up-to-date density guard cells
  gcmask(:) = .false.
  gcmask(DENS_VAR) = .true.
  call Grid_fillGuardCells(CENTER,ALLDIR,masksize=NUNK_VARS,mask=gcmask,makeMaskConsistent=.true., &
                           minLayers=1,selectBlockType=LEAF,doLogMask=.not.gcMaskLogged)
  gcMaskLogged = .true.

  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)
    call Grid_getDeltas(blockID,del)
    dx_min = minval(del(1:NDIM))
    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

          grad_rho = ((solnData(DENS_VAR,i+1,j,k) - solnData(DENS_VAR,i-1,j,k))/(2*del(IAXIS)))**2
#if NDIM>1
          grad_rho = grad_rho + ((solnData(DENS_VAR,i,j+1,k) - solnData(DENS_VAR,i,j-1,k))/(2*del(JAXIS)))**2
#if NDIM>2
          grad_rho = grad_rho + ((solnData(DENS_VAR,i,j,k+1) - solnData(DENS_VAR,i,j,k-1))/(2*del(KAXIS)))**2
#endif
#endif
          grad_rho = sqrt(grad_rho)
          lshield = get_shield_length_gradient(solnData(DENS_VAR,i,j,k), grad_rho, dx_min, lw_lshield_max)

          !Number densities of the shielding species
          nH2 = 0.0
          nH = 0.0
          nC = 0.0
          nCO = 0.0
#ifdef H2_SPEC
          nH2 = solnData(DENS_VAR,i,j,k) * solnData(H2_SPEC,i,j,k)/h2A
#endif
#ifdef H_SPEC
          nH = solnData(DENS_VAR,i,j,k) * solnData(H_SPEC,i,j,k)/hA
#endif
#ifdef C_SPEC
          nC = solnData(DENS_VAR,i,j,k) * solnData(C_SPEC,i,j,k)/cA
#endif
#ifdef CO_SPEC
          nCO = solnData(DENS_VAR,i,j,k) * solnData(CO_SPEC,i,j,k)/coA
#endif
          call get_local_fshield(nH2, nH, nC, nCO, lshield, bfive, fsh_H2, fsh_H, fsh_C, fsh_CO)

          solnData(FSHM_VAR,i,j,k) = fsh_H2
          solnData(FSHA_VAR,i,j,k) = fsh_H
          solnData(IC_VAR,i,j,k) = fsh_C
          solnData(DICO_VAR,i,j,k) = fsh_CO

        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID,solnData)
  end do

END SUBROUTINE Compute_LocalShield


SUBROUTINE Compute_Rates()
  use rt_lwmodule
  implicit none
  real :: fshield, totalflux
  !Unattenuated ISRF rates (s^-1) and J_LW_ISRF match the values used for the diffuse background
  !(RayTrace_LW_diffuse) and in KromeChemistry, so sink and background LW rates are normalised consistently
  real, parameter :: H2dissrate_ISRF = 5.6e-11, COdissrate_ISRF = 2.592e-10, Cdissrate_ISRF = 3.1e-10
  real, parameter :: J_LW_ISRF = 3.7e-5 !Mean intensity of ISRF in LW band (Draine; erg s^-1 cm^-2 sr^-1)

  do b = 1, blockCount_
    blockID = blockList_(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)
    do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
      do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)

          !Now set the self-shielded flux = totalflux * fshield
          solnData(DIH2_VAR,i,j,k) = solnData(LWFL_VAR,i,j,k) * solnData(FSHM_VAR,i,j,k)
          if(useHIshield) solnData(DIH2_VAR,i,j,k) = solnData(DIH2_VAR,i,j,k) * solnData(FSHA_VAR,i,j,k)
          solnData(IC_VAR,i,j,k) = solnData(LWFL_VAR,i,j,k) * solnData(IC_VAR,i,j,k)
          solnData(DICO_VAR,i,j,k) = solnData(LWFL_VAR,i,j,k) * solnData(DICO_VAR,i,j,k)

          !Now convert this to dissociation rates by i) converting to mean intensity, ii)scaling with the ISRF value, and iii) multiplying by the optically-thin rate
          solnData(DIH2_VAR,i,j,k) = solnData(DIH2_VAR,i,j,k) / (4*PI * J_LW_ISRF) * H2dissrate_ISRF
          solnData(IC_VAR,i,j,k) = solnData(IC_VAR,i,j,k) / (4*PI * J_LW_ISRF) * Cdissrate_ISRF
          solnData(DICO_VAR,i,j,k) = solnData(DICO_VAR,i,j,k) / (4*PI * J_LW_ISRF) * COdissrate_ISRF
          
        end do
      end do
    end do
    call Grid_releaseBlkPtr(blockID,solnData)
  end do



END SUBROUTINE Compute_Rates

END SUBROUTINE rt_lw
