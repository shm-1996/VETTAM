!!****f* source/physics/sourceTerms/RadTrans/RadTransMain/VETTAM/rt_dustTerms
!!
!! NAME
!!  
!!  rt_dustTerms
!!
!!
!! SYNOPSIS
!! 
!!  rt_dustTerms()
!!  
!! DESCRIPTION
!!
!!	Imparts the dust radiation pressure and any reprocessing of UV to IR wavelengths
!!
!!***
!!***

SUBROUTINE rt_dustTerms(dt)
  use RadTrans_data
  use Grid_interface, ONLY: Grid_getBlkIndexLimits,Grid_getBlkPtr, Grid_releaseBlkPtr
  use Eos_interface, ONLY:Eos_wrapped
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none


#include "Flash.h"
#include "constants.h"

    
  real, intent(in) :: dt
  real             :: rho,temp,vx,vy,vz
  real             :: dE,dMomx,dMomy,dMomz, ekin_old, ekin_new, d_ekin
  real             :: cellsize(MDIM)
  integer          :: b,blockID
  integer          :: i,j,k,l 
  real             :: d1x,d1y,d1z,d2,opac_planck,opac_rosseland
  real             :: prad_ratio
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(:,:,:,:), pointer :: solnData

  do b = 1, blockCount
    blockID = blockList(b)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    call Grid_getBlkPtr(blockID,solnData)
    call Grid_getDeltas(blockID,cellsize)
    do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
      do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

          rho    = solnData(DENS_VAR,i,j,k)
#ifdef TEMP_VAR
          temp   = solnData(TEMP_VAR,i,j,k)
#endif
#ifdef VELX_VAR 
          vx = solnData(VELX_VAR,i,j,k)
#endif          
#ifdef VELY_VAR
          vy = solnData(VELY_VAR,i,j,k)
#endif
#ifdef VELZ_VAR
          vz = solnData(VELZ_VAR,i,j,k)
#endif          
          d1x = dt/(2*cellsize(IAXIS))
          d1y = dt/(2*cellsize(JAXIS))
          d1z = dt/(2*cellsize(KAXIS))
          d2 = dt * rt_speedlt
          !Planck and Rosseland opacity of dust; this is obtained by subtracting out the opacity due to hydrogen
          opac_planck = solnData(TAUP_VAR,i,j,k)
          opac_rosseland = solnData(TAUR_VAR,i,j,k)
          !Subtract hydrogen contribution for EUV band (since this is added by the photoionization module)
#ifdef TAUH_VAR
          if(current_band .eq. 'EUV' .or. current_band .eq. 'EUV_13P6_15P2' .or. &
          & current_band .eq. 'EUV_15P2_INFTY') then
            opac_planck = opac_planck - solnData(TAUH_VAR,i,j,k)
            opac_rosseland = opac_rosseland - solnData(TAUH_VAR,i,j,k)
          endif
#endif

#ifdef REIR_VAR
          !Accumulate the rate per volume (erg cm^-3 s^-1) at which dust absorbs energy in all bands bluer than the IR,
          !assumed to be re-emitted in the IR completely (rt_fillMatrix multiplies by dt; KROME heats dust with it)
          if(current_band .ne. 'IR') solnData(REIR_VAR,i,j,k) = solnData(REIR_VAR,i,j,k) + opac_planck * rt_speedlt * solnData(ERAD_VAR,i,j,k)
#endif
          !PE and LW bands: store the radiation energy density (erg cm^-3) of the band, which KROME converts to a
          !mean intensity J = U c/(4pi) (chi_FUV for PE heating, and the stellar FUV field). The diffuse background
          !ray-trace (RayTrace_LW_diffuse) adds its energy density to the same containers afterwards.
#ifdef UPED_VAR
          if(current_band .eq. 'PE') then
            solnData(UPED_VAR,i,j,k) = solnData(ERAD_VAR,i,j,k)
          endif
#endif

#ifdef ULWD_VAR
          if(current_band .eq. 'LW' .or. current_band .eq. 'LYMAN_WERNER') then
            solnData(ULWD_VAR,i,j,k) = solnData(ERAD_VAR,i,j,k)
          endif
#endif

          !Reset to zero
          dMomx = 0.0
          dMomy = 0.0
          dMomz = 0.0
          ekin_old = 0.0
          ekin_new = 0.0

          ekin_old = 0.5*(solnData(VELX_VAR,i,j,k)**2)
#if NDIM > 1
          ekin_old = ekin_old+ 0.5*(solnData(VELY_VAR,i,j,k)**2)

#if NDIM >2
          ekin_old = ekin_old + 0.5*(solnData(VELZ_VAR,i,j,k)**2)
#endif
#endif
                
          !Computing explicit terms to add to the gas momentum

          !Gas-Radiation momentum exchange term
          dMomx = dMomx + dt * opac_rosseland * solnData(MOHX_VAR,i,j,k)/rt_speedlt
#if NDIM>1
          dMomy = dMomy + dt * opac_rosseland * solnData(MOHY_VAR,i,j,k)/rt_speedlt
#if NDIM>2
          dMomz = dMomz + dt * opac_rosseland * solnData(MOHZ_VAR,i,j,k)/rt_speedlt
#endif
#endif

          !Only add these terms if VET switched on
          if(rt_ovcterms) then
          
            dMomx = dMomx - &
              d2*opac_rosseland *vx*(1+solnData(XXED_VAR,i,j,k))*solnData(ERAD_VAR,i,j,k)/rt_speedlt**2
#if NDIM>1

            dMomx = dMomx - d2*opac_rosseland*vy*solnData(XYED_VAR,i,j,k)*solnData(ERAD_VAR,i,j,k)/rt_speedlt**2 

            dMomy = dMomy - &
              d2*opac_rosseland * (vy*(1+solnData(YYED_VAR,i,j,k)) + &
                vx*solnData(XYED_VAR,i,j,k)) * solnData(ERAD_VAR,i,j,k)/rt_speedlt**2
#if NDIM>2


            dMomx = dMomx - d2*opac_rosseland*vz*solnData(XZED_VAR,i,j,k)*solnData(ERAD_VAR,i,j,k)/rt_speedlt**2
            dMomy = dMomy - d2*opac_rosseland*vz*solnData(YZED_VAR,i,j,k)*solnData(ERAD_VAR,i,j,k)/rt_speedlt**2

            dMomz = dMomz - &
              d2*opac_rosseland * (vz*(1+solnData(ZZED_VAR,i,j,k)) + vx*solnData(XZED_VAR,i,j,k) + &
              vy*solnData(YZED_VAR,i,j,k)) * solnData(ERAD_VAR,i,j,k)/rt_speedlt**2
#endif
#endif
          endif

          !Non-thermal sputtering (Kim et al. 2023, App. D): grains drift at v_d with
          !s = v_d/sqrt(2kT/m_H) ~ sqrt(P_rad/P_gas); above s^2 = rt_sputter_pratio
          !(v_d ~ 150 km/s) they are destroyed, so the dust radiation force is omitted.
          !P_rad = |F|/c of the current band only, a lower bound on the total.
          !Not applied to the IR band: the criterion assumes <Q_pr> ~ 1, which holds for UV
          !but not IR, so IR cannot drive grains to sputtering drift velocities.
#ifdef PRES_VAR
          if(rt_nonthermal_sputter .and. current_band .ne. 'IR') then
            prad_ratio = solnData(MOHX_VAR,i,j,k)**2
#if NDIM>1
            prad_ratio = prad_ratio + solnData(MOHY_VAR,i,j,k)**2
#if NDIM>2
            prad_ratio = prad_ratio + solnData(MOHZ_VAR,i,j,k)**2
#endif
#endif
            prad_ratio = sqrt(prad_ratio)/rt_speedlt/max(solnData(PRES_VAR,i,j,k), tiny(1.0))
            if(prad_ratio .gt. rt_sputter_pratio) then
              dMomx = 0.0
              dMomy = 0.0
              dMomz = 0.0
            endif
          endif
#endif

        ! Update Gas velocities now.
        
        solnData(VELX_VAR,i,j,k) = solnData(VELX_VAR,i,j,k) + dMomx/solnData(DENS_VAR,i,j,k)
#if NDIM>1
        solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) + dMomy/solnData(DENS_VAR,i,j,k)
#if NDIM>2
        solnData(VELZ_VAR,i,j,k) = solnData(VELZ_VAR,i,j,k) + dMomz/solnData(DENS_VAR,i,j,k)  
#endif 
#endif

        
        ekin_new = 0.5*(solnData(VELX_VAR,i,j,k)**2)
#if NDIM > 1
        ekin_new = ekin_new+ 0.5*(solnData(VELY_VAR,i,j,k)**2)

#if NDIM >2
        ekin_new = ekin_new + 0.5*(solnData(VELZ_VAR,i,j,k)**2)
#endif
#endif

        ! compute energy injection
        d_ekin = ekin_new - ekin_old

#ifdef ENER_VAR
        ! update the total energy to be consistent with new velocities
        solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) + d_ekin
#endif
       
        end do
      end do
    end do 
    call Eos_wrapped(MODE_DENS_EI,blkLimits, blockID)
    call Grid_releaseBlkPtr(blockID,solnData)
  end do 


END SUBROUTINE rt_dustTerms