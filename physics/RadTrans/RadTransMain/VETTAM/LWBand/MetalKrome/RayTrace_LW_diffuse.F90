!===============================================================================
!
! Subroutine: RayTrace_LW
! Path: source/physics/RadTrans/RadTransMain/VETTAM/LWBand/MetalKrome/RayTrace_LW.F90
!
!===============================================================================
!
! Description: Ray-Trace step for computing the dissociation of H2, CO and the ionization of C via LW photons. 
! Author: Shyam Menon (2024)
! Email: shyam.menon@rutgers.edu / smenon@flatironsinstitute.org
!


#undef RT_CORRECT
#undef DEBUG_RT


#include "Flash.h"
#include "constants.h"
subroutine RayTrace_LW_diffuse()
  use RadTrans_RayTrace_3DRT
  use rt_lwdata
  use rt_lwmodule
  use Opacity_data, ONLY: dust_lwopac_value, dust_opacity_fuv, dust_peopac_value, dusttoGasRatio
  use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs
  use calc_local_mod, ONLY : calc_local_block_contributions_3DRT
  use create_cut_block_mod, ONLY : create_cut_block_list_3DRT
  use Grid_data, ONLY: gr_oneBlock, gr_nBlockX, gr_nBlockY, gr_nBlockZ
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp
  use rt_data_raytrace_3drt, ONLY : domainSizeX, domainSizeY, domainSizeZ, &
       nrOfAnglesPerGroup, nrOfAngles, dOmega
  use Particles_sinkData, ONLY : particles_global, localnp, localnpf
  use pt_sinkSort,      ONLY: NewQsort_IN
  use pt_sinkInterface, ONLY: pt_sinkGatherGlobal
  use RadTrans_hybridCharModule
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use RadTrans_data, ONLY: rt_useRadTrans
#ifdef CO_SPEC
  use KromeChemistry_data, ONLY: chi_LW_bg, chi_PE_bg
#endif

  implicit none

  !The number of optical depth/column density calculations done in this ray-trace
  !5 species: Dust taud, and H2, CO, C, H column densities for the shielding factors
  integer, parameter  :: NLW_SPEC = 5
  integer             :: iAngle,iAb,iAe,iAindex,iAnglemax
  real                :: dirFact(NDIM)
  real                :: phi,theta
  integer             :: iAngleGroup
  real                :: lambda_akku
  ! The regular Factors are used to convert coordinates from
  ! physical values [cm] to cell units on the maximum level of
  ! refinement. The unit of regFact is then [cells cm^-1].
  ! Regular coordinates mean, that they are given in units
  ! of cell numbers on the maximum refinement level.
  !
  real,    save       :: regFact(NDIM)
  real,    save       :: PointSource(NDIM)

  ! The domain size in units of cells at the maximum
  ! refinement level.
  !
  integer, save       :: regTotI,regTotJ,regTotK

  !
  ! block indices
  !
  !  b: the index of blk in leafBlockList
  !  blk: the index of a block in the flash database
  !
  integer             :: b, blk

  !
  ! Some general flash loop iteration counter and parameters.
  !
  integer             :: idx(NDIM), i, j, k, n, p, q
  real, dimension(MAXCELLS,NDIM) :: xr
  real :: ColDens, ColDens_local, ColDens_global, taud
  real :: intensity, intensity_local, intensity_global
  real :: dest(NDIM)
  real, parameter :: J_FUV_ISRF = 2.1e-4, J_PE_ISRF = 1.8e-4, J_LW_ISRF = 3.7e-5 ! in ergs cm^-2 s^-1 sr^-1
  real :: c_light, mH
  !Containers to store the column density of H2 and the dust optical depth, needed for other species
  !Containers would store these for every cell in the block, for every particle -- i.e. shape (nAllParticles, NXB, NYB, NZB)
  real, dimension(:,:,:,:,:), allocatable :: NH2_star, taud_star
  !Flag for which species the shielding is computed; 0: H2, 1: H
  character(len=MAX_STRING_LENGTH) :: fsh_spec
  logical, parameter :: debug = .true.



  !-------------------------------------------------------------------------------
  !  
  !
  !======================================================================================
  !

  !Return if no angles defined
  if(nrOfAngles .le. 0) return

  if(rt_nrOfAngleGroups .gt. 1) call Driver_abortFlash("RayTrace_LW_diffuse currently only supports a single angle group.")

  call Timers_start("rad_raytrace_bg")

  if(dr_globalMe.eq.MASTER_PE) call Logfile_stamp("Starting Raytracing for diffuse Lyman Werner contribution", "[RT_SINK] ")

  call PhysicalConstants_get("speed of light",c_light)
  call PhysicalConstants_get("proton mass",mH)


  ! Get the current resolution on the highest refinement level
  ! and the conversion factors
  regTotI = gr_nBlockX * 2**(maxLevel-1) * NXB
  regTotJ = gr_nBlockY * 2**(maxLevel-1) * NYB
  regTotK = gr_nBlockZ * 2**(maxLevel-1) * NZB
  regFact(1) = real(regTotI) / domainSizeX
  regFact(2) = real(regTotJ) / domainSizeY
  regFact(3) = real(regTotK) / domainSizeZ

  ! allocate containers for NH2 and taud (stored as C/CO shielding depends on these, and they aren't available simultaneously)
  allocate(NH2_star(nrOfLeafBlocks,NXB,NYB,NZB,nrOfAngles), taud_star(nrOfLeafBlocks,NXB,NYB,NZB,nrOfAngles))
  
  fsh_spec = "dust" !Indicate that dust optical depth to be calculated
  call raytrace_species()
  

  fsh_spec = "H2" !Indicate that shielding function for H2 is to be done
  call raytrace_species()  

  !Now do the same steps for H
  fsh_spec = "H" !Indicate that shielding function for H is to be done
  call raytrace_species()

#ifdef C_SPEC
  !Now do the same steps for C
  fsh_spec = "C" !Indicate that shielding function for C is to be done
  call raytrace_species()
#endif

#ifdef CO_SPEC
  !Now do the same steps for CO
  fsh_spec = "CO" !Indicate that shielding function for CO is to be done
  call raytrace_species()
#endif

  deallocate(NH2_star, taud_star)

  if(dr_globalMe.eq.MASTER_PE) call Logfile_stamp("Raytracing done for diffuse Lyman Werner contribution", "[RT_SINK] ")

  call Timers_stop("rad_raytrace_bg")

  CONTAINS

  subroutine raytrace_species()
    !This subroutine is used to perform the ray-tracing for each species column density from the boundary
    !Steps: i) Set OPAC_VAR to the appropriate quantity for the species (noDens for column calculations, or kappa for taud)
    !       iii) Compute face contribution, communicate face values, and combine the SC solve using face values with LC solve to get the total column density

    !Reset the containers if not using RadTrans as it would reset it otherwise
    if(.not. rt_useRadTrans) then
      call rad_resetLWContainers()
    end if

    !Set OPAC_VAR
    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call SetOPAC_LW(blk)
    end do

    ! Initialize Angle Indices
    iAb=0
    iAe=0

    

    do iAngleGroup=1,rt_nrOfAngleGroups
      iAb = iAe+1 ! Start where we stopped 
      iAe = iAngleGroup*nrOfAnglesPerGroup

      call Timers_start("rad_calc_face_values")
      call rad_facevalues()
      call Timers_stop("rad_calc_face_values")      

      call Timers_start("rad_MPI_ALLGATHER_faceValues")
      call gather_faces(iAb,iAe)
      call Timers_stop("rad_MPI_ALLGATHER_faceValues")


      call Timers_start("rad_calcContrib")
      call rad_calcContrib()
      call Timers_stop("rad_calcContrib")
    end do


  end subroutine raytrace_species

  subroutine rad_facevalues()
    real :: dColDens, dtaucell, dlcell

    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call Grid_getBlkPtr(blk, solnData)
      do n=IAXIS,KAXIS
        call Grid_getCellCoords(n, blk, RIGHT_EDGE, .TRUE., xr(:,n), MAXCELLS)
      end do

      do iAngle=iAb,iAe
        iAindex = (iAngle-iAb)+1
        
        !Get angles to integrate rays over
        call getAngles(iAngle,rand_angles(1),rand_angles(2),rand_angles(3),theta,phi,dirFact)

        ! Compute optical depths to faces
        ! The optical depth should be zero, if either:
        ! - the destination is equal the source
        ! - the face can not be reached in that direction
        ! (in all blocks but the on including the point source,
        !  only three faces have valid ColDens values.)
        do n=1,NDIM
          !Direction of face required. For x -y,z, for y: z,x, for z: x,y
          p = MOD(n,3) + 1
          q = MOD(n+1,3) + 1

          IF(dirFact(n).lt.0.) THEN
           i = ib(n)
          ELSE
           i = ie(n)
          END IF
          do j = ib(q),ie(q)
            do k = ib(p),ie(p)
              idx(n) = i
              idx(q) = j
              idx(p) = k
              ! Here is space for optimization
              dest = (/ xr(idx(1),1), xr(idx(2),2), xr(idx(3),3) /)

              call calc_local_block_contributions_3DRT(    &
                  regFact, blk, &
                  PointSource, dest, dirFact,   &
                  maxLevel,         &
                  ColDens, intensity, lambda_akku, &
                  n, &
                  .true.,dtaucell,dlcell)
              faceValueAll(k,j,n,1,b,iAindex,dr_globalMe) = ColDens
              faceValueAll(k,j,n,2,b,iAindex,dr_globalMe) = intensity
            end do ! k
          end do !j

        end do !n
      end do !iAngle
      call Grid_releaseBlkPtr(blk, solnData)
    enddo ! b

  END SUBROUTINE rad_facevalues

  SUBROUTINE rad_calcContrib
    real :: luminosity, localColDens, rho, del(NDIM), radius_sink, diss_rate, fshield
    real :: dColDens, dtaucell, dlcell, NH2
    logical :: printval
    REAL,DIMENSION(NDIM) :: target_point
    real, parameter :: dissh2_draine = 5.6e-11, ionC_draine = 3.1e-10,dissCO_draine = 2.592e-10, cr_sol = 2e-16, &
                      cr_shield_NH = 9.35e20 !Rates in s^-1; Draine rates scaled to chi=1; H2 dissoc rate from Draine & Bertoldi 1996; C ion rate from van Dishoeck & Black 1988; CO dissoc rate from Visser et al 2009

    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call Grid_getBlkPtr(blk, solnData)
      call Grid_getDeltas(blk, del) 
  
      do iAngle=iAb,iAe
        iAindex = (iAngle-iAb)+1
        call getAngles(iAngle,rand_angles(1),rand_angles(2),rand_angles(3),&
          theta,phi,dirFact)      

        do k=ib(3)+K3D,ie(3)
          do j=ib(2)+K2D,ie(2)
            do i=ib(1)+1,ie(1)
                dest(1) = gr_oneBlock(blk)%firstAxisCoords(CENTER, i)
                dest(2) = gr_oneBlock(blk)%secondAxisCoords(CENTER, j)
                dest(3) = gr_oneBlock(blk)%thirdAxisCoords(CENTER, k)

                call calc_local_block_contributions_3DRT(    &
                  regFact, blk, &
                  PointSource, dest, dirFact, &
                  maxLevel,                               &
                  ColDens_local, intensity_local, lambda_akku, &
                  0, .true., dtaucell, dlcell)

                call create_cut_block_list_3DRT(&
                        dr_globalMe, blk,ColDens_global, intensity_global, PointSource, dest, &
                        regFact, dirFact, iAindex, maxNrOfLeafBlocks, dr_globalNumProcs, maxLevel, .false.)

                ! Now add up global and local column densities
                ColDens = ColDens_global + ColDens_local

                if(fsh_spec == "dust") then
                  !For dust, no shielding factor calculation. Just store the optical depth for later use. 
                  taud_star(b,i-ib(1),j-ib(2),k-ib(3),iAindex) = ColDens
                  !return because no shielding calculations required here
                else
                  !Use the stored value to account for attenuation due to dust
                  taud = taud_star(b,i-ib(1),j-ib(2),k-ib(3),iAindex)
                endif

                if(fsh_spec == "H2") then
                  !For dust, no shielding factor calculation. Just store the optical depth for later use. 
                  NH2_star(b,i-ib(1),j-ib(2),k-ib(3),iAindex) = ColDens
                else if(fsh_spec == "C" .or. fsh_spec == "CO" .or. fsh_spec == 'H') then
                  !Use the stored value to account for attenuation due to dust
                  NH2 = NH2_star(b,i-ib(1),j-ib(2),k-ib(3),iAindex)
                endif

                
                if(fsh_spec == "dust") then
                  solnData(ULWD_VAR,i,j,k) = solnData(ULWD_VAR,i,j,k) + (chi_LW_bg * exp(-ColDens) * J_LW_ISRF)*dOmega/c_light
                  !PE heating rate obtained by adjusting for the different opacities in PE and LW bands (assuming constant opacity)
                  if(dust_opacity_fuv .ne. 'Fixed') then !Safety check
                    call Driver_abortFlash("Diffuse LW ray-tracing assumes fixed dust opacities for the PE band background field.")
                  end if
                  solnData(UPED_VAR,i,j,k) = solnData(UPED_VAR,i,j,k) + (chi_PE_bg * exp(-ColDens*dust_peopac_value/dust_lwopac_value) &
                                              * J_LW_ISRF)*dOmega/c_light

#ifdef REIR_VAR
                  !Add the dust heating from the background LW field
                  solnData(REIR_VAR,i,j,k) = solnData(REIR_VAR,i,j,k) + (chi_LW_bg * exp(-ColDens) * J_LW_ISRF)*dOmega* &
                                                    dust_lwopac_value*mH*dusttoGasRatio !in erg s^-1
                  !Add the dust heating from the background PE field
                  solnData(REIR_VAR,i,j,k) = solnData(REIR_VAR,i,j,k) + (chi_PE_bg * exp(-ColDens*dust_peopac_value/ &
                                                dust_lwopac_value) * J_LW_ISRF)*dOmega* dust_peopac_value*mH*dusttoGasRatio !in erg s^-1
                                                    
#endif
                endif
                !Now compute the contribution to dissociation/ionization rate from the background
                !Note: the H2 self-shielding is included in the step for H below since both H and H2 can decrease DIH2
                if(fsh_spec == "H") then
                  diss_rate = chi_LW_bg * dissh2_draine * exp(-taud) !Dust-attenuated dissociation rate
                  fshield = get_fshield_H2(NH2,bfive)
                  if(useHIShield) fshield = fshield * get_fshield_H(ColDens,bfive)
                  !Use the fshield container to cumulatively add the total unattenuated flux
                  solnData(DIH2_VAR,i,j,k) = solnData(DIH2_VAR,i,j,k) + diss_rate * fshield
                else if(fsh_spec == "C") then
                  diss_rate = chi_LW_bg * ionC_draine * exp(-taud) !Dust-attenuated dissociation rate
                  !Obtain the shielding due to H2 and C for C
                  fshield = get_fshield_C(NH2, ColDens, bfive)
                  solnData(IC_VAR,i,j,k) = solnData(IC_VAR,i,j,k) + diss_rate * fshield

                else if(fsh_spec == "CO") then
                  diss_rate = chi_LW_bg * dissCO_draine * exp(-taud) !Dust-attenuated dissociation rate
                  !Obtain the shielding due to H2 and CO for CO
                  fshield = get_fshield_CO(NH2, ColDens, bfive)
                  solnData(DICO_VAR,i,j,k) = solnData(DICO_VAR,i,j,k) + diss_rate * fshield
                endif
            enddo ! i
          enddo !j
        enddo !k
      enddo !iAngle
      call Grid_releaseBlkPtr(blk, solnData)
    enddo ! blocks

  END SUBROUTINE rad_calcContrib

  SUBROUTINE gather_faces(iAb,iAe)
    implicit none
    integer, intent(in) :: iAb,iAe
    integer :: angles
    integer, dimension(dr_globalNumProcs) :: displs,recvcounts
    !
    !  Gather all face values from all processors and store
    !  the results in faceValueAll.
    !
#ifdef DEBUG_RT
    write(*,*) 'gathering face-values...'
#endif
    angles = iAe-iAb+1
    IF(angles.eq.nrOfAnglesPerGroup) THEN
      call MPI_ALLGATHER( &
           MPI_IN_PLACE, &
           0, &
           MPI_DATATYPE_NULL, &
           faceValueAll,&
           (nb(1)+1)*(nb(2)+1)*3*2*maxNrOfLeafBlocks*nrOfAnglesPerGroup, &
           FLASH_REAL,&
           MPI_COMM_WORLD,ierr)
    ELSE IF(iAe-iAb+1.ge.1) THEN
      do i=1,dr_globalNumProcs
        displs(i) = (nb(1)+1)*(nb(2)+1)*3*2*maxNrOfLeafBlocks*nrOfAnglesPerGroup*(i-1)
      end do
      recvcounts(:) = (nb(1)+1)*(nb(2)+1)*3*2*maxNrOfLeafBlocks*angles

      call MPI_ALLGATHERV( &
           MPI_IN_PLACE, &
           0, &
           MPI_DATATYPE_NULL, &
           faceValueAll, &
           recvcounts, &
           displs, &
           FLASH_REAL, &
           MPI_COMM_WORLD,ierr)
    END IF
#ifdef DEBUG_RT
    write(*,*) 'gathering face-values... done'
#endif
    !
  END SUBROUTINE gather_faces

  SUBROUTINE SetOPAC_LW(blk)
    use rt_ionisedata, ONLY: h2A, hA
    IMPLICIT NONE
    integer, intent(in) :: blk
    real :: nH2, nH, nC, nCO
    integer :: i,j,k
    real, dimension(NDIM) :: position, sink_pos

    if(fsh_spec == "dust" .and. dust_opacity_fuv .ne. 'Fixed') then
      call Driver_abortFlash("Currently, only fixed dust opacity models are supported for LW ray-tracing.")
    endif

    call Grid_getBlkPtr(blk, solnData)
    do k=GRID_KLO_GC,GRID_KHI_GC
      do j=GRID_JLO_GC,GRID_JHI_GC
        do i=GRID_ILO_GC,GRID_IHI_GC
            
            if(fsh_spec == "H2") then
              !Obtain the number density of nH2
#ifdef H2_SPEC
              nH2 = solnData(DENS_VAR,i,j,k) * solnData(H2_SPEC,i,j,k)/h2A
#else
              call Driver_abortFlash("The Multispecies container H2_SPEC has to be present for LW dissociation.")
#endif
              !Now set the container OPAC_VAR which is the quantity integrated in the ray-tracer to the number density;
              !this ensures the column density is returned
              solnData(OPAC_VAR,i,j,k) = nH2
            else if(fsh_spec == "H") then
              !Obtain the number density of nH
#ifdef H_SPEC
              nH = solnData(DENS_VAR,i,j,k) * solnData(H_SPEC,i,j,k)/hA
#else
              call Driver_abortFlash("The Multispecies container H_SPEC has to be present for LW dissociation.")
#endif
              !Now set the container OPAC_VAR which is the quantity integrated in the ray-tracer to the number density;
              !this ensures the column density is returned
              solnData(OPAC_VAR,i,j,k) = nH
            else if(fsh_spec == "C") then
            !Obtain the number density of C
#ifdef C_SPEC
              nC = solnData(DENS_VAR,i,j,k) * solnData(C_SPEC,i,j,k)/cA
            
#else
              call Driver_abortFlash("The Multispecies container C_SPEC has to be present for C ionization.")
#endif
              solnData(OPAC_VAR,i,j,k) = nC
            else if(fsh_spec == "CO") then
            !Obtain the number density of CO
#ifdef CO_SPEC
              nCO = solnData(DENS_VAR,i,j,k) * solnData(CO_SPEC,i,j,k)/coA
#else
              call Driver_abortFlash("The Multispecies container CO_SPEC has to be present for CO dissociation.")
#endif
              solnData(OPAC_VAR,i,j,k) = nCO

            else if(fsh_spec == "dust") then
              !Can't use TAUP_VAR as no guarantee that VETTAM is being used and/or the last band was LW
              !Instead, use the constant dust opacity per unit length to recompute the dust opacity here
              solnData(OPAC_VAR,i,j,k) = dust_lwopac_value * dusttoGasRatio * solnData(DENS_VAR,i,j,k)
            endif

        end do 
      end do
    end do
    call Grid_releaseBlkPtr(blk, solnData)

  END SUBROUTINE SetOPAC_LW

  SUBROUTINE rad_resetLWContainers
    implicit none
    INTEGER :: b, blk

    call Timers_start("rad_resetData")

    do b=1,nrOfLeafBlocks
      !
      blk = leafList(b)

      call Grid_getBlkPtr(blk, solnData)
      !Molecular
      if(fsh_spec == "H2") solnData(DIH2_VAR,:,:,:) = 0.0
      if(fsh_spec == "C") solnData(IC_VAR,:,:,:) = 0.0
      if(fsh_spec == "CO") solnData(DICO_VAR,:,:,:) = 0.0
      if(fsh_spec == "dust") solnData(ULWD_VAR,:,:,:) = 0.0
      if(fsh_spec == "dust") solnData(UPED_VAR,:,:,:) = 0.0
      solnData(LWFL_VAR,:,:,:) = 0.0
      call Grid_releaseBlkPtr(blk, solnData)

    enddo ! b
    !
    call Timers_stop("rad_resetData")
    !
  end subroutine rad_resetLWContainers

END SUBROUTINE RayTrace_LW_diffuse


