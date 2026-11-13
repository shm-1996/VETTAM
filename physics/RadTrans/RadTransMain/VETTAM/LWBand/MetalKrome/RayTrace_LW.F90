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
subroutine RayTrace_LW()
  use RadTrans_RayTrace_3DRT
  use rt_lwdata
  use rt_lwmodule
  use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs
  use calc_local_sink_mod, ONLY: calc_local_block_contributions_sink
  use create_cut_block_mod, ONLY : create_cut_block_list_3DRT
  use Grid_data, ONLY: gr_oneBlock, gr_nBlockX, gr_nBlockY, gr_nBlockZ
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Logfile_interface, ONLY : Logfile_stamp
  use rt_data_raytrace_3drt, ONLY : domainSizeX, domainSizeY, domainSizeZ, &
       nrOfAnglesPerGroup, nrOfAngles
  use Particles_sinkData, ONLY : particles_global, localnp, localnpf
  use pt_sinkSort,      ONLY: NewQsort_IN
  use pt_sinkInterface, ONLY: pt_sinkGatherGlobal
  use RadTrans_hybridCharModule
  use Grid_interface, ONLY : Grid_getMinCellSize
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use RadTrans_data, ONLY : rt_useRadTrans, rt_speedlt
  use rt_faceComm, ONLY : rt_faceLayout, rt_faceExchange, RT_GEOM_POINT, RT_GEOM_DIR
  use rt_checksum, ONLY : rt_checksumLog, rt_checksumAdd

  implicit none

  !The number of optical depth/column density calculations done in this ray-trace
  !5 species: Dust taud, and H2, CO, C, H column densities for the shielding factors
  integer, parameter  :: NLW_SPEC = 5
  integer             :: iAngle,iAb,iAe,iAindex,iAnglemax
  real                :: dirFact(NDIM)
  integer             :: iAngleGroup
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
  real                :: distSqr_sink
  !Near-sink floor on r^2: the Gaussian width used for JSTR and UVF0 in rt_sinkInject
  real, save          :: sigma_star
  real                :: sigma_sink
  logical, save       :: firstCall = .true.
  !lwdiss_type = "raytrace" skips the LW moment solve, so this ray-trace also supplies what it would have given the
  !rest of VETTAM/KROME: E = sum_s F_s/c in ERAD and F = sum_s F_s rhat_s in MOH (for rt_dustTerms: dust heating
  !via REIR, ULWD, radiation pressure), and the unattenuated stellar LW flux in UVF0 (for KROME's chi0 and Av).
  !Only when the moment solver runs, which resets UVF0/ULWD/REIR each update.
  logical             :: fillMoments
  real, dimension(MAXCELLS,NDIM) :: xr
  real :: ColDens, ColDens_local, ColDens_global, taud
  real :: intensity_global
  real :: dest(NDIM)
  ! sink positions / ray directions of the current angle group (face-value need sets)
  real :: groupGeom(MDIM, nrOfAnglesPerGroup)

  integer, dimension(:), allocatable :: id_sorted, QSindex
  !Containers to store the column density of H2 and the dust optical depth, needed for other species
  !Containers would store these for every cell in the block, for every particle -- i.e. shape (nAllParticles, NXB, NYB, NZB)
  real, dimension(:,:,:,:,:), allocatable :: NH2_star, taud_star
  ! Particle stuff
  integer              :: nAllParticles
  !Flag for which species the shielding is computed; 0: H2, 1: H
  character(len=MAX_STRING_LENGTH) :: fsh_spec



  !-------------------------------------------------------------------------------
  !  
  !
  !======================================================================================
  !

  ! gather all particle positions and luminosities
  !
  !Sinks without an LW/FUV luminosity property (no stellar evolution model compiled in) emit no LW radiation,
  !so they are treated like having no sources
#if defined(SINK_PART_TYPE) && (defined(LUMLW_PART_PROP) || defined(LUMFUV_PART_PROP) || defined(LUMINOSITY_PART_PROP))
  call pt_sinkGatherGlobal()
  nAllParticles = localnpf
#else
  nAllParticles = 0
#endif

  !Return if no sink particles present
  if(nAllParticles .eq. 0) return

  call Timers_start("rad_raytrace_sink")

  if(dr_globalMe.eq.MASTER_PE) call Logfile_stamp("Starting Raytracing for sink Lyman Werner contribution", "[RT_SINK] ")

  if(firstCall) then
    call RuntimeParameters_get("sigma_star", sigma_star)
    firstCall = .false.
  endif
  call Grid_getMinCellSize(sigma_sink)
  sigma_sink = sigma_star*sigma_sink

  fillMoments = (lwdiss_type == "raytrace") .and. rt_useRadTrans
  if(fillMoments) then
    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call Grid_getBlkPtr(blk, solnData)
      solnData(ERAD_VAR,:,:,:) = 0.0
      solnData(MOHX_VAR,:,:,:) = 0.0
#if NDIM>1
      solnData(MOHY_VAR,:,:,:) = 0.0
#if NDIM>2
      solnData(MOHZ_VAR,:,:,:) = 0.0
#endif
#endif
      call Grid_releaseBlkPtr(blk, solnData)
    end do
  endif

  ! Initialize Angle Indices
  iAb=0
  iAe=0

  ! Get the current resolution on the highest refinement level
  ! and the conversion factors
  regTotI = gr_nBlockX * 2**(maxLevel-1) * NXB
  regTotJ = gr_nBlockY * 2**(maxLevel-1) * NYB
  regTotK = gr_nBlockZ * 2**(maxLevel-1) * NZB
  regFact(1) = real(regTotI) / domainSizeX
  regFact(2) = real(regTotJ) / domainSizeY
  regFact(3) = real(regTotK) / domainSizeZ

  ! particles_global is not in the same order on all procs. Order it!
  allocate(id_sorted(localnpf),QSindex(localnpf))
  do iAngle=1,localnpf
    id_sorted(iAngle) = int(particles_global(TAG_PART_PROP,iAngle))
  end do
  call NewQsort_IN(id_sorted, QSindex)

  ! allocate containers for NH2 and taud
  allocate(NH2_star(nrOfLeafBlocks,NXB,NYB,NZB,nAllParticles), taud_star(nrOfLeafBlocks,NXB,NYB,NZB,nAllParticles))

  fsh_spec = "dust" !Indicate that dust optical depth to be calculated
  call raytrace_species()

  fsh_spec = "H2" !Indicate that shielding function for H2 is to be done
  call raytrace_species()  

  !Compute the LW-flux weighted average shielding factor for H2
  do b=1, nrOfLeafBlocks
    blk=leafList(b)
    call average_fshield(blk)
  end do

  if(useHIshield) then
    !Now do the same steps for H
    fsh_spec = "H" !Indicate that shielding function for H is to be done
    call raytrace_species()
    !Compute the LW-flux weighted average shielding factor for H
    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call average_fshield(blk)
    end do
  endif

#ifdef C_SPEC
  !Now do the same steps for C
  fsh_spec = "C" !Indicate that shielding function for C is to be done
  call raytrace_species()
  !Compute the LW-flux weighted average shielding factor for C
  do b=1, nrOfLeafBlocks
    blk=leafList(b)
    call average_fshield(blk)
  end do
#endif

#ifdef CO_SPEC
  !Now do the same steps for CO
  fsh_spec = "CO" !Indicate that shielding function for CO is to be done
  call raytrace_species()
  !Compute the LW-flux weighted average shielding factor for CO
  do b=1, nrOfLeafBlocks
    blk=leafList(b)
    call average_fshield(blk)
  end do
#endif

  !Compute the LW dissociation rate, C ionization rate and CO dissociation
  ! if(lwdiss_type .eq. "raytrace") then
  !   fsh_spec = "flux" !This is done to ensure only LWFL_VAR is reset in next step
  !   call rad_resetLWContainers()
  !   !Now combine above shielding factors to get the LW flux
  !   call calculate_LWFlux()
  ! endif

  !E = sum_s F_s/c for the streaming stellar field (LWFL now holds the dust-attenuated flux summed over sinks)
  if(fillMoments) then
    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call Grid_getBlkPtr(blk, solnData)
      solnData(ERAD_VAR,:,:,:) = solnData(LWFL_VAR,:,:,:)/rt_speedlt
      call Grid_releaseBlkPtr(blk, solnData)
    end do
  endif

  deallocate(id_sorted,QSindex)
  deallocate(NH2_star, taud_star)
    

  if(dr_globalMe.eq.MASTER_PE) call Logfile_stamp("Raytracing for Lyman Werner Radiation completed", "[RT_SINK] ")

  if(rt_debugChecksum) call checksum_output()

  call Timers_stop("rad_raytrace_sink")

  CONTAINS

  ! rt_debugChecksum: exact checksum of every field written by this tracer
  subroutine checksum_output()
    integer :: vars(10), ncs
    character(len=4) :: names(10)
    character(len=32) :: ident
    ncs = 0
    call rt_checksumAdd(LWFL_VAR,'LWFL',ncs,vars,names)
    call rt_checksumAdd(FSHM_VAR,'FSHM',ncs,vars,names)
    if(useHIshield) call rt_checksumAdd(FSHA_VAR,'FSHA',ncs,vars,names)
#ifdef C_SPEC
    call rt_checksumAdd(IC_VAR,'IC  ',ncs,vars,names)
#endif
#ifdef CO_SPEC
    call rt_checksumAdd(DICO_VAR,'DICO',ncs,vars,names)
#endif
    ! lwdiss_type = "raytrace": this tracer also fills the LW moments
    if(fillMoments) then
      call rt_checksumAdd(ERAD_VAR,'ERAD',ncs,vars,names)
      call rt_checksumAdd(MOHX_VAR,'MOHX',ncs,vars,names)
#if NDIM>1
      call rt_checksumAdd(MOHY_VAR,'MOHY',ncs,vars,names)
#endif
#if NDIM>2
      call rt_checksumAdd(MOHZ_VAR,'MOHZ',ncs,vars,names)
#endif
#ifdef UVF0_VAR
      call rt_checksumAdd(UVF0_VAR,'UVF0',ncs,vars,names)
#endif
    endif
    write(ident,'(a,i0,a,i0)') 'nsinks=', nAllParticles, ' groups=', &
         CEILING(nAllParticles/REAL(nrOfAnglesPerGroup))
    call rt_checksumLog('LW_sink', ident, ncs, vars(1:ncs), names(1:ncs))
  end subroutine checksum_output

  subroutine raytrace_species()
    !This subroutine is used to perform the ray-tracing for each species column density
    !Steps: i) reset containers where flux-weighted sums are stored, ii) Set OPAC_VAR to the appropriate quantity for the species (noDens for column calculations, or kappa for taud)
    !       iii) Compute face contribution, communicate face values, and combine the SC solve using face values with LC solve to get the total column density

    !Reset the containers
    call rad_resetLWContainers()

    !Set OPAC_VAR
    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call SetOPAC_LW(blk)
    end do

    !Make sure the communication is limited to the memory-related constraint
    iAnglemax = CEILING(nAllParticles/REAL(nrOfAnglesPerGroup))
    do iAngleGroup=1,iAnglemax
      iAb = (iAngleGroup-1)*nrOfAnglesPerGroup + 1
      iAe = MIN(iAngleGroup*nrOfAnglesPerGroup,nAllParticles)

      do iAngle=iAb,iAe
        groupGeom(:,iAngle-iAb+1) = particles_global(POSX_PART_PROP:POSZ_PART_PROP,QSindex(iAngle))
      end do
      call rt_faceLayout(iAe-iAb+1, RT_GEOM_POINT, groupGeom)

      call Timers_start("rad_calc_face_values")
      call rad_facevalues()
      call Timers_stop("rad_calc_face_values")

      call rt_faceExchange(iAe-iAb+1)

      call Timers_start("rad_calcContrib")
      call rad_calcContrib()
      call Timers_stop("rad_calcContrib")
    end do


  end subroutine raytrace_species

  subroutine rad_facevalues()
    real :: dColDens, dlcell

    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call Grid_getBlkPtr(blk, solnData)
      do n=IAXIS,KAXIS
        call Grid_getCellCoords(n, blk, RIGHT_EDGE, .TRUE., xr(:,n), MAXCELLS)
      end do

      do iAngle=iAb,iAe
        iAindex = (iAngle-iAb)+1
        
        !Position of point source
        PointSource = particles_global(POSX_PART_PROP:POSZ_PART_PROP,QSindex(iAngle))

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

          i=ib(n)
          do j = ib(q),ie(q)
            do k = ib(p),ie(p)
              idx(n) = i
              idx(q) = j
              idx(p) = k
              ! Here is space for optimization
              dest = (/ xr(idx(1),1), xr(idx(2),2), xr(idx(3),3) /)
              dirFact = dest-PointSource
              
              if(dirFact(n) .ge. 0) then
                ColDens = 0.0
              else
                call calc_local_block_contributions_sink(&
                        regFact, blk, PointSource, dest, dirFact, maxLevel, ColDens, &
                        n,dColDens,dlcell,.false.)
              end if
              faceValueAll(k,j,n,1,faceSlot(b,iAindex,dr_globalMe)) = ColDens
            end do ! k
          end do !j


          i=ie(n)
          do j = ib(q),ie(q)
            do k = ib(p),ie(p)
              idx(n) = i
              idx(q) = j
              idx(p) = k
              ! Here is space for optimization
              dest = (/ xr(idx(1),1), xr(idx(2),2), xr(idx(3),3) /)
              dirFact = dest-PointSource
              if(dirFact(n) .lt. 0) then
                ColDens = 0.0
              else
                call calc_local_block_contributions_sink(&
                        regFact, blk, PointSource, dest, dirFact, maxLevel, ColDens, &
                        n,dColDens,dlcell,.false.)
              end if
              faceValueAll(k,j,n,2,faceSlot(b,iAindex,dr_globalMe)) = ColDens
            end do ! k
          end do !j

        end do !n
      end do !iAngle
    enddo ! b

  END SUBROUTINE rad_facevalues

  SUBROUTINE calculate_LWFlux
    real :: luminosity, del(NDIM), radius_sink, flux, fshield_eff

    !Loop over all blocks and sinks
    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call Grid_getBlkPtr(blk, solnData)
      call Grid_getDeltas(blk, del)
      do n=IAXIS,KAXIS
        call Grid_getCellCoords(n, blk, RIGHT_EDGE, .TRUE., xr(:,n), MAXCELLS)
      end do      
      do iAngle=iAb,iAe
        iAindex = (iAngle-iAb)+1
        
        !Position of point source
        PointSource = particles_global(POSX_PART_PROP:POSZ_PART_PROP,QSindex(iAngle))

        do k=ib(3)+K3D,ie(3)
          do j=ib(2)+K2D,ie(2)
            do i=ib(1)+1,ie(1)
                dest(1) = gr_oneBlock(blk)%firstAxisCoords(CENTER, i)
                dest(2) = gr_oneBlock(blk)%secondAxisCoords(CENTER, j)
                dest(3) = gr_oneBlock(blk)%thirdAxisCoords(CENTER, k)

                dirFact = dest-PointSource
                dirFact = dirFact/SQRT(SUM(dirFact**2))

                distSqr_sink = MAX(SUM((dest-PointSource)**2),(0.5/regFact(1))**2)

                !Prevent any radiation within the radius of the star
#ifdef STELLAR_RADIUS_PART_PROP
                radius_sink = particles_global(STELLAR_RADIUS_PART_PROP,QSindex(iAngle))
#else
                radius_sink = 0.0
#endif
                if (distSqr_sink .ge. radius_sink**2) then
#ifdef LUMLW_PART_PROP
                    luminosity = particles_global(LUMLW_PART_PROP,QSindex(iAngle))
#elif LUMFUV_PART_PROP
                    luminosity = particles_global(LUMFUV_PART_PROP,QSindex(iAngle))
#elif defined(LUMINOSITY_PART_PROP)
                    luminosity = particles_global(LUMINOSITY_PART_PROP,QSindex(iAngle))
#else
                    luminosity = 0.0
#endif
                
                    !Unattenuated radiation flux
                    flux = luminosity / (4.0*PI*distSqr_sink)

                else 
                    flux = 0.0

                end if

                fshield_eff = solnData(FSHM_VAR,i,j,k)
                !Include shielding due to atomic H if useHIShield included
                if(useHIshield) fshield_eff = fshield_eff * solnData(FSHA_VAR,i,j,k) 

                solnData(LWFL_VAR,i,j,k) = solnData(LWFL_VAR,i,j,k) + flux * fshield_eff
              enddo ! i
          enddo !j
        enddo !k
      enddo !iAngle
      call Grid_releaseBlkPtr(blk, solnData)
    enddo ! blocks
  END SUBROUTINE calculate_LWFlux



  SUBROUTINE rad_calcContrib
    !The approach is to accumulate the total (dust-attenuated) flux in variable LWFL_VAR
    !The (dust-attenuated) flux-weighted sum over all sources for the shielding factors are stored in FSHM, FSHA, IC, DICO
    !The combination of the two is used to compute the effective shielfing factor for each species (this is done later)
    real :: luminosity, localColDens, rho, del(NDIM), radius_sink, flux, fshield
    real :: dColDens, dlcell, NH2
    logical :: printval
    REAL,DIMENSION(NDIM) :: target_point

    do b=1, nrOfLeafBlocks
      blk=leafList(b)
      call Grid_getBlkPtr(blk, solnData)
      call Grid_getDeltas(blk, del)
      do n=IAXIS,KAXIS
        call Grid_getCellCoords(n, blk, RIGHT_EDGE, .TRUE., xr(:,n), MAXCELLS)
      end do      
      do iAngle=iAb,iAe
        iAindex = (iAngle-iAb)+1
        
        !Position of point source
        PointSource = particles_global(POSX_PART_PROP:POSZ_PART_PROP,QSindex(iAngle))

        do k=ib(3)+K3D,ie(3)
          do j=ib(2)+K2D,ie(2)
            do i=ib(1)+1,ie(1)
                dest(1) = gr_oneBlock(blk)%firstAxisCoords(CENTER, i)
                dest(2) = gr_oneBlock(blk)%secondAxisCoords(CENTER, j)
                dest(3) = gr_oneBlock(blk)%thirdAxisCoords(CENTER, k)

                dirFact = dest-PointSource
                dirFact = dirFact/SQRT(SUM(dirFact**2))

                printval = .false.

                call calc_local_block_contributions_sink(&
                        regFact, blk, PointSource, dest, dirFact, maxLevel, ColDens_local, &
                        0,dColDens,dlcell,printval)

                call create_cut_block_list_3DRT(&
                        dr_globalMe, blk,ColDens_global, intensity_global, PointSource, dest, &
                        regFact, dirFact, iAindex, maxNrOfLeafBlocks, dr_globalNumProcs, maxLevel, .true.)

                ! Now add up global and local column densities
                ColDens = ColDens_global + ColDens_local
                !calc_local_block_contributions_sink returns the column to the destination cell's entry face
                !(photon-conserving convention); dColDens is the destination cell's own segment. The LW rates
                !and shielding factors are for the cell, so evaluate both at its centre: dust and species columns
                !alike (at the entry face the flux is ~2x too high for tau_cell ~ 1).
                !In the sink's own cell the segment starts at the sink, so take the sink-to-centre part of it.
                if(dlcell > 0.0 .and. all(abs(dest-PointSource) .le. 0.5*del)) then
                  ColDens = ColDens + dColDens*min(1.0, sqrt(sum((dest-PointSource)**2))/dlcell)
                else
                  ColDens = ColDens + 0.5*dColDens
                endif

                if(fsh_spec == "dust") then
                  !For dust, no shielding factor calculation. Just store the optical depth for later use. 
                  taud_star(b,i-ib(1),j-ib(2),k-ib(3),QSindex(iAngle)) = ColDens
                  !return because no shielding calculations required here
                  cycle
                else
                  !Use the stored value to account for attenuation due to dust
                  taud = taud_star(b,i-ib(1),j-ib(2),k-ib(3),QSindex(iAngle))
                endif

                if(fsh_spec == "H2") then
                  !For dust, no shielding factor calculation. Just store the optical depth for later use. 
                  NH2_star(b,i-ib(1),j-ib(2),k-ib(3),QSindex(iAngle)) = ColDens
                else if(fsh_spec == "C" .or. fsh_spec == "CO") then
                  !Use the stored value to account for attenuation due to dust
                  NH2 = NH2_star(b,i-ib(1),j-ib(2),k-ib(3),QSindex(iAngle))
                endif

                !Floor r^2 at sigma_sink^2 (geometric dilution only; the columns above are the
                !actual path), as for JSTR and UVF0 in rt_sinkInject, so LWFL = UVF0_LW*exp(-taud).
                distSqr_sink = MAX(SUM((dest-PointSource)**2),sigma_sink**2)

                !Prevent any radiation within the radius of the star
#ifdef STELLAR_RADIUS_PART_PROP
                radius_sink = particles_global(STELLAR_RADIUS_PART_PROP,QSindex(iAngle))
#else
                radius_sink = 0.0
#endif
                if (distSqr_sink .ge. radius_sink**2) then
#ifdef LUMLW_PART_PROP
                    luminosity = particles_global(LUMLW_PART_PROP,QSindex(iAngle))
#elif LUMFUV_PART_PROP
                    luminosity = particles_global(LUMFUV_PART_PROP,QSindex(iAngle))
#elif defined(LUMINOSITY_PART_PROP)
                    luminosity = particles_global(LUMINOSITY_PART_PROP,QSindex(iAngle))
#else
                    luminosity = 0.0
#endif
                
                    !Unattenuated radiation flux
                    flux = luminosity / (4.0*PI*distSqr_sink)

                else 
                    flux = 0.0

                end if


#ifdef UVF0_VAR
                !Unattenuated stellar LW flux for chi0 (the PE part comes from rt_sinkInject); same floor and
                !stellar-radius cut as LWFL, so its ratio to LWFL is exactly exp(-taud)
                if(fillMoments .and. fsh_spec == "H2") solnData(UVF0_VAR,i,j,k) = solnData(UVF0_VAR,i,j,k) + flux
#endif

                !Attenuate flux due to dust absorption
                flux = flux * exp(-taud)

                !Flux vector of the streaming stellar field, for the radiation pressure in rt_dustTerms
                if(fillMoments .and. fsh_spec == "H2") then
                  solnData(MOHX_VAR,i,j,k) = solnData(MOHX_VAR,i,j,k) + flux*dirFact(1)
#if NDIM>1
                  solnData(MOHY_VAR,i,j,k) = solnData(MOHY_VAR,i,j,k) + flux*dirFact(2)
#if NDIM>2
                  solnData(MOHZ_VAR,i,j,k) = solnData(MOHZ_VAR,i,j,k) + flux*dirFact(3)
#endif
#endif
                endif

                !Get the shielding factor for the obtained column density
                if(fsh_spec == "H2") then
                  fshield = get_fshield_H2(ColDens,bfive)
                  !Use the fshield container to cumulatively add the total unattenuated flux
                  solnData(FSHM_VAR,i,j,k) = solnData(FSHM_VAR,i,j,k) + flux * fshield
                else if(fsh_spec == "H") then
                  fshield = get_fshield_H(ColDens,bfive)
                  !Use the fshield container to cumulatively add the total unattenuated flux
                  solnData(FSHA_VAR,i,j,k) = solnData(FSHA_VAR,i,j,k) + flux * fshield

                else if(fsh_spec == "C") then
                  !Obtain the shielding due to H2 and C for C
                  fshield = get_fshield_C(NH2, ColDens, bfive)
                  solnData(IC_VAR,i,j,k) = solnData(IC_VAR,i,j,k) + flux * fshield

                else if(fsh_spec == "CO") then
                  !Obtain the shielding due to H2 and CO for CO
                  fshield = get_fshield_CO(NH2, ColDens, bfive)
                  solnData(DICO_VAR,i,j,k) = solnData(DICO_VAR,i,j,k) + flux * fshield
                endif

                !Flux-weighted shielding: this is not the final value of the LW flux used
                !This counter is used to compute the average flux-weighted fshield over all species; will be replaced later
                solnData(LWFL_VAR,i,j,k) = solnData(LWFL_VAR,i,j,k) + flux 
            enddo ! i
          enddo !j
        enddo !k
      enddo !iAngle
      call Grid_releaseBlkPtr(blk, solnData)
    enddo ! blocks

  END SUBROUTINE rad_calcContrib


  SUBROUTINE SetOPAC_LW(blk)
    use rt_ionisedata, ONLY: h2A, hA
    use Opacity_data, ONLY: dust_lwopac_value, dust_opacity_fuv, dusttoGasRatio
    IMPLICIT NONE
    integer, intent(in) :: blk
    real :: nH2, nH, nC, nCO
    integer :: i,j,k
    real, dimension(NDIM) :: position, sink_pos

    if(fsh_spec == "dust" .and. lwdiss_type == "raytrace" .and. dust_opacity_fuv .ne. 'Fixed') then
      call Driver_abortFlash("lwdiss_type = raytrace currently only supports dust_opacity_fuv = Fixed.")
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
              if(lwdiss_type == "raytrace") then
                !TAUP_VAR can't be used here: RadTrans skips rt_setOpacity for the LW band in this mode, so it holds the
                !previous band's opacity (e.g. EUV dust + HI photoionization). Recompute the LW dust opacity as rt_setOpacity does.
                solnData(OPAC_VAR,i,j,k) = max(dust_lwopac_value * dusttoGasRatio * solnData(DENS_VAR,i,j,k), &
                                               1.e-10 * solnData(DENS_VAR,i,j,k))
#ifdef TEMP_VAR
                !Thermal sputtering of dust
                if(solnData(TEMP_VAR,i,j,k) .gt. 1.e6) solnData(OPAC_VAR,i,j,k) = 1.e-10 * solnData(DENS_VAR,i,j,k)
#endif
              else
                !hybrid: called after the LW moment solve, so TAUP_VAR holds the LW-band dust opacity
                solnData(OPAC_VAR,i,j,k) = solnData(TAUP_VAR,i,j,k)
              endif
            endif

            
          

        end do 
      end do
    end do
    call Grid_releaseBlkPtr(blk, solnData)

  END SUBROUTINE SetOPAC_LW

  SUBROUTINE average_fshield(blk)
    !This subroutine computes the dust-attenuate flux-weighted shielding factor for H2, H, C and CO
    IMPLICIT NONE
    integer, intent(in) :: blk
    integer :: i,j,k
    real :: totalflux, fshield

    call Grid_getBlkPtr(blk, solnData)
    do k=GRID_KLO_GC,GRID_KHI_GC
        do j=GRID_JLO_GC,GRID_JHI_GC
            do i=GRID_ILO_GC,GRID_IHI_GC
                totalflux = solnData(LWFL_VAR,i,j,k)
                !The total cumulative flux from sinks is stored in FSHM_VAR or FSHA_VAR (depending on which is being updated now)
                if(fsh_spec == "H2") then
                  fshield = solnData(FSHM_VAR,i,j,k)
                else if(fsh_spec == "H") then
                  fshield = solnData(FSHA_VAR,i,j,k)
                else if(fsh_spec == "C") then
                  fshield = solnData(IC_VAR,i,j,k)
                else if(fsh_spec == "CO") then
                  fshield = solnData(DICO_VAR,i,j,k)
                endif
                !Now replace this var with the appropriate data, i.e. the avg shielding factor in FSHL_VAR
                if(totalflux .gt. 0) then
                    fshield = fshield/totalflux
                else
                    fshield = 0.0
                endif
                
                !Now store the shielding quantities in their containers
                if(fsh_spec == "H2") then
                  solnData(FSHM_VAR,i,j,k) = fshield
                else if(fsh_spec == "H") then
                  solnData(FSHA_VAR,i,j,k) = fshield
                else if(fsh_spec == "C") then
                  solnData(IC_VAR,i,j,k) = fshield
                else if(fsh_spec == "CO") then
                  solnData(DICO_VAR,i,j,k) = fshield
                endif
                    
            end do
        end do
    end do
    call Grid_releaseBlkPtr(blk, solnData)

  END SUBROUTINE average_fshield

  SUBROUTINE rad_resetLWContainers
    implicit none
    INTEGER :: b, blk

    call Timers_start("rad_resetData")

    do b=1,nrOfLeafBlocks
      !
      blk = leafList(b)

      call Grid_getBlkPtr(blk, solnData)
      !Molecular
      if(fsh_spec == "H2") solnData(FSHM_VAR,:,:,:) = 0.0
      if(fsh_spec == "H") solnData(FSHA_VAR,:,:,:) = 0.0
      if(fsh_spec == "C") solnData(IC_VAR,:,:,:) = 0.0
      if(fsh_spec == "CO") solnData(DICO_VAR,:,:,:) = 0.0
      solnData(LWFL_VAR,:,:,:) = 0.0
      call Grid_releaseBlkPtr(blk, solnData)

    enddo ! b
    !
    call Timers_stop("rad_resetData")
    !
  end subroutine rad_resetLWContainers




END SUBROUTINE RayTrace_LW


