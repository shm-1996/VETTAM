!!****if* source/Simulation/SimulationMain/StirFromFile/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!  Simulation_initBlock(integer, intent(IN)  :: blockid)
!!
!! DESCRIPTION
!!  Initializes data (density, pressure, velocity, etc.) for
!!  a specified block for the Marshak Wave Test. 
!!
!! ARGUMENTS
!!   blockid : ID of block in current processor
!!
!! AUTHOR
!!   Shyam Harimohan Menon, 2020
!!
!!***


subroutine Simulation_initBlock (blockID)

#include "Flash.h"
#include "constants.h"

  use Driver_interface, ONLY : Driver_getSimTime, Driver_getMype
  use Driver_data, ONLY : dr_dtInit
  use Simulation_data
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Grid_interface, ONLY   : Grid_getBlkIndexLimits, Grid_getCellCoords, &
                               Grid_getBlkPtr, Grid_releaseBlkPtr,Grid_putPointData, &
                               Grid_getDeltas
  use Eos_interface, ONLY    : Eos_wrapped
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use RadTrans_data, ONLY: rt_speedlt, rt_radconst, rt_dtFactor, &
                           rt_boundary_T, rt_boundary_T_face, rt_cv_alpha
  use Eos_data, ONLY : eos_gasConstant, eos_gamma, eos_singleSpeciesA
  use eos_idealGammaData, ONLY          : eos_gammam1
  use Opacity_interface, ONLY: Opacity

#ifdef RAYTRACE_3DRT
    use RadTrans_data, ONLY : rt_useRadTrans
#endif

  implicit none 

  integer, intent(IN) :: blockID

  integer, save                :: myPE
  integer                      :: sizeZ, sizeY, sizeX, i, j, k
  integer, dimension(2,MDIM)   :: blkLimits, blkLimitsGC
  logical, parameter           :: gcell = .true.
  real                         :: del(MDIM), dvol, delx
  real                         :: rho, trad, temp, eint, pres, ener,vx,vy,vz,erad,frad,F_inc
  real                         :: ggprod, abar, abarinv, gam1inv
  real                         :: sim_kB, sim_mH
  real, DIMENSION(:,:,:,:), POINTER :: solnData
#if NFACE_VARS > 0
  real, DIMENSION(:,:,:,:), POINTER :: facexData, faceyData, facezData
#endif
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC) :: xCoord
  real, dimension(GRID_JHI_GC) :: yCoord
  real, dimension(GRID_KHI_GC) :: zCoord
#else
  real, allocatable, dimension(:) :: xCoord, yCoord, zCoord
  integer :: istat
#endif
  integer, dimension(MDIM) :: axis
  logical, parameter :: Debug = .false.


  call Driver_getMype(GLOBAL_COMM, myPE)

  if (Debug) print *, '[', myPE, '] Simulation_initBlock entering.'

  !get the index limits of the block
  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

  ! get a pointer to the current block of data
  call Grid_getBlkPtr(blockID, solnData,CENTER)
  if (Debug) print *, '[', myPE, '] Simulation_initBlock: got solnData pointer.'

  !getting the dx's
  call Grid_getDeltas(blockID, del)
  call Grid_getMinCellSize(delx)
#if NDIM == 1
   dvol = del(IAXIS)
#endif
#if NDIM == 2
   dvol = del(IAXIS) * del(JAXIS)
#endif
#if NDIM == 3
   dvol = del(IAXIS) * del(JAXIS) * del(KAXIS)
#endif

   !Set initial timestep to light crossing time
   dr_dtInit = rt_dtFactor * delx/rt_speedlt

   sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
   sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
   sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE
   allocate(xCoord(sizeX),stat=istat)
   if (istat .ne. 0) call Driver_abortFlash("could not allocate xCoord in Simulation_initBlock.F90")
   allocate(yCoord(sizeY),stat=istat)
   if (istat .ne. 0) call Driver_abortFlash("could not allocate yCoord in Simulation_initBlock.F90")
   allocate(zCoord(sizeZ),stat=istat)
   if (istat .ne. 0) call Driver_abortFlash("could not allocate zCoord in Simulation_initBlock.F90")
#endif

   ! x coordinates
   call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell,xCoord,sizeX)
#if NDIM > 1
   ! y coordinates
   call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell,yCoord,sizeY)
#if NDIM > 2
   ! z coordinates
   call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell,zCoord,sizeZ)
#endif
#endif

  ! Initialise some thermodynamic quantities
  abar      = eos_singleSpeciesA
  abarinv   = 1. / eos_singleSpeciesA
  ggprod    = eos_gammam1 * eos_gasConstant
  gam1inv   = 1. / eos_gammam1

  !Set simulation alpha value: alpha = 4*a_R/epsilon. The solver needs it for
  !the rt_temp_type=5 closed-form update, so hand it over; RadTrans_init has
  !already run by this point.
  sim_alpha = 4*rt_radconst/sim_epsilon
  rt_cv_alpha = sim_alpha

  !Check if rt_boundary_T consistent with the incoming flux
  F_inc = rt_radconst*rt_speedlt/4 *rt_boundary_T**4 
  if(F_inc .ne. sim_Finc) then 
    if(myPE .eq. 0) then 
      print *, "rt_boundary_T not consistent with sim_Finc. Changing rt_boundary_T."
    endif
    rt_boundary_T = ((4*sim_finc)/(rt_speedlt*rt_radconst))**0.25
    !The Marshak boundary reads the per-face array, which RadTrans_init already
    !resolved, so it has to be corrected alongside rt_boundary_T.
    rt_boundary_T_face(:,:) = rt_boundary_T
  endif

  call PhysicalConstants_get("Boltzmann",sim_kB)
  call PhysicalConstants_get("proton mass",sim_mH)

   ! loop over cells in block
   do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
       do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              rho = sim_rhoAmbient
              vx = 0.
              vy = 0.
              vz = 0.
              erad = sim_erAmbient
              frad = 0.0
              temp = sim_tempAmbient
              solnData(TEMP_VAR,i,j,k) = temp
              solnData(ERAD_VAR,i,j,k) = erad
              solnData(MOHX_VAR,i,j,k) = 0.0
#if NDIM>1              
              solnData(MOHY_VAR,i,j,k) = 0.0
#if NDIM>2              
              solnData(MOHZ_VAR,i,j,k) = 0.0
#endif
#endif              
           !Assuming thermal equilibrium 
#ifdef PRES_VAR
              pres = rho*sim_kB*temp/(sim_mH*abar)
#endif
#ifdef EINT_VAR
              eint = pres/(sim_gamma-1.0) * rho
#endif
#ifdef ENER_VAR
              ener = eint + 0.5*(vx**2 + vy**2 + vz**2)
#endif


           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           call Grid_putPointData(blockID,CENTER,DENS_VAR,EXTERIOR,axis,rho)
           call Grid_putPointData(blockID,CENTER,VELX_VAR,EXTERIOR,axis,vx)
#if NDIM>1           
           call Grid_putPointData(blockID,CENTER,VELY_VAR,EXTERIOR,axis,vy)
#if NDIM>2           
           call Grid_putPointData(blockID,CENTER,VELZ_VAR,EXTERIOR,axis,vz)
#endif
#endif           
           call Grid_putPointData(blockID,CENTER,TEMP_VAR,EXTERIOR,axis,temp)
           call Grid_putPointData(blockID,CENTER,PRES_VAR,EXTERIOR,axis,pres)
           call Grid_putPointData(blockID,CENTER,EINT_VAR,EXTERIOR,axis,eint)
           call Grid_putPointData(blockID,CENTER,ENER_VAR,EXTERIOR,axis,ener)

           solnData(XXED_VAR,i,j,k) = sim_eddvalue
#if NDIM>1
           solnData(YYED_VAR,i,j,k) = 0.33
           solnData(XYED_VAR,i,j,k) = 0.0

#if NDIM>2
           solnData(ZZED_VAR,i,j,k) = 0.33
           solnData(XZED_VAR,i,j,k) = 0.0
           solnData(YZED_VAR,i,j,k) = 0.0
#endif
#endif



           enddo ! i
       enddo ! j
   enddo ! k

   call Grid_releaseBlkPtr(blockID, solnData,CENTER)

#ifndef FIXEDBLOCKSIZE
   deallocate(xCoord)
#if NDIM>1   
   deallocate(yCoord)
#if NDIM>2   
   deallocate(zCoord)
#endif
#endif   
#endif

  ! update eint and temp by calling EOS
  call Eos_wrapped(MODE_DENS_TEMP, blkLimits, blockID)

  if (Debug) print *, '[', myPE, '] Simulation_initBlock exiting.'

end subroutine Simulation_initBlock
