!!****if* source/physics/RadTrans/RadTransMain/VETTAM/RayTrace/raytrace_data
!!
!!  NAME 
!!    raytrace_data
!!
!!  SYNOPSIS
!!    use raytrace_data
!!
!!  DESCRIPTION 
!!    Stores data for HybridChar3DRT Ray Tracer
!!
!!***

#include "Flash.h"

module raytrace_data
  use gr_interfaceTypeDecl, ONLY: AllBlockRegions_t
  implicit none

  logical, save :: rt_useRayTrace
  integer, save :: rt_nPhi
  integer, save :: rt_nTheta

  integer, save :: rt_nrOfAngleGroups
  integer, save :: rt_maxNrOfBoundIter
  integer, save :: rt_ALI
  real, save    :: rt_epsilon

  integer, save :: rt_healpix_nSide
  integer, save :: rt_healpix_randomize
  
  real, save :: rt_dirX
  real, save :: rt_dirY
  real, save :: rt_dirZ

  real, save :: irradiation
  !
  ! Internal reference at which dr_simGeneration the memory allocation
  ! has been done. If dr_simGeneration is raised, a refinement or
  ! derefinement has happend and we have to reallocate memory as well
  ! as initialize it again.
  INTEGER, save :: allocGeneration

  ! rt_needOnlyFaces = .true. : every rank stores its own face values plus only
  ! the remote face values its rays will read (see 3DRT/rt_faceComm.F90).
  ! rt_needOnlyFaces = .false.: legacy behaviour, full MPI_ALLGATHER of all faces.
  logical, save :: rt_needOnlyFaces = .true.

  ! rt_debugChecksum = .true.: log an exact checksum of the output fields at the end of
  ! every ray-tracer call (3DRT/rt_checksum.F90)
  logical, save :: rt_debugChecksum = .false.

  ! Face values (tau per length / source fn, or column densities for sink rays).
  ! faceValueAll(ib(1):ie(1), ib(2):ie(2), face(3), var(2), slot)
  ! The slot of (leaf index, angle index in group, proc) is faceSlot(leaf, iAindex, proc);
  ! 0 means "not stored on this rank". In legacy mode the slot layout is identical to the
  ! old 7-D array faceValueAll(...,maxNrOfLeafBlocks,nrOfAnglesPerGroup,0:nProcs-1).
  real, dimension(:,:,:,:,:), ALLOCATABLE :: faceValueAll
  integer, dimension(:,:,:), ALLOCATABLE    :: faceSlot
  integer, dimension(:), ALLOCATABLE        :: nrOfLeafBlocksList
  integer, dimension(:,:), ALLOCATABLE      :: leafListAll

  real, dimension(:,:,:,:), pointer :: solnData

  integer, dimension(:,:,:,:), allocatable :: SurrBlkSumArray
  integer, dimension(:,:), allocatable :: SurrBlkSumNumNegh
  type (AllBlockRegions_t), allocatable :: SurrBlkSum(:)
  integer, DIMENSION(:,:), allocatable :: rblockListAll

end module raytrace_data
