!===============================================================================
!
! Module: rt_faceComm
! Path: source/physics/RadTrans/RadTransMain/VETTAM/RayTrace/3DRT/rt_faceComm.F90
!
!===============================================================================
!
! Storage layout and MPI exchange of block face values for the hybrid
! characteristics ray tracers (diffuse VETTAM rays, LW sink rays, LW diffuse rays).
!
! Usage per angle group (all ranks must call both routines collectively):
!
!   call rt_faceLayout(nAng, geomType, geom)   ! sets faceSlot, sizes faceValueAll
!   ... write own faces: faceValueAll(:,:,n,v,faceSlot(b,iAindex,dr_globalMe)) ...
!   call rt_faceExchange(nAng)                 ! fills the remote slots
!   ... create_cut_block_list_3DRT reads faceValueAll(:,:,:,:,faceSlot(leaf,iA,proc))
!
! rt_needOnlyFaces = .false. (legacy): faceSlot is the dense layout of the old 7-D array
!   faceValueAll(...,maxNrOfLeafBlocks,nrOfAnglesPerGroup,0:nProcs-1) and the exchange is
!   the original MPI_ALLGATHER(V), so results are identical to the old code.
!
! rt_needOnlyFaces = .true.: every rank computes its "need set", i.e. the remote leaf
!   blocks (per angle) whose face values its backward-traced rays can read, requests
!   exactly those (MPI_Alltoall of counts + MPI_Alltoallv of request lists) and receives
!   them with one MPI_Alltoallv into a compact buffer behind its own face values.
!
!   Need set (geometric, conservative, AMR-safe): rays start at cell centres of a local
!   block b, i.e. inside box(b) shrunk by half a cell of b's own level, and run
!     geomType = RT_GEOM_DIR  : backwards along -dir to the domain boundary
!                               (points P - s*dir, s >= 0)
!     geomType = RT_GEOM_POINT: towards the point source c (sink rays)
!                               (points (1-s)*P + s*c, 0 <= s <= 1)
!   A remote leaf block B is needed iff its closed box, enlarged by 1/4 cell at the
!   finest level, intersects that swept set for some local block. For each axis the
!   overlap condition is linear in s, so the test is exact for boxes (interval of s).
!   The traversal in create_cut_block_list_3DRT only reads the face values of leaf
!   blocks the (fully resolved) ray enters, so this is a superset of what is read for
!   any mix of refinement levels. Need sets are recomputed every angle group, so
!   randomised HEALPix orientations are handled. A miss is detected in
!   create_cut_block_list_3DRT (faceSlot == 0 -> Driver_abortFlash), always on.
!   Rays wrapping around periodic boundaries are not followed: rt_init_raytrace_3drt
!   forces rt_needOnlyFaces = .false. (with a warning) if any boundary is periodic.
!
!===============================================================================

#include "Flash.h"
#include "constants.h"

module rt_faceComm

  implicit none
  private
#include "Flash_mpi.h"

  integer, parameter, public :: RT_GEOM_DIR = 1, RT_GEOM_POINT = 2

  public :: rt_faceLayout, rt_faceExchange

  ! number of reals per slot
  integer, parameter :: nFace = (NXB+1)*(NYB+1)*3*2

  integer, save :: faceType = MPI_DATATYPE_NULL
  logical, save :: layoutNeedOnly = .false.
  integer, save :: nOwnSlots = 0, nRecvSlots = 0, nSendSlots = 0

  integer, allocatable, save :: reqCnt(:), reqDsp(:), sendCnt(:), sendDsp(:)
  integer, allocatable, save :: sendList(:)
  real,    allocatable, save :: sendBuf(:,:,:,:,:)

contains

  !-----------------------------------------------------------------------------
  subroutine rt_faceLayout(nAng, geomType, geom)
    use raytrace_data, ONLY : rt_needOnlyFaces, faceValueAll, faceSlot
    use RadTrans_RayTrace_3DRT, ONLY : maxNrOfLeafBlocks
    use rt_data_raytrace_3drt, ONLY : nrOfAnglesPerGroup
    use Driver_data, ONLY : dr_globalNumProcs
    use Driver_interface, ONLY : Driver_abortFlash
    use Timers_interface, ONLY : Timers_start, Timers_stop
    implicit none
    integer, intent(in) :: nAng, geomType
    real, intent(in)    :: geom(MDIM, nAng)
    integer :: p, a, b, nSlots, ierr

    if(nAng > nrOfAnglesPerGroup .or. nAng < 1) &
      call Driver_abortFlash("rt_faceLayout: nAng out of range (must be 1..nrOfAnglesPerGroup)")

    if(faceType == MPI_DATATYPE_NULL) then
      call MPI_Type_contiguous(nFace, FLASH_REAL, faceType, ierr)
      call MPI_Type_commit(faceType, ierr)
    end if

    if(.not. rt_needOnlyFaces) then
      ! ---------------- legacy dense layout ----------------
      layoutNeedOnly = .false.
      nSlots = maxNrOfLeafBlocks * nrOfAnglesPerGroup * dr_globalNumProcs
      call ensureCapacity(nSlots, .true.)
      do p = 0, dr_globalNumProcs-1
        do a = 1, nrOfAnglesPerGroup
          do b = 1, maxNrOfLeafBlocks
            faceSlot(b,a,p) = b + maxNrOfLeafBlocks*((a-1) + nrOfAnglesPerGroup*p)
          end do
        end do
      end do
      return
    end if

    ! ---------------- need-only layout ----------------
    layoutNeedOnly = .true.
    call Timers_start("rad_faceNeedSet")
    call buildNeedSet(nAng, geomType, geom)
    call Timers_stop("rad_faceNeedSet")

  end subroutine rt_faceLayout

  !-----------------------------------------------------------------------------
  subroutine ensureCapacity(nSlots, exact)
    use raytrace_data, ONLY : faceValueAll
    use RadTrans_RayTrace_3DRT, ONLY : ib, ie
    use Driver_interface, ONLY : Driver_abortFlash
    implicit none
    integer, intent(in) :: nSlots
    logical, intent(in) :: exact
    integer :: newCap, ierr
    if(allocated(faceValueAll)) then
      if(exact) then
        if(size(faceValueAll,5) == nSlots) return
      else
        ! grow if too small; shrink if far too large (e.g. after a direction change)
        if(size(faceValueAll,5) >= nSlots .and. size(faceValueAll,5) <= 2*nSlots+16) return
      end if
      deallocate(faceValueAll)
    end if
    if(exact) then
      newCap = max(nSlots,1)
    else
      newCap = max(nSlots + nSlots/8 + 1, 1)
    end if
    allocate(faceValueAll(ib(1):ie(1), ib(2):ie(2), 3, 2, newCap), stat=ierr)
    if(ierr /= 0) call Driver_abortFlash("rt_faceComm: could not allocate faceValueAll")
  end subroutine ensureCapacity

  !-----------------------------------------------------------------------------
  ! Does the set { (1-s)*Pbox + s*c } (point mode, s in [0,1]) or
  ! { Pbox - s*dir } (direction mode, s >= 0) intersect box [lo2,hi2] ?
  ! In both cases the lower/upper faces of the moving box are lo + s*w1 and hi + s*w2.
  pure logical function sweepHits(lo, hi, w1, w2, smax, lo2, hi2)
    implicit none
    real, intent(in) :: lo(MDIM), hi(MDIM), w1(MDIM), w2(MDIM), smax, lo2(MDIM), hi2(MDIM)
    real :: smin, sup, aa, bb
    integer :: i, c
    smin = 0.0
    sup  = smax
    sweepHits = .false.
    do i = 1, NDIM
      do c = 1, 2
        ! constraint aa + s*bb <= 0
        if(c == 1) then
          aa = lo(i) - hi2(i);  bb = w1(i)        ! lo + s*w1 <= hi2
        else
          aa = lo2(i) - hi(i);  bb = -w2(i)       ! hi + s*w2 >= lo2
        end if
        if(bb == 0.0) then
          if(aa > 0.0) return
        else if(bb > 0.0) then
          sup = min(sup, -aa/bb)
        else
          smin = max(smin, -aa/bb)
        end if
      end do
    end do
    sweepHits = (smin <= sup)
  end function sweepHits

  !-----------------------------------------------------------------------------
  subroutine buildNeedSet(nAng, geomType, geom)
    use raytrace_data, ONLY : faceSlot, nrOfLeafBlocksList, leafListAll
    use RadTrans_RayTrace_3DRT, ONLY : leafList, nrOfLeafBlocks, maxNrOfLeafBlocks, maxLevel
    use RadTrans_TreeCommData, ONLY : g_bndbox
    use RadTrans_hybridCharModule, ONLY : io
    use Grid_data, ONLY : gr_globalDomain, gr_nBlockX, gr_nBlockY, gr_nBlockZ
    use Driver_data, ONLY : dr_globalMe, dr_globalNumProcs
    use Driver_interface, ONLY : Driver_abortFlash
    use Timers_interface, ONLY : Timers_start, Timers_stop
    implicit none
    integer, intent(in) :: nAng, geomType
    real, intent(in)    :: geom(MDIM, nAng)

    real :: lbLo(MDIM, max(nrOfLeafBlocks,1)), lbHi(MDIM, max(nrOfLeafBlocks,1))
    real :: uLo(MDIM), uHi(MDIM), lo2(MDIM), hi2(MDIM), tol(MDIM), half(MDIM)
    real :: w1(MDIM), w2(MDIM), smax, nbCells(MDIM)
    integer, allocatable :: reqList(:), recvReq(:)
    integer :: p, a, b, i, k, blk, nReqMax, nTot, ierr, me, np, code
    integer :: stats(2), statsSum(2), statsMax(1), statsMin(1)
    logical :: hit

    me = dr_globalMe
    np = dr_globalNumProcs
    nOwnSlots = nrOfLeafBlocks * nAng

    if(.not. allocated(reqCnt)) then
      allocate(reqCnt(0:np-1), reqDsp(0:np-1), sendCnt(0:np-1), sendDsp(0:np-1))
    end if

    ! ---- own slots ----
    faceSlot = 0
    do a = 1, nAng
      do b = 1, nrOfLeafBlocks
        faceSlot(b,a,me) = b + nrOfLeafBlocks*(a-1)
      end do
    end do

    call Timers_start("rad_faceNeedGeom")
    ! ---- geometry of local ray origins (cell centres) ----
    nbCells = (/ real(NXB), real(NYB), real(NZB) /)
    uLo =  huge(1.0)
    uHi = -huge(1.0)
    do b = 1, nrOfLeafBlocks
      blk = leafList(b)
      half = 0.5 * (g_bndbox(2,:,blk,me) - g_bndbox(1,:,blk,me)) / nbCells
      lbLo(:,b) = g_bndbox(1,:,blk,me) + half
      lbHi(:,b) = g_bndbox(2,:,blk,me) - half
      uLo = min(uLo, lbLo(:,b))
      uHi = max(uHi, lbHi(:,b))
    end do
    ! tolerance: 1/4 of a cell at the finest level (>> the 1e-9 cell fudges of the tracer)
    tol(1) = 0.25 * (gr_globalDomain(HIGH,IAXIS)-gr_globalDomain(LOW,IAXIS)) &
                  / (real(gr_nBlockX) * 2.0**(maxLevel-1) * real(NXB))
    tol(2) = 0.25 * (gr_globalDomain(HIGH,JAXIS)-gr_globalDomain(LOW,JAXIS)) &
                  / (real(gr_nBlockY) * 2.0**(maxLevel-1) * real(NYB))
    tol(3) = 0.25 * (gr_globalDomain(HIGH,KAXIS)-gr_globalDomain(LOW,KAXIS)) &
                  / (real(gr_nBlockZ) * 2.0**(maxLevel-1) * real(NZB))

    nTot = SUM(nrOfLeafBlocksList)
    nReqMax = max((nTot - nrOfLeafBlocks) * nAng, 1)
    allocate(reqList(nReqMax))

    ! ---- need set: loop over remote procs (ascending), angles, remote leaves ----
    k = 0
    reqCnt = 0
    if(nrOfLeafBlocks > 0) then
      do p = 0, np-1
        if(p == me) cycle
        do a = 1, nAng
          if(geomType /= RT_GEOM_POINT) then
            ! backward ray: P - s*dir, s in [0,inf)
            w1 = -geom(:,a)
            w2 = -geom(:,a)
            smax = huge(1.0)
          end if
          do i = 1, nrOfLeafBlocksList(p)
            blk = leafListAll(i,p)
            lo2 = g_bndbox(1,:,blk,p) - tol
            hi2 = g_bndbox(2,:,blk,p) + tol
            if(geomType == RT_GEOM_POINT) then
              ! segment towards the point source: (1-s)*P + s*c, s in [0,1]
              w1 = geom(:,a) - uLo
              w2 = geom(:,a) - uHi
              smax = 1.0
            end if
            if(.not. sweepHits(uLo, uHi, w1, w2, smax, lo2, hi2)) cycle
            hit = .false.
            do b = 1, nrOfLeafBlocks
              if(geomType == RT_GEOM_POINT) then
                w1 = geom(:,a) - lbLo(:,b)
                w2 = geom(:,a) - lbHi(:,b)
              end if
              if(sweepHits(lbLo(:,b), lbHi(:,b), w1, w2, smax, lo2, hi2)) then
                hit = .true.
                exit
              end if
            end do
            if(hit) then
              k = k + 1
              reqList(k) = i + maxNrOfLeafBlocks*(a-1)
              reqCnt(p) = reqCnt(p) + 1
              faceSlot(i,a,p) = nOwnSlots + k
            end if
          end do
        end do
      end do
    end if
    nRecvSlots = k
    call Timers_stop("rad_faceNeedGeom")

    ! ---- tell owners what we need ----
    ! (the first collective of an angle group: includes waiting for the slowest rank)
    call MPI_Alltoall(reqCnt, 1, MPI_INTEGER, sendCnt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    reqDsp(0) = 0
    sendDsp(0) = 0
    do p = 1, np-1
      reqDsp(p)  = reqDsp(p-1)  + reqCnt(p-1)
      sendDsp(p) = sendDsp(p-1) + sendCnt(p-1)
    end do
    nSendSlots = sendDsp(np-1) + sendCnt(np-1)

    if(allocated(sendList)) then
      if(size(sendList) < nSendSlots) deallocate(sendList)
    end if
    if(.not. allocated(sendList)) allocate(sendList(max(nSendSlots,1)))
    allocate(recvReq(max(nSendSlots,1)))
    call MPI_Alltoallv(reqList, reqCnt, reqDsp, MPI_INTEGER, &
                       recvReq, sendCnt, sendDsp, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    do k = 1, nSendSlots
      code = recvReq(k) - 1
      i = MOD(code, maxNrOfLeafBlocks) + 1
      a = code / maxNrOfLeafBlocks + 1
      if(i > nrOfLeafBlocks .or. a > nAng) &
        call Driver_abortFlash("rt_faceComm: invalid face request received")
      sendList(k) = i + nrOfLeafBlocks*(a-1)   ! own slot
    end do
    deallocate(reqList, recvReq)

    ! ---- storage: own + needed remote ----
    ! (+1 so that the receive address faceValueAll(...,nOwnSlots+1) is always valid)
    call ensureCapacity(nOwnSlots + nRecvSlots + 1, .false.)

    ! ---- statistics (master writes to the RT log) ----
    stats(1) = nRecvSlots
    stats(2) = (nTot - nrOfLeafBlocks) * nAng      ! remote (block,angle) pairs in all-gather
    call MPI_Reduce(stats, statsSum, 2, MPI_INTEGER, MPI_SUM, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(stats(1:1), statsMax, 1, MPI_INTEGER, MPI_MAX, MASTER_PE, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(stats(1:1), statsMin, 1, MPI_INTEGER, MPI_MIN, MASTER_PE, MPI_COMM_WORLD, ierr)
    if(me == MASTER_PE) then
      write(io,'(a,i2,a,i5,a,i10,a,i8,a,i8,a,f8.5,a,f8.5)') ' [faceComm] geom=', geomType, &
        ' nAng=', nAng, ' recvSlots(sum)=', statsSum(1), ' min=', statsMin(1), ' max=', statsMax(1), &
        ' frac(vs remote all-gather)=', real(statsSum(1))/real(max(statsSum(2),1)), &
        ' frac(vs padded all-gather)=', real(statsSum(1)) / &
           real(max(maxNrOfLeafBlocks*nAng*(np-1)*np,1))
    end if
  end subroutine buildNeedSet

  !-----------------------------------------------------------------------------
  subroutine rt_faceExchange(nAng)
    use raytrace_data, ONLY : faceValueAll, rt_needOnlyFaces
    use RadTrans_RayTrace_3DRT, ONLY : ib, ie, maxNrOfLeafBlocks
    use rt_data_raytrace_3drt, ONLY : nrOfAnglesPerGroup
    use Driver_data, ONLY : dr_globalNumProcs
    use Timers_interface, ONLY : Timers_start, Timers_stop
    implicit none
    integer, intent(in) :: nAng
    integer :: i, k, ierr
    integer, dimension(dr_globalNumProcs) :: displs, recvcounts

    if(.not. layoutNeedOnly) then
      ! legacy all-gather (identical to the former gather_faces routines)
      call Timers_start("rad_MPI_ALLGATHER_faceValues")
      if(nAng .eq. nrOfAnglesPerGroup) then
        call MPI_ALLGATHER(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
             faceValueAll, maxNrOfLeafBlocks*nrOfAnglesPerGroup, faceType, &
             MPI_COMM_WORLD, ierr)
      else if(nAng .ge. 1) then
        do i = 1, dr_globalNumProcs
          displs(i) = maxNrOfLeafBlocks*nrOfAnglesPerGroup*(i-1)
        end do
        recvcounts(:) = maxNrOfLeafBlocks*nAng
        call MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
             faceValueAll, recvcounts, displs, faceType, &
             MPI_COMM_WORLD, ierr)
      end if
      call Timers_stop("rad_MPI_ALLGATHER_faceValues")
      return
    end if

    call Timers_start("rad_faceExchange")
    if(allocated(sendBuf)) then
      if(size(sendBuf,5) < nSendSlots .or. size(sendBuf,5) > 2*nSendSlots+16) deallocate(sendBuf)
    end if
    if(.not. allocated(sendBuf)) &
      allocate(sendBuf(ib(1):ie(1), ib(2):ie(2), 3, 2, max(nSendSlots + nSendSlots/8 + 1, 1)))
    do k = 1, nSendSlots
      sendBuf(:,:,:,:,k) = faceValueAll(:,:,:,:,sendList(k))
    end do
    call MPI_Alltoallv(sendBuf, sendCnt, sendDsp, faceType, &
         faceValueAll(ib(1),ib(2),1,1,nOwnSlots+1), reqCnt, reqDsp, faceType, &
         MPI_COMM_WORLD, ierr)
    call Timers_stop("rad_faceExchange")
  end subroutine rt_faceExchange

end module rt_faceComm
