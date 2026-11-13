!===============================================================================
!
! Module: rt_checksum
! Path: source/physics/RadTrans/RadTransMain/VETTAM/RayTrace/3DRT/rt_checksum.F90
!
!===============================================================================
!
! Debug diagnostic (runtime parameter rt_debugChecksum): exact, order-independent
! checksum of the unk fields written by a ray tracer. Called at the end of each
! ray-tracer call; the master rank writes one line to the ray-trace log (unit io).
!
! For every interior cell of every local leaf block the IEEE bit pattern of the
! double is split into its low and high 32-bit halves (as non-negative integers)
! which are summed separately in integer*8 (no overflow below 2**31 cells), then
! reduced over all ranks with MPI_SUM. The result does not depend on the order of
! summation or on the block distribution, and any single-bit difference changes it.
!
!===============================================================================

#include "Flash.h"
#include "constants.h"

module rt_checksum

  implicit none
  private
  public :: rt_checksumLog, rt_checksumAdd

  ! running number of checksum lines (same on all ranks; makes lines easy to align)
  integer, save :: callCount = 0

contains

  ! append unk index v with 4-character label nm to the list (vars, names, n)
  subroutine rt_checksumAdd(v, nm, n, vars, names)
    implicit none
    integer, intent(in)             :: v
    character(len=*), intent(in)    :: nm
    integer, intent(inout)          :: n
    integer, intent(inout)          :: vars(:)
    character(len=*), intent(inout) :: names(:)
    n = n + 1
    vars(n)  = v
    names(n) = nm
  end subroutine rt_checksumAdd

  subroutine rt_checksumLog(tracer, ident, nv, vars, names)
    use Grid_interface, ONLY : Grid_getListOfBlocks, Grid_getBlkPtr, Grid_releaseBlkPtr, &
         Grid_getBlkIndexLimits
    use Driver_data, ONLY : dr_globalMe, dr_nStep
    use RadTrans_hybridCharModule, ONLY : io
    implicit none
#include "Flash_mpi.h"
    character(len=*), intent(in) :: tracer, ident
    integer, intent(in)          :: nv
    integer, intent(in)          :: vars(nv)
    character(len=*), intent(in) :: names(nv)

    integer(8), parameter :: LO32 = 4294967295_8   ! 2**32-1
    integer(8) :: loc(2*nv+1), glob(2*nv+1), bits
    integer :: blockList(MAXBLOCKS), blockCount, b, v, i, j, k, ierr
    integer :: blkLimits(2,MDIM), blkLimitsGC(2,MDIM)
    real, dimension(:,:,:,:), pointer :: solnData
    character(len=2048) :: line
    character(len=64)   :: word

    callCount = callCount + 1
    loc = 0_8

    call Grid_getListOfBlocks(LEAF, blockList, blockCount)
    do b = 1, blockCount
      call Grid_getBlkIndexLimits(blockList(b), blkLimits, blkLimitsGC)
      call Grid_getBlkPtr(blockList(b), solnData)
      do v = 1, nv
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
          do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
            do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              bits = transfer(solnData(vars(v),i,j,k), bits)
              loc(2*v-1) = loc(2*v-1) + iand(bits, LO32)
              loc(2*v)   = loc(2*v)   + ishft(bits, -32)
            end do
          end do
        end do
      end do
      call Grid_releaseBlkPtr(blockList(b), solnData)
      loc(2*nv+1) = loc(2*nv+1) + int((blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1) * &
                                       (blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1) * &
                                       (blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1), 8)
    end do

    call MPI_Reduce(loc, glob, 2*nv+1, MPI_INTEGER8, MPI_SUM, MASTER_PE, MPI_COMM_WORLD, ierr)

    if(dr_globalMe == MASTER_PE) then
      write(line,'(a,i6,a,i8,5a,i0)') ' [checksum] call=', callCount, ' step=', dr_nStep, &
           ' tracer=', trim(tracer), ' ', trim(ident), ' ncells=', glob(2*nv+1)
      do v = 1, nv
        ! <name>=<sum of high halves>:<sum of low halves>
        write(word,'(1x,a,a,i0,a,i0)') trim(names(v)), '=', glob(2*v), ':', glob(2*v-1)
        line = trim(line)//trim(word)
      end do
      write(io,'(a)') trim(line)
      flush(io)
    end if
  end subroutine rt_checksumLog

end module rt_checksum
