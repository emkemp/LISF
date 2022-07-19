
!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: LIS_avoid_null_dereference
!  \label{LIS_avoid_null_dereference}
!
! !DESCRIPTION: Gracefully aborts to avoid null pointer dereference.
!
! !REVISION HISTORY:
!
! 19 Jul 2022: Eric Kemp; Initial implementation
!EOP

subroutine LIS_avoid_null_deref(stringname, length)

  ! Imports
  use LIS_coreMod, only: LIS_masterproc
  use LIS_logMod, only: LIS_abort, LIS_alert, LIS_logunit
  use LIS_mpiMod, only: LIS_MPI_COMM

  ! Defaults
  implicit none

  ! Arguments
  character(*), intent(in) :: stringname
  integer, intent(in) :: length

  ! Locals
  character(100) :: message(20)
  integer :: ierr

  ! Write to standard log
  write(LIS_logunit,*)'[ERR] Cannot find ', stringname(1:length), &
       ', aborting...'
  flush(LIS_logunit)

  ! Catch all MPI processes other than the master process here.  Only the
  ! master process should initiate the abort.
  if (.not. LIS_masterproc) then
#if (defined SPMD)
     call MPI_Barrier(LIS_MPI_COMM, ierr)
#endif
  endif

  message = ''
  message(1) = '[ERR] Program: LIS'
  message(2) = '  Cannot find '//stringname(1:length)//', aborting...'
  call LIS_alert('LIS.avoid_null_dereference', 1, message)
  call LIS_abort(message)
end subroutine LIS_avoid_null_deref
