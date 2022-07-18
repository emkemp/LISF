
subroutine avoid_null_dereference()
  use LIS_coreMod, only: LIS_masterproc
  use LIS_logMod, only: LIS_abort, LIS_alert, LIS_logunit
  use LIS_mpiMod, only: LIS_MPI_COMM
  implicit none
  character(100) :: message(20)
  integer :: ierr
  write(LIS_logunit,*)'[ERR] Avoiding null pointer dereferencing, aborting...'
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
  message(2) = '  Avoiding null pointer deferencing, aborting...'
  call LIS_alert('LIS.avoid_null_dereference', 1, message)
  call LIS_abort(message)
end subroutine avoid_null_dereference
