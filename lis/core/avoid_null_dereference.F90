
subroutine avoid_null_dereference()
  use LIS_coreMod, only: LIS_masterproc
  use LIS_logMod, only: LIS_abort, LIS_alert, LIS_logunit
  use LIS_mpiMod, only: LIS_MPI_COMM
  implicit none
  character(100) :: message(20)
  integer :: ierr
  write(LIS_logunit,*)'[ERR] Avoiding null pointer dereferencing, aborting...'
  flush(LIS_logunit)
  message = ''
  message(1) = '[ERR] Program: LIS'
  message(2) = '  Avoiding null pointer deferencing, aborting...'
  if (LIS_masterproc) then
     call LIS_alert('LIS.avoid_null_dereference', 1, message)
     call LIS_abort(message)
  end if
#if (defined SPMD)
  call MPI_Barrier(LIS_MPI_COMM, ierr)
#endif
end subroutine avoid_null_dereference
