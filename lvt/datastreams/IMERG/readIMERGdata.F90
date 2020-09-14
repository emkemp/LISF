!-----------------------BEGIN NOTICE -- DO NOT EDIT----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V7.3
!-------------------------END NOTICE -- DO NOT EDIT----------------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readIMERGdata
! \label{readIMERGdata}
!
! !INTERFACE: 
subroutine readIMERGdata(source)
! 
! !USES:   
   use ESMF
#if (defined USE_HDF5)
   use HDF5
#endif
   use IMERG_dataMod
   use LVT_coreMod
   use LVT_logMod
   use LVT_timeMgrMod
   use LVT_histDataMod

   implicit none

!
! !INPUT PARAMETERS: 
   integer, intent(in)      :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  
!  This subroutine provides the data reader for 30-min, 0.1 deg IMERG HDF5 
!  files.  The routine reads in the precipitation estimates, and spatially 
!  interpolates to the LVT target grid.
! 
! !REVISION HISTORY: 
!  14 Dec 2018: Eric Kemp, Initial Specification
! 
!EOP

   integer :: nc,nr
   character*200                :: filename
   logical                      :: file_exists
   integer                      :: nunit,ufn,iret,ierr
   integer                      :: c,r
   integer                      :: ftn
   real                         :: &
        precip_f(imergdata(source)%nc*imergdata(source)%nr)
   integer :: &
        nprecip_f(imergdata(source)%nc*imergdata(source)%nr)
   real                         :: varfield(LVT_rc%lnc,LVT_rc%lnr)
   integer                      :: yr1, mo1, da1, hr1, mn1, ss1
   integer                      :: yr2, mo2, da2, hr2, mn2, ss2
   type(ESMF_Time)              :: time1
   type(ESMF_Time)              :: time2
   type(ESMF_TimeInterval)      :: lis_ts
   integer :: hdferr   
#if (defined USE_HDF5)
   integer(HID_T) :: file_id, dataset_id, datatype_id 
   integer(HSIZE_T) :: dims(3)
#endif
   real, allocatable :: tmp_precip_cal(:,:,:)
   logical :: fail
   integer :: status

   ! Before we proceed, make sure LVT was compiled with HDF5
#if (defined USE_HDF5)
#else
   write(LVT_logunit,*) &
        '[ERR] LVT was not compiled with HDF5 support for IMERG!'
   write(LVT_logunit,*) &
        'Reconfigure to use HDF5, recompile, and try again!'
   write(LVT_logunit,*) &
        'Or, avoid use of IMERG!'
   call LVT_endrun()
#endif

   ! Initialize variables
   nc = imergdata(source)%nc
   nr = imergdata(source)%nr
   precip_f(:)  = LVT_rc%udef
   nprecip_f(:) = 0
   varfield(:,:) = LVT_rc%udef

   ! Set the start time.  Note that IMERG data indicate starting time of
   ! precipitation, so we need to look back 30 minutes to get "present"
   ! precipitation.
   yr1 = LVT_rc%dyr(source)
   mo1 = LVT_rc%dmo(source)
   da1 = LVT_rc%dda(source)
   hr1 = LVT_rc%dhr(source)
   mn1 = LVT_rc%dmn(source)
   ss1 = 0
   call ESMF_TimeSet(time1,yy=yr1, mm=mo1, dd=da1, &
       h=hr1,m=mn1,s=ss1,calendar=LVT_calendar, rc=status)
   call LVT_verify(status)

   call ESMF_TimeIntervalSet(lis_ts, s = LVT_rc%ts, &
        rc=status)
   call LVT_verify(status)  
   time2 = time1 - lis_ts

   call ESMF_TimeGet(time2,yy=yr2, mm=mo2, dd=da2, &
        h=hr2,m=mn2,s=ss2,calendar=LVT_calendar, rc=status)
   call LVT_verify(status)

   call create_imergdata_filename(source, &
        yr2, mo2, da2, hr2, mn2, filename)

   inquire(file=trim(filename),exist=file_exists)
   if (.not. file_exists) then
      write(LVT_logunit,*)'[INFO] Cannot find file : ',trim(filename)
   end if

   if (file_exists) then
      ! Initialize IDs.  Useful later for error handling.
      file_id = -1
      dataset_id = -1
      datatype_id = -1
   
      ! Initialize HDF5 Fortran interface.
      call open_hdf5_f_interface(fail)
      if (fail) then 
         call LVT_endrun()
      end if

      ! Open the file
      call open_imerg_file(trim(filename),file_id,fail)
      if (fail) then
         call LVT_endrun()
      end if

      ! Open the calibrated precipitation
      call open_imerg_dataset(file_id,"/Grid/precipitationCal",dataset_id,&
           fail)
      if (fail) then
         call LVT_endrun()
      end if
       
      ! Get the datatype
      call get_imerg_datatype(dataset_id,datatype_id,fail)
      if (fail) then
         call LVT_endrun()
      end if

      ! Check the data type
      call check_imerg_type(datatype_id,H5T_IEEE_F32LE,fail)
      if (fail) then
         call LVT_endrun()
      end if

      ! Check the IMERG dimensions
      ! NOTE:  IMERG V05B dimensions are r,c
      ! NOTE:  IMERG V06B dimensions are r,c,1
      dims(1) = imergdata(source)%nr
      dims(2) = imergdata(source)%nc      
      dims(3) = 1
      call check_imerg_dims(dataset_id,3,dims,fail)
      if (fail) then
         call LVT_endrun()
      end if

      ! Check the units
      call check_imerg_units(dataset_id,"mm/hr",fail)
      if (fail) then
         call LVT_endrun()
      end if

      ! We are ready for the data.
      allocate(tmp_precip_cal(dims(1),dims(2),1))
      tmp_precip_cal = 0
      call h5dread_f(dataset_id, H5T_IEEE_F32LE, tmp_precip_cal, dims, hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot read IMERG precipitation!'
         call LVT_endrun()
      end if

      ! Clean up
      call close_imerg_datatype(datatype_id,fail)
      call close_imerg_dataset(dataset_id,fail)
      call close_imerg_file(file_id,fail)
      call close_hdf5_f_interface(fail)

      ! Convert mm/hr to kg/m2s.  Note we also change the r,c convention
      ! to equivalent c,r in 1D
      do c = 1,nc
         do r = 1,nr
            if (.not. tmp_precip_cal(r,c,1) .lt. 0) then
               precip_f(c + (r-1)*nc) = &
                    tmp_precip_cal(r,c,1)/3600
               nprecip_f(c + (r-1)*nc) = 1
            end if
         end do
      end do

      deallocate(tmp_precip_cal)
   end if


   ! Now interpolate the data.  First as rate and then as depth
   call interp_imergvar(source, nc, nr, precip_f, nprecip_f, varfield)
   call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip, source, varfield, &
        vlevel=1, units="kg/m2s")

   do r = 1,nr
      do c = 1,nc
         if (nprecip_f(c+(r-1)*nc) .gt. 0) then
            precip_f(c+(r-1)*nc) = precip_f(c+(r-1)*nc) * 1800 ! 30 min
         end if
      end do
   end do
   
   call interp_imergvar(source, nc, nr, precip_f, nprecip_f, varfield)
   call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip, source, varfield, &
        vlevel=1, units="kg/m2")

contains

#if (defined USE_HDF5)
   ! Internal subroutine.  Open the HDF5 Fortran interface
   subroutine open_hdf5_f_interface(fail)
      use HDF5
      use LVT_logMod, only: LVT_logunit
      implicit none
      logical,intent(out) :: fail
      integer :: hdferr
      fail = .false.
      call h5open_f(hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot initialize HDF5 ', &
              'Fortran interface!'
         fail = .true.
      end if
   end subroutine open_hdf5_f_interface

   ! Internal subroutine.  Open the IMERG HDF5 file
   subroutine open_imerg_file(filename,file_id,fail)
      use HDF5
      use LVT_logMod, only: LVT_logunit
      implicit none
      character(len=*), intent(in) :: filename
      integer(HID_T), intent(out) :: file_id
      logical, intent(out) :: fail
      integer :: hdferr
      fail = .false.
      call h5fopen_f(trim(filename),H5F_ACC_RDONLY_F,file_id,hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*)&
              '[ERR] Cannot open IMERG HDF5 file ', &
              trim(filename)
         fail = .true.
      else
         write(LVT_logunit,*) &
              '[INFO] Opened IMERG file ',trim(filename)
      end if
   end subroutine open_imerg_file
   
   ! Internal subroutine.  Open HDF5 dataset.
   subroutine open_imerg_dataset(file_id,dataset_name,dataset_id,fail)
      use HDF5
      use LVT_logMod, only: LVT_logunit
      implicit none
      integer(HID_T),intent(in) :: file_id
      character(len=*),intent(in) :: dataset_name
      integer(HID_T),intent(out) :: dataset_id
      logical, intent(out) :: fail
      fail = .false.
      call h5dopen_f(file_id,trim(dataset_name),dataset_id, hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*)&
              '[ERR] Cannot open IMERG HDF5 dataset ', &
              trim(dataset_name)
         fail = .true.
      end if
   end subroutine open_imerg_dataset
   
   ! Internal subroutine.  Sanity check IMERG precipitation units.
   subroutine check_imerg_units(dataset_id,units,fail)
      
      ! Imports
      use HDF5
      use ISO_C_BINDING
      use LVT_logMod, only: LVT_logunit
      
      ! Defaults
      implicit none
      
      ! Arguments
      integer(HID_T),intent(in) :: dataset_id
      character(len=*),intent(in) :: units
      logical,intent(out) :: fail
      
      ! Local variables
      integer(HID_T) :: attr_id, type_id, space_id, memtype_id
      integer :: hdferr
      integer(size_t) :: size
      integer(SIZE_T), parameter :: sdim = 5
      integer(HSIZE_T), dimension(1:1) :: dims = (/1/)
      integer(HSIZE_T), dimension(1:1) :: maxdims
      character(len=sdim), dimension(:), allocatable, target :: rdata
      type(C_PTR) :: f_ptr
      integer :: i
      
      fail = .false.
      
      ! Open the attribute
      call h5aopen_f(dataset_id,'Units',attr_id,hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot open IMERG HDF5 attribute'
         fail = .true.
         return
      end if
      
      ! Get the attribute datatype
      call h5aget_type_f(attr_id, type_id, hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot get IMERG HDF5 attribute datatype'
         call h5aclose_f(attr_id,hdferr)
         fail = .true.
         return
      end if
      
      ! Get the size of the attribute datatype, and sanity check.
      call h5tget_size_f(type_id,size,hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot get IMERG HDF5 attribute ', &
              'datatype size'
         call h5tclose_f(type_id,hdferr)
         call h5aclose_f(attr_id,hdferr)
         fail = .true.
         return
      end if
      if (size .gt. sdim+1) then
         write(LVT_logunit,*) &
              '[ERR] Expected smaller IMERG HDF5 attribute',&
              'datatype size'
         write(LVT_logunit,*)'Expected ',sdim+1
         write(LVT_logunit,*)'Found ',size
         call h5tclose_f(type_id,hdferr)
         call h5aclose_f(attr_id,hdferr)
         fail = .true.
         return
      end if
      
      ! Get the attribute dataspace
      call h5aget_space_f(attr_id, space_id, hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot get IMERG HDF5 attribute', &
              'dataspace'
         call h5tclose_f(type_id,hdferr)
         call h5aclose_f(attr_id,hdferr)
         fail = .true.
         return
      end if
      
      ! Get the dimensions of the dataspace
      call h5sget_simple_extent_dims_f(space_id,dims,maxdims,hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot get IMERG HDF5 attribute ', &
              'dataspace dimensions'
         call h5sclose_f(space_id,hdferr)
         call h5tclose_f(type_id,hdferr)
         call h5aclose_f(attr_id,hdferr)
         fail = .true.
         return
      end if
      ! Create the memory datatype
      call h5tcopy_f(H5T_FORTRAN_S1, memtype_id, hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot copy IMERG HDF5 attribute ', &
              'memory datatype.'
         call h5sclose_f(space_id,hdferr)
         call h5tclose_f(type_id,hdferr)
         call h5aclose_f(attr_id,hdferr)
         fail = .true.
         return
      end if
      call h5tset_size_f(memtype_id, sdim, hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot set IMERG HDF5 attribute ', &
              'memory datatype size.'
         call h5tclose_f(memtype_id,hdferr)
         call h5sclose_f(space_id,hdferr)
         call h5tclose_f(type_id,hdferr)
         call h5aclose_f(attr_id,hdferr)
         fail = .true.
         return
      end if
      
      ! Read the attribute
      allocate(rdata(1:dims(1)))
      f_ptr = C_LOC(rdata(1)(1:1))
      call h5aread_f(attr_id, memtype_id, f_ptr, hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot read IMERG HDF5 attribute.'
         deallocate(rdata)
         call h5tclose_f(memtype_id,hdferr)
         call h5sclose_f(space_id,hdferr)
         call h5tclose_f(type_id,hdferr)
         call h5aclose_f(attr_id,hdferr)
         fail = .true.
         return
      end if
      
      ! Check the units
      if (trim(rdata(1)) .ne. trim(units)) then
         write(LVT_logunit,*) &
              '[ERR] Found wrong IMERG precipitation', &
              'units'
         write(LVT_logunit,*) 'Expected mm/hr'
         write(LVT_logunit,*) 'Found ',trim(rdata(1))
         deallocate(rdata)
         call h5tclose_f(memtype_id,hdferr)
         call h5sclose_f(space_id,hdferr)
         call h5tclose_f(type_id,hdferr)
         call h5aclose_f(attr_id,hdferr)
         fail = .true.
         return
      end if

      ! Clean up
      deallocate(rdata)
      call h5tclose_f(memtype_id,hdferr)
      call h5sclose_f(space_id,hdferr)
      call h5tclose_f(type_id,hdferr)
      call h5aclose_f(attr_id,hdferr)
      
   end subroutine check_imerg_units

   ! Internal subroutine.  Get HDF5 datatype
   subroutine get_imerg_datatype(dataset_id,datatype_id,fail)
      use HDF5
      use LVT_logMod, only: LVT_logunit
      implicit none
      integer(HID_T),intent(in) :: dataset_id
      integer(HID_T),intent(out) :: datatype_id
      logical, intent(out) :: fail
      integer :: hdferr
      fail = .false.
      call h5dget_type_f(dataset_id, datatype_id, hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*)&
              '[ERR] Cannot determine IMERG HDF5 datatype'
         fail = .true.
      end if
   end subroutine get_imerg_datatype

   ! Internal function.  Check datatype 
   subroutine check_imerg_type(datatype_id, datatype, fail)
      use HDF5
      use LVT_logMod, only: LVT_logunit
      integer(HID_T), intent(in) :: datatype_id
      integer(HID_T), intent(in) :: datatype
      logical, intent(out) :: fail
      logical :: flag
      integer :: hdferr
      fail = .false.
      call h5tequal_f(datatype_id, datatype, flag, hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot confirm IMERG HDF5 datatype!'
         fail = .true.
         return
      end if
      if (.not. flag) then
         write(LVT_logunit,*)&
              '[ERR] IMERG HDF5 datatype is wrong type!'
         fail = .true.
         return
      end if
   end subroutine check_imerg_type
   
   ! Internal subroutine.  Check the rank/dimensions of dataset.
   subroutine check_imerg_dims(dataset_id,rank,dims,fail)
      
      ! Imports
      use HDF5
      use LVT_logMod, only: LVT_logunit
      
      ! Defaults
      implicit none
      
      ! Arguments
      integer(HID_T), intent(in) :: dataset_id
      integer, intent(in) :: rank
      integer(HSIZE_T), intent(in) :: dims(rank)
      logical, intent(out) :: fail
      
      ! Local variables
      integer(HID_T) :: dataspace_id
      integer :: hdferr, dataspace_rank
      integer(HSIZE_T), allocatable :: dataspace_dims(:)
      integer(HSIZE_T), allocatable :: dataspace_maxdims(:)
      logical :: flag
      integer :: i
      
      ! First, get the dataspace for the dataset
      call h5dget_space_f(dataset_id, dataspace_id, hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*)&
              '[ERR] Could not get IMERG HDF5 dataspace'
         fail = .true.
         return
      end if
      
      ! Sanity check:  Make sure this dataspace is "simple"
      call h5sis_simple_f(dataspace_id, flag, hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot determine if IMERG HDF5 ', &
              'dataspace is simple'
         fail = .true.
         return
      end if
      if (.not. flag) then
         write(LVT_logunit,*) &
              '[ERR] IMERG HDF5 dataspace is not simple'
         fail = .true.
         return
      end if
      
      ! Check the rank (number of dimensions)
      call h5sget_simple_extent_ndims_f(dataspace_id,dataspace_rank,hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*)&
              '[ERR] Cannot get rank of IMERG HDF5 dataspace '
         fail = .true.
         return
      end if
      if (dataspace_rank .ne. rank) then
         write(LVT_logunit,*)&
              '[ERR] Expected IMERG HDF5 precipitation rank ', rank
         write(LVT_logunit,*)'But found rank ',dataspace_rank
         fail = .true.
         return
      end if
      
      ! Check the dimensions
      allocate(dataspace_dims(rank))
      allocate(dataspace_maxdims(rank))
      call h5sget_simple_extent_dims_f(dataspace_id,dataspace_dims, &
           dataspace_maxdims,hdferr)
      if (hdferr .ne. rank) then
         write(LVT_logunit,*) &
              '[ERR] Cannot get dims for IMERG HDF5 dataspace'
         deallocate(dataspace_dims)
         deallocate(dataspace_maxdims)
         fail = .true.
         return
      end if
      
      do i = 1, rank
         if (dataspace_dims(i) .ne. dims(i)) then
            fail = .true.
            exit
         end if
      end do
      if (fail) then
         write(LVT_logunit,*) &
              '[ERR] Found bad dimensions for IMERG HDF5 ', &
              'dataspace'
         write(LVT_logunit,*)'Expected ',dims(:)
         write(LVT_logunit,*)'Found ',dataspace_dims(:)
         deallocate(dataspace_dims)
         deallocate(dataspace_maxdims)
         return
      end if

      ! Clean up
      deallocate(dataspace_dims)
      deallocate(dataspace_maxdims)
      call h5sclose_f(dataspace_id,hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot close IMERG HDF5 dataspace'
         fail = .true.
         return
      end if

   end subroutine check_imerg_dims

   ! Internal subroutine.  Close datatype
   subroutine close_imerg_datatype(datatype_id,fail)
      use HDF5
      use LVT_logMod, only: LVT_logunit
      integer(HID_T), intent(inout) :: datatype_id
      logical, intent(out) :: fail
      integer :: hdferr
      fail = .false.
      call h5tclose_f(datatype_id,hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot close IMERG HDF5 datatype '
         fail = .true.
      end if
      datatype_id = -1
   end subroutine close_imerg_datatype

   ! Internal function.  Close the dataset
   subroutine close_imerg_dataset(dataset_id,fail)
      use HDF5
      use LVT_logMod, only: LVT_logunit
      integer(HID_T), intent(inout) :: dataset_id
      logical, intent(out) :: fail
      integer :: hdferr
      fail = .false.
      call h5dclose_f(dataset_id,hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot close IMERG HDF5 dataset '
         fail = .true.
      end if
      dataset_id = -1
   end subroutine close_imerg_dataset

   ! Internal subroutine.  Close the IMERG HDF5 file.
   subroutine close_imerg_file(file_id,fail)
      use HDF5
      use LVT_logMod, only: LVT_logunit
      implicit none
      integer(HID_T), intent(inout) :: file_id
      logical, intent(out) :: fail
      integer :: hdferr
      fail = .false.
      call h5fclose_f(file_id,hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot close IMERG HDF5 file ', &
              trim(filename)
         fail = .true.
      end if
      file_id = -1
   end subroutine close_imerg_file

   ! Internal subroutine.  Close the HDF5 Fortran interface.
   subroutine close_hdf5_f_interface(fail)
      use HDF5
      use LVT_logMod, only: LVT_logunit
      implicit none
      logical, intent(out) :: fail
      integer :: hdferr
      fail = .false.
      call h5close_f(hdferr)
      if (hdferr .ne. 0) then
         write(LVT_logunit,*) &
              '[ERR] Cannot close HDF5 Fortran ',&
              'interface!'
         fail = .true.
      end if
   end subroutine close_hdf5_f_interface
   
#endif

end subroutine readIMERGdata

!------------------------------------------------------------------------------

subroutine create_imergdata_filename(source, yr, mo, da, hr, mn, filename)

   ! Import
   use ESMF
   use IMERG_dataMod
   use LVT_logMod
   use LVT_timeMgrMod, only: LVT_calendar

   ! Defaults
   implicit none

   ! Arguments
   integer,intent(in) :: source
   integer,intent(in) :: yr
   integer,intent(in) :: mo
   integer,intent(in) :: da
   integer,intent(in) :: hr
   integer,intent(in) :: mn
   character(len=*),intent(inout) :: filename
   
   ! Local variables
   integer :: tmp_yr, tmp_mo, tmp_da, tmp_hr, tmp_mn, tmp_ss
   integer :: tmp_minutes_in_day
   character(len=4) :: syr,sminutes_in_day
   character(len=2) :: smo,sda,shr,smn,sss
   type(ESMF_TIME) :: start_time, end_time, start_of_day
   type(ESMF_TIMEINTERVAL) :: half_hour
   type(ESMF_TIMEINTERVAL) :: time_diff
   integer :: rc

   ! Construct start of filename.
   write(unit=syr, fmt='(i4.4)') yr
   write(unit=smo, fmt='(i2.2)') mo
   filename = &
        trim(imergdata(source)%odir) &
        //"/"//trim(syr)//trim(smo)//"/"
   filename = trim(filename)//trim(imergdata(source)%product)
   filename = trim(filename)//".MS.MRG.3IMERG."

   ! Construct filename up through start date/time.
   write(unit=syr, fmt='(i4.4)') yr
   write(unit=smo, fmt='(i2.2)') mo
   write(unit=sda, fmt='(i2.2)') da
   write(unit=shr, fmt='(i2.2)') hr
   write(unit=smn, fmt='(i2.2)') mn
   write(unit=sss, fmt='(i2.2)') 0
   filename = trim(filename)//syr//smo//sda
   filename = trim(filename)//"-S"//shr//smn//sss

   ! Determine end date/time.
   call esmf_timeset(start_time, &
        yy=yr, &
        mm=mo, &
        dd=da, &
        h=hr, &
        m=mn, &
        s=0, &
        calendar = LVT_calendar, &
        rc = rc)
   call esmf_timeintervalset(half_hour,m=29,s=59,rc=rc)
   end_time = start_time + half_hour
   call esmf_timeget(end_time, &
        yy = tmp_yr, &
        mm = tmp_mo, &
        dd = tmp_da, &
        h  = tmp_hr, &
        m  = tmp_mn, &
        s  = tmp_ss, &
        rc = rc)

   ! Add end date/time to filename
   write(unit=syr, fmt='(i4.4)') tmp_yr
   write(unit=smo, fmt='(i2.2)') tmp_mo
   write(unit=sda, fmt='(i2.2)') tmp_da
   write(unit=shr, fmt='(i2.2)') tmp_hr
   write(unit=smn, fmt='(i2.2)') tmp_mn
   write(unit=sss, fmt='(i2.2)') tmp_ss
   filename = trim(filename)//"-E"//shr//smn//sss

   ! Determine number of minutes between start_of_day and start_time
   call esmf_timeset(start_of_day, &
        yy=yr, &
        mm=mo, &
        dd=da, &
        h=0, &
        m=0, &
        s=0, &
        calendar = LVT_calendar, &
        rc = rc)
   time_diff = start_time - start_of_day
   call esmf_timeintervalget(time_diff, m = tmp_minutes_in_day)

   ! Append minutes from start of month to filename
   write(unit=sminutes_in_day, fmt='(i4.4)') tmp_minutes_in_day
   filename = trim(filename)//"."//sminutes_in_day
   
   ! Finish filename construction
   if (trim(imergdata(source)%product) == "3B-HHR") then
      filename = trim(filename)//"."//trim(imergdata(source)%version)//".HDF5"
   else
      filename = trim(filename)//"."//trim(imergdata(source)%version)//".RT-H5"
   end if
end subroutine create_imergdata_filename

!------------------------------------------------------------------------------

subroutine interp_imergvar(source, nc, nr, var_input, nvar_input, var_output)

   ! Imports
   use LVT_coreMod,   only : LVT_rc, LVT_isAtAFinerResolution
   use IMERG_dataMod, only : imergdata

   ! Defaults
   implicit none

   ! Arguments
   integer,intent(in)          :: source
   integer,intent(in)          :: nc
   integer,intent(in)          :: nr
   real,intent(inout)          :: var_input(nc*nr)
   integer,intent(in)          :: nvar_input(nc*nr)
   real,intent(out)            :: var_output(LVT_rc%lnc, LVT_rc%lnr)

   ! Local variables
   logical*1          :: lb(nc*nr)
   integer            :: iret
   integer            :: c,r
   logical*1          :: lo(LVT_rc%lnc*LVT_rc%lnr)
   real               :: go(LVT_rc%lnc*LVT_rc%lnr)

   ! Initialize arrays
   var_output(:,:) = LVT_rc%udef
   lb(:) = .false.

  ! Average values and construct logical bit map
   do r=1,nr
     do c=1,nc
        if(nvar_input(c+(r-1)*nc).gt.0) then 
           var_input(c+(r-1)*nc) = &
                var_input(c+(r-1)*nc) / nvar_input(c+(r-1)*nc)
           lb(c+(r-1)*nc) = .true.
        else
           var_input(c+(r-1)*nc) = LVT_rc%udef
        endif
     enddo
  enddo

  ! Use budget interpolation if IMERG is at a coarser reoslution than
  ! the evaluation grid; otherwise, use upscale averaging.
  if (LVT_isAtAFinerResolution(imergdata(source)%datares)) then
     call conserv_interp(LVT_rc%gridDesc,lb,var_input, &
          lo,go, imergdata(source)%nc*imergdata(source)%nr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          imergdata(source)%rlat, &
          imergdata(source)%rlon, &
          imergdata(source)%w112, imergdata(source)%w122, &
          imergdata(source)%w212, imergdata(source)%w222, &
          imergdata(source)%n112, imergdata(source)%n122, &
          imergdata(source)%n212, imergdata(source)%n222, &
          LVT_rc%udef, iret)     
  else
     call upscaleByAveraging(&
          imergdata(source)%nc*imergdata(source)%nr, &
          LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
          imergdata(source)%n11, lb, &
          var_input, lo, go)
  end if

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        var_output(c,r) = go(c+(r-1)*LVT_rc%lnc)
     enddo
  enddo

end subroutine interp_imergvar



   
   
